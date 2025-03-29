function [ERR_EST,Ared,pars] = App_SEVES(A,theta,thetap,bounds,options)
%% For EIGOPT credits goes to: [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim, SIMAX 2014
%% The rest of the function:   [1] M. Manucci, B. Stamm, and Z. Zeng, Certified Model Order Reduction for parametric Hermitian eigenproblems, 2025
%% Some useful inizializations
if isfield(options,'num_init_inter')
    num_init_inter = options.num_init_inter;
else
    num_init_inter = 1;
end
if isfield(options,'EigOptMaxIt')
    pars.itertol = options.EigOptMaxIt; 
else
    pars.itertol = 2000;
end
if isfield(options,'Rel_Error')
    pars.Rel_Error = options.Rel_Error;
else
    pars.Rel_Error = 1;
end
if isfield(options,'Rel_Error')
    pars.Rel_Error = options.Rel_Error;
else
    pars.Rel_Error = 1;
end
if isfield(options,'tol')
    tol = options.tol;
else
    tol = 1e-4;
end
if isfield(options,'RSG_tol')
    options.RSG_tol = options.RSG_tol;
else
    options.RSG_tol = 1e-6;
end
if isfield(options,'gamma')
    pars.gamma = options.gamma;
else
    pars.gamma = -4e5; 
end
if isfield(options,'Vin')
    V_in = options.Vin;
    V_in = orth(V_in);
    flag_In_SUB=1;
else
    flag_In_SUB = 0; 
end



ff=[];
RES=[];
EEV=[];
EIG_ERR=[];
ff=[];
kappa = length(A);
n=size(A{1},1);
pars.RSG_tol = options.RSG_tol;
dim = length(bounds.lb);
sp = issparse(A{1});

pars.theta = theta;
pars.thetap = thetap;

pars.tol = tol*10^-1;
pars.minmax = 1;

curerror = 10000;
opts.maxit=900000;

V_GAP=options.Space_Gap;

%% Computing the box-constraints for the LP
if sp==1
    pars.lambounds = [];
    for j = 1:kappa
        Bmax = eigs(A{j},1,'largestreal',opts);
        Bmin = eigs(A{j},1,'smallestreal',opts);
        pars.lambounds = [pars.lambounds; Bmin Bmax];
    end
else    
    pars.lambounds = [];
    for j = 1:kappa
        D = eig(A{j});        
        pars.lambounds = [pars.lambounds; min(D) max(D)];
    end
end
pars.options = optimoptions('linprog','Display','none');
%% Loop to compute the starting subspace
if flag_In_SUB==0
    seed=123; rng(seed);
    h = bounds.ub - bounds.lb;
    for j = 1:num_init_inter
        mulist(:,j) = bounds.lb + rand.*h;
    end
    pars.mu=mulist;
    P = [];
    eiglist = [];

    for j = 1:num_init_inter

        ne(j)=1;
        mu = mulist(:,j);
        thetanew = theta(mu);
        thetalist(j,:) = thetanew;
        Amu = thetanew(1)*A{1};

        for k = 2:kappa
            Amu = Amu + thetanew(k)*A{k};
        end

        if sp
            [V,D] = eigs(Amu,3*ne(j)+1,'smallestreal',opts);
            while (abs(D(1,1)-D(ne(j),ne(j))))<options.RSG_tol
                ne(j)=ne(j)+1;
                [V,D] = eigs(Amu,3*ne(j)+1,'smallestreal',opts);
            end
            ne(j)=ne(j)-1;
            Pext = V(:,1:ne(j));
            eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];

        else
            [V,D] = eig(Amu);
            [eigAj,inds] = sort(diag(D));
            for i=1:n
                if abs((eigAj(i)- eigAj(i+1)))/abs(eigAj(i))>options.RSG_tol %Check if expression is corrected
                    ne(j)=i;
                    break
                end
            end
            Pext = V(:,inds(1:ne(j))); eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];
        end
        P = [P Pext];
        pars.eigvecs{j} = Pext;
    end

    [P,~] = qr(P,0);
    for j=1:numel(pars.eigvecs)
        pars.premult{j}=pars.eigvecs{j}'*P;
    end
else 
    thetalist= options.thetalist; 
    mulist= options.mu_in(:,1:end-1);
    pars.mu=mulist;
    eiglist=options.eiglist;
    P=V_in;
    q(1)= options.ne(1);
    pars.premult{1}=options.Vin(:,1:q(1))'*P;
    ne(1)=options.ne(1);
    pars.eigvecs{1}=options.Vin(:,1:q(1));
    for num_init_inter=2:(numel(options.ne)-1)
        q(num_init_inter) = q(num_init_inter-1)+options.ne(num_init_inter);
        pars.premult{num_init_inter}=options.Vin(:,q(num_init_inter-1)+1:q(num_init_inter))'*P;
        pars.eigvecs{num_init_inter}=options.Vin(:,q(num_init_inter-1)+1:q(num_init_inter));
        ne(num_init_inter)=options.ne(num_init_inter);
    end
end
Pext=P;
iter = num_init_inter+1;
Pold=[];
AP=cell(kappa,1); PA1=cell(kappa,1); PA2{j}=cell(kappa,1);
for j = 1:kappa
    AP{j}=[]; PA1{j}=[]; PA2{j}=[];
end
for j = 1:kappa
    for jj = 1:kappa
        AP2{kappa*(j-1)+jj} = [];
    end
end
%% Creating the orthonormal base for efficent and stable evaluation of the residual norm
tol_newnorm=1-1e-10;
Ares=A; Ares{end+1}=speye(n); 
Nres=size(P,2)*(kappa+1); vi=zeros(n,Nres);
for i=1:size(P,2)
    for j=1:(kappa+1)
       vi(:,(i-1)*(kappa+1)+j)=Ares{j}*P(:,i);
    end
end
vi(:,1)=vi(:,1)/norm(vi(:,1));
for j=2:Nres
    vi(:,j)=vi(:,j)/norm(vi(:,j)); newnorm=0;
    while newnorm<tol_newnorm
        for i=1:(j-1)
            vi(:,j) = vi(:,j) - (vi(:,j)'*vi(:,i))*vi(:,i);
        end
        newnorm = norm(vi(:,j));
        if newnorm<1e-14
            vi(:,j) = zeros(n,1);
            break
        else
            vi(:,j) = vi(:,j)/newnorm;
        end
    end
end
coef_Res_Base=zeros(Nres,Nres);
for i=1:size(P,2)
    for j=1:(kappa+1)
        coef_Res_Base(:,(i-1)*(kappa+1)+j) = vi'*(Ares{j}*P(:,i));
        coef_Res_Base(abs(coef_Res_Base)<1e-5)=0;
    end
end
pars.Pres=vi; pars.coef_ResB=coef_Res_Base; pars.Nres=Nres;
%% Matrix for the approximated gap
for j = 1:kappa
        Ar{j} = V_GAP'*(A{j}*V_GAP);
end
pars.AGAP = Ar;
%% MAIN LOOP FOR GREEDY SELECTION
while (curerror > tol)    
     
    % Project the problem
    ne(iter)=1;
    for j = 1:kappa
        AP{j} = P'*(A{j}*P);
    end
    % for j = 1:kappa
    %     PA1=A{j}*P;
    %     for jj=j:kappa
    %         PA2=P'*(A{jj})';
    %         AP2{kappa*(j-1)+jj} = PA2*PA1;
    %     end
    % end

    % for j = 2:kappa
    %     for jj=1:(j-1)
    %         AP2{kappa*(j-1)+jj} = AP2{kappa*(jj-1)+j};
    %     end
    % end


   
    for j = 1:kappa
        PA1{j}=A{j}*P;
        PA2{j}=P'*A{j};
    end
   
    for j=1:kappa
        for jj=j:kappa
            AP2{kappa*(j-1)+jj} = PA2{jj}*PA1{j};
        end
    end


    for j = 1:kappa
        for jj=1:(j-1)
            AP2{kappa*(j-1)+jj} = AP2{kappa*(jj-1)+j}';
        end
    end
    

    pars.A = AP;
    pars.Afull = AP2;
    pars.AP = PA1;
    pars.ne = ne;
    pars.thetalist = thetalist;
    pars.eiglist = eiglist;
    pars.P = P;
    [curerror,mu,parsout]=eigopt('EVALUATE_ERROR_ESTIMATE',bounds,pars); 
    fprintf('EigOpt at iteration %d required %d function evaluations\n',iter,parsout.nfevals);
    fprintf('Current surrogate error is %g \n',curerror);
    ff=[ff,curerror];
    pars.tol = min(curerror*1e-1,1e-2); %Dynamically adjust the exit tolerance of EigOpt
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    AREDmu = thetanew(1)*AP{1};
    Amur  = thetanew(1)*Ar{1};
    for k = 2:kappa
         AREDmu  =  AREDmu  + thetanew(k)*AP{k};
         Amur  =  Amur  + thetanew(k)*Ar{k};
    end
    [EV_EE(iter),RES(iter), GAP(iter)] = EVALUATE_ERROR_ESTIMATE_C(mu,pars);

    %% Updating the subspace
    if sp==1
        
        [~,D] = eigs(Amu,100*ne(iter)+1,'smallestreal',opts); %Check what to do when isnan
        while (abs(D(1,1)-D(ne(iter)+1,ne(iter)+1)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
        end
        [V,D] = eigs(Amu,5*ne(iter)+1,'smallestreal',opts);
        Pext = V(:,1:ne(iter));
        eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(1)- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(iter)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(iter))); eiglist{iter} =  diag(D(inds(1:ne(iter)+1),inds(1:ne(iter)+1)));
        
    end
    %% Compute error at that iteration
    [Vr,Dr] = eig(AREDmu);
    Dr=real(Dr);
    [~,inds_r] = sort(diag(Dr));
    nr=size(Pext,2);
    EIG_E(iter) = norm(P*Vr(:,inds_r(1:nr))-Pext*(Pext'*P*Vr(:,inds_r(1:nr))));
    %% ------------------------------------------------------------------------------------
    pars.eigvecs{iter} =  Pext;
    Pold=P;
    
    Pext2=Pext;
    for jj = 1:ne(iter)
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        if norm(Pext2(:,jj))>eps
             Pext2(:,jj)=Pext2(:,jj)/norm(Pext2(:,jj));
             P = [P Pext2(:,jj)];
        end
    end
    

    mulist = [mulist mu];
    thetalist = [thetalist; thetanew];
    pars.mu = [pars.mu, mu];
    for j=1:(numel(pars.eigvecs)-1)
        pars.premult{j}=[pars.premult{j},pars.eigvecs{j}'*Pext];
        %pars.premult{j}=pars.eigvecs{j}'*P;
    end
    pars.premult{numel(pars.eigvecs)}=pars.eigvecs{numel(pars.eigvecs)}'*P;
    iter = iter+1;

    %% Updating Residual Base
    Nres_new=size(P,2)*(kappa+1); nRESnew=size(Pext2,2);
    for i=1:nRESnew
        for j=1:(kappa+1)
            vi(:,Nres+j+(kappa+1)*(i-1))=Ares{j}*Pext2(:,i);
        end
    end
    for j=(Nres+1):Nres_new
        vi(:,j)=vi(:,j)/norm(vi(:,j)); newnorm=0;
        while newnorm<tol_newnorm
            for i=1:(j-1)
                vi(:,j) = vi(:,j) - (vi(:,j)'*vi(:,i))*vi(:,i);
            end
            newnorm = norm(vi(:,j));
            if newnorm<1e-14
                vi(:,j) = zeros(n,1);
                newnorm=1;
            else
                vi(:,j) = vi(:,j)/newnorm;
            end
        end
    end
    % Here one could also implement this in a more efficent way
    coef_Res_Base=zeros(Nres_new,Nres_new);
    for i=1:size(P,2)
        for j=1:(kappa+1)
           if numel(coef_Res_Base(:,j+(kappa+1)*(i-1)))==numel(vi'*(Ares{j}*P(:,i)))
               coef_Res_Base(:,j+(kappa+1)*(i-1)) = vi'*(Ares{j}*P(:,i));
           else
               disp(coef_Res_Base(:,j+(kappa+1)*(i-1)))
           end
           coef_Res_Base(abs(coef_Res_Base)<1e-15)=0;
        end
    end
    Nres=Nres_new;
    pars.Pres=vi; pars.coef_ResB=coef_Res_Base; pars.Nres=Nres;
end

Ared = AP;
ERR_EST.ff=ff;
ERR_EST.EV_EE=EV_EE;
ERR_EST.RES=RES;
ERR_EST.GAP=GAP;
ERR_EST.EIG_E=EIG_E;

return
