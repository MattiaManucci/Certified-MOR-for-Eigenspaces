function [ERR_EST,Ared,pars] = subspace_SCMM(A,theta,thetap,bounds,options)
%% Reference: [1] M. Manucci, B. Stamm, and Z. Zeng, Certified Model Order Reduction for parametric Hermitian eigenproblems, 2025
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
if isfield(options,'tol_trunc')
    tol_trunc = options.tol_trunc;
else
    tol_trunc = tol*1e-2; 
end
if isfield(options,'flag_PC')
    flag_PC = options.flag_PC;
else
    flag_PC = 0; 
end

ff=[];
RES=[];
EEV=[];
kappa = length(A);
n=size(A{1},1);

dim = length(bounds.lb);
sp = issparse(A{1});
Ntrain=options.Nt;
pars.theta = theta;
pars.thetap = thetap;
pars.tol = tol*10^-1;
pars.RSG_tol = options.RSG_tol;

curerror = 10000;
opts.maxit=30000;

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
seed=123; rng(seed);
h = bounds.ub - bounds.lb;
mulist=zeros(dim,num_init_inter);
for j = 1:num_init_inter
    mulist(:,j) = bounds.lb + rand.*h;
end
pars.mu=mulist;
P = [];

ne = zeros(1,num_init_inter);
eiglist = cell(1,num_init_inter);
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

iter = num_init_inter+1;
AP=cell(kappa,1); PA1=cell(kappa,1); PA2{j}=cell(kappa,1);
for j = 1:kappa
    AP{j}=[]; PA1{j}=[]; PA2{j}=[];
end
for j = 1:kappa
    for jj = 1:kappa
        AP2{kappa*(j-1)+jj} = [];
    end
end
Ar = cell(kappa,1);
for j = 1:kappa
        Ar{j} = V_GAP'*(A{j}*V_GAP);
end
pars.AGAP = Ar;
%% Discreet Domain definition
% Cheb Points
nn=Ntrain-3;
i = 0:1:nn;
cx = [1,cos(((2*i + 1)/(2*(nn+1)))*pi),-1];
cx=(cx+1)/2;
mu_t=bounds.lb+cx.*h;

if dim>1
    mu_1=kron(mu_t(1,:),ones(1,Ntrain));
    mu_2=kron(ones(1,Ntrain),mu_t(2,:));
    mu_t=[mu_1;mu_2]; Ntrain=Ntrain^2;
end
Ntrain_bis= Ntrain; mu_t_bis=mu_t;
%% Evaluation of the spectral GAP in the grid points
GAP=zeros(Ntrain);
parfor ii=1:Ntrain
    mu=mu_t(:,ii);
    thetanew = theta(mu);
    PAP_GAP = thetanew(1)*pars.AGAP{1};
    for k = 2:kappa
        PAP_GAP = PAP_GAP + thetanew(k)*pars.AGAP{k};
    end
    [~,D] = eig(PAP_GAP); D=real(diag(D));
    [~,inds] = sort(D);
    for i=2:numel(inds)
        if abs(D(inds(i))-D(inds(1)))>options.RSG_tol
            GAP(ii) = D(inds(i))-D(inds(1));
            break
        end
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
        coef_Res_Base(abs(coef_Res_Base)<1e-14)=0;
    end
end
pars.Pres=vi; pars.coef_ResB=coef_Res_Base; pars.Nres=Nres;
%% MAIN LOOP FOR GREEDY SELECTION
while (curerror > tol)    
     
    % Project the problem
    ne(iter)=1;
    for j = 1:kappa
        if iter==(num_init_inter+1)
            AP{j} = P'*A{j}*P;
        else
            AP_off_diag=Pold'*(A{j}*Pext);
            AP{j} = [AP{j}, AP_off_diag;  AP_off_diag', Pext'*A{j}*Pext ];
        end
    end
    for j = 1:kappa
        if iter==(num_init_inter+1)
            PA1{j}=[PA1{j},A{j}*P];
            PA2{j}=[PA2{j};P'*A{j}];
        else
            PA1{j}=[PA1{j},A{j}*Pext];
            PA2{j}=[PA2{j};Pext'*A{j}];
        end
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
    app_err_EIGV=zeros(Ntrain,1);
    %% Evaluate the error estimate over Xi
    if flag_PC==1
        parfor ii=1:Ntrain
            [app_err_EIGV(ii)] = Error_Estimate_EigVec_Disc(mu_t(:,ii),GAP(ii),pars);
        end
    end
    if flag_PC==0
        for ii=1:Ntrain
            [app_err_EIGV(ii)] = Error_Estimate_EigVec_Disc(mu_t(:,ii),GAP(ii),pars);
        end
    end
    %% Determine the parameter that maximize the error estimate
    [f,ind]=max(app_err_EIGV); mu=mu_t(:,ind); mu_t(:,ind)=[]; app_err_EIGV(ind)=[];
    GAP_IT(iter)=GAP(ind); GAP(ind)=[]; Ntrain=Ntrain-1;
    curerror=f;
    fprintf('Iteration %g \nCurrent surrogate error is %g \n',iter,curerror);
    ff=[ff,curerror]; %Store the error estimate decay
    pars.tol = min(curerror*1e-1,1e-2); %Dynamically adjust the exit tolerance of EigOpt
    %% Compute the error estimate components
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
    [EEV(iter),RES(iter)] = EVALUATE_ERROR_ESTIMATE_C(mu,pars);
    %% Delate all parameters that verifies the exit condition
    ind_2=find(app_err_EIGV<tol_trunc); 
    mu_t(:,ind_2)=[]; GAP(ind_2)=[]; Ntrain=Ntrain-numel(ind_2);
    %% Updating the subspace
    if sp==1
        
        [~,D] = eigs(Amu,40*ne(iter)+1,'smallestreal',opts); 
        while (abs(D(1,1)-D(ne(iter)+1,ne(iter)+1)))<options.RSG_tol
            ne(iter)=ne(iter)+1;
        end
        [V,D] = eigs(Amu,20*ne(iter)+1,'smallestreal',opts);
        Pext = V(:,1:ne(iter));
        eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(1)- eigAj(i+1)))>options.RSG_tol 
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
    %% ----------------------------------------------------------------------------------------
    pars.eigvecs{iter} =  Pext;
    Pold=P;
    Pext2=Pext;
    for jj = 1:ne(iter)
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));  Pext2(:,jj) = Pext2(:,jj)/norm(Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));  Pext2(:,jj) = Pext2(:,jj)/norm(Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        if norm(Pext2(:,jj))>1e-15
             Pext2(:,jj)=Pext2(:,jj)/norm(Pext2(:,jj));
             P = [P Pext2(:,jj)];
        end
    end
    Pext=Pext2;
    mulist = [mulist mu];
    thetalist = [thetalist; thetanew];
    pars.mu = [pars.mu, mu];
    for j=1:(numel(pars.eigvecs)-1)
         pars.premult{j}=[pars.premult{j},pars.eigvecs{j}'*Pext];
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
    % This way of storing the residual coefficent is not optimal 
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

ERR_EST.ff=ff;
ERR_EST.EEV=EEV;
ERR_EST.RES=RES;
ERR_EST.GAP_IT=GAP_IT;
ERR_EST.EIG_E=EIG_E;

%% Check if eigenspace dimension is correctly captured at convergence
pars.lambda_1_red=zeros(Ntrain_bis,1); pars.eta=zeros(Ntrain_bis,1); pars.eta_epsilon=zeros(Ntrain_bis,1);
for ii=1:Ntrain_bis
    % Estimate EIG Error
    [output] = lambda_eta_eps_eig(mu_t_bis(:,ii),pars);
    pars.lambda_1_red(ii) = output.eigr;
    pars.eta(ii) = output.eta;
    pars.eta_epsilon(ii) = output.eta_epsilon;
end

cond_cho_eig = pars.eta-pars.lambda_1_red-pars.eta_epsilon;
l = find(cond_cho_eig<0); iii=1;
while isempty(l)==0
    mu=mu_t_bis(:,l(iii)); ne(iter)=1;
    %% Update the base, interpolation points and constraints
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    if sp==1
        
        [~,D] = eigs(Amu,20*(ne(iter)+1),'smallestreal',opts); 
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
           if abs((eigAj(1)- eigAj(i+1)))>options.RSG_tol 
               ne(iter)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(iter))); eiglist{iter} =  diag(D(inds(1:ne(iter)+1),inds(1:ne(iter)+1)));
        
    end

    pars.eigvecs{iter} =  Pext;
    Pext2=Pext;
    for jj = 1:ne(iter)
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));  Pext2(:,jj) = Pext2(:,jj)/norm(Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));  Pext2(:,jj) = Pext2(:,jj)/norm(Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        if norm(Pext2(:,jj))>1e-15
             Pext2(:,jj)=Pext2(:,jj)/norm(Pext2(:,jj));
             P = [P Pext2(:,jj)];
        end
    end
    Pext=Pext2;
    mulist = [mulist mu];
    thetalist = [thetalist; thetanew];
    pars.mu = [pars.mu, mu];
    for j=1:(numel(pars.eigvecs)-1)
         pars.premult{j}=[pars.premult{j},pars.eigvecs{j}'*Pext];
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
    % This way of storing the residual coefficent is not optimal 
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
    %% Update Areduced
    for j = 1:kappa
        AP{j} = P'*A{j}*P;
    end
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
    %% Check again over l
    for j = 1:numel(l)
        % Estimate EIG Error
        [output] = lambda_eta_eps_eig(mu_t_bis(:,l(j)),pars);
        pars.lambda_1_red(l(j)) = output.eigr;
        pars.eta(l(j)) = output.eta;
        pars.eta_epsilon(l(j)) = output.eta_epsilon;
    end
    cond_cho_eig = pars.eta-pars.lambda_1_red-pars.eta_epsilon;
    l = find(cond_cho_eig<0);
end
Ared = AP;
return
