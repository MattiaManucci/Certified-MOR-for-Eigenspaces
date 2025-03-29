function [ff,Ared,pars] = approx_sGAP_Dset(A,theta,thetap,bounds,options)
%% Some useful inizializations
if isfield(options,'num_init_inter')
    num_init_inter = options.num_init_inter;
else
    num_init_inter = 2;
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
if isfield(options,'Nt')
    Ntrain = options.Nt;
else
    Ntrain = 50; 
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
kappa = length(A);
n=size(A{1},1);

dim = length(bounds.lb);
sp = issparse(A{1});

pars.RSG_tol = options.RSG_tol;
pars.theta = theta;
pars.thetap = thetap;

pars.tol = tol*10^-1;
pars.minmax = 1;

curerror = 10000;
opts.maxit=30000;

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
    
    ne(j)=2;
    mu = mulist(:,j);
    thetanew = theta(mu);    
    thetalist(j,:) = thetanew;
    Amu = thetanew(1)*A{1};

    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};    
    end

    if sp
        [V,D] = eigs(Amu,10*ne(j)+1,'smallestreal',opts);
        while (abs(D(1,1)-D(ne(j),ne(j))))<options.RSG_tol
            ne(j)=ne(j)+1;
            [V,D] = eigs(Amu,10*ne(j)+1,'smallestreal',opts);
        end
        E2=ne(j)+1;
        while (abs(D(E2,E2)-D(ne(j)+1,ne(j)+1)))<options.RSG_tol
            ne(j)=ne(j)+1;
            [V,D] = eigs(Amu,10*ne(j)+1,'smallestreal',opts);
        end

        Pext = V(:,1:ne(j));
        eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];
       
    else
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(1)- eigAj(ne(j))))>options.RSG_tol 
               ne(j)=i+1;
               break
           end
        end
        E2=ne(j)+1;
        for i=1:n
           if abs((eigAj(E2)- eigAj(ne(j)+1)))>options.RSG_tol 
               ne(j)=i+1;
               break
           end
        end
        Pext = V(:,inds(1:ne(j)));  
        eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];
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
%% Discreet grid
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
%% MAIN LOOP FOR GREEDY SELECTION
ff2=[]; mulist2=[];
while (curerror > tol)    
    ne(iter)=2;
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
    pars.ne = ne;
    pars.thetalist = thetalist;
    pars.eiglist = eiglist;
    pars.P = P;
    %% Loop over the discreet set
    app_err_sGAP=zeros(Ntrain,1);
    if flag_PC==1
        parfor ii=1:Ntrain
            [app_err_sGAP(ii)] = approximation_SG_d(mu_t(:,ii),pars);
        end
    end
    if flag_PC==0
        for ii=1:Ntrain
            [app_err_sGAP(ii)] = approximation_SG_d(mu_t(:,ii),pars);
        end
    end
    %% Compute the maximal value of the surrogate error
    [f,ind]=max(app_err_sGAP); mu=mu_t(:,ind);
    curerror=f;
    thetanew = theta(mu_t(:,ind));
    mu_t(:,ind)=[]; app_err_sGAP(ind)=[]; Ntrain=Ntrain-1;
    APmu = thetanew(1)*AP{1};
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu  = Amu  + thetanew(k)*A{k};
        APmu = APmu + thetanew(k)*AP{k};
    end
    % Compute the reduced spectral GAP
    [~,D] = eig(APmu); D=diag(D); [~,inds] = sort(D);
    for i=2:numel(inds)
        if abs(D(inds(i))-D(inds(1)))> options.RSG_tol 
            GAP_RED = D(inds(i))-D(inds(1));
            break
        end
    end
    fprintf('Itaration: %g\n Current surrogate error is %g \n',iter,curerror);
    ff=[ff,curerror];
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    %% Delate all parameters that verifies the exit condition
    ind_2=find(app_err_sGAP<tol_trunc); 
    mu_t(:,ind_2)=[]; Ntrain=Ntrain-numel(ind_2);
    %% Updating the subspace
    if sp==1
        
        [~,D] = eigs(Amu,20*ne(iter)+1,'smallestreal',opts);
        for jj=2:size(D,2)
            if abs(D(jj,jj)-D(1,1))> options.RSG_tol
                ne(iter)=jj-1;
                GAP = D(jj,jj)-D(1,1);
                break
            end
        end

        for ii=(ne(iter)+2):size(D,2)
            if abs(D((ii),(ii))-D((ne(iter)+1),(ne(iter)+1)))> options.RSG_tol
                ne(iter)=ii-1;
                break
            end
        end
        ne(iter)=ne(iter);
        [V,D] = eigs(Amu,20*(ne(iter)+1),'smallestreal',opts);
        Pext = V(:,1:ne(iter));
        eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))>options.RSG_tol 
               ne(iter)=i+1;
               GAP = eigAj(i)- eigAj(i+1);
               break
           end
        end
        for i=ne(iter):n
           if abs((eigAj(ne(iter))- eigAj(i+1)))>options.RSG_tol
               ne(iter)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(iter))); eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
        
    end

    ERR_GAP(iter)= abs(GAP-GAP_RED)/GAP_RED;
    
    pars.eigvecs{iter} =  Pext;
    Pold=P;
    Pext2=Pext;
    for jj = 1:ne(iter)
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj)); Pext2(:,jj)=Pext2(:,jj)/norm(Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj)); Pext2(:,jj)=Pext2(:,jj)/norm(Pext2(:,jj));
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
end
pars.ERR=ERR_GAP;
pars.ff2=ff2; pars.mulist2=mulist2;
Ared = AP;
%% Check subspaces condtion
pars.lambda_1_red=zeros(Ntrain_bis,1); pars.eta=zeros(Ntrain_bis,1); pars.eta_epsilon=zeros(Ntrain_bis,1);
for ii=1:Ntrain_bis
    % Estimate EIG Error
    [output] = lambda_eta_eps(mu_t_bis(:,ii),pars);
    pars.lambda_1_red(ii) = output.eigr;
    pars.eta(ii) = output.eta;
    pars.eta_epsilon(ii) = output.eta_epsilon;
end
cond_cho_GAP = pars.eta-pars.lambda_1_red-pars.eta_epsilon;
l = find(cond_cho_GAP<0);
% Lines of code to enforce the condition needs to be added here (for the
% Eigenspace those lines are added), if l is an empty variable, as for this test problems,
% then those lines are not necessary
return
