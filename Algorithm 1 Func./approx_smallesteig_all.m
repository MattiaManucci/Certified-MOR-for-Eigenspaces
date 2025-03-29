function [ff,Ared,pars] = approx_smallesteig_all(A,theta,thetap,bounds,options)
%% For EIGOPT credits goes to: [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim, SIMAX 2014
%% The rest of the function:   [1] M. Manucci, B. Stamm, and Z. Zeng, Certified Model Order Reduction for parametric Hermitian eigenproblems, 2025
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
opts.maxit = 30000;

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
for j = 1:num_init_inter
    mulist(:,j) = bounds.lb + rand.*h;
end
pars.mu=mulist;
P = [];
eiglist = [];
P_1 = [];

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
        i=1;
        [V,D] = eigs(Amu,ne(j)+1,'smallestreal',opts);
        while (abs(D(1,1)-D(ne(j),ne(j))))<options.RSG_tol
            ne(j)=ne(j)+1;
            [V,D] = eigs(Amu,ne(j)+1,'smallestreal',opts);
        end
        Pext = V(:,1:ne(j)); 
        eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(j)=i+1;
               break
           end
        end
        Pext = V(:,inds(1:ne(j)));  eiglist{j} = [diag(D(1:ne(j)+1,1:ne(j)+1))];
    end
    P = [P Pext]; P_1 = [P_1, Pext ]; ne1(j)=ne(j);
    pars.eigvecs{j} = Pext;
end

[P,~] = qr(P,0);
for j=1:numel(pars.eigvecs)
    pars.premult{j}=pars.eigvecs{j}'*P;
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
%% MAIN LOOP FOR GREEDY SELECTION
Pold=P;
while (curerror > tol)    
    ne(iter)=2;
    for j = 1:kappa
        % if iter==(num_init_inter+1)
        %     AP_off_diag=[];
        % else
        %     AP_off_diag=Pold'*(A{j}*Pext);
        % end
        % AP{j} = [AP{j}, AP_off_diag;  AP_off_diag', Pext'*A{j}*Pext ];
        AP{j} = P'*A{j}*P;
    end
    for j = 1:kappa
        % PA1{j}=[PA1{j},A{j}*Pext];
        % PA2{j}=[PA2{j};Pext'*A{j}];
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
    pars.ne = ne;
    pars.thetalist = thetalist;
    pars.eiglist = eiglist;
    pars.P = P;
    pars.lb = bounds.lb;
    pars.ub = bounds.ub;
    [curerror,mu,parsout]=eigopt('approximation_SG',bounds,pars); 
    fprintf('EigOpt at iteration %d required %d function evaluations\n',iter,parsout.nfevals);
    fprintf('Current surrogate error is %g \n',curerror);
       
    ff=[ff,curerror];
    pars.tol = min(curerror*1e-1,1e-2); %Dynamically adjust the exit tolerance of EigOpt
    thetanew = theta(mu);
    Amu = thetanew(1)*A{1};
    for k = 2:kappa
        Amu = Amu + thetanew(k)*A{k};
    end
    %% Updating the subspace
    if sp==1
        
        [~,D] = eigs(Amu,20*ne(iter)+1,'smallestreal',opts);
        for jj=2:size(D,2)
            if abs(D(jj,jj)-D(1,1))> options.RSG_tol
                ne(iter)=jj-1;
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
        [V,D] = eigs(Amu,5*(ne(iter)+1),'smallestreal',opts);
        Pext = V(:,1:ne(iter));
        P_1 = [P_1,V(:,1:(jj-1))];
        ne1(iter) = jj-1;
        eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
       
    else
        
        [V,D] = eig(Amu);
        [eigAj,inds] = sort(diag(D));
        for i=1:n
           if abs((eigAj(i)- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(iter)=i+1;
               break
           end
        end
        for i=ne(iter):n
           if abs((eigAj(ne(iter))- eigAj(i+1)))>options.RSG_tol %Check if expression is corrected
               ne(iter)=i;
               break
           end
        end
        Pext = V(:,inds(1:ne(iter))); 
        P_1 = [P_1, V(:,inds(1:i+1))]; 
        ne1(iter) = i+1;
        eiglist{iter} = diag(D(1:ne(iter)+1,1:ne(iter)+1));
        
    end
    
    pars.eigvecs{iter} =  Pext;
    
    Pext2=Pext;
    for jj = 1:ne(iter)
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        Pext2(:,jj) = Pext2(:,jj) - P*(P'*Pext2(:,jj));
        if norm(Pext2(:,jj))>1e-14
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
end

pars.ne1=ne1;
pars.P_1=P_1;
Ared = AP;

return
