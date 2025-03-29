function [mu,allCombinations, OUTPUT] = plot_lambdamin_MD(A,Ared,P,theta,bounds,opts,pars)


if isfield(opts,'RSG_tol')
    RSG_tol = opts.RSG_tol;
else
    RSG_tol = 0;
end
if isfield(opts,'mu')
    mu_given = opts.mu;
else
     mu_given = [];
end
if isfield(opts,'flag_PC')
    flag_PC = opts.flag_PC;
else
    flag_PC = 0;
end
Nmu=opts.Nmu;
kappa=numel(A);
p=numel(bounds.ub);
mu=cell(p,1);
opt.maxit=3000;
I=speye(size(A{1},1));

for i=1:p
    mu{i}=zeros(p,Nmu(i));

    if opts.log_space(i)==1
        mu{i}=logspace(log10(bounds.lb(i)-bounds.lb(i)+1),log10(bounds.ub(i)-bounds.lb(i)+1),Nmu(i));
        mu{i}=mu{i}+(bounds.lb(i))-1;
    else
        nn=Nmu(i)-3;
        j = 0:1:nn;
        cx = [1,cos(((2*j + 1)/(2*(nn+1)))*pi),-1];
        cx=(cx+1)/2;
        mu{i}=bounds.lb(i)+cx.*(bounds.ub(i)-bounds.lb(i));
    end
end

indices = cell(1, p);
[indices{:}] = ndgrid(1:Nmu(i));
allCombinations = cell2mat(cellfun(@(x) x(:), indices, 'UniformOutput', false));
numCombinations = size(allCombinations, 1);
Eig_Val = zeros(numCombinations,1);
EV_EE = zeros(numCombinations,1);

if flag_PC==1
    parfor i=1:numCombinations
        currentCombination = allCombinations(i, :);
        v=[];
        for j=1:p
            v=[v;mu{j}(currentCombination(j))];
        end

        thetanew = theta(v);
        Amu=thetanew(1)*A{1};
        Amu_RED=thetanew(1)*Ared{1};
        for k=2:kappa
            Amu=Amu+thetanew(k)*A{k};
            Amu_RED=Amu_RED+thetanew(k)*Ared{k};
        end
        [Vred,Dred] = eig(Amu_RED); Dred=real(Dred);
        [~,ind_red]=sort(diag(Dred));
        for jj=2:size(Dred,2)
            if (Dred(ind_red(jj),ind_red(jj))-Dred(ind_red(1),ind_red(1)))> RSG_tol
                ind=jj-1;
                break
            end
        end
        if jj==size(Dred,2)
            ind=size(Amu,2);
            fprintf('For parameter %d we have the matrix to have all the same eigenvalue, therefore the ground state is all the space\n Error is set to 0\n',v)
            Eig_Val(i)=eps;
            ind=[];
        else
            [V,~] = eigs(Amu,20*ind,'smallestreal',opt);
            g = V(:,1:ind)'*P*Vred(:,ind_red(1:ind));
            g = V(:,1:ind)*g;
            Eig_Val(i) = norm(I*P*Vred(:,ind_red(1:ind))-g);
            ind=[];
        end
     % Evaluation of the error estimate (need for efficency index comparison)
     [EV_EE(i),RES, GAP] = EVALUATE_ERROR_ESTIMATE_C(v,pars);
     EV_EE(i) = (RES+EV_EE(i))/GAP;
    end
end
if flag_PC==0
    for i=1:numCombinations
        currentCombination = allCombinations(i, :);
        v=[];
        for j=1:p
            v=[v;mu{j}(currentCombination(j))];
        end

        thetanew = theta(v);
        Amu=thetanew(1)*A{1};
        Amu_RED=thetanew(1)*Ared{1};
        for k=2:kappa
            Amu=Amu+thetanew(k)*A{k};
            Amu_RED=Amu_RED+thetanew(k)*Ared{k};
        end
        tic
        [Vred,Dred] = eig(Amu_RED); Dred=real(Dred);
        time_red=toc
        [~,ind_red]=sort(diag(Dred));
        for jj=2:size(Dred,2)
            if (Dred(ind_red(jj),ind_red(jj))-Dred(ind_red(1),ind_red(1)))> RSG_tol
                ind=jj-1;
                break
            end
        end
        if jj==size(Dred,2)
            ind=size(Amu,2);
            fprintf('For parameter %d we have the matrix to have all the same eigenvalue, therefore the ground state is all the space\n Error is set to 0\n',v)
            Eig_Val(i)=eps;
        else
            tic
            [V,~] = eigs(Amu,20*ind,'smallestreal',opt);
            time_FOM=toc
            g = V(:,1:ind)'*P*Vred(:,ind_red(1:ind));
            g = V(:,1:ind)*g;
            Eig_Val(i) = norm(I*P*Vred(:,ind_red(1:ind))-g);
        end
        % Evaluation of the error estimate (need for efficency index comparison)
        [EV_EE(i),RES, GAP] = EVALUATE_ERROR_ESTIMATE_C(v,pars);
        EV_EE(i) = (RES+EV_EE(i))/GAP;

    end
end

for i=(numCombinations+1):(numCombinations+numel(mu_given))

    v = mu_given(:,i-numCombinations);

    thetanew = theta(v);
    Amu=thetanew(1)*A{1};
    Amu_RED=thetanew(1)*Ared{1};
    for k=2:kappa
        Amu=Amu+thetanew(k)*A{k};
        Amu_RED=Amu_RED+thetanew(k)*Ared{k};
    end
    [Vred,Dred] = eig(Amu_RED); Dred=real(Dred);
    [~,ind_red]=sort(diag(Dred));
    for jj=2:size(Dred,2)
        if (Dred(ind_red(jj),ind_red(jj))-Dred(ind_red(1),ind_red(1)))> RSG_tol
            ind=jj-1;
            break
        end
    end
    if jj==size(Dred,2)
        ind=size(Amu,2);
        fprintf('For parameter %d we have the matrix to have all the same eigenvalue, therefore the ground state is all the space\n Error is set to 0\n',v)
        Eig_Val(i)=eps;
    else
        [V,~] = eigs(Amu,10*ind,'smallestreal',opt);
        Eig_Val(i) = norm((I-V(:,1:ind)*V(:,1:ind)')*P*Vred(:,ind_red(1:ind)));
    end

end

for i=1:p
    if isempty(mu_given)==0
        mu{i}= [mu{i},mu_given(i,:)];
    end
end

OUTPUT.Eig_Val=Eig_Val;
OUTPUT.EV_EE=EV_EE;

end

