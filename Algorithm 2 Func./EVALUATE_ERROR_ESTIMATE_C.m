function [f_eig,RES,GAP] = EVALUATE_ERROR_ESTIMATE_C(mu,pars)
%% Output
% -f_eig:  eigenvalue surrogate error evaluated at mu;
% -RES:    residul evaluated at mu;
% -GAP:    approximated spectral gap evaluated at mu.
%% ------------------------------------------------------------------------
RSG_tol=pars.RSG_tol;
ne = pars.ne;
alength = length(pars.A);
thetanew = pars.theta(mu);
APmu = thetanew(1)*pars.A{1};

for k = 2:alength
    APmu = APmu + thetanew(k)*pars.A{k};
end
[V,D] = eig(APmu);
D=real(D); 
[~,inds] = sort(diag(D));

% Check if the smallest eigenvalue is not simple up to the given threshold
for jj=2:size(D,2)
    if (D(inds(jj),inds(jj))-D(inds(1),inds(1)))> RSG_tol
        ind_2=jj-1;
        break
    end
end
%% Check Ared is not similar to an identity matrix, if yes terminate
if exist('ind_2', 'var')==0
    if norm(thetanew)>1-15
        f_eig=10; RES=0; GAP=eps;
        return
    else
        f_eig=0; RES=0; GAP=eps;
        return
    end
end
nr=ind_2;
%% LB from SirK16---->Algorithm 1 in [1]
options = pars.options;
[~,kappa] = size(pars.eiglist);
%% Computing beta_j(mu)
beta=zeros(kappa,1); b=zeros(kappa+2*alength,1);
for j = 1:kappa
    Li = diag(pars.eiglist{j}(1:ne(j)));
    li = pars.eiglist{j}(1);
    lip = pars.eiglist{j}(ne(j)+1);
    beta(j,1) = min(eig((Li-li*eye(ne(j))) - pars.premult{j}*pars.premult{j}'*(Li-lip*eye(ne(j)))));
    b(j,1)=-pars.eiglist{j}(1)-beta(j,1);
end
%% Setting the LP constraints
M = -pars.thetalist;
for j = 1:alength
    vecadd = zeros(1,alength);
    vecadd(j) = 1;
    M = [M; vecadd; -vecadd];
    b(kappa+2*j-1,1)=pars.lambounds(j,2);
    b(kappa+2*j,1)= -pars.lambounds(j,1);
   % b = [b; pars.lambounds(j,2); -pars.lambounds(j,1)];
end
%% Linear programming for lower bound
x = linprog(thetanew',M,b,[],[],[],[],options);
eta = thetanew*x;
%% Computation of rho(mu)
AAmu=zeros(size(pars.A{1},1));
for j=1:alength
    for jj=1:alength
        AAmu=AAmu+thetanew(j)*thetanew(jj)*pars.Afull{(alength)*(j-1)+jj};
    end
end
temp = (V(:,inds(1:nr)))'*AAmu*V(:,inds(1:nr));
rhosq = abs(max(eig(temp - V(:,inds(1:nr))'*APmu*APmu*V(:,inds(1:nr)))));
%% Evaluating f
slb = min(D(inds(1),inds(1)),eta) - ...
    (2*rhosq)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*rhosq));
f_eig = (D(inds(1),inds(1)) - slb);
%% Stable evaluation of the residual
nn=nr; RES = zeros(1,nn); rd = size(D,1); Nres = pars.Nres; alpha_Res = zeros(Nres,1);
thetaRES = [thetanew';-D(inds(1),inds(1))];
for k = 1: nn
    for i=1:rd
        for j=1:(alength+1)
            alpha_Res((i-1)*(alength+1)+j) = thetaRES(j)*V(i,inds(k));
        end
    end
    for i=1:Nres
        RES_k=0;
        for j=1:Nres
            RES_k = RES_k+alpha_Res(j)*pars.coef_ResB(i,j);
        end
        RES(k) = RES(k)+RES_k^2;
    end
    RES(k)=sqrt(abs(RES(k)));
end
RES=max(RES);
slb_2 = min(D(inds(1),inds(1)),eta) - ...
    (2*RES^2)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*RES^2));
if slb_2>slb
    f_eig = (D(inds(1),inds(1)) - slb_2);
end
% Evaluate spectral gap 
PAP_GAP = thetanew(1)*pars.AGAP{1};
for k = 2:alength
    PAP_GAP = PAP_GAP + thetanew(k)*pars.AGAP{k};
end
[~,D] = eig(PAP_GAP); D=real(D);
[~,inds] = sort(diag(D));
for i=2:numel(inds)
    if abs(D(inds(i),inds(i))-D(inds(1),inds(1)))>RSG_tol
        GAP = D(inds(i),inds(i))-D(inds(1),inds(1));
        break
    end
end

end


