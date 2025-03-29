function [f] = approximation_SG_d(mu,pars)
% f: spectral gap surrogate error evaluated at mu

RSG_tol=pars.RSG_tol;
% Option for Relative or Absolute Errors
RE=pars.Rel_Error;

ne = pars.ne;
alength = length(pars.A);
thetanew = pars.theta(mu);
APmu = thetanew(1)*pars.A{1};

for k = 2:alength
    APmu = APmu + thetanew(k)*pars.A{k};
end
[V,D] = eig(APmu); D=real(D);
[~,inds] = sort(diag(D));

% Need to ensure the right separation
for jj=2:size(D,2)
    if (D(inds(jj),inds(jj))-D(inds(1),inds(1)))> RSG_tol
        ind_2=jj-1;
        break
    end
end

for ii=(ind_2+2):size(D,2)
    if (D(inds(ii),inds(ii))-D(inds(ind_2+1),inds(ind_2+1)))> RSG_tol
        ind_2=ii-1;
        break
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
    b(j,1) = -pars.eiglist{j}(1)-beta(j,1);
end
%% Setting the LP constraints
M = -pars.thetalist; 
for j = 1:alength
    vecadd = zeros(1,alength);
    vecadd(j) = 1;
    M = [M; vecadd; -vecadd];
    b(kappa+2*j-1,1)=pars.lambounds(j,2);
    b(kappa+2*j,1)= -pars.lambounds(j,1);
end
if max(abs(imag(b)))>0
   b=real(b);
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
%LB for lambda_1
slb_1 = min(D(inds(1),inds(1)),eta) - ...
    (2*rhosq)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*rhosq));
%LB for lambda_2
g = min(abs(eta-D(inds(1),inds(1))),abs(eta-D(inds(nr),inds(nr))));
slb_2 = min(D(inds(nr),inds(nr)),eta) - ...
    (2*rhosq)/(g + sqrt(g^2 + 4*rhosq));

if RE
    f = (D(inds(nr),inds(nr)) - slb_1 - slb_2 + D(inds(1),inds(1)))/(D(inds(nr),inds(nr))-D(inds(1),inds(1)));
else
    f = (D(inds(nr),inds(nr)) - slb_1 - slb_2 + D(inds(1),inds(1)));
end

end

