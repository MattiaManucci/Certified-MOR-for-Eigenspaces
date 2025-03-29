function [output] = lambda_eta_eps_eig(mu,pars)

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

% Need to ensure the right separation
ind_2=1;
for jj=2:size(D,2)
    if (D(inds(jj),inds(jj))-D(inds(1),inds(1)))> RSG_tol
        ind_2=jj-1;
        break
    end
end
nr=ind_2;
%% LB from SirK16---->Algorithm 1 in [1]
options = pars.options;
[~,kappa] = size(pars.eiglist);
%% Computing beta_j(mu)
beta=zeros(kappa,1); b=zeros(kappa,1);
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
    b = [b; pars.lambounds(j,2); -pars.lambounds(j,1)];
end
%% Linear programming for lower bound
x = linprog(thetanew',M,b,[],[],[],[],options);
eta = thetanew*x;
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
%% Evaluating output
output.eigr = D(inds(1),inds(1));
output.eta = eta;
output.eta_epsilon = (2*RES^2)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*RES^2));

end

