function [f,fd] = EVALUATE_ERROR_ESTIMATE(mu,pars)
% f:  target function of EigOPt evaluated at mu
% fd: derivative of target function of EigOPt evaluated at mu
if isfield(pars,'opt_method')
    opt_method = pars.opt_method;
else
    opt_method = 1;
end
if abs(imag(mu))>0
    mu=real(mu);
end

% Option for Relative or Absolute Errors
RSG_tol=pars.RSG_tol;

ne = pars.ne;
alength = length(pars.A);
dim = length(mu);
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
        if ind_2>1
            not_simplicity_lambda1=1;
        else
            not_simplicity_lambda1=0;
        end
        break
    end
end
%% Check Ared is not similar to an identity matrix, if yes terminate
if exist('not_simplicity_lambda1', 'var')==0
    if norm(thetanew)>RSG_tol
        f=10; fd=0;
        return
    else
        f=0; fd=5000;
        return
    end
end
%% -------------------------------------------------------------------
nr=ind_2;
%% LB from SirK16---->Algorithm 1 in [1]
options = pars.options;
[~,kappa] = size(pars.eiglist);
%% Computing beta_j(mu)
beta=zeros(kappa,1);
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
f = (D(inds(1),inds(1)) - slb);
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
RES = max(RES);

slb_2 = min(D(inds(1),inds(1)),eta) - ...
    (2*RES^2)/(abs(D(inds(1),inds(1)) - eta) + sqrt((D(inds(1),inds(1))-eta)^2 + 4*RES^2));
if slb_2>slb
    rhosq = RES^2;
    f = (D(inds(1),inds(1)) - slb_2);
end
%% Evaluating df
h = 10^-8;
slbh = min(D(inds(1),inds(1)),eta+h) - ...
    (2*rhosq)/(abs(D(inds(1),inds(1)) - eta - h) + sqrt((D(inds(1),inds(1)) - eta - h)^2 + 4*rhosq));
% One should consider that eta constrain depends on mu, this would require
% to solve another LP problem to evaluate the discrete derivative.
% For computational sake we avoid this considering eta independnt from (mu)
thetanewp = pars.thetap(mu); fd=zeros(dim,1);
for j = 1:dim

    APmup = thetanewp(1,j)*pars.A{1};
    for k = 2:alength
        APmup = APmup + thetanewp(k,j)*pars.A{k};
    end
    fd(j,1) =  real(V(:,inds(1))'*(APmup*V(:,inds(1)))) - ((slbh - slb)/h)*thetanewp(:,j)'*x;
        if not_simplicity_lambda1==1
            % mu_h=mu;
            % mu_h(j)=mu_h(j)+h;
            % 
            % thetanew_h=pars.theta(mu_h);
            % APmu_h = thetanew_h(1)*pars.A{1};
            % 
            % for k = 2:alength
            %     APmu_h = APmu_h + thetanew_h(k)*pars.A{k};
            % end
            % [~,D_h] = eig(APmu_h); D_h=real(D_h);
            % [~,inds_h] = sort(diag(D_h));
            der=zeros(nr,1);
            for ii=1:nr
               der(ii) = real(V(:,inds(ii))'*(APmup*V(:,inds(ii))));
            end
            der_m = min(der);
            fd(j,1) =  der_m - ((slbh - slb)/h)*thetanewp(:,j)'*x;
        end

end
f_eig = f;
f_eig_d = fd;
%% Evaluating Residual derivative
for j = 1:dim
    mu_h=mu;
    mu_h(j) = mu(j)+h;
    thetanew_h = pars.theta(mu_h);
    Amu_h = thetanew_h(1)*pars.A{1};
    for k = 2:alength
        Amu_h = Amu_h + thetanew_h(k)*pars.A{k};
    end
    [V_h,D_h] = eig(Amu_h); D_h=real(D_h);
    [~,inds_h] = sort(diag(D_h));
    %% Efficent Residual Evaluation
    nn=nr; RES_h = zeros(1,nn); alpha_Res = zeros(Nres,1);
    thetaRES_h = [thetanew_h';-D_h(inds(1),inds(1))];
    for k = 1: nn
        for i=1:rd
            for j=1:(alength+1)
                alpha_Res((i-1)*(alength+1)+j) = thetaRES_h(j)*V_h(i,inds_h(k));
            end
        end
        for i=1:Nres
            RES_k=0;
            for j=1:Nres
                RES_k = RES_k+alpha_Res(j)*pars.coef_ResB(i,j);
            end
            RES_h(k) = RES_h(k)+RES_k^2;
        end
        RES_h(k)=sqrt(abs(RES_h(k)));
    end
    RES_h = max(RES_h);

    RES_d(j,1) = (RES-RES_h)/h;

end
%% Evaluate spectral gap and derivative
PAP_GAP = thetanew(1)*pars.AGAP{1};
for k = 2:alength
    PAP_GAP = PAP_GAP + thetanew(k)*pars.AGAP{k};
end
[V,D] = eig(PAP_GAP); D=real(D);
[~,inds] = sort(diag(D));
for i=2:numel(inds)
    if abs(D(inds(i),inds(i))-D(inds(1),inds(1)))>RSG_tol
        GAP = D(inds(i),inds(i))-D(inds(1),inds(1));
        ind_2r=i;
        if i>2
            mult_GAP=1;
        else
            mult_GAP=0;
        end
        break
    end
end
for j = 1:dim

    APmup = thetanewp(1,j)*pars.AGAP{1};
    for k = 2:alength
        APmup = APmup + thetanewp(k,j)*pars.AGAP{k};
    end

    GAP_d(j,1) =  real(V(:,inds(2))'*(APmup*V(:,inds(2)))) -...
        real(V(:,inds(1))'*(APmup*V(:,inds(1))));

    if mult_GAP==1
        mu_h=mu; GAP_der=zeros(2,1);
        for i=1:2
           
            mu_h(j)=mu(j)+((-1)^i)*h;
            thetanew_h=pars.theta(mu_h);
            PAP_GAP_h = thetanew_h(1)*pars.AGAP{1};
            for k = 2:alength
                PAP_GAP_h = PAP_GAP_h + thetanew_h(k)*pars.AGAP{k};
            end
            [~,D_h] = eig(PAP_GAP_h); D_h=real(D_h);
            [~,inds_h] = sort(diag(D_h));

            GAP_der(i) = ((D_h(inds_h(ind_2r),inds_h(ind_2r))-D(inds(ind_2r),inds(ind_2r)))/h)-((D_h(inds_h(1),inds_h(1))-D(inds(1),inds(1)))/h);
        
        end
        GAP_d(j,1) = min(GAP_der);
    end

end
% Final evaluations
f = (1/GAP)*(RES+f_eig);

if opt_method ~= 1
    f = -f;
end

if abs(imag(f))>0
    display(f)
    fprintf('Complex f should not apper in the Hermitian contest, if imaginary part larger than machine precision then something is going wrong')
    f=real(f);
end

for j = 1:dim

    fd(j,1) =  ((1/GAP)^2)*(GAP*(f_eig_d(j,1)+RES_d(j,1))+GAP_d(j,1)*(RES+f_eig));

    if opt_method ~= 1
        fd(j,1) = -fd(j,1);
    end

end

if opt_method ~= 1
        fd(j,1) = -fd(j,1);
end

end

