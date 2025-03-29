function [f,fd] = approximation_SG(mu,pars)
% f:  target function of EigOPt evaluated at mu
% fd: derivative of target function of EigOPt evaluated at mu
if isfield(pars,'opt_method')
    opt_method = pars.opt_method;
else
    opt_method = 1;
end

if abs(imag(mu))>0
    display(mu);
    mu=real(mu);
end

RSG_tol=pars.RSG_tol;
% Option for Relative or Absolute Errors
RE=pars.Rel_Error;

ne = pars.ne;
alength = length(pars.A);
dim = length(mu);
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
        if ind_2>1
            not_simplicity_lambda1=1;
        else
            not_simplicity_lambda1=0;
        end
        break
    end
end

for ii=(ind_2+2):size(D,2)
    if (D(inds(ii),inds(ii))-D(inds(ind_2+1),inds(ind_2+1)))> RSG_tol
        ind_2=ii-1;
        if ii>ind_2+2
             not_simplicity_lambda2=1;
        else
             not_simplicity_lambda2=0;
        end
        break
    end
end

nr=ind_2;


%% LB from SirK16---->Algorithm 1 in [1]
Lu = D(inds(1:nr),inds(1:nr));
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
if opt_method ~= 1
    f = -f;
end
if abs(imag(f))>0
    display(f)
    fprintf('Complex f should not apper in the Hermitian contest, if imaginary part larger than machine precision then something is going wrong')
    f=real(f);
end
%% Evaluating df
h = 10^-8;
slbh_1 = min(D(inds(1),inds(1)),eta+h) - ...
    (2*rhosq)/(abs(D(inds(1),inds(1)) - eta - h) + sqrt((D(inds(1),inds(1)) - eta - h)^2 + 4*rhosq));
g_h = min(abs(eta +h -D(inds(1),inds(1))),abs(eta+ h-D(inds(nr),inds(nr))));
slbh_2 = min(D(inds(nr),inds(nr)),eta+h) - ...
    (2*rhosq)/(g_h + sqrt(g_h^2 + 4*rhosq));
% One should consider that eta constrain depends on mu, this would require
% to solve another LP problem to evaluate the discrete derivative.
% For computational sake we avoid this considering eta independnt from (mu)
thetanewp = pars.thetap(mu); fd=zeros(dim,1);
for j = 1:dim

    APmup = thetanewp(1,j)*pars.A{1};
    for k = 2:alength
        APmup = APmup + thetanewp(k,j)*pars.A{k};
    end

    % Compute the derivative for the absolute error
    fd(j,1) =  real(V(:,inds(nr))'*(APmup*V(:,inds(nr)))/norm(V(:,inds(nr)))^2) - ((slbh_1 - slb_1)/h)*thetanewp(:,j)'*x...
            +  real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2) - ((slbh_2 - slb_2)/h)*thetanewp(:,j)'*x;

    if  ((not_simplicity_lambda1==1)||(not_simplicity_lambda2==1))
        mu_h=mu;
        mu_h(j)=mu_h(j)+h;

        thetanew_h=pars.theta(mu_h);
        APmu_h = thetanew_h(1)*pars.A{1};

        for k = 2:alength
            APmu_h = APmu_h + thetanew_h(k)*pars.A{k};
        end
        [~,D_h] = eig(APmu_h); D_h=real(D_h);
        [~,inds_h] = sort(diag(D_h));

        if not_simplicity_lambda1==1

            fd(j,1) =  real(V(:,inds(nr))'*(APmup*V(:,inds(nr)))/norm(V(:,inds(nr)))^2) - ((slbh_1 - slb_1)/h)*thetanewp(:,j)'*x...
                +  ((D_h(inds_h(1),inds_h(1))- D(inds(1),inds(1)))/h) - ((slbh_2 - slb_2)/h)*thetanewp(:,j)'*x;


            if not_simplicity_lambda2==1
                fd(j,1) =  ((D_h(inds_h(nr),inds_h(nr))- D(inds(nr),inds(nr)))/h) - ((slbh_1 - slb_1)/h)*thetanewp(:,j)'*x...
                    +  ((D_h(inds_h(1),inds_h(1))- D(inds(1),inds(1)))/h) - ((slbh_2 - slb_2)/h)*thetanewp(:,j)'*x;

            end
        end

        if (not_simplicity_lambda2==1)&&(not_simplicity_lambda1==0)

            fd(j,1) =  ((D_h(inds_h(nr),inds_h(nr))- D(inds(nr),inds(nr)))/h) - ((slbh_1 - slb_1)/h)*thetanewp(:,j)'*x...
                +  real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2) - ((slbh_2 - slb_2)/h)*thetanewp(:,j)'*x;

        end

    end

    % Compute the derivative for the relative error
    if RE
        fd(j,1) =  (1/((D(inds(nr),inds(nr))-D(inds(1),inds(1))))^2)*(fd(j,1)*(D(inds(nr),inds(nr))-D(inds(1),inds(1)))...
            -(D(inds(nr),inds(nr)) - slb_1 - slb_2 + D(inds(1),inds(1)))*(real(V(:,inds(nr))'*(APmup*V(:,inds(nr)))/norm(V(:,inds(nr)))^2)...
            -real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2)));

        if  ((not_simplicity_lambda1==1)&&(not_simplicity_lambda2==1))
            fd(j,1) =  (1/((D(inds(nr),inds(nr))-D(inds(1),inds(1))))^2)*(fd(j,1)*(D(inds(nr),inds(nr))-D(inds(1),inds(1)))...
                -(D(inds(nr),inds(nr)) - slb_1 - slb_2 + D(inds(1),inds(1)))*(((D_h(inds_h(nr),inds_h(nr))- D(inds(nr),inds(nr)))/h)...
                -((D_h(inds_h(1),inds_h(1))- D(inds(1),inds(1)))/h)));

        end
        if  ((not_simplicity_lambda1==0)&&(not_simplicity_lambda2==1))
            fd(j,1) =  (1/((D(inds(nr),inds(nr))-D(inds(1),inds(1))))^2)*(fd(j,1)*(D(inds(nr),inds(nr))-D(inds(1),inds(1)))...
                -(D(inds(nr),inds(nr)) - slb_1 - slb_2 + D(inds(1),inds(1)))*(((D_h(inds_h(nr),inds_h(nr))- D(inds(nr),inds(nr)))/h)...
                -real(V(:,inds(1))'*(APmup*V(:,inds(1)))/norm(V(:,inds(1)))^2)));


        end
        if  ((not_simplicity_lambda1==1)&&(not_simplicity_lambda2==0))
            fd(j,1) =  (1/((D(inds(nr),inds(nr))-D(inds(1),inds(1))))^2)*(fd(j,1)*(D(inds(nr),inds(nr))-D(inds(1),inds(1)))...
                -(D(inds(nr),inds(nr)) - slb_1 - slb_2 + D(inds(1),inds(1)))*(real(V(:,inds(nr))'*(APmup*V(:,inds(nr)))/norm(V(:,inds(nr)))^2)...
                -((D_h(inds_h(1),inds_h(1))- D(inds(1),inds(1)))/h)));

        end

    end

    if abs(rhosq)<5e-14
        fd(j,1)=0;
    end


    if opt_method ~= 1
        fd(j,1) = -fd(j,1);
    end

end
end

