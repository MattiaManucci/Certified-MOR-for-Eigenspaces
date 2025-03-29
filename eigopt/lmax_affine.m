function [f,g] = lmax_affine(y,pars)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% Auxiliary routine called by Hstructured_radii_e, Hstructured_radii_e_invQs
%
% TASK:
% Computes lambda_min(pars.H0 + t pars.H1) and its derivative at a given t,
% returns the computed values inside f and g.


[eigvecs,eigvals] = eig(pars.C - y*pars.A);

[f,indx] = max(real(diag(eigvals)));


f = ((pars.b)/2)*f + pars.b*y;
g = -(pars.b/2)*real(eigvecs(:,indx)'*pars.A*eigvecs(:,indx)) + pars.b;


return