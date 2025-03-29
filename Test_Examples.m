%% ------------------ Test Problems related to [1]----------------------------------------------------------
% [1] M. Manucci, B. Stamm, and Z. Zeng, 
%     Certified Model Order Reduction for parametric Hermitian
%     eigenproblems, 2025
% [2] Mustafa Kilic, Emre Mengi and E. Alper Yildirim, SIMAX 2014
%% -----------------------------------------------------------------
clc 
clearvars
close all
%% LEGEND OF NUMERICAL EXAMPLES
fprintf('Digit the corresponding integer number to run the test problem:\n\n');
fprintf(['1 for full Hermitian matrix with n=2000 and 1 parameter, ' ...
        'please make sure to have downloaded the Data2.mat file (see readme)\n']);
fprintf('2 for the QSS Example xxz chain model\n');
fprintf('3 for the QSS Example blbq chain model\n');
flag=input('Please insert the corresponding number for the problem you want to run\n');
%% Folders inclusions for functions and plot setting 
addpath(genpath('./eigopt/'))
addpath(genpath('./Data_Matrices/'))
addpath(genpath('./export_fig/'))
addpath(genpath('./Algorithm 2 Func./'))
addpath(genpath('./Algorithm 1 Func./'))
addpath(genpath('./Plot_Functions/'))
addpath(genpath('./github_repo/src'))
FS = 15;       % Fontsize
FN = 'times';  % Fontname
LW = 1.6;      % Linewidth
MS = 7.8;      % Markersize
%% Example 1----->Large Full Matrix Example
if flag==1
    load('Data2') % The data matrices needs to be dowloanded, please have a look to the readme.
    % Settle the proiblem
    theta = @(x)[x.^2, x]; % theta(x) are the analytic functions for the affine decomposition
    theta_d = @(x)[2*x;1]; % thetap(x) is the Jacobian of theta(x)
    bounds.lb = 1e-1; bounds.ub = 10; % Parametric domain
    A{1} = A1; A{2} = A2; % Matrices of the affine decomposition
    % Options for the greedy algorithm
    opts.num_init_inter=2; %Number of intial parameters to interpolate
    opts.Rel_Error = 1; %Run Alg. 1 to construct $V_{\gamma}$ such that (4.1) from [1] smaller than $\varepsilon_{\gamma}$
    opts.tol = 1e-8; %$\varepsilon_{\gamma}$ used for Alg. 1 from [1]
    opts.RSG_tol = 1e-7; %Relative threshold to include eigenvalues coalescence. Default is 1e-6.
    %% Spectral GAP greedy-algorithm (Algorithm 1 from [1] with EigOpt from [2])
    [f,Ared_GAP,pars_GAP] = approx_smallesteig_all(A,theta,theta_d,bounds,opts); 
    %% EIGENSPACE greedy-algorithm (Algorithm 2 from [1] with EigOpt from [2])
    opts.Space_Gap=pars_GAP.P; %Subspace constructed for the spectral gap approximation
    opts.tol = 1e-8; %$\varepsilon_{W}$ used for Alg. 2 of [1]
    [ERR_EST,Ared,pars] = App_SEVES(A,theta,theta_d,bounds,opts); 
    %% Post processing and plots
    ff =ERR_EST.ff;
    RES=ERR_EST.RES(opts.num_init_inter+1:end); 
    GAP=ERR_EST.GAP(opts.num_init_inter+1:end); 
    EEV=ERR_EST.EV_EE(opts.num_init_inter+1:end); 
    EIG=ERR_EST.EIG_E(opts.num_init_inter+1:end); 
    % Plot options
    plot_opts.Nmu=200; %Number of parameters for the discrete grid where the error is evaluated
    plot_opts.sparse=issparse(A{1}); %Check if matrices are sparse or full
    plot_opts.log_space=1; %Set 1 to have a logarithmic distribution in the discrete grid
    plot_opts.Hermitian=1; %set 0 if the matrix is not Hermitian
    plot_opts.mu=[]; %To include specific points in the plot
    plot_opts.Rel_Error = 1; %Set 1 to run the problem for relative accuracy
    plot_opts.flag_PC = 0; %Set 1 to use parallel computing, 0 otherwise 
    % Spectral gap computations
    [GAP_DATA] = plot_Spectral_GAP(A,Ared_GAP,theta,bounds,plot_opts);
    % Eigenspace error computation
    [mu_d,~, OUTPUT_ES_PLOT] = plot_lambdamin_MD(A,Ared,pars.P,theta,bounds,plot_opts,pars);
    eigvec_err=OUTPUT_ES_PLOT.Eig_Val;
    eigvec_err_est=OUTPUT_ES_PLOT.EV_EE;
    %% GAP Plots
    %Interpolation points for GAP
    GAP_int=zeros(numel(pars_GAP.eiglist),1);
    for i=1:numel(pars_GAP.eiglist)
       GAP_int(i)= pars_GAP.eiglist{i}(2)-pars_GAP.eiglist{i}(1);
    end
    %Figure 1a from [1] 
    [GAP_DATA_mu,ind_GAP] = sort(GAP_DATA.mu{1});
    figure
    loglog(GAP_DATA_mu,GAP_DATA.GAP(ind_GAP),'-b','LineWidth',LW)
    hold on
    loglog(GAP_DATA_mu,GAP_DATA.GAP_red(ind_GAP),'--r','LineWidth',LW)
    loglog(pars_GAP.mu(1:end-1),GAP_int,'*k','Markersize',MS)
    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
    %Figure 1b from [1] 
    figure
    loglog(GAP_DATA_mu,GAP_DATA.GAP_err(ind_GAP),'-b','LineWidth',LW)
    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
    %Figure 1c from [1] 
    figure
    semilogy(f,'-ob','LineWidth',LW)
    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
    %% Eigenspace Plots
    %Figure 2a from [1]
    [mu_eig,ind_EIGS]=sort(mu_d{1});
    EEV(EEV==0)=eps;
    figure
    loglog(mu_eig,eigvec_err(ind_EIGS),'-b','LineWidth',LW)
    hold on
    loglog(mu_eig,eigvec_err_est(ind_EIGS),'--r','LineWidth',LW)

    lgd=legend('$\|\mathcal{P}^{\mathcal{W}_1^{\bot}(\mu)}\mathbf{W}_1(\mu)\|$','$\Delta(\mu)$','Location','best');
    set(lgd,'Interpreter','Latex');

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
    %Figure 2b from [1]
    figure
    semilogy(ff,'-ob','LineWidth',LW)
    hold on
    semilogy(EIG,'-og','LineWidth',LW)
    semilogy(RES./GAP,'-*r','LineWidth',LW)
    semilogy(EEV./GAP,'-+k','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')

    lgd=legend('$\Delta(\mu^*)$','$P^{\mathcal{W}_{1}(\mu^*)}-P^{\mathcal{W}^{\mathcal{V}}_{1}(\mu^*)}$','$\|R(\mu^*)\|_2/\gamma(\mu)$','$(\lambda^{{UB}}(\mu^*)-\lambda^{{LB}}(\mu^*))/\gamma(\mu)$','Location','best');
    set(lgd,'Interpreter','Latex');

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

end
if flag==2
    % QSS Example xxz chain model
    load('DataQSS.mat')
    theta = @(x)[1, x(1), -x(2)];
    theta_d = @(x)[0,0; 1,0; 0, -1]; % thetap(x) is the Jacobian of theta(x)
    bounds.lb = [-1; 0];
    bounds.ub = [2.5; 3.5];
    A{1} = A1;
    A{2} = A2; 
    A{3} = A3;
    %% Options to run the problem
    opts.Rel_Error = 1; %Set 0 if you want to run Algorithm 2 of [1] for H(\mu) or 1 for H_r(\mu). Defaoult is 1.
    opts.tol = 1e-8; %\varepsilon in Algorithm 2 of [1]. Default is 1e-4.
    opts.RSG_tol = 1e-8; %Relative threshold to include eigenvalues coalescence. Default is 1e-6.
    opts.gamma=-4e4; %Lower bound of the curvature for the target function, see [2], default -4e5.
    opts.Nt = 35;
    opts.num_init_inter=2;
    opts.flag_PC=1; %Set to 1 to use parallel comupting toolbox, set to 0 to not use it
    opts.opt_method = 1; % To use eigopt for subproblems
    %% Spectral GAP greedy-algorithm (Algorithm 1 from [1])
    [f,Ared_GAP,pars_GAP] = approx_sGAP_Dset(A,theta,theta_d,bounds,opts);
    %% Eigenspace approximation (Algorithm 2 from [1])
    opts.Space_Gap=pars_GAP.P;
    opts.tol = 1e-8; 
    opts.num_init_inter=2; %Number of intial parameters to interpolate
    opts.flag_PC=1; %Set to 1 to use parallel comupting toolbox, set to 0 to not use it
    [ERR_EST,Ared,pars] = subspace_SCMM(A,theta,theta_d,bounds,opts);
    
    ff =ERR_EST.ff;
    RES=ERR_EST.RES(opts.num_init_inter+1:end);
    GAP=ERR_EST.GAP_IT(opts.num_init_inter+1:end);
    EEV=ERR_EST.EEV(opts.num_init_inter+1:end);
    EIG=ERR_EST.EIG_E(opts.num_init_inter+1:end); 
    %% Plots for QSS Example
    plot_opts.Nmu=[opts.Nt, opts.Nt]; %Number of parameters for the discrete grid
    plot_opts.sparse=issparse(A{1}); %Check if problem is sparse or full
    plot_opts.log_space=[0,0]; %Set 1 to have a logarithmic distribution in the discrete grid
    plot_opts.Hermitian=1; %set 0 if the matrix is not Hermitian
    plot_opts.mu=[]; %To include specific points in the plot
    plot_opts.RSG_tol=opts.RSG_tol;
    plot_opts.flag_PC = 1; %Set 1 to use parallel computing, 0 otherwise 
    % Spectral gap computations
    [GAP_DATA] = plot_Spectral_GAP(A,Ared_GAP,theta,bounds,plot_opts);
    % Eigenspace error computation
    [mu_d,~, OUTPUT_ES_PLOT] = plot_lambdamin_MD(A,Ared,pars.P,theta,bounds,plot_opts,pars);
    eigvec_err=OUTPUT_ES_PLOT.Eig_Val;
    eigvec_err_est=OUTPUT_ES_PLOT.EV_EE;
  
    GAP_int=zeros(numel(pars_GAP.eiglist),1);
    for i=1:numel(pars_GAP.eiglist)
       GAP_int(i)= pars_GAP.eiglist{i}(pars_GAP.ne(i))-pars_GAP.eiglist{i}(1);
    end
    
    %Figure 3a from [1]
    GAP_full = reshape(GAP_DATA.GAP, plot_opts.Nmu(1)+size(plot_opts.mu,2), plot_opts.Nmu(2)+size(plot_opts.mu,2));
    GAP_red  = reshape(GAP_DATA.GAP_red,  plot_opts.Nmu(1)+size(plot_opts.mu,2), plot_opts.Nmu(2)+size(plot_opts.mu,2));
    GAP_full = GAP_full'; GAP_red = GAP_red';
    figure
    h=gca;
    surf(GAP_DATA.mu{1}, GAP_DATA.mu{2}, GAP_red)
    set(h,'zscale','log')
     
    %Figure 3b from [1]
    GAP_DATA_err = reshape(GAP_DATA.GAP_err, plot_opts.Nmu(1)+size(plot_opts.mu,2), plot_opts.Nmu(2)+size(plot_opts.mu,2));
    GAP_DATA_err(GAP_DATA_err==0)=eps;
    GAP_DATA_err = GAP_DATA_err';
    figure
    Z_log = log10(GAP_DATA_err);
    surf(GAP_DATA.mu{1}, GAP_DATA.mu{2}, Z_log)
    shading interp; % Smooth colors across surface
    view(2); % Top-down view
    c = colorbar;
    % Customize color bar to show original scale
    c.Ticks = log10([1e-15, 1e-14, 1e-13, 1e-12, 1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3]);  % Adjust based on data range
    c.TickLabels = {'1e-15', '1e-14', '1e-13', '1e-12', '1e-11','1e-10','1e-9','1e-8','1e-7','1e-6','1e-5','1e-4','1e-3'};
    hold on
    plot(pars_GAP.mu(1,:),pars_GAP.mu(2,:),'rx')
    %% Eigenspaces plot
    %Figure 4a from [1]
    eigvec_err = reshape(eigvec_err,opts.Nt,opts.Nt);
    eigvec_err = eigvec_err';
    figure
    Z_log = log10(eigvec_err);
    h=gca;
    surf(mu_d{1},mu_d{2}, Z_log)
    shading interp; % Smooth colors across surface
    view(2); % Top-down view
    c = colorbar;
    % Customize color bar to show original scale
    c.Ticks = log10([1e-15, 1e-14, 1e-13, 1e-12, 1e-11,1e-10,1e-9,1e-8]);  % Adjust based on data range
    c.TickLabels = {'1e-15', '1e-14', '1e-13', '1e-12', '1e-11','1e-10','1e-9','1e-8'};
    hold on
    plot(pars.mu(1,:),pars.mu(2,:),'rx')

    %Figure 4b from [1]
    EEV(EEV==0)=eps; %To have it in log plots
    figure
    semilogy(ff,'-ob','LineWidth',LW)
    hold on
    semilogy(EIG,'-.xg','LineWidth',LW)
    semilogy(RES./GAP,'-*r','LineWidth',LW)
    semilogy(EEV./GAP,'-+k','LineWidth',LW)
    xlabel('$j$','Interpreter','Latex')

    lgd=legend('$\Delta(\mu^*)$','$\|P^{\mathcal{W}_{1}(\mu^*)}-P^{\mathcal{W}^{\mathcal{V}}_{1}(\mu^*)}\|$','$\|R(\mu^*)\|_2/\gamma(\mu)$','$(\lambda^{{UB}}(\mu^*)-\lambda^{{LB}}(\mu^*))/\gamma(\mu)$','Location','best');
    set(lgd,'Interpreter','Latex');
    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
end

if flag==3
    % QSS Example blbq chain model
    % There is a problem with the approximation of the eigenvalue GAP
    load('datablbq_OBC_chain.mat')
    theta = @(x)[cos(x(1)), sin(x(1)), x(2)];
    theta_d = @(x)[-sin(x(1)),0; cos(x(1)),0; 0, 1]; % thetap(x) is the Jacobian of theta(x)
    bounds.lb = [-pi; -2];
    bounds.ub = [pi; 3];
    A{1} = matrices{1};
    A{2} = matrices{2}; % Matrix-Valued FuncGon : exp(x(1))*A{1} + x(2)*A(2)
    A{3} = matrices{3};
    %% Options to run the problem
    opts.Rel_Error = 1; %Set 0 if you want to run Algorithm 2 of [1] for H(\mu) or 1 for H_r(\mu). Defaoult is 1.
    opts.tol = 1e-2; %\varepsilon in Algorithm 2 of [1]. Default is 1e-4.
    opts.RSG_tol = 1e-8; %Relative threshold to include eigenvalues coalescence. Default is 1e-6.
    opts.Nt = 50;
    opts.num_init_inter=2;
    opts.flag_PC=1; %Set to 1 to use parallel comupting toolbox, set to 0 to not use it
    %% Spectral GAP greedy-algorithm (Algorithm 1 [1] with finite grid)
    [f,Ared_GAP,pars_GAP] = approx_sGAP_Dset(A,theta,theta_d,bounds,opts);
    %% Eigenspace approximation (Algorithm 2 [1] with finite grid)
    opts.Space_Gap=pars_GAP.P;
    opts.tol = 1e-4; 
    opts.flag_PC=0; %Set to 1 to use parallel comupting toolbox, set to 0 to not use it
    [ERR_EST,Ared,pars] = subspace_SCMM(A,theta,theta_d,bounds,opts);
    ff =ERR_EST.ff;
    RES=ERR_EST.RES(opts.num_init_inter+1:end); 
    GAP=ERR_EST.GAP_IT(opts.num_init_inter+1:end); 
    EEV=ERR_EST.EEV(opts.num_init_inter+1:end); 
    %% Plots for QSS Example
    plot_opts.Nmu=[opts.Nt, opts.Nt]; %Number of parameters for the discrete grid
    plot_opts.sparse=issparse(A{1}); %Check if problem is sparse or full
    plot_opts.log_space=[0,0]; %Set 1 to have a logarithmic distribution in the discrete grid
    plot_opts.Hermitian=1; %set 0 if the matrix is not Hermitian
    plot_opts.mu=[]; %To include specific points in the plot
    plot_opts.RSG_tol=opts.RSG_tol;
    plot_opts.flag_PC = 1; %Set 1 to use parallel computing, 0 otherwise 
    % Spectral gap computations
    [GAP_DATA] = plot_Spectral_GAP(A,Ared_GAP,theta,bounds,plot_opts);
    % Eigenspace error computation
    [mu_d,~, OUTPUT_ES_PLOT] = plot_lambdamin_MD(A,Ared,pars.P,theta,bounds,pars);
    eigvec_err=OUTPUT_ES_PLOT.Eig_Val;
    eigvec_err_est=OUTPUT_ES_PLOT.EV_EE;

    GAP_int=zeros(numel(pars_GAP.eiglist),1);
    for i=1:numel(pars_GAP.eiglist)
        GAP_int(i)= pars_GAP.eiglist{i}(pars_GAP.ne(i))-pars_GAP.eiglist{i}(1);
    end
    
    % Figure 5a from [1]
    GAP_full = reshape(GAP_DATA.GAP, plot_opts.Nmu(1)+size(plot_opts.mu,2), plot_opts.Nmu(2)+size(plot_opts.mu,2));
    GAP_red  = reshape(GAP_DATA.GAP_red,  plot_opts.Nmu(1)+size(plot_opts.mu,2), plot_opts.Nmu(2)+size(plot_opts.mu,2));
    GAP_full=GAP_full'; GAP_red=GAP_red';
    figure
    h=gca;
    surf(GAP_DATA.mu{1}, GAP_DATA.mu{2}((1:2:end)), GAP_red((1:2:end),:))
    set(h,'zscale','log')
    
    % Figure 5b from [1]
    GAP_DATA_err = reshape(GAP_DATA.GAP_err, plot_opts.Nmu(1)+size(plot_opts.mu,2), plot_opts.Nmu(2)+size(plot_opts.mu,2));
    GAP_DATA_err(GAP_DATA_err<=0)=eps;
    GAP_DATA_err(isnan(GAP_DATA_err))=eps; %This may be needed if eigs does not converge
    GAP_DATA_err=GAP_DATA_err';
    figure
    Z_log = log10(GAP_DATA_err);
    surf(GAP_DATA.mu{1}, GAP_DATA.mu{2}, Z_log)
    shading interp; % Smooth colors across surface
    view(2); % Top-down view
    c = colorbar;
    % Customize color bar to show original scale
    c.Ticks = log10([1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3]);  % Adjust based on data range
    c.TickLabels = {'1e-16','1e-15', '1e-14', '1e-13', '1e-12', '1e-11','1e-10','1e-9','1e-8','1e-7','1e-6','1e-5','1e-4','1e-3'};
    hold on
    plot(pars_GAP.mu(1,:),pars_GAP.mu(2,:),'rx')

    % Figure 7a from [1]
    figure
    plot(pars_GAP.ne,'-b','LineWidth',LW)
    hold on
    plot(pars.ne,'-r','LineWidth',LW)
    xlabel('$i$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    % Figure 7b from [1]
    cond_cho = pars_GAP.eta-pars_GAP.lambda_1_red-pars_GAP.eta_epsilon;
    cond_cho_eig = pars.eta-pars.lambda_1_red-pars.eta_epsilon;
    figure
    semilogy(cond_cho,'-b','LineWidth',LW)
    hold on
    semilogy(cond_cho_eig,'-r','LineWidth',LW)
    xlabel('$mu_i$','Interpreter','Latex')

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');

    %% Eigenspaces approximation plots
    % Plot 6a from [1]
    eigvec_err = reshape(eigvec_err,opts.Nt,opts.Nt);
    eigvec_err = eigvec_err';
    figure
    Z_log = log10(eigvec_err);
    h=gca;
    surf(mu_d{1},mu_d{2}, Z_log)
    shading interp; % Smooth colors across surface
    view(2); % Top-down view
    c = colorbar;
    % Customize color bar to show original scale
    c.Ticks = log10([1e-15, 1e-14, 1e-13, 1e-12, 1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4]);  % Adjust based on data range
    c.TickLabels = {'1e-15', '1e-14', '1e-13', '1e-12', '1e-11','1e-10','1e-9','1e-8','1e-7','1e-6','1e-5','1e-4'};
    hold on
    plot(pars.mu(1,:),pars.mu(2,:),'rx')

    % Plot 6b from [1]
    figure
    semilogy(ff(1:2:end),'-ob','LineWidth',LW)
    hold on
    semilogy(RES(1:2:end)./GAP(1:2:end),'-*r','LineWidth',LW)
    semilogy(EEV(1:2:end)./GAP(1:2:end),'-+k','LineWidth',LW)
    xlabel('$r$','Interpreter','Latex')

    lgd=legend('$\Delta(\mu^*)$','$\|R(\mu^*)\|_2$','$\lambda^{{UB}}(\mu^*)-\lambda^{{LB}}(\mu^*)$','$\|v_1(\mu^*)-\tilde{v}_1(\mu^*)\|_2$','Location','best');
    set(lgd,'Interpreter','Latex');

    set(gca,'Fontname',FN,'Fontsize',FS);
    set(gcf, 'Color', 'w');
end
