function [figures, plots_data] = BGMG_cpp_loglike_plot(params)
    figures = [];
    plots_data = [];

	logspace_43 = logspace(log10(3/4), log10(4/3), 11);
    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));

    pi1u = sum(params.pi_vec([1 3]));
    pi2u = sum(params.pi_vec([2 3]));
    pi12 = params.pi_vec(3);
    dat_trait1.sigsq = params.sig2_beta(1, 3);
    dat_trait2.sigsq = params.sig2_beta(2, 3);

    plot_name = 'pi12_frac';
    plots_data.(plot_name).x = 0:0.05:1; plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [(1-x), (1-x), x] * pi1u;
        sig2_beta = [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero; rho_beta = params.rho_beta(end); rho_zero = params.rho_zero;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;
    if 0
    plot_name = 'pivec_scale';
    plots_data.(plot_name).x = logspace_43; plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [x, x, x] .* [pi1u-pi12, pi2u-pi12, pi12];
        sig2_beta = [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero; rho_beta = params.rho_beta(end); rho_zero = params.rho_zero;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;

    plot_name = 'sig2beta_scale';
    plots_data.(plot_name).x = logspace_43; plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
        sig2_beta = [x, x] .* [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero; rho_beta = params.rho_beta(end); rho_zero = params.rho_zero;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;

    plot_name = 'pivec_scale_const_h2';
    plots_data.(plot_name).x = logspace_43; plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [x, x, x] .* [pi1u-pi12, pi2u-pi12, pi12];
        sig2_beta = [1./x, 1./x] .* [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero; rho_beta = params.rho_beta(end); rho_zero = params.rho_zero;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;

    plot_name = 'sig2zero';
    plots_data.(plot_name).x = logspace(log10(0.5), log10(2), 11); plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
        sig2_beta = [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero * x; rho_beta = params.rho_beta(end); rho_zero = params.rho_zero;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;

    plot_name = 'rhozero';
    plots_data.(plot_name).x = -0.1:0.01:0.1; plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
        sig2_beta = [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero; rho_beta = params.rho_beta(end); rho_zero = x;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;

    plot_name = 'rhobeta';
    plots_data.(plot_name).x = -0.2:0.02:0.2; plots_data.(plot_name).y = [];
    for x = plots_data.(plot_name).x, 
        pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
        sig2_beta = [dat_trait1.sigsq, dat_trait2.sigsq];
        sig2_zero = params.sig2_zero; rho_beta = x; rho_zero = params.rho_zero;
        plots_data.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
    end;
    end
    fnames = fieldnames(plots_data);
    for fni=1:length(fnames)
        subplot(3,3,fni);hold on; 
        plot_name=fnames{fni};
        y = plots_data.(plot_name).y; y(y>1e99)=nan;
        plot(plots_data.(plot_name).x, y-min(y), '-*');
        plot_name(plot_name=='_')=' ';title(plot_name);
        %if fig==3
        %    legend({'complete', 'partial25', 'random'})
        %end
    end
end