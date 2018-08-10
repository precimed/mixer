function params = BGMG_cpp_fit_bivariate_fast(params1, params2)
    % params1   - univariate parameters for the first trait
    % params2   - univariate parameters for the second trait    
    % options   - options; see the below for a full list of configurable options

    fminsearch_options = struct('Display', 'on');
    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), fminsearch_options));

    bgmglib = BGMG_cpp();
    bgmglib.set_option('fast_cost', 1);
    zmat = [bgmglib.zvec1, bgmglib.zvec2];
    corrmat = corr(zmat(all(~isnan(zmat), 2), :));

    BGMG_cpp.log('Trait 1,2: fit infinitesimal model to find initial rho_zero and rho_beta (fast cost function)\n');
    params_inft0 = struct('rho_zero', corrmat(1,2), 'rho_beta', corrmat(1,2));
    options_inft = struct('sig2_zero', [params1.sig2_zero; params2.sig2_zero], ...
                          'sig2_beta', [params1.pi_vec * params1.sig2_beta; params2.pi_vec * params2.sig2_beta], ...
                          'pi_vec', 1);
    params_inft  = fit(params_inft0, @(x)BGMG_util.BGMG_mapparams1_rho(x, options_inft));

    % Step2. Final unconstrained optimization (jointly on all parameters), using fast cost function
    func_map_params = @(x)BGMG_util.BGMG_mapparams3(x, struct('sig2_zero', [params1.sig2_zero; params2.sig2_zero], 'sig2_beta', [params1.sig2_beta; params2.sig2_beta], 'rho_zero', params_inft.rho_zero));

    BGMG_cpp.log('Trait 1,2: final unconstrained optimization (fast cost function)\n');
    params_final0 = params_inft;
    params_final0.sig2_beta = [params1.sig2_beta 0 params1.sig2_beta; 0 params2.sig2_beta params2.sig2_beta];
    init_pi12 = min(params1.pi_vec, params2.pi_vec)/exp(1);
    params_final0.pi_vec = [params1.pi_vec-init_pi12 params2.pi_vec-init_pi12 init_pi12];
    params_final0.rho_beta = [0 0 params_inft.rho_beta];

    params = fit(params_final0, func_map_params);
end



