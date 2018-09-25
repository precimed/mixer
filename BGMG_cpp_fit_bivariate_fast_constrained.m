function params = BGMG_cpp_fit_bivariate_fast_constrained(params)
    % params1   - univariate parameters for the first trait
    % params2   - univariate parameters for the second trait    
    % options   - options; see the below for a full list of configurable options

    fminsearch_options = struct('Display', 'on');
    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), fminsearch_options));

    bgmglib = BGMG_cpp();
    bgmglib.set_option('fast_cost', 1);

    BGMG_cpp.log('Trait 1,2: Fit sig2_zero and rho_zero constrained to all other params (fast cost function)\n');
    params0 = struct('sig2_zero', params.sig2_zero, 'rho_zero', params.rho_zero);
    constrain = struct('sig2_beta', params.sig2_beta(:, 3), 'pi_vec', params.pi_vec, 'rho_beta', params.rho_beta);
    params = fit(params0, @(x)BGMG_util.BGMG_mapparams3(x, constrain));
end



