function params = BGMG_cpp_fit_univariate_fast_constrained(trait_index, params)
    bgmglib = BGMG_cpp();
    bgmglib.set_option('fast_cost', 1);

    fminsearch_options = struct('Display', 'on');
    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), fminsearch_options));

    BGMG_cpp.log('Fit sig2_zero constrained to all other params (fast cost function)\n');
    params = fit(struct('sig2_zero', params.sig2_zero), @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', params.pi_vec, 'sig2_beta', params.sig2_beta)));
end
