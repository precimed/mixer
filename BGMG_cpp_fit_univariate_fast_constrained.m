function params = BGMG_cpp_fit_univariate_fast_constrained(trait_index, params)
    bgmglib = BGMG_cpp();
    bgmglib.set_option('fast_cost', 1);

    BGMG_cpp.log('Fit sig2_zero constrained to all other params (fast cost function)\n');
    params = BGMG_util.fit_mixer_model( ...
        struct('sig2_zero', params.sig2_zero), ... % x0
        @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', params.pi_vec, 'sig2_beta', params.sig2_beta)), ... % mapparams
        @(x)BGMG_util.UGMG_fminsearch_cost(x, trait_index), ... % costfunc
        struct('Display', 'on') ... % fminsearch options
    );
end
