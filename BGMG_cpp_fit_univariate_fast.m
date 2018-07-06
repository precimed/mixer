function params = BGMG_cpp_fit_univariate_fast(trait_index)
    bgmglib = BGMG_cpp();
    bgmglib.set_option('fast_cost', 1);
    zvec = bgmglib.get_zvec(trait_index);
    nvec = bgmglib.get_nvec(trait_index);
    
    fminsearch_options = struct('Display', 'on');
    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), fminsearch_options));

    fprintf('Fit infinitesimal model to find initial sig2_zero (fast cost function)\n');
    params_inft0 = struct('sig2_zero', var(zvec(~isnan(zvec))), 'sig2_beta', 1 ./ mean(nvec(~isnan(nvec))));
    params_inft  = fit(params_inft0, @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', 1.0)));

    fprintf('Fit pi_vec and sig2_beta, constrained on sig2_zero and h2 (fast cost function)\n');
    params_mix0  = fit(struct('pi_vec', 0.01), @(x)UGMG_mapparams1_h2(x, params_inft.sig2_beta, params_inft.sig2_zero));

    fprintf('Final unconstrained optimization (fast cost function)\n');
    params = fit(params_mix0, @(x)BGMG_util.UGMG_mapparams1(x));
end

function ov = UGMG_mapparams1_h2(iv, h2, sig2_zero)
    % Optimization with constrained product h2=pi_vec*sig2_beta and fixed sig2_zero
    
    if isstruct(iv),
        ov = BGMG_util.logit_amd(iv.pi_vec,0);
    else
        ov = struct();
        ov.pi_vec=BGMG_util.logit_amd(iv,1);
        ov.sig2_zero = sig2_zero;
        ov.sig2_beta = h2 / ov.pi_vec;
    end
end
