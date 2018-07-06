function result = BGMG_cpp_fit_univariate(trait_index, options)
    % options   - options; see the below for a full list of configurable options
    %

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.05; end;
    if ~isfield(options, 'ci_sample'), options.ci_sample = 10000; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    if ~isfield(options, 'fit_infinitesimal'), options.fit_infinitesimal = false; end;  % infinitesimal model (implemented only for univariate analysis)
    if ~isfield(options, 'fit_full_model'), options.fit_full_model = true; end;  % whether to fit full model
    if ~isfield(options, 'params0'), options.params0 = []; end; % initial approximation
    result = [];

    bgmglib = BGMG_cpp();
    zvec = bgmglib.get_zvec(trait_index);
    nvec = bgmglib.get_nvec(trait_index);
    
    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), fminsearch_options));

    % Try to find initial approximation
    if isempty(options.params0)
        fprintf('Fit infinitesimal model to find initial sig2_zero (fast cost function)\n');
        bgmglib.set_option('fast_cost', 1);

        params_inft0 = struct('sig2_zero', var(zvec(~isnan(zvec))), 'sig2_beta', 1 ./ mean(nvec(~isnan(nvec))));
        params_inft  = fit(params_inft0, @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', 1.0)));
        
        % Stop at infinitesimal model (if requested by user)
        if options.fit_infinitesimal
            result.params = params_inft;
            return;
        end

        fprintf('Fit pi_vec and sig2_beta, constrained on sig2_zero and h2 (fast cost function)\n');
        params_mix0  = fit(struct('pi_vec', 0.01), @(x)UGMG_mapparams1_h2(x, params_inft.sig2_beta, params_inft.sig2_zero));

        fprintf('Final unconstrained optimization (fast cost function)\n');
        fast_params = fit(params_mix0, @(x)BGMG_util.UGMG_mapparams1(x));

        % Stop after fitting model with constrained cost function (if requested by user)
        if ~options.fit_full_model
            result.params = fast_params;
            return
        end

        options.params0 = fast_params;
    end
    
    fprintf('Final unconstrained optimization (full cost function)\n');
    bgmglib.set_option('fast_cost', 0);
    bgmglib.clear_loglike_cache();
    result.params = fit(options.params0, @(x)BGMG_util.UGMG_mapparams1(x));
    result.loglike_fit_trajectory = bgmglib.extract_univariate_loglike_trajectory();

    bgmglib.clear_loglike_cache();
    if ~isnan(options.ci_alpha)  % not implemented
        fprintf('Uncertainty estimation\n');
        %ws=warning; warning('off', 'all'); 
        [ci_hess, ci_hess_err] = hessian(@(x)BGMG_util.UGMG_fminsearch_cost(BGMG_util.UGMG_mapparams1(x), trait_index), BGMG_util.UGMG_mapparams1(result.params)); 
        result.loglike_ci_trajectory = bgmglib.extract_univariate_loglike_trajectory();
        result.ci_hess = ci_hess;
        result.ci_hess_err = ci_hess_err;
        %warning(ws);
        result.ci_params = [];
        try
            result.ci_hess(diag(any(abs(inv(result.ci_hess)) > 1e10))) = +Inf;
            ci_sample = mvnrnd(BGMG_util.UGMG_mapparams1(result.params), inv(result.ci_hess), options.ci_sample);
            result.ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, result.ci_params{i} = BGMG_util.UGMG_mapparams1(ci_sample(i, :)); end;
        catch err
            fprintf('Error, %s\n', err.message);
        end

        [ci_univariate_funcs, ~] = BGMG_util.find_extract_funcs(options);
        result.ci = BGMG_util.extract_ci_funcs(result.ci_params, ci_univariate_funcs, result.params, options.ci_alpha);
    end

    if options.verbose
       fprintf('Done, results are:\n');
       disp(BGMG_util.struct_to_display(result.params));
    end
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
