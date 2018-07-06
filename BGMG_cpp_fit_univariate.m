function result = BGMG_cpp_fit_univariate(trait_index, params0, options)
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.05; end;
    if ~isfield(options, 'ci_sample'), options.ci_sample = 10000; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    result = [];

    bgmglib = BGMG_cpp();
    
    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), fminsearch_options));

    fprintf('Final unconstrained optimization\n');
    bgmglib.clear_loglike_cache();
    result.params = fit(params0, @(x)BGMG_util.UGMG_mapparams1(x));
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
end
