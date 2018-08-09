function result = BGMG_cpp_fit_bivariate(params, options)
    % params    - bivariate params (struct with sig2_zero, sig2_beta, pi_vec, rho_zero, rho_beta)
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'TolX'), options.TolX = NaN; end;
    if ~isfield(options, 'TolFun'), options.TolFun = NaN; end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.05; end;
    if ~isfield(options, 'ci_sample'), options.ci_sample = 10000; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    
    bgmglib = BGMG_cpp();
    bgmglib.clear_loglike_cache();
 
    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;
    if ~isnan(options.TolX), fminsearch_options.MaxFunEvals=options.TolX; end;
    if ~isnan(options.TolFun), fminsearch_options.MaxFunEvals=options.TolFun; end;

    mapparams = @(x)BGMG_util.BGMG_mapparams3_rho_and_pifrac(x, struct(...
        'sig2_zero', params.sig2_zero, ...
        'sig2_beta', params.sig2_beta(:, end), ...
        'pi_vec', [sum(params.pi_vec([1, 3])), sum(params.pi_vec([2, 3]))]));

    fit = @(x0)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), fminsearch_options));

    % Step2. Final unconstrained optimization (jointly on all parameters), using fast cost function
    fprintf('Trait 1,2: final optimization (full cost function)\n');
    result.params = fit(params);

    % Step3. Uncertainty estimation. 
    if ~isnan(options.ci_alpha)
        fprintf('Trait 1,2: uncertainty estimation\n');
        %ws=warning; warning('off', 'all');
        [ci_hess, ci_hess_err] = hessian(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(result.params));
        result.ci_hess = ci_hess;
        result.ci_hess_err = ci_hess_err;
        %warning(ws);

        ci_params = [];
        try
            % Some elements in hessian might be misscalculated
            % due to various reasons:
            % - derivest library may return some hessian elements as NaN,
            %   for example due to very small values in the the pivec
            % - derivest library may return too high error bars
            % Below we detect which variables in hessian appear to be
            % miscalculated, save the list in nanidx, and ignore
            % such hessian elements in the code below.
            ci_hess(diag(any(abs(inv(ci_hess)) > 1e10))) = +Inf;

            ci_sample = mvnrnd(mapparams(result.params), inv(ci_hess), options.ci_sample);

            ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, ci_params{i} = mapparams(ci_sample(i, :)); end;
        catch err
            fprintf('Error, %s\n', err.message);
        end

        result.ci_params = ci_params;
        [~, ci_bivariate_funcs] = BGMG_util.find_extract_funcs(options);
        result.ci = BGMG_util.extract_ci_funcs(ci_params, ci_bivariate_funcs, result.params, options.ci_alpha);
    end

    result.loglike_fit_trajectory = bgmglib.extract_bivariate_loglike_trajectory();
    
    if options.verbose
       fprintf('Done, results are:\n');
       disp(BGMG_util.struct_to_display(result.params));
    end
end



