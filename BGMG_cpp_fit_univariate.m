function [result, options] = BGMG_cpp_fit_univariate(trait_index, params0, options)
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'TolX'), options.TolX = NaN; end;
    if ~isfield(options, 'TolFun'), options.TolFun = NaN; end;
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.05; end;
    if ~isfield(options, 'ci_sample'), options.ci_sample = 10000; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    if ~isfield(options, 'mapparams'), options.mapparams = @(x)BGMG_util.UGMG_mapparams1_decorrelated_parametrization(x); end;
    result = [];

    bgmglib = BGMG_cpp();
    
    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;
    if ~isnan(options.TolX), fminsearch_options.MaxFunEvals=options.TolX; end;
    if ~isnan(options.TolFun), fminsearch_options.MaxFunEvals=options.TolFun; end;

    fitfunc = @(x0)options.mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(options.mapparams(x), trait_index), options.mapparams(x0), fminsearch_options));

    fprintf('Final unconstrained optimization\n');
    bgmglib.clear_loglike_cache();
    result.params = fitfunc(params0);
    result.loglike_fit_trajectory = bgmglib.extract_univariate_loglike_trajectory();

    % For the last parameter, fit quadradic curve and take the minimum
    % Decide what range to include (filter out points that are too far to fit the curve)
    if 1
        bgmglib.clear_loglike_cache();
        x0 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(result.params);
        arg_index = 3;
        xgrid = linspace(x0(arg_index) - 1, x0(arg_index) + 1, 30)';
        for xi = 1:length(xgrid)
            x=x0; x(arg_index)=xgrid(xi);
            BGMG_util.UGMG_fminsearch_cost(options.mapparams(x), trait_index);
        end
        result.loglike_adj_trajectory = bgmglib.extract_univariate_loglike_trajectory();
        def = isfinite(xgrid+result.loglike_adj_trajectory.cost);
        x = xgrid(def);
        y = result.loglike_adj_trajectory.cost(def);
        curve3 = fit(x, y, 'poly3');
        x0opt = fminsearch(@(x)curve3(x), x0(arg_index));
        
        x0(arg_index) = x0opt;
        result.params = options.mapparams(x0);
        
        %figure;
        %subplot(1,2,1); plot(result.loglike_adj_trajectory.pivec, [result.loglike_adj_trajectory.cost, curve3(xgrid)], '.-');
        %subplot(1,2,2); plot(xgrid, [result.loglike_adj_trajectory.cost, curve3(xgrid)], '.-'); hold on; plot(x0opt, curve3(x0opt), '*k');
    end
    
    bgmglib.clear_loglike_cache();
    if ~isnan(options.ci_alpha)  % not implemented
        fprintf('Uncertainty estimation\n');
        %ws=warning; warning('off', 'all'); 
        [ci_hess, ci_hess_err] = hessian(@(x)BGMG_util.UGMG_fminsearch_cost(options.mapparams(x), trait_index), options.mapparams(result.params)); 
        result.loglike_ci_trajectory = bgmglib.extract_univariate_loglike_trajectory();
        result.ci_hess = ci_hess;
        result.ci_hess_err = ci_hess_err;
        %warning(ws);
        result.ci_params = [];
        try
            result.ci_hess(diag(any(abs(inv(result.ci_hess)) > 1e10))) = +Inf;
            ci_sample = mvnrnd(options.mapparams(result.params), inv(result.ci_hess), options.ci_sample);
            result.ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, result.ci_params{i} = options.mapparams(ci_sample(i, :)); end;
        catch err
            fprintf('Error, %s\n', err.message);
        end

        [ci_univariate_funcs, ~] = BGMG_util.find_extract_funcs(options);
        result.ci = BGMG_util.extract_ci_funcs(result.ci_params, ci_univariate_funcs, result.params, options.ci_alpha);
    end
end
