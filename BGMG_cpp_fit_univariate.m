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

    BGMG_cpp.log('Final unconstrained optimization\n');
    bgmglib.clear_loglike_cache();
    result.params = fitfunc(params0);
    result.loglike_fit_trajectory = bgmglib.extract_univariate_loglike_trajectory();

    % For the last parameter, fit quadradic curve and take the minimum
    % Decide what range to include (filter out points that are too far to fit the curve)
    if options.fit_full_model
        bgmglib.clear_loglike_cache();
        x0 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(result.params);
        arg_index = 3;
        xgrid = linspace(x0(arg_index) - 1, x0(arg_index) + 1, 30)';
        for xi = 1:length(xgrid)
            x=x0; x(arg_index)=xgrid(xi);
            BGMG_util.UGMG_fminsearch_cost(options.mapparams(x), trait_index);
        end
        result.loglike_adj_trajectory = bgmglib.extract_univariate_loglike_trajectory();
        result.loglike_adj_trajectory.xgrid = xgrid;
        def = isfinite(xgrid+result.loglike_adj_trajectory.cost);
        x = xgrid(def);
        y = result.loglike_adj_trajectory.cost(def);
        try
            curve3 = fit(x, y, 'poly3');
            x0opt = fminsearch(@(x)curve3(x), x0(arg_index));
        catch
            x0opt = nan;
        end

        if isfinite(x0opt) && (x0opt > quantile(xgrid(def), 0.2)) && (x0opt < quantile(xgrid(def), 0.8))
            BGMG_cpp.log('Change UGMG solution (%.3f -> %.3f) based on smooth curve3 fit\n', x0(arg_index), x0opt);
            x0(arg_index) = x0opt;
            result.params = options.mapparams(x0);
        else
            BGMG_cpp.log('Refuse to change BGMG solution (%.3f -> %.3f) based on smooth curve3 fit\n', x0(arg_index), x0opt);
        end

        %figure;
        %subplot(1,2,1); plot(result.loglike_adj_trajectory.pivec, [result.loglike_adj_trajectory.cost, curve3(xgrid)], '.-');
        %subplot(1,2,2); plot(xgrid, [result.loglike_adj_trajectory.cost, curve3(xgrid)], '.-'); hold on; plot(x0opt, curve3(x0opt), '*k');
    end
    
    if ~isnan(options.ci_alpha)
        BGMG_cpp.log('Uncertainty estimation\n');
        bgmglib.clear_loglike_cache();
        bgmglib.set_option('fast_cost', 1);
        result.ci_params = [];
        try
            [ci_hess, ci_hess_info] = BGMG_util.hessian_robust(@(x)BGMG_util.UGMG_fminsearch_cost(options.mapparams(x), trait_index), options.mapparams(result.params)); 
            result.ci_hess = ci_hess; result.ci_hess_info = ci_hess_info;
            result.ci_sample = mvnrnd(options.mapparams(result.params), inv(ci_hess), options.ci_sample);
            result.ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, result.ci_params{i} = options.mapparams(result.ci_sample(i, :)); end;
        catch err
			result.ci_err = err;
            BGMG_cpp.log('Error, %s\n', err.message);
			for i=1:length(err.stack), BGMG_cpp.log('%s:%i %s\n', err.stack(i).file, err.stack(i).line, err.stack(i).name); end
        end

        [ci_univariate_funcs, ~] = BGMG_util.find_extract_funcs(options);
        result.ci = BGMG_util.extract_ci_funcs(result.ci_params, ci_univariate_funcs, result.params, options.ci_alpha);
    end
end
