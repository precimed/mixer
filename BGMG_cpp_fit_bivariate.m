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

    mapparams_generic = @(x, params0)BGMG_util.BGMG_mapparams3_rho_and_pifrac(x, struct(...
        'sig2_zero', params0.sig2_zero, ...
        'sig2_beta', params0.sig2_beta(:, end), ...
        'pi_vec', [sum(params0.pi_vec([1, 3])), sum(params0.pi_vec([2, 3]))]));
    mapparams = @(x)mapparams_generic(x, params);
    
    fitfunc = @(x0)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), fminsearch_options));

    % Step2. Final unconstrained optimization (jointly on all parameters), using fast cost function
    BGMG_cpp.log('Trait 1,2: final optimization (full cost function)\n');
    result.params = fitfunc(params);
    result.loglike_fit_trajectory = bgmglib.extract_bivariate_loglike_trajectory();

    % Adjust based on smooth approximation
    if options.fit_full_model
        bgmglib.clear_loglike_cache();
        x0 = BGMG_util.BGMG_mapparams3_decorrelated_parametrization(result.params);
        arg_index=9;
        xgrid = linspace(x0(arg_index) - 2, x0(arg_index) + 2, 30)';
        for xi = 1:length(xgrid)
            x=x0; x(arg_index)=xgrid(xi);
            BGMG_util.BGMG_fminsearch_cost(BGMG_util.BGMG_mapparams3_decorrelated_parametrization(x));
        end
        result.loglike_adj_trajectory = bgmglib.extract_bivariate_loglike_trajectory();
        result.loglike_adj_trajectory.xgrid = xgrid;
        def = isfinite(xgrid+result.loglike_adj_trajectory.cost);
        x = xgrid(def);
        y = result.loglike_adj_trajectory.cost(def);
        try
            curve3 = fit(x, y, 'poly3');
            x0opt = fminsearch(@(x)curve3(x), x0(arg_index));
        catch
            x0opt=nan;
        end
        
        if isfinite(x0opt) && (x0opt > quantile(xgrid(def), 0.2)) && (x0opt < quantile(xgrid(def), 0.8))
            BGMG_cpp.log('Change BGMG solution (%.3f -> %.3f) based on smooth curve3 fit\n', x0(arg_index), x0opt);
            x0(arg_index) = x0opt;
            result.params = BGMG_util.BGMG_mapparams3_decorrelated_parametrization(x0);
        else
            BGMG_cpp.log('Refuse to change BGMG solution (%.3f -> %.3f) based on smooth curve3 fit\n', x0(arg_index), x0opt);
        end

        if 0
         figure(20);
         bt = result.loglike_adj_trajectory;
         subplot(1,2,1); plot(bt.pivec(:, 3) ./ min(sum(bt.pivec(:, [1 3]), 2), sum(bt.pivec(:, [2 3]), 2)), [bt.cost, curve3(xgrid)])
         subplot(1,2,2); plot(xgrid, [bt.cost, curve3(xgrid)])
        end
    end
    
    % Step3. Uncertainty estimation. 
    if ~isnan(options.ci_alpha)
        BGMG_cpp.log('Trait 1,2: uncertainty estimation\n');
        
        % Compined ci sample: 3 columns for trait 1 params, 3 columns for trait3 params, finally 3 columns for bivariate params
        result.ci_sample = nan(options.ci_sample, 9);
        ci_params = [];
        bgmglib.set_option('fast_cost', 1);

        try
            % find uncertainty from univariate optimization
            ugmg_mapparams = @(x)BGMG_util.UGMG_mapparams1_decorrelated_parametrization(x);
            for trait_index=1:2
                params_ugmg = struct('pi_vec', sum(result.params.pi_vec([trait_index 3])), 'sig2_zero',params.sig2_zero(trait_index), 'sig2_beta', params.sig2_beta(trait_index, 3));
                [ci_hess, ci_hess_err] = hessian(@(x)BGMG_util.UGMG_fminsearch_cost(ugmg_mapparams(x), trait_index), ugmg_mapparams(params_ugmg)); 
                result.ci_hess_ugmg{trait_index} = ci_hess;
                result.ci_hess_err_ugmg{trait_index} = ci_hess_err;
                ci_hess(diag(any(abs(inv(ci_hess)) > 1e10))) = +Inf;
                if any(~isfinite(ci_hess(:))), BGMG_cpp.log('Warning: unable to estimate hessian for confidence intervals'); end;
                result.ci_sample(:, (1:3) + 3*(trait_index-1)) = mvnrnd(ugmg_mapparams(params_ugmg), inv(ci_hess), options.ci_sample);
            end

            % find uncertainty from bivariate optimization
            mapparams = @(x)mapparams_generic(x, result.params);
            [ci_hess, ci_hess_err] = hessian(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(result.params));
            result.ci_hess = ci_hess;
            result.ci_hess_err = ci_hess_err;
            ci_hess(diag(any(abs(inv(ci_hess)) > 1e10))) = +Inf;
            if any(~isfinite(ci_hess(:))), BGMG_cpp.log('Warning: unable to estimate hessian for confidence intervals'); end;
            result.ci_sample(:, 7:9) = mvnrnd(mapparams(result.params), inv(ci_hess), options.ci_sample);

            ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, ci_params{i} = BGMG_util.BGMG_mapparams3_decorrelated_parametrization_9arguments(result.ci_sample(i, :)); end;
        catch err
            result.ci_err = err;
            BGMG_cpp.log('Error in BGMG uncertainty estimation, %s\n', err.message);
			for i=1:length(err.stack), BGMG_cpp.log('%s:%i %s\n', err.stack(i).file, err.stack(i).line, err.stack(i).name); end
        end

        result.ci_params = ci_params;
        [~, ci_bivariate_funcs] = BGMG_util.find_extract_funcs(options);
        result.ci = BGMG_util.extract_ci_funcs(ci_params, ci_bivariate_funcs, result.params, options.ci_alpha);
    end

    if options.verbose
       fprintf('Done, results are:\n');
       disp(BGMG_util.struct_to_display(result.params));
    end
end



