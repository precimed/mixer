function result = BGMG_cpp_fit_bivariate(zmat, Nmat, params1, params2, options)
    % zmat      - matrix SNP x 2, z scores
    % Nmat      - number of subjects genotyped per SNP
    % params1   - univariate parameters for the first trait
    % params2   - univariate parameters for the second trait    
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.05; end;
    if ~isfield(options, 'ci_sample'), options.ci_sample = 10000; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    if ~isfield(options, 'fit_infinitesimal'), options.fit_infinitesimal = false; end;  % infinitesimal model (implemented only for univariate analysis)
    if ~isfield(options, 'fit_full_model'), options.fit_full_model = true; end;  % whether to fit full model
    if ~isfield(options, 'fit_with_constrains'), options.fit_with_constrains = true; end;  % whether to keep constrains from univariate model
    
    if any(Nmat(:) <= 0), error('Nmat values must be positive'); end;

    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
    calllib('bgmg', 'bgmg_clear_loglike_cache', 0); check();
 
    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    if size(zmat, 2) ~= 2, error('BGMG_cpp_fit_bivariate supports 2 traits'); end;

    fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), fminsearch_options));

    corrmat = corr(zmat(all(~isnan(zmat), 2), :));

    fprintf('Trait 1,2: fit infinitesimal model to find initial rho_zero and rho_beta (fast cost function)\n');
    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', 1); check();

    params_inft0 = struct('rho_zero', corrmat(1,2), 'rho_beta', corrmat(1,2));
    options_inft = struct('sig2_zero', [params1.sig2_zero; params2.sig2_zero], ...
                          'sig2_beta', [params1.pi_vec * params1.sig2_beta; params2.pi_vec * params2.sig2_beta], ...
                          'pi_vec', 1);
    params_inft  = fit(params_inft0, @(x)BGMG_util.BGMG_mapparams1_rho(x, options_inft));

    % Step2. Final unconstrained optimization (jointly on all parameters), using fast cost function
    func_map_params = @(x)BGMG_util.BGMG_mapparams3(x, struct('sig2_zero', [params1.sig2_zero; params2.sig2_zero], 'sig2_beta', [params1.sig2_beta; params2.sig2_beta], 'rho_zero', params_inft.rho_zero));
    if ~options.fit_with_constrains,  func_map_params = @(x)BGMG_util.BGMG_mapparams3(x); end

    fprintf('Trait 1,2: final unconstrained optimization (fast cost function)\n');

    params_final0 = params_inft;
    params_final0.sig2_beta = [params1.sig2_beta 0 params1.sig2_beta; 0 params2.sig2_beta params2.sig2_beta];
    init_pi12 = min(params1.pi_vec, params2.pi_vec)/exp(1);
    params_final0.pi_vec = [params1.pi_vec-init_pi12 params2.pi_vec-init_pi12 init_pi12];
    params_final0.rho_beta = [0 0 params_inft.rho_beta];

    params_fast = fit(params_final0, func_map_params);

    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', 0); check();
    if options.fit_full_model
        fprintf('Trait 1,2: final unconstrained optimization (full cost function)\n');
        result.params = fit(params_fast, func_map_params);
    else
        result.params = params_fast;
    end

    % Step3. Uncertainty estimation. 
    if ~isnan(options.ci_alpha)
        fprintf('Trait 1,2: uncertainty estimation\n');
        %ws=warning; warning('off', 'all');
        [ci_hess, ci_hess_err] = hessian(@(x)BGMG_util.BGMG_fminsearch_cost(BGMG_util.BGMG_mapparams3(x)), BGMG_util.BGMG_mapparams3(result.params));
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

            ci_sample = mvnrnd(BGMG_util.BGMG_mapparams3(result.params), inv(ci_hess), options.ci_sample);

            ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, ci_params{i} = BGMG_util.BGMG_mapparams3(ci_sample(i, :)); end;
        catch err
            fprintf('Error, %s\n', err.message);
        end

        result.ci_params = ci_params;
        [~, ci_bivariate_funcs] = BGMG_util.find_extract_funcs(options);
        result.ci = BGMG_util.extract_ci_funcs(ci_params, ci_bivariate_funcs, result.params, options.ci_alpha);
    end

    result.loglike_fit_trajectory = BGMG_util.extract_bivariate_loglike_trajectory();
    
    if options.verbose
       fprintf('Done, results are:\n');
       disp(BGMG_util.struct_to_display(result.params));
    end
end



