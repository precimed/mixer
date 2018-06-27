function result = BGMG_cpp_fit_univariate(zvec, Nvec, options)
    % zvat      - vector of length #SNP, z scores
    % Nmat      - number of subjects genotyped per SNP
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
    result = [];

    if any(Nvec(:) <= 0), error('Nmat values must be positive'); end;

    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));

    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    ntraits = size(zvec, 2);
    if ntraits ~= 1, error('BGMG_cpp_fit_univariate supports only 1 trait'); end;

    fit = @(x0, mapparams)mapparams(fminsearch(@(x)UGMG_fminsearch_cost(mapparams(x)), mapparams(x0), fminsearch_options));

    fprintf('Fit infinitesimal model to find initial sig2_zero (fast cost function)\n');
    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', 1); check();

    params_inft0 = struct('sig2_zero', var(zvec(~isnan(zvec))), 'sig2_beta', 1 ./ mean(Nvec(~isnan(Nvec))));
    params_inft  = fit(params_inft0, @(x)UGMG_mapparams1(x, struct('pi_vec', 1.0)));

    if options.fit_infinitesimal
        result.params = params_inft;
    else
        fprintf('Fit pi_vec and sig2_beta, constrained on sig2_zero and h2 (fast cost function)\n');
        params_mix0  = fit(struct('pi_vec', 0.01), @(x)UGMG_mapparams1_h2(x, params_inft.sig2_beta, params_inft.sig2_zero));

        fprintf('Final unconstrained optimization (fast cost function)\n');
        fast_params = fit(params_mix0, @(x)UGMG_mapparams1(x));

        calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', 0); check();
        if options.fit_full_model
            fprintf('Final unconstrained optimization (full cost function)\n');
            result.params = fit(fast_params, @(x)UGMG_mapparams1(x));
        else
            result.params = fast_params;
        end
    end

    if ~isnan(options.ci_alpha)  % not implemented
        fprintf('Uncertainty estimation\n');
        %ws=warning; warning('off', 'all'); 
        [ci_hess, ci_hess_err] = hessian(@(x)UGMG_fminsearch_cost(UGMG_mapparams1(x)), UGMG_mapparams1(result.params)); 
        result.ci_hess = ci_hess;
        result.ci_hess_err = ci_hess_err;
        %warning(ws);
        result.ci_params = [];
        try
            result.ci_hess(diag(any(abs(inv(result.ci_hess)) > 1e10))) = +Inf;
            ci_sample = mvnrnd(UGMG_mapparams1(result.params), inv(result.ci_hess), options.ci_sample);
            result.ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, result.ci_params{i} = UGMG_mapparams1(ci_sample(i, :)); end;
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

function cost = UGMG_fminsearch_cost(ov)
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, ov.pi_vec, ov.sig2_zero, ov.sig2_beta);
    fprintf('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost);
end
 
function ov = UGMG_mapparams1(iv, options)
    % mapparams for univariate mixture with a single causal component
    if ~exist('options', 'var'), options=[]; end;
    if ~isfield(options, 'pi_vec'), options.pi_vec = nan; end;
    if ~isfield(options, 'sig2_zero'), options.sig2_zero = nan; end;
    if ~isfield(options, 'sig2_beta'), options.sig2_beta = nan; end;

    is_packing = isstruct(iv); cnti = 1;
    if is_packing, ov = []; else ov = struct(); end;

    [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.logit_amd, 'pi_vec');
    [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_zero');
    [ov, ~] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
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
