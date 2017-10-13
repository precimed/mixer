function results = BGMG_fit2(zmat, Hvec, Nmat, w_ld, ref_ld, options)
    % zmat      - matrix SNP x 2, z scores
    % Hvec      - vector SNP x 1, heterozigosity per SNP
    % Nmat      - number of subjects genotyped per SNP
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.05; end;
    if ~isfield(options, 'ci_sample'), options.ci_sample = 10000; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    
    if any(Nmat(:) <= 0), error('Nmat values must be positive'); end;

    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    ntraits = size(zmat, 2);
    if ntraits > 2, error('Support only 1 or 2 traits'); end;

    % Fit univariate model for each trait
    % (user may skip this by providing options.params1 and
    % options.params2 with univariate parameters for each trait)
    for itrait=1:ntraits
        zvec = zmat(:, itrait);
        Nvec = Nmat(:, itrait);

        if isfield(options, sprintf('params%i', itrait)),
            results.univariate{itrait}.params = options.(sprintf('params%i', itrait));
        else
            fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_univariate_cost(mapparams(x), zvec, Hvec, Nvec, w_ld, ref_ld, options), mapparams(x0), fminsearch_options));

            % Step1. Fit infinitesimal model (pivec=1) to initialize sig2_zero
            params_inft0 = struct('sig2_zero', var(zvec(~isnan(zvec))), 'sig2_beta', 1 ./ mean(Nvec(~isnan(Nvec))));
            params_inft  = fit(params_inft0, @(x)UGMG_mapparams1(x, struct('pi_vec', 1.0)));

            % Step2. Fit pi_vec and sig2_beta, constrained on sig2_zero and product h2=pi_vec*sig2_beta from x_inft
            params_mix0  = fit(struct('pi_vec', 0.5), @(x)UGMG_mapparams1_h2(x, params_inft.sig2_beta, params_inft.sig2_zero));

            % Step3. Final unconstrained optimization (jointly on all parameters)
            results.univariate{itrait}.params = fit(params_mix0, @(x)UGMG_mapparams1(x));
 
            % Step4. Uncertainty estimation.
            if ~isnan(options.ci_alpha)
                ws=warning; warning('off', 'all'); hess = hessian(@(x)BGMG_univariate_cost(UGMG_mapparams1(x), zvec, Hvec, Nvec, w_ld, ref_ld, options), UGMG_mapparams1(results.univariate{itrait}.params)); warning(ws);
                ci_sample = mvnrnd(UGMG_mapparams1(results.univariate{itrait}.params), inv(hess), options.ci_sample);

                ci_params = cell(options.ci_sample, 1);
                for i=1:options.ci_sample, ci_params{i} = UGMG_mapparams1(ci_sample(i, :)); end;

                [ci_univariate_funcs, ~] = find_extract_funcs(options);
                results.univariate{itrait}.ci = extract_ci_funcs(ci_params, ci_univariate_funcs, results.univariate{itrait}.params, options.ci_alpha);
            end
        end
    end

    % Fit bivariate model
    if ntraits == 2
        p1 = results.univariate{1}.params;
        p2 = results.univariate{2}.params;
        
        fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_bivariate_cost(mapparams(x), zmat, Hvec, Nmat, w_ld, ref_ld, options), mapparams(x0), fminsearch_options));

        corrmat = corr(zmat(all(~isnan(zmat), 2), :));
        
        % Step1. Fit infinitesimal model (pivec=1) to initialize rho_zero and rho_beta
        params_inft0 = struct('rho_zero', corrmat(1,2), 'rho_beta', corrmat(1,2));
        options_inft = struct('sig2_zero', [p1.sig2_zero; p2.sig2_zero], ...
                              'sig2_beta', [p1.pi_vec * p1.sig2_beta; p2.pi_vec * p2.sig2_beta], ...
                              'pi_vec', 1);
        params_inft  = fit(params_inft0, @(x)BGMG_mapparams1_rho(x, options_inft));
        
        % Step2. Final unconstrained optimization (jointly on all parameters)
        params_final0 = params_inft;
        params_final0.sig2_beta = [p1.sig2_beta 0 p1.sig2_beta; 0 p2.sig2_beta p2.sig2_beta];
        params_final0.pi_vec = [p1.pi_vec p2.pi_vec min(p1.pi_vec, p2.pi_vec)];
        params_final0.rho_beta = [0 0 params_inft.rho_beta];
        results.bivariate.params = fit(params_final0, @(x)BGMG_mapparams3(x));

        % Step3. Uncertainty estimation. 
        if ~isnan(options.ci_alpha)
            ws=warning; warning('off', 'all'); hess = hessian(@(x)BGMG_bivariate_cost(BGMG_mapparams3(x), zmat, Hvec, Nmat, w_ld, ref_ld, options), BGMG_mapparams3(results.bivariate.params)); warning(ws);
            ci_sample = mvnrnd(BGMG_mapparams3(results.bivariate.params), inv(hess), options.ci_sample);

            ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, ci_params{i} = BGMG_mapparams3(ci_sample(i, :)); end;

            [~, ci_bivariate_funcs] = find_extract_funcs(options);
            results.bivariate.ci = extract_ci_funcs(ci_params, ci_bivariate_funcs, results.bivariate.params, options.ci_alpha);
        end
    end

    if options.verbose
       fprintf('Done\n');
    end
end

function [univariate_ci_funcs, bivariate_ci_funcs] = find_extract_funcs(options)
    univariate_ci_funcs.sig2_zero        = @(params)(params.sig2_zero);
    univariate_ci_funcs.sig2_zero_minus1 = @(params)(params.sig2_zero - 1);
    univariate_ci_funcs.sig2_beta        = @(params)(params.sig2_beta);
    univariate_ci_funcs.pi_vec           = @(params)(params.pi_vec);
    univariate_ci_funcs.h2               = @(params)((params.sig2_beta*params.pi_vec')*options.total_het);

    bivariate_ci_funcs = univariate_ci_funcs;
    bivariate_ci_funcs.rho_zero          = @(params)(params.rho_zero);
    bivariate_ci_funcs.rho_beta          = @(params)(params.rho_beta);
    bivariate_ci_funcs.pi1u              = @(params)(sum(params.pi_vec([1 3])));
    bivariate_ci_funcs.pi2u              = @(params)(sum(params.pi_vec([2 3])));
    bivariate_ci_funcs.rg                = @(params)(params.rho_beta * params.pi_vec(3) / sqrt(sum(params.pi_vec([1 3])) * sum(params.pi_vec([2 3]))));
end

function ci = extract_ci_funcs(ci_params, ci_funcs, params, ci_alpha)
    ci_func_names = fieldnames(ci_funcs);
    for i=1:length(ci_func_names)
        ci_func_name = ci_func_names{i};
        ci_func      = ci_funcs.(ci_func_name);

        pe = ci_func(params);  % pe = point estimate
        dist = nan(length(ci_params), numel(pe));
        for j=1:length(ci_params), dist(j, :) = BGMG_util.rowvec(ci_func(ci_params{j})); end

        ci_result.point_estimate = pe;
        ci_result.mean = reshape(mean(dist), size(pe));
        ci_result.median = reshape(median(dist), size(pe));
        ci_result.lower = reshape(quantile(dist,     ci_alpha/2), size(pe));
        ci_result.upper = reshape(quantile(dist, 1 - ci_alpha/2), size(pe));
        ci_result.se = reshape(std(dist), size(pe));
        ci_result.pval = reshape(2*normcdf(-abs(ci_result.mean ./ ci_result.se)), size(pe));

        ci.(ci_func_name) = ci_result;
    end
end

function [ov, cnti] = mapparams(iv, ov, cnti, options, transform, field)
    idx_pack = isnan(options.(field));
    transform_forward = 0;
    transform_backward = 1;
    if isstruct(iv)
        % transform from struct to vector
        if any(idx_pack)
            ov = cat(2,ov,BGMG_util.rowvec(transform(iv.(field)(idx_pack),transform_forward)));
        end
    else
        % transform from vector to struct
        ov.(field) = options.(field);
        if any(idx_pack)
            ov.(field)(idx_pack) = transform(iv(cnti : (cnti + sum(idx_pack) - 1)), transform_backward);
            cnti=cnti+sum(idx_pack);
        end
    end
end

function ov = UGMG_mapparams1(iv, options)
    % mapparams for univariate mixture with a single causal component
    if ~exist('options', 'var'), options=[]; end;
    if ~isfield(options, 'pi_vec'), options.pi_vec = nan; end;
    if ~isfield(options, 'sig2_zero'), options.sig2_zero = nan; end;
    if ~isfield(options, 'sig2_beta'), options.sig2_beta = nan; end;

    is_packing = isstruct(iv); cnti = 1;
    if is_packing, ov = []; else ov = struct(); end;
    
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.logit_amd, 'pi_vec');
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_zero');
    [ov, ~] = mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
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

function ov = BGMG_mapparams1_rho(iv, options)
    % mapparams for bivaraite mixture with a one causal pleiotropic component
    % all params except rho_zero and rho_beta must be fixed by options
    if ~isfield(options, 'rho_zero'), options.rho_zero = nan; end;
    if ~isfield(options, 'rho_beta'), options.rho_beta = nan; end;

    is_packing = isstruct(iv); cnti = 1;
    if is_packing, ov = []; else ov = options; end;
    
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.erf_amd, 'rho_zero');
    [ov, ~] = mapparams(iv, ov, cnti, options, @BGMG_util.erf_amd, 'rho_beta');
end

function ov = BGMG_mapparams3(iv, options)
    % mapparams for saturated bivaraite mixture with a three causal component
    % (trait1-specific, trait2-specific, and pleiotropic components)
    % pleiotropic component re-use the same sig2_beta as trait-specific
    % components.

    if ~exist('options', 'var'), options=[]; end;
    if ~isfield(options, 'pi_vec'), options.pi_vec = [nan nan nan]; end;
    if ~isfield(options, 'sig2_zero'), options.sig2_zero = [nan; nan]; end;
    if ~isfield(options, 'sig2_beta'), options.sig2_beta = [nan; nan]; end;
    if ~isfield(options, 'rho_zero'), options.rho_zero = nan; end;
    if ~isfield(options, 'rho_beta'), options.rho_beta = [0 0 nan]; end;

    is_packing = isstruct(iv); cnti = 1;
    if is_packing, ov = []; else ov = struct(); end;
    
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.logit_amd, 'pi_vec');
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_zero');
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.erf_amd, 'rho_zero');
    [ov, cnti] = mapparams(iv, ov, cnti, options, @BGMG_util.erf_amd, 'rho_beta');
    
    % Tricks to map sig2_beta from two-vector into 2x3 matrix with two zero
    % elements and equal variances in trait-specific and pleiotropic components.
    if is_packing
        iv.sig2_beta = iv.sig2_beta(:,3);
        [ov, ~] = mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
    else
        [ov, ~] = mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
        s1 = ov.sig2_beta(1); s2 = ov.sig2_beta(2);
        ov.sig2_beta = [s1 0 s1; 0 s2 s2];
    end
end
