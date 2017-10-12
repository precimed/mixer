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
            params_mix0  = fit(struct('pi_vec', 0.5), @(x)UGMG_mapparams1(x, struct('sig2_zero', params_inft.sig2_zero, 'h2', params_inft.sig2_beta)));

            % Step3. Final unconstrained optimization (jointly on all parameters)
            results.univariate{itrait}.params = fit(params_mix0, @(x)UGMG_mapparams1(x));
        end

        % Step4. Uncertainty estimation.
        if ~isnan(options.ci_alpha)
            ws=warning; warning('off', 'all'); hess = hessian(@(x)BGMG_univariate_cost(UGMG_mapparams1(x), zvec, Hvec, Nvec, w_ld, ref_ld, options), UGMG_mapparams1(results.univariate{itrait}.params)); warning(ws);
            ci_sample = mvnrnd(UGMG_mapparams1(results.univariate{itrait}.params), inv(hess), options.ci_sample);

            ci_params = cell(options.ci_sample, 1);
            for i=1:options.ci_sample, ci_params{i} = UGMG_mapparams1(ci_sample(i, :)); end;

            [ci_funcs, ~] = find_extract_funcs(options);
            results.univariate{itrait}.ci = extract_ci_funcs(ci_params, ci_funcs, results.univariate{itrait}.params, options.ci_alpha);
        end
    end

    return;

    % Fit bivariate model
    if ntraits == 2
        params1 = results.univariate{1}.params;
        params2 = results.univariate{2}.params;

        sigma0      = [params1.sigma0, params2.sigma0];
        sigma_beta  = [params1.sigma_beta, params2.sigma_beta];
        pivec       = [params1.pivec, params2.pivec, sqrt(params1.pivec * params2.pivec)];
        rho0        = corr(zmat(isfinite(sum(zmat, 2)), :)); rho0 = rho0(1,2);
        rho_beta    = rho0;

        if 0
            % unconstrained optimization, adjusting all parameters
            s0          = struct('sigma_beta', sigma_beta, 'rho_beta', rho_beta, 'sigma0', sigma0, 'rho0', rho0, 'pivec', pivec);
            mapparams   = @BGMG_bivariate_mapparams;
            costfuncBVT = @(x)BGMG_bivariate_cost(x, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options);
        else
            % constrain sigma_beta and sigma0 to the univariate estimates;
            % optimize the remaining 5 parameters (3 in pivec, rho and rho0).
            s0          = struct('rho_beta', rho_beta, 'rho0', rho0, 'pivec', pivec);
            mapparams   = @(x)BGMG_bivariate_mapparams(x, struct('sigma_beta', sigma_beta, 'sigma0', sigma0));
            costfuncBVT = @(x)BGMG_bivariate_cost(x, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options);
        end

        s_FMS = mapparams(fminsearch(costfuncBVT, mapparams(s0), fminsearch_options));

        if ~isnan(options.alpha)
            costfuncBVT2 = @(x)BGMG_bivariate_cost(struct('sigma0', s_FMS.sigma0, 'sigma_beta', s_FMS.sigma_beta, 'pivec', x(1:3), 'rho0', x(4), 'rho_beta', x(5)), zmat, Hvec, Nmat, w_ld, ref_ld, [], options);
            x = [s_FMS.pivec, s_FMS.rho0, s_FMS.rho_beta];
            values = {'pi1', 'pi2', 'pi3', 'pi1_plus_pi3', 'pi2_plus_pi3', 'rho0', 'rho_beta'};
            tests  = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 1 0 1 0 0; 0 1 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
            ws=warning; warning('off', 'all'); [hess, err] = hessian_rescaled(costfuncBVT2, x); warning(ws);

            % Some elements in hessian might be misscalculated
            % due to various reasons:
            % - derivest library may return some hessian elements as NaN,
            %   for example due to very small values in the the pivec
            % - derivest library may return too high error bars
            % Below we detect which variables in hessian appear to be
            % miscalculated, save the list in nanidx, and ignore
            % such hessian elements in the code below.
            nanidx = false(size(hess, 1), 1)';
            for ii=1:length(nanidx),
                nanmat = isnan(hess) | abs(err ./ hess) > 1 | hess == 0;
                [m, i] = max(sum(nanmat) + sum(nanmat'));
                if m == 0, break; end;
                if (m > 0), nanidx(i) = true; hess(i, :) = +Inf; hess(:, i) = +Inf; end;
            end
            hess = hess(~nanidx, ~nanidx);

            c = -norminv(options.alpha/2);  % 1.96 for 95-confidence intervals
            for j = 1:length(values)
                if any(nanidx & tests(j, :))
                    se = nan;               % se = standard error
                else
                    se2 = tests(j, ~nanidx) * inv(hess) * tests(j, ~nanidx)';
                    if se2 >= 0, se = sqrt(se2); else se = nan; end;
                end

                value = sum(x(logical(tests(j, :))));
                s_FMS.uncertainty.(values{j}).value  = value;
                s_FMS.uncertainty.(values{j}).se     = se;
                s_FMS.uncertainty.(values{j}).ci     = [value - c*se, value + c*se];
                s_FMS.uncertainty.(values{j}).wald   = abs(value / se);
                s_FMS.uncertainty.(values{j}).pvalue = 2*normcdf(-abs(value) / se);
                s_FMS.uncertainty.(values{j}).H0     = sprintf('%s == 0', values{j});
            end
        end

        [~, result] = BGMG_bivariate_cost(s_FMS, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, discovery_options);
        results.bivariate = result;
        results.bivariate.params   = s_FMS;
        if options.plot_costlines, figure(2); qc_plot_bivariate_costline(s_FMS, zmat, Hvec, Nmat, w_ld, ref_ld, options); end;
    end

    if options.plot_costmaps, figure(1); saveas(gca, fullfile(options.outputdir, 'gmm_costmaps.png'), 'png'); end;
    if options.plot_costlines, figure(2); saveas(gca, fullfile(options.outputdir, 'gmm_costlines.png'), 'png'); end;

    if options.verbose
       fprintf('Done\n');
    end
end

function [univariate_ci_funcs, bivariate_ci_funcs] = find_extract_funcs(options)
    univariate_ci_funcs.sig2_zero        = @(params)(params.sig2_zero);
    univariate_ci_funcs.sig2_zero_minus1 = @(params)(params.sig2_zero - 1);
    univariate_ci_funcs.sig2_beta        = @(params)(params.sig2_beta);
    univariate_ci_funcs.pi_vec           = @(params)(params.pi_vec);
    univariate_ci_funcs.h2               = @(params)((params.pi_vec .* params.sig2_beta') * options.total_het);

    bivariate_ci_funcs = univariate_ci_funcs;
    bivariate_ci_funcs.rho_zero          = @(params)(params.rho_zero);
    bivariate_ci_funcs.rho_beta          = @(params)(params.rho_beta);
    bivariate_ci_funcs.pi1u              = @(params)(sum(params.pi_vec(1,3)));
    bivariate_ci_funcs.pi2u              = @(params)(sum(params.pi_vec(2,3)));
    bivariate_ci_funcs.rg                = @(params)(params.rho_beta * params.pi_vec(3) / sqrt(sum(params.pi_vec(1,3)) * sum(params.pi_vec(2,3))));
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
        ci_result.pval = reshape(2*normcdf(-abs(ci_result.mean / ci_result.se)), size(pe));

        ci.(ci_func_name) = ci_result;
    end
end

function ov = UGMG_mapparams1(iv, options)
    % mapparams for univariate mixture with a single causal component
    if ~exist('options', 'var'), options=[]; end;
    if ~isfield(options, 'pi_vec'), options.pi_vec = nan; end;
    if ~isfield(options, 'sig2_zero'), options.sig2_zero = nan; end;
    if ~isfield(options, 'sig2_beta'), options.sig2_beta = nan; end;

    % Hack-hack: setting "options.h2=pi_vec*sig2_beta" imply optimization
    % with constrained product pi_vec*sig2_beta.
    if ~isfield(options, 'h2'), options.h2 = nan; end;

    if isstruct(iv)
      ov = [];
      if isnan(options.sig2_zero), ov = cat(2,ov,BGMG_util.exp_amd(iv.sig2_zero,0)); end
      if isnan(options.sig2_beta) && isnan(options.h2), ov = cat(2,ov,BGMG_util.exp_amd(iv.sig2_beta,0)); end
      if isnan(options.pi_vec), ov = cat(2,ov,BGMG_util.logit_amd(iv.pi_vec,0)); end
    else
      ov = struct(); cnti = 1;
      if ~isnan(options.sig2_zero), ov.sig2_zero=options.sig2_zero; else ov.sig2_zero=BGMG_util.exp_amd(iv(cnti),1); cnti=cnti+1; end
      if ~isnan(options.sig2_beta) || ~isnan(options.h2), ov.sig2_beta=options.sig2_beta; else ov.sig2_beta=BGMG_util.exp_amd(iv(cnti),1); cnti=cnti+1; end
      if ~isnan(options.pi_vec), ov.pi_vec=options.pi_vec; else ov.pi_vec=BGMG_util.logit_amd(iv(cnti),1); cnti=cnti+1; end
      if ~isnan(options.h2), ov.sig2_beta = options.h2 ./ ov.pi_vec; end;
    end
end

% TBD
%   use BIC to prune redundant components; re-fit after pruning?
%   implement mapparams that constraint certain components (not everything)
%   posterior effect size estimation, for both bivariate and univariate
%   [DONE] show 9 figures with cost map along each component;
%   [????] show 9*8/2 = 36 figures with costmap along each pair of components
