function results = GMM2_fit(zmat, Hvec, Nmat, w_ld, ref_ld, options)
    % zmat      - matrix SNP x 2, z scores
    % Hvec      - vector SNP x 1, heterozigosity per SNP
    % Nmat      - number of subjects genotyped per SNP
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'zmax'), options.zmax = 12; end;           % throw away z scores larger than this value
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'fit_sbt_vec'), options.fit_sbt_vec = logspace(-3, 0, 13); end;
    if ~isfield(options, 'fit_pi_vec'),  options.fit_pi_vec  = logspace(-7, -0.5, 14); end;
    if ~isfield(options, 'fit_rho_vec'), options.fit_rho_vec = linspace(-0.9, 0.9, 20);  end;
    if ~isfield(options, 'fit_sigma0_vec'), options.fit_sigma0_vec = logspace(-0.1, 0.3, 25/2); end;    
    if ~isfield(options, 'alpha'), options.alpha = 0.05; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    if ~isfield(options, 'plot_costlines'), options.plot_costlines = false; end;
    if ~isfield(options, 'plot_costmaps'), options.plot_costmaps = false; end;
    if ~isfield(options, 'calculate_beta_hat'), options.calculate_beta_hat = false; end;
    if ~isfield(options, 'calculate_global_fdr'), options.calculate_global_fdr = false; end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'outputdir'), options.outputdir = '.'; end;

    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    discovery_options = options; discovery_options.zmax = +Inf;

    ntraits = size(zmat, 2);
    if ntraits > 2, error('Support only 1 or 2 traits'); end;
    
    % Fit univariate model for each trait
    % (user may skip this by providing options.params1 and options.params2 with univariate parameters for each trait)
    for itrait=1:ntraits
        if isfield(options, sprintf('params%i', itrait)),
            results.univariate{itrait}.params = options.(sprintf('params%i', itrait));
            continue;
        end;

        zvec = zmat(:, itrait);
        Nvec = Nmat(:, itrait);
        
        sbt_vec = options.fit_sbt_vec;
        pi_vec = options.fit_pi_vec;
        sbt_mat = repmat(sbt_vec, [length(pi_vec), 1]); sbt_mat=sbt_mat(:);
        pi_mat = repmat(pi_vec, [1, length(sbt_vec)]); pi_mat=pi_mat(:);

        costmap     = zeros(size(sbt_mat)); sigma0 = std(zvec(isfinite(zvec)));
        mapparams   = @GMM2_univariate_mapparams;
        costfuncUVT = @(x)GMM2_univariate_cost(x, zvec, Hvec, Nvec, w_ld, ref_ld, mapparams, options);
        for i=1:length(sbt_mat),
            costmap(i) = costfuncUVT(struct('sigma_beta', sbt_mat(i), 'sigma0', sigma0, 'pivec', pi_mat(i)));
        end;
        
        costmap = reshape(costmap, [length(pi_vec), length(sbt_vec)]);
        
        [i, j] = ind2sub(size(costmap), find(costmap(:) == min(costmap(:))));
        s0 = struct('sigma_beta', sbt_vec(j), 'sigma0', sigma0, 'pivec', pi_vec(i));
        if options.verbose, fprintf('pi=%.3e, sigma_beta^2=%.3e, sigma0^2=%.3f\n', s0.pivec, s0.sigma_beta.^2, s0.sigma0.^2); end
        sFMS = mapparams(fminsearch(costfuncUVT, mapparams(s0), fminsearch_options));

        if ~isnan(options.alpha)
            % Estimate p-value and confidence intervals
            % Note the following. The GMM2_univariate_cost accepts sigma0
            % and sigma_beta parameters, but we are interested in
            % confidence intervals for their squared values (sigma0^2 and
            % sigma_beta^2). This is already taken care below,
            costfuncUVT2 = @(x)GMM2_univariate_cost(struct('sigma0', sqrt(x(1)), 'sigma_beta', sqrt(x(2)), 'pivec', x(3)), zvec, Hvec, Nvec, w_ld, ref_ld, [], options);
            x = [sFMS.sigma0^2, sFMS.sigma_beta^2, sFMS.pivec];
            ws=warning; warning('off', 'all'); hess = hessian_rescaled(costfuncUVT2, x); warning(ws);

            se = sqrt(diag(inv(hess)));     % se = standard error
            c = -norminv(options.alpha/2);  % 1.96 for 95-confidence intervals

            sFMS.uncertainty.sigma0_sqr.value  = x(1);
            sFMS.uncertainty.sigma0_sqr.se     = se(1);
            sFMS.uncertainty.sigma0_sqr.ci     = [x(1) - c*se(1), x(1) + c*se(1)];
            sFMS.uncertainty.sigma0_sqr.wald   = abs((x(1)-1) / se(1));
            sFMS.uncertainty.sigma0_sqr.pvalue = 2*normcdf(-abs((x(1)-1) / se(1)));
            sFMS.uncertainty.sigma0_sqr.H0     = 'sigma0^2 == 1';

            sFMS.uncertainty.sigma_beta_sqr.value  = x(2);
            sFMS.uncertainty.sigma_beta_sqr.se     = se(2);
            sFMS.uncertainty.sigma_beta_sqr.ci     = [x(2) - c*se(2), x(2) + c*se(2)];
            sFMS.uncertainty.sigma_beta_sqr.wald   = abs(x(2) / se(2));
            sFMS.uncertainty.sigma_beta_sqr.pvalue = 2*normcdf(-abs(x(2) / se(2)));
            sFMS.uncertainty.sigma_beta_sqr.H0     = 'sigma_beta^2 == 0';

            sFMS.uncertainty.pi1.value  = x(3);
            sFMS.uncertainty.pi1.se     = se(3);
            sFMS.uncertainty.pi1.ci     = [x(3) - c*se(3), x(3) + c*se(3)];
            sFMS.uncertainty.pi1.wald   = abs(x(3) / se(3));
            sFMS.uncertainty.pi1.pvalue = 2*normcdf(-abs(x(3) / se(3)));
            sFMS.uncertainty.pi1.H0     = 'pi1 == 0';

            % Little magic to estimate heritability and its confidence interval
            if ~isnan(options.total_het)
                costfuncUVT3 = @(x)GMM2_univariate_cost(struct('sigma0', sqrt(x(1)), 'sigma_beta', sqrt(exp(x(2))), 'pivec', exp(x(3))), zvec, Hvec, Nvec, w_ld, ref_ld, [], options);
                x = [sFMS.sigma0^2, log(sFMS.sigma_beta^2), log(sFMS.pivec)];
                ws=warning; warning('off', 'all'); hess = hessian_rescaled(costfuncUVT3, x); warning(ws);
                h2_value = options.total_het * sFMS.sigma_beta^2 * sFMS.pivec;

                % Here, remember to multiply inv(hess) by h2 to get se(h2).
                % This is because inv(hess) estimates se of log(pi) and
                % log(sigma_beta^2). Let's denote t=log(pi)+log(sigma_beta^2),
                % and we are interested in exp(t). Now, by the definition
                % of deriviative,
                %   delta(exp(t)) = exp'(t) * delta(t)
                % and since exp'(t) = exp(t), we get
                %   delta(h2) = h2 * delta(t).
                se = h2_value * sqrt([0 1 1] * inv(hess) * [0 1 1]');

                sFMS.uncertainty.h2.value  = h2_value;
                sFMS.uncertainty.h2.se     = se;
                sFMS.uncertainty.h2.ci     = [h2_value - c * se, h2_value + c * se];
                sFMS.uncertainty.h2.wald   = abs(h2_value ./ se);
                sFMS.uncertainty.h2.pvalue = 2*normcdf(-abs(h2_value ./ se));
                sFMS.uncertainty.h2.H0     = 'h2 == 0';
            end
        end

        % Per-SNP estimates (posterior effect size, false discovery rates)
        [~, result] = GMM2_univariate_cost(sFMS,  zvec, Hvec, Nvec, w_ld, ref_ld, mapparams, discovery_options);
        results.univariate{itrait} = result;
        results.univariate{itrait}.costmap = costmap;
        results.univariate{itrait}.params  = sFMS;
        results.univariate{itrait}.grid_params = s0;
        if options.plot_costlines, figure(2); qc_plot_univariate_costline(sFMS, zvec, Hvec, Nvec, w_ld, ref_ld, options, itrait); end;
        
        if options.plot_costmaps
            figure(1); subplot(1,2,itrait); 
            imagesc(log(1 + costmap-min(costmap(:))), [0, 13]); colorbar;
            idx=find(floor(-log10(pi_vec)) == -log10(pi_vec)); tlabels={}; for i=idx, tlabels{end+1, 1}=sprintf('10^{-%i}',-log10(pi_vec(i))); end; yticks = idx; set(gca, 'YTick', yticks, 'YTickLabels', tlabels); ylabel('\pi_1');
            idx=find(floor(-log10(sbt_vec.^2)) == -log10(sbt_vec.^2)); tlabels={}; for i=idx, tlabels{end+1, 1}=sprintf('10^{-%i}',-log10(sbt_vec(i).^2)); end; xticks = idx; set(gca, 'XTick', xticks, 'XTickLabels', tlabels); xlabel('\sigma^2_{\beta}');
            title(sprintf('$$\\hat{\\pi_1}$$=%.2e, $$\\hat{\\sigma^2_{\\beta}}$$=%.2e', sFMS.pivec, sFMS.sigma_beta^2),'Interpreter','Latex')
        end        
    end
    
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
            mapparams   = @GMM2_bivariate_mapparams;
            costfuncBVT = @(x)GMM2_bivariate_cost(x, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options);
        else
            % constrain sigma_beta and sigma0 to the univariate estimates;
            % optimize the remaining 5 parameters (3 in pivec, rho and rho0).
            s0          = struct('rho_beta', rho_beta, 'rho0', rho0, 'pivec', pivec);
            mapparams   = @(x)GMM2_bivariate_mapparams(x, struct('sigma_beta', sigma_beta, 'sigma0', sigma0));
            costfuncBVT = @(x)GMM2_bivariate_cost(x, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options);
        end

        s_FMS = mapparams(fminsearch(costfuncBVT, mapparams(s0), fminsearch_options));

        if ~isnan(options.alpha)
            costfuncBVT2 = @(x)GMM2_bivariate_cost(struct('sigma0', s_FMS.sigma0, 'sigma_beta', s_FMS.sigma_beta, 'pivec', x(1:3), 'rho0', x(4), 'rho_beta', x(5)), zmat, Hvec, Nmat, w_ld, ref_ld, [], options);
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

        [~, result] = GMM2_bivariate_cost(s_FMS, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, discovery_options);
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

function [hess, err] = hessian_rescaled(f, x0)
    % hessian function from derivest suite has issues dealing with extreme
    % small parameter, like 1e-8 -- and that is typical for pivec.
    % hessian_rescaled makes a linear transformation to change the scale
    % of all input x variables, runs hessian from derives, and restore
    % the original scale. In cases when there are no small values in x0
    % hessian_rescaled(f, x0) should yield same result as hessian(f, x0).

    scale = x0;
    scale(x0 == 0) = 1;  % keep zeros as it is

    f_rescaled = @(x)(f(x .* scale));
    x0_rescaled = x0 ./ scale;

    [hess_rescaled, err_rescaled] = hessian(f_rescaled, x0_rescaled);
    hess = hess_rescaled ./ (scale(:) * scale(:)');
    err  = err_rescaled ./  (scale(:) * scale(:)');
end

function qc_plot_univariate_costline(s0, zvec, Hvec, Nvec, w_ld, ref_ld, options, plot_index)
    % Calculates cost along each parameter, starting with x0.
    costfuncUVT = @(x)GMM2_univariate_cost(x, zvec, Hvec, Nvec, w_ld, ref_ld, @GMM2_univariate_mapparams, options);
    cost        = costfuncUVT(s0);

    for i=1:length(options.fit_pi_vec), s=s0;     s.pivec=options.fit_pi_vec(i);       costline.pivec(i) = costfuncUVT(s); end;
    for i=1:length(options.fit_sbt_vec), s=s0;    s.sigma_beta=options.fit_sbt_vec(i); costline.sigma_beta(i) = costfuncUVT(s); end;
    for i=1:length(options.fit_sigma0_vec), s=s0; s.sigma0=options.fit_sigma0_vec(i);  costline.sigma0(i) = costfuncUVT(s); end;

    subplot(3,5,plot_index); plot(log10(options.fit_pi_vec), log10(costline.pivec), '-', log10(s0.pivec), log10(cost), '*'); xlabel(sprintf('pi=%.3e', s0.pivec));
    subplot(3,5,plot_index+5); plot(log10(options.fit_sbt_vec), log10(costline.sigma_beta), '-', log10(s0.sigma_beta), log10(cost), '*'); xlabel(sprintf('sigma^2_b_e_t_a=%.3e', s0.sigma_beta.^2));
    subplot(3,5,plot_index+10); plot(      options.fit_sigma0_vec, log10(costline.sigma0), '-', s0.sigma0, log10(cost), '*'); xlabel(sprintf('sigma^2_0=%.3f', s0.sigma0.^2));
end

function qc_plot_bivariate_costline(s0, zmat, Hvec, Nmat, w_ld, ref_ld, options)
    % Calculates cost along each parameter, starting with x0.
    costfuncBVT = @(s)GMM2_bivariate_cost(s, zmat, Hvec, Nmat, w_ld, ref_ld, @GMM2_bivariate_mapparams, options);
    cost        = costfuncBVT(s0);

    for itrait=1:2
        for i=1:length(options.fit_pi_vec), s=s0;     s.pivec(itrait)=options.fit_pi_vec(i);       costline.pivec(i) = costfuncBVT(s); end;
        for i=1:length(options.fit_sbt_vec), s=s0;    s.sigma_beta(itrait)=options.fit_sbt_vec(i); costline.sigma_beta(i) = costfuncBVT(s); end;
        for i=1:length(options.fit_sigma0_vec), s=s0; s.sigma0(itrait)=options.fit_sigma0_vec(i);  costline.sigma0(i) = costfuncBVT(s); end;
        subplot(3,5,itrait+2); plot(log10(options.fit_pi_vec), log10(costline.pivec), '-', log10(s0.pivec(itrait)), log10(cost), '*'); xlabel(sprintf('pi_%i=%.3e', itrait, s0.pivec(itrait)));
        subplot(3,5,itrait+2+5); plot(log10(options.fit_sbt_vec), log10(costline.sigma_beta), '-', log10(s0.sigma_beta(itrait)), log10(cost), '*'); xlabel(sprintf('sigma^2_b_e_t_a_,_%i=%.3e', itrait, s0.sigma_beta(itrait).^2));
        subplot(3,5,itrait+2+10); plot(      options.fit_sigma0_vec, log10(costline.sigma0), '-', s0.sigma0(itrait), log10(cost), '*'); xlabel(sprintf('sigma^2_0_,_%i=%.3f', itrait, s0.sigma0(itrait).^2));
    end

    for i=1:length(options.fit_pi_vec), s=s0;     s.pivec(3)=options.fit_pi_vec(i);  costline.pivec(i) = costfuncBVT(s); end;
    for i=1:length(options.fit_rho_vec), s=s0;    s.rho_beta=options.fit_rho_vec(i); costline.rho_beta(i) = costfuncBVT(s); end;
    for i=1:length(options.fit_rho_vec), s=s0;    s.rho0=options.fit_rho_vec(i); costline.rho0(i) = costfuncBVT(s); end;
    subplot(3,5,5); plot(log10(options.fit_pi_vec), log10(costline.pivec), '-', log10(s0.pivec(3)), log10(cost), '*'); xlabel(sprintf('pi_3=%.3e', s0.pivec(3)));
    subplot(3,5,10); plot(      options.fit_rho_vec, log10(costline.rho_beta), '-', s0.rho_beta, log10(cost), '*'); xlabel(sprintf('rho_b_e_t_a=%.3f', s0.rho_beta));
    subplot(3,5,15); plot(      options.fit_rho_vec, log10(costline.rho0), '-', s0.rho0, log10(cost), '*'); xlabel(sprintf('rho0=%.3f', s0.rho0));
end

% TBD
%   use BIC to prune redundant components; re-fit after pruning?
%   implement mapparams that constraint certain components (not everything)
%   posterior effect size estimation, for both bivariate and univariate
%   [DONE] show 9 figures with cost map along each component;
%   [????] show 9*8/2 = 36 figures with costmap along each pair of components

