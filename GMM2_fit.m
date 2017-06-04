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

