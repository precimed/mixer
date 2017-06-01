function results = GMM_ofrei_fit_observed(zmat, Hvec, Nmat, w_ld, ref_ld, options)
    % zmat      - matrix SNP x 2, z scores
    % Hvec      - vector SNP x 1, heterozigosity per SNP
    % Nmat      - number of subjects genotyped per SNP
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'zmax'), options.zmax = 12; end;           % throw away z scores larger than this value
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'fit_sbt_vec'), options.fit_sbt_vec = logspace(-3, 0, 25/2); end;
    if ~isfield(options, 'fit_pi_vec'),  options.fit_pi_vec  = logspace(-7, -0.5, 29/2); end;
    if ~isfield(options, 'fit_rho_vec'), options.fit_rho_vec = linspace(-0.9, 0.9, 20);  end;
    if ~isfield(options, 'fit_sigma0_vec'), options.fit_sigma0_vec = logspace(-0.1, 0.3, 25/2); end;    
    if ~isfield(options, 'plot_costlines'), options.plot_costlines = false; end;
    if ~isfield(options, 'plot_costmaps'), options.plot_costmaps = false; end;
    if ~isfield(options, 'calculate_beta_hat'), options.calculate_beta_hat = false; end;

    ntraits = size(zmat, 2);
    if ntraits > 2, error('Support only 1 or 2 traits'); end;
    
    % Fit univariate model for each trait
    for itrait=1:ntraits
        zvec = zmat(:, itrait);
        Nvec = Nmat(:, itrait);
        
        sbt_vec = options.fit_sbt_vec;
        pi_vec = options.fit_pi_vec;
        sbt_mat = repmat(sbt_vec, [length(pi_vec), 1]); sbt_mat=sbt_mat(:);
        pi_mat = repmat(pi_vec, [1, length(sbt_vec)]); pi_mat=pi_mat(:);

        costmap     = zeros(size(sbt_mat)); sigma0 = std(zvec(isfinite(zvec)));
        mapparams   = @GMM_ofrei_univariate_mapparams;
        costfuncUVT = @(x)GMM_ofrei_univariate_observed_cost(x, zvec, Hvec, Nvec, w_ld, ref_ld, mapparams, options);
        for i=1:length(sbt_mat),
            costmap(i) = costfuncUVT(struct('sigma_beta', sbt_mat(i), 'sigma0', sigma0, 'pivec', pi_mat(i)));
        end;
        
        costmap = reshape(costmap, [length(pi_vec), length(sbt_vec)]);
        if options.plot_costmaps, subplot(1,2,itrait); imagesc(log(1 + costmap(1:end, :)-min(costmap(:)))); colorbar; end;
        
        [i, j] = ind2sub(size(costmap), find(costmap(:) == min(costmap(:))));
        s0 = struct('sigma_beta', sbt_vec(j), 'sigma0', sigma0, 'pivec', pi_vec(i));
        if options.verbose, fprintf('pi=%.3e, sigma_beta^2=%.3e, sigma0^2=%.3f\n', s0.pivec, s0.sigma_beta.^2, s0.sigma0.^2); end
        sFMS = mapparams(fminsearch(costfuncUVT, mapparams(s0), struct('Display', 'on')));

        [~, result] = costfuncUVT(sFMS);
        results.univariate{itrait} = result;
        results.univariate{itrait}.costmap = costmap;
        results.univariate{itrait}.params  = sFMS;
        results.univariate{itrait}.grid_params = s0;
        if options.plot_costlines, plot_univariate_costline(sFMS, zvec, Hvec, Nvec, w_ld, ref_ld, options, itrait); end;
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
            mapparams   = @GMM_ofrei_bivariate_mapparams;
            costfuncBVT = @(x)GMM_ofrei_bivariate_observed_cost(x, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options);
        else
            % constrain sigma_beta and sigma0 to the univariate estimates;
            % optimize the remaining 5 parameters (3 in pivec, rho and rho0).
            s0          = struct('rho_beta', rho_beta, 'rho0', rho0, 'pivec', pivec);
            mapparams   = @(x)GMM_ofrei_bivariate_mapparams(x, struct('sigma_beta', sigma_beta, 'sigma0', sigma0));
            costfuncBVT = @(x)GMM_ofrei_bivariate_observed_cost(x, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options);
        end

        s_FMS = mapparams(fminsearch(costfuncBVT, mapparams(s0), struct('Display', 'on')));
        [~, result] = costfuncBVT(s_FMS);
        results.bivariate = result;
        results.bivariate.params   = s_FMS;
        if options.plot_costlines, plot_bivariate_costline(s_FMS, zmat, Hvec, Nmat, w_ld, ref_ld, options); end;
    end

    if options.verbose
       fprintf('Done\n');
    end
end

function plot_univariate_costline(s0, zvec, Hvec, Nvec, w_ld, ref_ld, options, plot_index)
    % Calculates cost along each parameter, starting with x0.
    costfuncUVT = @(x)GMM_ofrei_univariate_observed_cost(x, zvec, Hvec, Nvec, w_ld, ref_ld, @GMM_ofrei_univariate_mapparams, options);
    cost        = costfuncUVT(s0);

    for i=1:length(options.fit_pi_vec), s=s0;     s.pivec=options.fit_pi_vec(i);       costline.pivec(i) = costfuncUVT(s); end;
    for i=1:length(options.fit_sbt_vec), s=s0;    s.sigma_beta=options.fit_sbt_vec(i); costline.sigma_beta(i) = costfuncUVT(s); end;
    for i=1:length(options.fit_sigma0_vec), s=s0; s.sigma0=options.fit_sigma0_vec(i);  costline.sigma0(i) = costfuncUVT(s); end;

    subplot(3,5,plot_index); plot(log10(options.fit_pi_vec), log10(costline.pivec), '-', log10(s0.pivec), log10(cost), '*'); xlabel(sprintf('pi=%.3e', s0.pivec));
    subplot(3,5,plot_index+5); plot(log10(options.fit_sbt_vec), log10(costline.sigma_beta), '-', log10(s0.sigma_beta), log10(cost), '*'); xlabel(sprintf('sigma^2_b_e_t_a=%.3e', s0.sigma_beta.^2));
    subplot(3,5,plot_index+10); plot(      options.fit_sigma0_vec, log10(costline.sigma0), '-', s0.sigma0, log10(cost), '*'); xlabel(sprintf('sigma^2_0=%.3f', s0.sigma0.^2));
end

function plot_bivariate_costline(s0, zmat, Hvec, Nmat, w_ld, ref_ld, options)
    % Calculates cost along each parameter, starting with x0.
    costfuncBVT = @(s)GMM_ofrei_bivariate_observed_cost(s, zmat, Hvec, Nmat, w_ld, ref_ld, @GMM_ofrei_bivariate_mapparams, options);
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

