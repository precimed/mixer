function plot_data = UGMG_power_plot(params, zvec, Hvec, Nvec, ref_ld, options)
    % Power plot. Proportion of narrow-sense chip heritability, S(N),
    % captured by genome-wide significant SNPs as a function of sample size.
    plot_data = [];
    if ~isfield(options, 'thin'), options.thin = 100; end; % speedup by random thining SNPs; setting this to 50 will take one-out-of-50 SNPs
    if ~isfield(options, 'nvalues'), options.nvalues = logspace(3, 8, 100); end;
    if ~isfield(options, 'zthresh'), options.zthresh = 5.45; end;  % genome-wide significance threshold 5e-8
    if ~isfield(options, 'title'), options.title = 'UNKNOWN TRAIT'; end;

    fprintf('Perform random thinning, exclude %i variants out of each %i variants\n', options.thin-1, options.thin); 
    cdf_idx     = find(isfinite(zvec + Nvec + Hvec));
    cdf_idx     = cdf_idx(1:options.thin:end);
    cdf_mask    = false(size(zvec));
    cdf_mask(cdf_idx) = true;

    zvec(~cdf_mask) = nan;

    options.verbose = true;
    options.calculate_z_cdf = true;
    options.calculate_z_cdf_limit = 300;  % must be really large, so that we don't lose any variance beyond the grid
    options.calculate_z_cdf_step = 0.2;
    options.calculate_z_cdf_weights = ones(size(zvec));  % no weighting (because of random thinning)
    options.calculate_z_cdf_nscale = options.nvalues ./ nanmedian(Nvec);
    options.calculate_z_cdf_func = @normpdf;
    options.poisson_sigma_grid_chunks = 20;
    options.use_poisson = true;
    w_ld_dummy = ones(size(zvec)); % dummy parameter --- it is used in cost calculation, which we are not interested in.
    [~, ugmg_cdf] = BGMG_univariate_cost(params, zvec, Hvec, Nvec, w_ld_dummy, ref_ld, options);
    
    expected_shape = [length(options.calculate_z_cdf_nscale), length(ugmg_cdf.cdf_z_grid) / length(options.calculate_z_cdf_nscale)];
    ugmg_cdf.cdf_z_grid = reshape(ugmg_cdf.cdf_z_grid, expected_shape);
    ugmg_cdf.cdf_nscale = reshape(ugmg_cdf.cdf_nscale, expected_shape);
    ugmg_cdf.cdf        = reshape(ugmg_cdf.cdf, expected_shape);

    if 0
        % normalize CDF to compensate for pdf lost in the tail (both options
        % are now disabled - best solution is to extend z_cdf_limit)
        for i=1:size(ugmg_cdf.cdf, 1)
            % Option 1: 
            % ugmg_cdf.cdf(i, end) = ugmg_cdf.cdf(i, end) +  (1/options.calculate_z_cdf_step - sum(ugmg_cdf.cdf(i, :)));

            % Option 2:
            % ugmg_cdf.cdf(i, :) = ugmg_cdf.cdf(i, :) / (sum(ugmg_cdf.cdf(i, :)) * options.calculate_z_cdf_step);
        end
    end

    delta_cond_z = @(z2, t)(max(0, z2 - 1));  % Empirical formula for faster computation.

    cdf = ugmg_cdf.cdf; ugmg_cdf.power_denominator = sum(delta_cond_z(ugmg_cdf.cdf_z_grid.^2, params.sig2_zero) .* cdf, 2);
    cdf = ugmg_cdf.cdf; cdf(abs(ugmg_cdf.cdf_z_grid) < options.zthresh) = 0; ugmg_cdf.power_numerator = sum(delta_cond_z(ugmg_cdf.cdf_z_grid.^2, params.sig2_zero) .* cdf, 2);

    plot_data.nvalues = options.nvalues;
    plot_data.prop_h2 = ugmg_cdf.power_numerator./ugmg_cdf.power_denominator; plot_data.prop_h2 = plot_data.prop_h2(:)';
    plot_data.nvalue_current = nanmedian(Nvec);
    plot_data.prop_h2_current = interp1(log10(plot_data.nvalues), plot_data.prop_h2, log10(plot_data.nvalue_current));
    plot_data.title = options.title;
end