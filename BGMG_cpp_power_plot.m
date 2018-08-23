function plot_data = BGMG_cpp_power_plot(params, trait_index, options)
    % QQ plot for data and model

    if ~isfield(options, 'title'), options.title = 'UNKNOWN TRAIT'; end;
    if ~isfield(options, 'downscale'), options.downscale = 10; end;
    if ~isfield(options, 'power_nvec'), options.power_nvec = 10.^(3:0.1:8); end;
    plot_data = []; hold on;

    bgmglib = BGMG_cpp();
    weights_bgmg = bgmglib.weights;
    zvec = bgmglib.get_zvec(trait_index);
    
    % Check that weights and zvec are all defined
    if any(~isfinite(zvec)), error('all values must be defined'); end;
    if any(~isfinite(weights_bgmg)), error('all values must be defined'); end;

    % We may limit QQ plots to certain group of SNPs (let's say a bin of mafvec and TLD)
    % Such subset should be defined via mask (same convention as defvec, 1 in mask means "use for QQ plot)
    if ~isfield(options, 'mask'), options.mask = true(size(zvec)); end;

    % To speedup calculations we set temporary weights where many elements
    % are zeroed. For data QQ plots all elements are set to 1.
    % At the end of this function we restore weights_bgmg.
    downscale_mask_indices = false(sum(options.mask), 1); 
    downscale_mask_indices(1:options.downscale:end)=true;
    mask_indices = find(options.mask); downscale_indices = mask_indices(downscale_mask_indices);
    
    model_weights = zeros(size(weights_bgmg)); 
    model_weights(downscale_indices) = weights_bgmg(downscale_indices); 
    model_weights = model_weights ./ sum(model_weights);
    bgmglib.weights = model_weights;
    
    % Calculate data for power plot
    power_svec = bgmglib.calc_univariate_power(trait_index, params.pi_vec, params.sig2_zero, params.sig2_beta, 5.45, options.power_nvec);
    plot_data.power_nvec = options.power_nvec;
    plot_data.power_svec = power_svec;

    % Restore original weights
    bgmglib.weights = weights_bgmg;

    plot(log10(options.power_nvec), power_svec, '-', 'LineWidth',1); hold on;
end
