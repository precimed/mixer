function result = BGMG_cpp_fit_kl_univariate(zvec, Nvec, options)
    % zvat      - vector of length #SNP, z scores
    % Nmat      - number of subjects genotyped per SNP
    % options   - options; see the below for a full list of configurable options
    %

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'MaxFunEvals'), options.MaxFunEvals = NaN; end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'params0'), options.params0 = []; end; % initial approximation
    result = [];

    if any(Nvec(:) <= 0), error('Nmat values must be positive'); end;

    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));

    fminsearch_options = struct('Display', 'on');
    if ~isnan(options.MaxFunEvals), fminsearch_options.MaxFunEvals=options.MaxFunEvals; end;

    if isempty(options.params0), error('options.param0 is required'); end;

    % Retrieve weights from c++ library
    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
    pBuffer = libpointer('singlePtr', zeros(length(zvec), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', 0, length(zvec), pBuffer);  check(); 
    weights_bgmg = pBuffer.Value;
    clear pBuffer
    
    data_weights = weights_bgmg; data_weights = data_weights/sum(data_weights);
    model_weights = weights_bgmg;

    hv_z = single(-15:0.05:0); 
    hv_logp = -log10(2*normcdf(hv_z));

    % Calculate data_logpvec
    [data_y, si] = sort(-log10(2*normcdf(-abs(zvec))));
    data_x=-log10(cumsum(data_weights(si),1,'reverse'));
    data_idx = ([data_y(2:end); +Inf] ~= data_y);
    data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

    % Calculate model_logpvec
   % model_logpvec = calc_model_logpvec(options.params0, hv_z, model_weights);
    
% clf; hold on
  %  hData  = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
 %  hModel = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;

    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', 0); check();
    kl_cost_func = @(params)klcost(params, data_logpvec, model_weights, hv_z);
    fit = @(x0, mapparams)mapparams(fminsearch(@(x)kl_cost_func(mapparams(x)), mapparams(x0), fminsearch_options));
    result.params = fit(options.params0, @(x)BGMG_util.UGMG_mapparams1(x));
    
    if options.verbose
       fprintf('Done, results are:\n');
       disp(BGMG_util.struct_to_display(result.params));
    end
end

function cost = UGMG_fminsearch_cost(ov)
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, ov.pi_vec, ov.sig2_zero, ov.sig2_beta);
    fprintf('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost);
end

function cost = klcost(ov, data_logpvec, model_weights, hv_z)
   model_logpvec = calc_model_logpvec(ov, hv_z, model_weights);
   % cost = nansum((10.^data_logpvec) .* abs(data_logpvec - model_logpvec));
   diff = model_logpvec - data_logpvec;
   cost = nansum(diff.*diff);git 
    fprintf('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost);
    if 1 
      clf; hold on
        hv_logp = -log10(2*normcdf(hv_z));

      hData  = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
      hModel = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;
      hNull = plot(hv_logp,hv_logp, 'k--', 'LineWidth',1); hold on;
      xlim([0 7])
      ylim([0 50])
        drawnow
    end
end

function model_logpvec = calc_model_logpvec(params, hv_z, model_weights)
    % Calculate model_logpvec
    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));

    pBuffer = libpointer('singlePtr', zeros(length(hv_z), 1, 'single'));
    calllib('bgmg', 'bgmg_calc_univariate_pdf', 0, params.pi_vec, params.sig2_zero, params.sig2_beta, length(hv_z), hv_z, pBuffer);  check(); 
    pdf = pBuffer.Value'; clear pBuffer
    pdf = pdf / sum(model_weights);
    model_cdf = cumsum(pdf)  * (hv_z(2) - hv_z(1)) ;
    X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
    model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
    model_logpvec = -log10(2*model_cdf);
end