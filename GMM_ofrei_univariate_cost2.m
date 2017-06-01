function [cost, zprobvec] = GMM_ofrei_univariate_cost2(x, zvec, Hvec, Nvec, ref_ld_hist, ref_ld_bins, w_ld, mapparams, options)
    % x         - vectors of params
    % mapparams - function that maps x into a struct with fields pi1, sigma_beta, sigma0
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld_hist - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs,
    % ref_ld_bins - calculated as a histogram binned by levels of ld r2
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options
    error('This code is not finished.');
    params = mapparams(x);
    if params.sigma0 == 0, error('sigma0 can not be zero'); end;
    if params.sigma_beta == 0, error('sigma_beta can not be zero'); end;

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'x_limit'), options.x_limit = 20; end;             % tabulate probability distributions for [-xlim : xstep : xlim]
    if ~isfield(options, 'x_step'), options.x_step = 0.01; end;
    if ~isfield(options, 'ld_bins'), options.ld_bins = 30; end;             % number of LD bins to use

    defvec = isfinite(zvec + Hvec + Nvec + sum(ref_ld_hist, 2) + w_ld) & (Hvec > 0);
    defvec(abs(zvec) > 12) = false;
    %defvec(10000:end) = false;
    zvec = zvec(defvec); Hvec = Hvec(defvec); Nvec = Nvec(defvec); ref_ld_hist = ref_ld_hist(defvec, :); w_ld = w_ld(defvec);

    pivec = params.pivec;              % proportion of causal SNPs
    sbt = params.sigma_beta;           % 'discoverability', aka variance of the effect sizes
    sigma0 = params.sigma0;            % variance from non-heritable (noise) component
    nsnp = length(zvec);               % number of SNPs

    % preliminary stuff - analyze ref_ld_hist
    ref_ld_hist_index = nan(nsnp, size(tld_hist, 2)); ld_table_hist = nan(size(tld_hist, 2), options.ld_bins); mean_ld_r2_hist = [];
    for tld_hist_bini = 1:size(tld_hist, 2)
        mean_ld_r2 = mean(ref_ld_bins(tld_hist_bini, :));
        ref_ld = ref_ld_hist(:, tld_hist_bini) / mean_ld_r2;
        ld_table = unique(ceil(logspace(0, log10(ceil(max(ref_ld)) + 1), options.ld_bins)));
        for ld_index=1:length(ld_table)
            index = ref_ld <= ld_table(ld_index);
            if (ld_index > 1), index = index & ref_ld > ld_table(ld_index - 1); end;
            ref_ld_hist_index(index, tld_hist_bini) = ld_index;
        end
        ld_table_hist(tld_hist_bini, 1:length(ld_table)) = ld_table;
        mean_ld_r2_hist(tld_hist_bini) = mean_ld_r2;
    end
    if any(any(isnan(ref_ld_hist_index))), error('failed to find LD bins for some SNPs'); end;

    conv_n = @(x, n)op_power(x, @(a,b)conv(a, b, 'same'), n);
    sbt = sbt .* mean(sqrt(Nvec .* Hvec));  % Hack-hack, use the mean here.

    % prior distribution 
    x = (-options.x_limit:options.x_step:options.x_limit); % * sbt;   % <------ scale x by sbt.
    middle_index = (length(x)+1)/2;  % assume length(x) is odd
    if (mod(length(x), 2) ~= 1), error('length(x) must be odd'); end;
    
    beta_conv = nan(size(tld_hist, 2), options.ld_bins, length(x));
    for tld_hist_bini = 1:size(tld_hist, 2)
        mean_ld_r2 = mean_ld_r2_hist(tld_hist_bini);
        ld_table = ld_table_hist(tld_hist_bini, :); ld_table = ld_table(~isnan(ld_table));
        
        beta = normpdf(x, 0, sqrt(mean_ld_r2) * sbt); beta = beta * pivec / sum(beta);
        beta(middle_index) = beta(middle_index) + (1 - pivec);

        % noise component
        noise = normpdf(x, 0, sigma0); noise = noise / sum(noise);

        % calculate probability density for $z = sum_{i=1}^ld beta_i + noise$
        % across exponentially-spaced set of ld values (see ld_table)
        for ld_index=1:length(ld_table)
            % (conv_n(beta, ld_table(ld_index)), noise, 'same');
            beta_conv(tld_hist_bini, ld_index, :) = conv_n(beta, ld_table(ld_index));
        end
    end

    % adjusted up to here.
    % next steps: manually perform convolution for each SNP.
    
    % perform lookup of z probabilities
    zprobvec = nan(size(zvec)); 
    for ld_index=1:length(ld_table)
        index = ref_ld <= ld_table(ld_index);
        if (ld_index > 1), index = index & ref_ld > ld_table(ld_index - 1); end;
        zprobvec(index) = interp1(x, beta_conv(ld_index, :), zvec(index)) / (x(2) - x(1));
    end
    if any(isnan(zprobvec)), error('failed to calculate p(z) for some SNPs'); end;
    
    % Likelihood term, weighted by inverse TLD
    weights = 1 ./ w_ld;
    cost = sum(weights .* -log(zprobvec));
    if ~isfinite(cost), cost = NaN; end;

    % project zprobvec back into the original SNPs indexing (which may include undefined values)
    tmp = zprobvec; zprobvec = nan(size(defvec)); zprobvec(defvec) = tmp;
    fprintf('UVT: pi1 = %.3e, sigma_beta^2 = %.3e (eff. %.3e), sigma_0 = %.3e, cost = %.3e\n', pivec, params.sigma_beta.^2, sbt.^2, sigma0, cost);
end


% TBD
%   - bin convolution calculations by sample size and heterozigosity (today just use mean)
%   - fix how we interpolate; use interp2 instead of interp1, because actually we lookup in 2D plane (x,LD)