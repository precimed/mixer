function [cost, result] = GMM2_univariate_cost(params, zvec, Hvec, Nvec, w_ld, ref_ld, mapparams, options)
    % INPUT:
    % params    - struct with fields pi1, sigma_beta, sigma0
    % mapparams - function that maps params into a vector
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld    - a struct with two fields
    %   .sum_r2 - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs
    %             ref_ld.sum_r2(j) = sum_i r^2_{ij}
    %   .sum_r4 - Similar to ref_ld.sum_r2, but contains a sum of r_{ij}^4
    %             ref_ld.sum_r4(j) = sum_i r^4_{ij}
    %             WARNING: the sum_r2 and sum_r4 must be consistent with each other
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options
    %
    % OUTPUT
    % cost      - log-likelihood
    % result    - struct with fields:
    %   .pdf             - probability density function p(z|params)
    %   .cdf             - cumulated distribution function of p(z|params)
    %                      this is a matrix of size nsnps X len(zgrid),
    %                      where zgrid is defined by options.calculate_z_cdf_limit and
    %                      options.calculate_z_cdf_step)
    %   .fdr_local       - local false discovery rate (fdr)
    %   .fdr_global      - false discovery rate (FDR)
    %   .beta_hat_mean   - posterior effect size estimate (mean)
    %   .beta_hat_se     - standard error of the mean posterior effect size estimate
    %   .beta_hat_median - posterior effect size estimate (median)
    %                      (use options.calculate_beta_hat=false to disable
    %                      beta_hat calculation)

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'zmax'), options.zmax = 12; end;                % throw away z scores larger than this value
    if ~isfield(options, 'verbose'), options.verbose = false; end;       % enable or disable verbose logging
    if ~isfield(options, 'use_ref_ld'), options.use_ref_ld = true; end;  % enable approximation that uses ref_ld
    if ~isfield(options, 'use_convolution'), options.use_convolution = false; end;  % experimental option to calculate pdf via convolution

    % beta_hat_std_limit and beta_hat_std_step are used in posterior effect size
    % estimation. They express the grid to calculate posterior beta
    % distribution before taking mean or median statistic.
    % The values are expressed in 'normalized' form, e.g. in units of the
    % params.sigma_beta.
    if ~isfield(options, 'beta_hat_std_limit'), options.beta_hat_std_limit = 20; end;
    if ~isfield(options, 'beta_hat_std_step'),  options.beta_hat_std_step  = 0.01; end;
    if ~isfield(options, 'calculate_beta_hat'), options.calculate_beta_hat = true; end;

    % z_cdf_limit and z_cdf_step are used in cdf estimation
    if ~isfield(options, 'calculate_z_cdf_limit'), options.calculate_z_cdf_limit = 15; end;
    if ~isfield(options, 'calculate_z_cdf_step'), options.calculate_z_cdf_step = 0.25; end;
    if ~isfield(options, 'calculate_z_cdf'), options.calculate_z_cdf = false; end;

    if ~isstruct(params), params = mapparams(params); end;
    % params.sigma0     - scalar in (0, +inf]
    % params.sigma_beta - scalar in (0, +inf]
    % params.pivec      - scalar in [0, 1]

    % Validate input params
    if any(params.sigma0 <= 0), warning('sigma0 can not be zero'); cost = nan; return; end;
    if any(params.sigma_beta) <= 0, warning('sigma_beta can not be zero'); cost = nan; return; end;
    if any(params.pivec < 0), warning('pivec must be from 0 to 1'); cost = nan; return; end;
    if sum(params.pivec) > 1, warning('sum(pivec) can not exceed 1'); cost = nan; return; end;
    if isempty(ref_ld), error('ref_ld argument is required'); cost = nan; return; end;
    
    % Symmetrization trick - disable, doesn't seems to make a difference
    % zvec = [zvec; -zvec];
    % Hvec = [Hvec; Hvec]; Nvec = [Nvec; Nvec]; w_ld = 2 * [w_ld; w_ld]; % twice w_ld to preserve sum(1./w_ld).
    % ref_ld = struct('sum_r2', [ref_ld.sum_r2; ref_ld.sum_r2], 'sum_r4', [ref_ld.sum_r4; ref_ld.sum_r4]);

    ref_ld_r2 = ref_ld.sum_r2; ref_ld_r4 = ref_ld.sum_r4;
    r2eff   = ref_ld_r4 ./ ref_ld_r2;

    defvec = isfinite(zvec + Hvec + Nvec + w_ld + ref_ld_r2 + ref_ld_r4) & (Hvec > 0);
    defvec(any(abs(zvec) > options.zmax, 2)) = false;
    defvec((r2eff < 0) | (r2eff > 1) | (ref_ld_r4 < 0) | (ref_ld_r2 < 0)) = false;

    zvec = zvec(defvec, :); Hvec = Hvec(defvec); Nvec = Nvec(defvec, :);
    w_ld = w_ld(defvec); ref_ld_r2 = ref_ld_r2(defvec); ref_ld_r4 = ref_ld_r4(defvec); r2eff = r2eff(defvec);

    sigma0_sqr = params.sigma0 .^ 2;
    sbt_sqr    = repmat(params.sigma_beta .^ 2, size(Hvec));

    % p1 - proportion of causal SNPs, and p0 - proportion of null SNPs (p0) 
    pi1 = repmat(params.pivec, size(Hvec)); pi0 = 1 - pi1;
    
    if options.use_ref_ld
        %pi0 = pi0 .^ ref_ld_r2; pi1 = 1 - pi0; % previous approximation from GMM code

        % Approximation that preserves variance and kurtosis
        sbt_sqr = sbt_sqr .* (pi1 .* ref_ld_r2 + pi0 .* r2eff);
        pi1     = (pi1 .* ref_ld_r2) ./ (pi1 .* ref_ld_r2 + pi0 .* r2eff); 
        pi0     = 1 - pi1;
    end

    % Components
    pdf0 = pi0 .* fast_normpdf1(zvec,                          sigma0_sqr);
    pdf1 = pi1 .* fast_normpdf1(zvec, sbt_sqr .* Hvec .* Nvec + sigma0_sqr);
    pdf  = pdf0 + pdf1;

    % Likelihood term, weighted by inverse TLD
    weights = 1 ./ w_ld;
    weights(w_ld < 1) = 1;
    cost = sum(weights .* -log(pdf));
    if ~isfinite(cost), cost = NaN; end;
    if ~isreal(cost), cost = NaN; end;

    if options.use_convolution
        distr   = @(x)(x ./ sum(x));
        conv_n  = @(x, n)op_power(x, @(a,b)conv(a, b, 'same'), n);
        dirac   = @(z)([zeros(1, (length(z)-1)/2), 1, zeros(1,(length(z)-1)/2)]);
        mixture = @(z, p, s)(p*distr(normpdf(z, 0, s)) + (1-p) * dirac(z));

        pdf_by_conv = nan(size(pdf));

        ld_block = round(ref_ld_r2./r2eff);
        ld_sigma = sqrt(Nvec.*Hvec.*r2eff.*params.sigma_beta.^2);
        
        num_ld_block_bins = 20;
        num_ld_sigma_bins = 30;
        
        ld_block_bins = [-Inf unique(floor(logspace(0, log10(max(ld_block)+1), num_ld_block_bins)))];
        ld_sigma_bins = [-Inf quantile(ld_sigma, num_ld_sigma_bins) Inf];
         
        for ld_block_bini = 2:length(ld_block_bins)
        for ld_sigma_bini = 2:length(ld_sigma_bins)
            idx1 = (ld_block > ld_block_bins(ld_block_bini - 1)) & (ld_block <= ld_block_bins(ld_block_bini));
            idx2 = (ld_sigma > ld_sigma_bins(ld_sigma_bini - 1)) & (ld_sigma <= ld_sigma_bins(ld_sigma_bini));
            idx = idx1 & idx2;
            if (sum(idx) == 0), continue; end;

            mean_ld_block = round(mean(ld_block(idx)));
            mean_ld_sigma = mean(ld_sigma(idx));
            
            z_max = max(100, round(max(abs(zvec(idx))) / mean_ld_sigma) + 1);
            if z_max <= 100, step = 0.1; 
            elseif z_max <= 200, step = 0.2;
            elseif z_max <= 500, step = 0.5;
            elseif z_max <= 1000, step = 1.0;
            else warning('too large z_max'); step = 1.0;
            end
            z = -z_max:step:z_max; z = z * mean_ld_sigma; 
            
            table_to_pdf_factor = sum(normpdf(z)); % or, 1./(step*mean_ld_sigma) - this is the same
            pdf_bini = table_to_pdf_factor * conv(conv_n(mixture(z, params.pivec, mean_ld_sigma), mean_ld_block), distr(normpdf(z, 0, params.sigma0)), 'same');
            
            % Hack-hack to avoid concentrating the entire distribution at zero.
            pdf_bini((length(z)+1)/2) = pdf_bini((length(z)-1)/2);

            pdf_by_conv(idx) = interp1(z, pdf_bini, zvec(idx));
            if any(isnan(pdf_by_conv(idx))),
                fprintf('error %i %i\n', ld_block_bini, ld_sigma_bini)
            end
        end
        end
        
        cost_by_conv = sum(weights .* -log(pdf_by_conv));
        if ~isfinite(cost_by_conv), cost_by_conv = NaN; end;
        if ~isreal(cost_by_conv), cost_by_conv = NaN; end;
        if isnan(cost_by_conv)
            fprintf('stop');
        end

        fprintf('cost: %.5e vs %.5e, delta = %.5e\n', cost, cost_by_conv, cost - cost_by_conv);
        pdf = pdf_by_conv;
        cost = cost_by_conv;

        %hold on; plot(-log10(pdf),-log10(pdf_by_conv),'.'); plot([0 7], [0 7]); hold off
    end

    if nargout > 1,
        % probability density function
        result.pdf = pdf;
        
        % local fdr
        result.fdr_local = pdf0 ./ pdf;

        % global FDR
        zvec_neg = min(zvec, -zvec);  % ensure negative values
        cdf0 = pi0 .* fast_normcdf1(zvec_neg,                          sigma0_sqr);
        cdf1 = pi1 .* fast_normcdf1(zvec_neg, sbt_sqr .* Hvec .* Nvec + sigma0_sqr);
        cdf  = cdf0 + cdf1;
        result.fdr_global = cdf0 ./ cdf;
    end 

    % Cumulated distribution function for Z scores
    if (nargout > 1) && options.calculate_z_cdf
        fprintf('Estimate cumulated distribution function for Z scores...\n');
        z_grid =  (-options.calculate_z_cdf_limit:options.calculate_z_cdf_step:options.calculate_z_cdf_limit);
        snps=length(zvec);
        result.cdf = nan(snps, length(z_grid));
        chunksize = floor(snps/length(z_grid));
        for snpi=1:chunksize:snps
            snpj = min(snpi+chunksize-1, snps);
            z_grid_data = struct('pi0', pi0, 'pi1', pi1, 'sbt_sqr', sbt_sqr, 'Hvec', Hvec, 'Nvec', Nvec);
            fn = fieldnames(z_grid_data); for fi = 1:length(fn), d = z_grid_data.(fn{fi}); z_grid_data.(fn{fi}) = repmat(d(snpi:snpj, :), [1, length(z_grid)]); end;
            z_grid_zvec = repmat(z_grid, [snpj-snpi+1, 1]);
            z_grid_pdf0 = z_grid_data.pi0 .* fast_normpdf1(z_grid_zvec, sigma0_sqr);
            z_grid_pdf1 = z_grid_data.pi1 .* fast_normpdf1(z_grid_zvec, sigma0_sqr + z_grid_data.sbt_sqr .* z_grid_data.Hvec .* z_grid_data.Nvec);
            z_grid_pdf  = z_grid_pdf0 + z_grid_pdf1;
            result.cdf(snpi:snpj, :)  = cumsum(z_grid_pdf, 2) ./ repmat(sum(z_grid_pdf, 2), [1, length(z_grid)]);
            fprintf('\tFinish %i SNPs out of %i\n', snpj, snps);
        end
    end
    
    % Posterior effect size estimates
    if (nargout > 1) && options.calculate_beta_hat
        fprintf('Estimate posterior effect sizes...\n');
        snps=length(zvec);
        result.beta_hat_mean = nan(size(zvec));
        result.beta_hat_se = nan(size(zvec));
        result.beta_hat_median = nan(size(zvec));
        for snpi=1:snps
            if (mod(snpi, 100000) == 0),  fprintf('\tFinish %i SNPs out of %i\n', snpi, snps); end;
            beta_grid = (-options.beta_hat_std_limit:options.beta_hat_std_step:options.beta_hat_std_limit) * sqrt(sbt_sqr(snpi));
            beta_prior = fast_normpdf1(beta_grid, sbt_sqr(snpi));
            beta_prior = beta_prior .* pi1(snpi) / sum(beta_prior);
            middle_index = (length(beta_grid)+1)/2;  % assume length(x) is odd
            if (mod(length(beta_grid), 2) ~= 1), error('length(x) must be odd'); end;
            beta_prior(middle_index) = beta_prior(middle_index) + pi0(snpi);

            like_func = normpdf(ones(size(beta_grid)) .* zvec(snpi), beta_grid .* sqrt(Nvec(snpi) * Hvec(snpi)), params.sigma0);
            beta_posterior = beta_prior .* like_func;
            beta_posterior = beta_posterior ./ sum(beta_posterior);  % normalize

            result.beta_hat_mean(snpi)   = sum(beta_posterior .* beta_grid);  % calc mean
            result.beta_hat_se(snpi)     = sqrt(sum(beta_posterior .* (beta_grid - result.beta_hat_mean(snpi)).^2)); % calc standard error
            result.beta_hat_median(snpi) = beta_grid(find(cumsum(beta_posterior)>0.5, 1, 'first'));  % calc median (works as a filter for noisy observations)
        end
    end

    % project all results back into the original SNPs indexing (which may include undefined values)
    if exist('result', 'var'), result = restore_original_indexing(result, defvec); end

    if options.verbose
        fprintf('Univariate: pi=%.3e, sigma_beta^2=%.3e, (eff. %.3f), sigma0^2=%.3e, cost=%.3e\n', ...
                 params.pivec, params.sigma_beta.^2, mean(sbt_sqr .* Hvec .* Nvec), params.sigma0.^2, cost);
    end
end

function pdf = fast_normpdf1(x, sigma_sqr)
    % Calculation of log-likelihood and pdf, specific to univariate normal
    % distribution with zero mean. In fact this works as fast as 
    % matlab's normpdf; I kept this to be consistent with bivariate code 
    % where fast_normpdf2 is indeed much faster.
    log_exp = -0.5 * x .* x ./ sigma_sqr;
    log_dt  = -0.5 * log(sigma_sqr);
    log_pi  = -0.5 * log(2*pi);
    pdf = exp(log_pi + log_dt + log_exp);
    % loglike = sum(log_exp+log_dt) + length(dt) * log_pi;
end

function cdf = fast_normcdf1(x, sigma_sqr)
    if any(x > 0), error('this should never happen; call this with x=-abs(x)'); end;
    cdf = normcdf(x, 0, sqrt(sigma_sqr));
end

function values = restore_original_indexing(values, defvec)
    if isstruct(values),
        % Apply to each field of the struct
        fn = fieldnames(values);
        for i = 1:length(fn),
            values.(fn{i}) = restore_original_indexing(values.(fn{i}), defvec);
        end;
    else
        tmp = values;
        values = nan(length(defvec), size(values, 2));
        values(defvec, :) = tmp;
    end
end
