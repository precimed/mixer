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

    % beta_hat_std_limit and beta_hat_std_step are used in posterior effect size
    % estimation. They express the grid to calculate posterior beta
    % distribution before taking mean or median statistic.
    % The values are expressed in 'normalized' form, e.g. in units of the
    % params.sigma_beta.
    if ~isfield(options, 'beta_hat_std_limit'), options.beta_hat_std_limit = 20; end;
    if ~isfield(options, 'beta_hat_std_step'),  options.beta_hat_std_step  = 0.01; end;
    if ~isfield(options, 'calculate_beta_hat'), options.calculate_beta_hat = true; end;

    if ~isstruct(params), params = mapparams(params); end;
    % params.sigma0     - scalar in (0, +inf]
    % params.sigma_beta - scalar in (0, +inf]
    % params.pivec      - scalar in [0, 1]

    % Validate input params
    if any(params.sigma0 <= 0), warning('sigma0 can not be zero'); cost = nan; return; end;
    if any(params.sigma_beta) <= 0, warning('sigma_beta can not be zero'); cost = nan; return; end;
    if any(params.pivec < 0), warning('pivec must be from 0 to 1'); cost = nan; return; end;
    if sum(params.pivec) > 1, warning('sum(pivec) can not exceed 1'); cost = nan; return; end;
    if isempty(ref_ld), ref_ld_r2 = zero(size(w_ld)); ref_ld_r4 = zero(size(w_ld));
    else ref_ld_r2 = ref_ld.sum_r2; ref_ld_r4 = ref_ld.sum_r4; end;
    
    defvec = isfinite(zvec + Hvec + Nvec + w_ld + ref_ld_r2 + ref_ld_r4) & (Hvec > 0);
    defvec(any(abs(zvec) > options.zmax, 2)) = false;
    zvec = zvec(defvec, :); Hvec = Hvec(defvec); Nvec = Nvec(defvec, :);
    w_ld = w_ld(defvec); ref_ld_r2 = ref_ld_r2(defvec); ref_ld_r4 = ref_ld_r4(defvec);

    sigma0_sqr = params.sigma0 .^ 2;
    sbt_sqr    = repmat(params.sigma_beta .^ 2, size(Hvec));

    % p1 - proportion of causal SNPs, and p0 - proportion of null SNPs (p0) 
    pi1 = repmat(params.pivec, size(Hvec)); pi0 = 1 - pi1;
    
    if options.use_ref_ld
        %pi0 = pi0 .^ ref_ld_r2; pi1 = 1 - pi0; % previous approximation from GMM code

        % Approximation that preserves variance and kurtosis
        r2eff   = ref_ld_r4 ./ ref_ld_r2;
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
    cost = sum(weights .* -log(pdf));
    if ~isfinite(cost), cost = NaN; end;

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
