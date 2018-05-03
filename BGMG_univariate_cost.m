function [cost, result] = BGMG_univariate_cost(params, zvec, Hvec, Nvec, w_ld, ref_ld, options)
    % INPUT:
    % params    - struct with fields pi_vec, sig2_zero, sig2_beta
    %             pi_vec --- mixture weights, one per mixture component (excluding null component)
    %             sig2_zero --- inflation from cryptic relatedness / sample stratification (sigma0^2)
    %             sig2_beta --- discoverability, one per mixture component (excluding null component)
    %
    %             Example1: standard two-component mixture model, null + causal
    %               struct('pi_vec', 0.1, 'sig2_zero', 1.05, 'sig2_beta', 1e-3)
    %
    %             Example2: Alzheimer disease with polygenic signal + APOE
    %               struct('pi_vec', [0.1 0.2], 'sig2_zero', 1.05, 'sig2_beta', [1e-3 1e-4])
    %
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld    - a struct with two fields
    %   .sum_r2 - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs
    %             ref_ld.sum_r2(j) = sum_i r^2_{ij}
    %   .chi_r4 - A fraction sum_r4 ./ sum_r2, where
    %             sum_r4(j) = sum_i r_{ij}^4
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
    %   .delta_hat_mean   - posterior effect size estimate (mean)
    %   .delta_hat_se     - standard error of the mean posterior effect size estimate
    %   .delta_hat_median - posterior effect size estimate (median)
    %                       (use options.calculate_delta_hat=false to disable
    %                       delta_hat calculation)

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;    % enable or disable verbose logging
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    if ~isfield(options, 'zmax'), options.zmax = +Inf; end;           % ignore z scores above zmax

    if ~isfield(options, 'use_poisson'), options.use_poisson = false; end;  % enable accurate pdf(z) estimation
    if ~isfield(options, 'poisson_sigma_grid_nodes'), options.poisson_sigma_grid_nodes = 25;  end;  % grid size of poisson projection
    if ~isfield(options, 'poisson_sigma_grid_scale'), options.poisson_sigma_grid_scale = 1;  end;
    if ~isfield(options, 'poisson_sigma_grid_chunks'), options.poisson_sigma_grid_chunks = 20; end;

    % z_cdf_limit and z_cdf_step are used in cdf estimation
    if ~isfield(options, 'calculate_z_cdf_limit'), options.calculate_z_cdf_limit = 15; end;
    if ~isfield(options, 'calculate_z_cdf_step'), options.calculate_z_cdf_step = 0.25; end;
    if ~isfield(options, 'calculate_z_cdf'), options.calculate_z_cdf = false; end;
    if ~isfield(options, 'calculate_z_cdf_weights'), options.calculate_z_cdf_weights = ones(size(zvec)) / length(zvec); end;
    if ~isfield(options, 'calculate_z_cdf_nscale'), options.calculate_z_cdf_nscale = 1; end;
    if ~isfield(options, 'calculate_z_cdf_func'), options.calculate_z_cdf_func = @normcdf; end;

    % delta_hat_std_limit and delta_hat_std_step are used in posterior effect size
    % estimation. They express the grid to calculate posterior delta
    % distribution before taking mean or median statistic.
    if ~isfield(options, 'delta_hat_std_limit'), options.delta_hat_std_limit = 20; end;
    if ~isfield(options, 'delta_hat_std_step'),  options.delta_hat_std_step  = 0.5; end;
    if ~isfield(options, 'calculate_delta_hat'), options.calculate_delta_hat = false; end;

    if ~isfield(options, 'calculate_fdr'), options.calculate_fdr = false; end;

    non_zero_mixi = (params.sig2_beta ~=0) & (params.pi_vec ~= 0);
    params.sig2_beta = params.sig2_beta(non_zero_mixi);
    params.pi_vec = params.pi_vec(non_zero_mixi);

    % Validate input params
    result=[];
    if ~BGMG_util.validate_params(params); cost = nan; return; end;
    if isempty(ref_ld), error('ref_ld argument is required'); cost = nan; return; end;
    if options.calculate_z_cdf && ~options.use_poisson, error('calculate_z_cdf is set to true; did you forget "options.use_poisson=true?"'); end;
    
    if ~options.use_poisson && (size(ref_ld.sum_r2, 2) ~= 1)
        ref_ld.sum_r2 = sum(ref_ld.sum_r2, 2);
        ref_ld.sum_r2_biased = sum(ref_ld.sum_r2_biased, 2);
        ref_ld.sum_r4_biased = sum(ref_ld.sum_r4_biased, 2);
    end

    ref_ld_sum_r2 = ref_ld.sum_r2;
    ref_ld_chi_r4 = ref_ld.sum_r4_biased ./ ref_ld.sum_r2_biased;
    ref_ld_chi_r4(ref_ld.sum_r4_biased == 0) = 0;

    % Symmetrization trick - disable, doesn't seems to make a difference
    % zvec = [zvec; -zvec];
    % Hvec = [Hvec; Hvec]; Nvec = [Nvec; Nvec]; w_ld = 2 * [w_ld; w_ld]; % twice w_ld to preserve sum(1./w_ld).
    % ref_ld = struct('sum_r2', [ref_ld.sum_r2; ref_ld.sum_r2], 'sum_r4', [ref_ld.sum_r4; ref_ld.sum_r4]);

    defvec = isfinite(zvec + Hvec + Nvec + w_ld + sum(ref_ld_sum_r2, 2) + sum(ref_ld_chi_r4, 2)) & (Hvec > 0);
    defvec(any(abs(zvec) > options.zmax, 2)) = false;
    defvec(any(ref_ld_chi_r4 < 0, 2) | any(ref_ld_chi_r4 > 1, 2) | any(ref_ld_sum_r2 < 0, 2)) = false;

    zvec = zvec(defvec, :); Hvec = Hvec(defvec); Nvec = Nvec(defvec, :);
    w_ld = w_ld(defvec); ref_ld_sum_r2 = ref_ld_sum_r2(defvec, :); ref_ld_chi_r4 = ref_ld_chi_r4(defvec, :);
    options.calculate_z_cdf_weights = options.calculate_z_cdf_weights(defvec, :)./repmat(nansum(options.calculate_z_cdf_weights(defvec, :)), [sum(defvec) 1]);
    
    % Weights as inverse TLD
    weights = 1 ./ w_ld;
    weights(w_ld < 1) = 1;

    num_traits = length(params.sig2_zero);  % number of traits in the anslysis
    num_mix    = length(params.pi_vec);  % number of components in the mixture
    num_snps   = length(Hvec);           % number of SNPs in the analysis
    num_ldr2bins = size(ref_ld_chi_r4, 2);  % number of bins in LDr2 histogram
    if num_traits ~= 1, error('Params define more than one trait, unable to calculate univariate cost.'); end;

    if ~options.use_poisson
        % Approximation that preserves variance and kurtosis
        pi_vec     = repmat(params.pi_vec, [num_snps 1]);
        eta_factor = (pi_vec .* repmat(ref_ld_sum_r2, [1, num_mix]) + (1-pi_vec) .* repmat(ref_ld_chi_r4, [1, num_mix]));
        pi_vec     = (pi_vec .* repmat(ref_ld_sum_r2, [1, num_mix])) ./ eta_factor;

        % Normalize pi vector, and compensate total variance by adjusting eta_factor
        [pi0, pi_vec_norm] = BGMG_util.normalize_pi_vec(pi_vec);
        eta_factor  = eta_factor .* (pi_vec ./ pi_vec_norm);
        eta_factor(pi_vec == 0) = 0;

        % Multiply sig2_beta by eta_factor
        sig2_beta  = repmat(params.sig2_beta, [num_snps 1]) .* eta_factor;

        % Calculate p.d.f. from the mixture model
        pdf0 = pi0 .* fast_normpdf1(zvec, params.sig2_zero);
        pdf = pi_vec_norm .* fast_normpdf1(repmat(zvec, [1 num_mix]), sig2_beta .* repmat(Hvec .* Nvec, [1 num_mix]) + params.sig2_zero);
        pdf  = pdf0 + sum(pdf, 2);
    else
        % Poisson approximation + "fixed sigma grid"
        poisson_lambda = cell(num_mix, 1); poisson_sigma = cell(num_mix, 1);
        snp_contribution = 0;
        for mixi = 1:num_mix
            pi1u = params.pi_vec(mixi);
            poisson_lambda{mixi} = pi1u * ref_ld_sum_r2 ./ ((1-pi1u) * ref_ld_chi_r4);
            poisson_lambda{mixi}(~isfinite(poisson_lambda{mixi})) = 0;
            poisson_sigma{mixi}  = repmat(Hvec .* Nvec, [1 num_ldr2bins]) .* ...
                             params.sig2_beta(mixi) .* (1-pi1u) .* ref_ld_chi_r4;
            snp_contribution = snp_contribution + sum(poisson_lambda{mixi} .* poisson_sigma{mixi}, 2);
        end
        % Compensate for low kmax (commented out because kmax is not known yet... this is some kind of "circular dependency")
        % poisson_sigma =poisson_sigma ./ poisscdf(repmat(options.poisson_kmax - 1, [size(poisson_lambda, 1), 1]), poisson_lambda);

        snp_chunk_size  = floor(num_snps/options.poisson_sigma_grid_chunks);
        [~, sorted_snps] = sort(snp_contribution);

        pdf = zeros(num_snps, 1);
        % Cumulated distribution function for Z scores
        if (nargout > 1) && options.calculate_z_cdf
            fprintf('Estimate cumulated distribution function for Z scores... \n');
            result_cdf_z_grid = (-options.calculate_z_cdf_limit:options.calculate_z_cdf_step:options.calculate_z_cdf_limit);
            [result_cdf_z_grid, result_cdf_nscale] = meshgrid(result_cdf_z_grid,options.calculate_z_cdf_nscale);
            result_cdf_z_grid = result_cdf_z_grid(:)';
            result_cdf_nscale = result_cdf_nscale(:)';
            num_cdf_weights = size(options.calculate_z_cdf_weights, 2);
            result_cdf = zeros(num_cdf_weights, length(result_cdf_z_grid));
            % power_numerator = zeros(size(options.calculate_z_cdf_nscale)); power_numerator=power_numerator(:);
            % power_denominator = zeros(size(options.calculate_z_cdf_nscale)); power_denominator=power_denominator(:);
        end

        % Make sure last chunk is as large as all previous
        snpi_vals = 1:snp_chunk_size:(num_snps-snp_chunk_size);
        snpj_vals = snpi_vals + (snp_chunk_size - 1); snpj_vals(end) = num_snps;
        
        has_speedup_info = isfield(options, 'speedup_info');
        if ~has_speedup_info
            options.speedup_info.kmax = nan(length(snpi_vals), num_ldr2bins);
            options.speedup_info.poisson_sigma_grid_limit = nan(length(snpi_vals), 1);
        end

        for snp_vals_index = 1:length(snpi_vals)
            snpi = snpi_vals(snp_vals_index);
            snpj = snpj_vals(snp_vals_index);

            num_snps_in_chunk = snpj - snpi + 1;
            snps_in_chunk = sorted_snps(snpi:snpj);

            if ~has_speedup_info
                % detect parameters of the grid
                quantile_value = 0.999;
                poisson_sigma_grid_limit = nan(num_mix, num_ldr2bins);
                for mixi = 1:num_mix
                    poisson_sigma_grid_limit(mixi, :) = quantile(max(1, poissinv(quantile_value, poisson_lambda{mixi}(snps_in_chunk, :))) .* poisson_sigma{mixi}(snps_in_chunk, :), quantile_value);
                end
                poisson_sigma_grid_limit = options.poisson_sigma_grid_scale * max(2.0, 1.05 * max(poisson_sigma_grid_limit(:)));
                options.speedup_info.poisson_sigma_grid_limit(snp_vals_index) = poisson_sigma_grid_limit;
            end
            poisson_sigma_grid_limit = options.speedup_info.poisson_sigma_grid_limit(snp_vals_index);
            poisson_sigma_grid_nodes = (options.poisson_sigma_grid_scale.^2) * options.poisson_sigma_grid_nodes;
            poisson_sigma_grid_delta = poisson_sigma_grid_limit / (poisson_sigma_grid_nodes - 1);
            poisson_sigma_grid       = 0:poisson_sigma_grid_delta:poisson_sigma_grid_limit;

            if (nargout > 1) && options.calculate_z_cdf
                cdfmat = zeros(length(result_cdf_z_grid), poisson_sigma_grid_nodes);
                for i=1:poisson_sigma_grid_nodes
                    cdfmat(:, i) = options.calculate_z_cdf_func(result_cdf_z_grid, ...
                                           0, sqrt(poisson_sigma_grid(i) * result_cdf_nscale + params.sig2_zero));
                end
            end

            sigma_grid_pi   = [];
            for mixi = 1:num_mix
            for ldr2_bini = 1:num_ldr2bins
                if ~has_speedup_info
                    % auto-detect reasonable kmax for poisson approximation
                    poisson_quantile_value = 0.999;
                    kmax      = poissinv(poisson_quantile_value, quantile(poisson_lambda{mixi}(:, ldr2_bini), poisson_quantile_value));
                    kmax      = options.poisson_sigma_grid_scale * max(min(kmax + 1, 40), 4);  assert(isfinite(kmax));
                    options.speedup_info.kmax(snp_vals_index, ldr2_bini) = kmax;
                end
                kmax      = options.speedup_info.kmax(snp_vals_index, ldr2_bini);
                k         = 0:kmax;
                p         = repmat(poisson_lambda{mixi}(snps_in_chunk, ldr2_bini), [1 kmax+1]) .^ repmat(k, [num_snps_in_chunk, 1]) .* ...
                            repmat(exp(-poisson_lambda{mixi}(snps_in_chunk, ldr2_bini)), [1 kmax+1]) ./ repmat(factorial(k), [num_snps_in_chunk, 1]);
                grid_idx  = repmat(k, [num_snps_in_chunk 1])  .* repmat(poisson_sigma{mixi}(snps_in_chunk, ldr2_bini), [1 1+kmax]) ./ poisson_sigma_grid_delta;
                grid_idx  = min(grid_idx, poisson_sigma_grid_nodes - 2);
                grid_idx_floor = floor(grid_idx);
                grid_idx_frac  = grid_idx - grid_idx_floor;

                grid_idx_floor = grid_idx_floor'; p = p'; grid_idx_frac = grid_idx_frac';
                grid_coef_tmp  = accummatrix(1+grid_idx_floor, p .* (1 - grid_idx_frac), [poisson_sigma_grid_nodes, size(p, 2)]) + ...
                                 accummatrix(2+grid_idx_floor, p .* (    grid_idx_frac), [poisson_sigma_grid_nodes, size(p, 2)]);
                grid_coef_tmp  = grid_coef_tmp';

                sigma_grid_pi = convmatrix(sigma_grid_pi, grid_coef_tmp);
            end
            end
            for k=1:poisson_sigma_grid_nodes
                pdf(snps_in_chunk) = pdf(snps_in_chunk) + sigma_grid_pi(:, k) .* normpdf(zvec(snps_in_chunk), 0, sqrt(poisson_sigma_grid(k) + params.sig2_zero));
            end
            
            if (nargout > 1) && options.calculate_z_cdf
                result_cdf = result_cdf + options.calculate_z_cdf_weights(snps_in_chunk, :)' * (sigma_grid_pi * cdfmat');
                fprintf('Finish %i SNPs out of %i - Poisson grid 0:%.2f:%.2f (%i nodes)\n', snpj, num_snps, poisson_sigma_grid_delta, poisson_sigma_grid_limit, poisson_sigma_grid_nodes);
            end
        end
    end

    % Likelihood term, weighted by inverse TLD
    cost = sum(weights .* -log(pdf));
    if ~isfinite(cost), cost = NaN; end;
    if ~isreal(cost), cost = NaN; end;

    % TBD: all options below this line are broken - must be fixed.
    if nargout > 1 && options.calculate_fdr
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
    if (nargout > 1) && options.calculate_delta_hat
        fprintf('Estimate posterior effect sizes...\n');
        snps=length(zvec);
        result.delta_hat_mean = nan(size(zvec));
        result.delta_hat_se = nan(size(zvec));
        result.delta_hat_median = nan(size(zvec));
        for snpi=1:snps
            if (mod(snpi, 100000) == 0),  fprintf('\tFinish %i SNPs out of %i\n', snpi, snps); end;
            delta_grid = (-options.delta_hat_std_limit:options.delta_hat_std_step:options.delta_hat_std_limit);
            delta_prior = fast_normpdf1(delta_grid, sbt_sqr(snpi) .* Hvec(snpi) .* Nvec(snpi) );
            delta_prior = delta_prior .* pi1(snpi) / sum(delta_prior);
            middle_index = (length(delta_grid)+1)/2;  % assume length(x) is odd
            if (mod(length(delta_grid), 2) ~= 1), error('length(x) must be odd'); end;
            delta_prior(middle_index) = delta_prior(middle_index) + pi0(snpi);

            like_func = normpdf(ones(size(delta_grid)) .* zvec(snpi), delta_grid, params.sigma0);
            delta_posterior = delta_prior .* like_func;
            delta_posterior = delta_posterior ./ sum(delta_posterior);  % normalize

            result.delta_hat_mean(snpi)   = sum(delta_posterior .* delta_grid);  % calc mean
            result.delta_hat_se(snpi)     = sqrt(sum(delta_posterior .* (delta_grid - result.delta_hat_mean(snpi)).^2)); % calc standard error
            result.delta_hat_median(snpi) = delta_grid(find(cumsum(delta_posterior)>0.5, 1, 'first'));  % calc median (works as a filter for noisy observations)
        end
    end

    % project all results back into the original SNPs indexing (which may include undefined values)
    if exist('result', 'var'), result = restore_original_indexing(result, defvec); end

    % save non-vectorized parameters
    if exist('result', 'var') && isfield(options, 'speedup_info'), result.speedup_info = options.speedup_info; end
    if exist('result_cdf', 'var')
        result.cdf = result_cdf;
        result.cdf_z_grid = result_cdf_z_grid;
        result.cdf_nscale = result_cdf_nscale;
        if exist('power_denominator', 'var'), result.power_denominator = power_denominator; end
        if exist('power_numerator', 'var'), result.power_numerator = power_numerator; end;
    end

    if num_mix == 1
        fprintf('Univariate: pi_vec=%.3e, sig2_beta^2=%.3e, sig2_zero^2=%.3f, h2=%.3f, cost=%.6e, nsnp=%i\n', ...
                 params.pi_vec, params.sig2_beta, params.sig2_zero, ...
                 (params.sig2_beta*params.pi_vec')*options.total_het, cost, sum(defvec));
    else
        fprintf('Univariate: pi_vec=[%s], sig2_beta^2=[%s], h2=[%s], sig2_zero^2=%.3f, cost=%.6e, nsnp=%i\n', ...
            sprintf('%.3e ', params.pi_vec), ...
            sprintf('%.3e ', params.sig2_beta), ...
            sprintf('%.3f ', params.sig2_beta.*params.pi_vec*options.total_het), ...
            params.sig2_zero, cost, sum(defvec));
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

function A = accummatrix(subs,val,sz)
    % as accumarray, but val can be a matrix
    % example:
    % val = [ 0.1 1; 0.2 2; 0.4 4];
    % subs= [ 1 1; 1 3; 3 3];
    % A   = accummatrix(subs, val, [3 2]);   % [0.3 1.0;  0 0; 0.4 6.0]
    subs = subs + repmat(sz(1) * (0:(size(val, 2) - 1)), [size(val, 1), 1]);
    A = accumarray(subs(:), val(:), [prod(sz) 1]);
    A = reshape(A, sz);
end

function A = convmatrix(B, C)
    % apply convolution to each row of B and C.
    if isempty(B), A=C; return; end;
    assert(isequal(size(B), size(C)));
    num_cols = size(B, 2);
    A=zeros(size(B));
    for i=1:num_cols
        for j=1:i
            A(:,i) = A(:, i) + B(:, j) .* C(:, i-j+1);
        end
    end
end

function make_power_plot_helper()
    % misc code to make power plot (not used; keep just for the history)
    options.power_zthresh = 4.75;
    XXX_pdf = (sigma_grid_pi * cdfmat');
    expected_shape = [length(options.calculate_z_cdf_nscale), size(XXX_pdf, 2) / length(options.calculate_z_cdf_nscale)];
    for snpi_in_chunk = 1:num_snps_in_chunk
        n = zeros(size(result_cdf_z_grid)); d = n;
        for i=1:size(sigma_grid_pi, 2)
            pdf_z_at_delta = 0.5 * normpdf(result_cdf_z_grid, sqrt(poisson_sigma_grid(i) * result_cdf_nscale), sqrt(params.sig2_zero)) + ...
                             0.5 * normpdf(result_cdf_z_grid, -sqrt(poisson_sigma_grid(i) * result_cdf_nscale), sqrt(params.sig2_zero));
            n = n + pdf_z_at_delta .* poisson_sigma_grid(i) .* result_cdf_nscale .* sigma_grid_pi(snpi_in_chunk, i);
            d = d + pdf_z_at_delta .*                                               sigma_grid_pi(snpi_in_chunk, i);
        end
        delta2_expected = n./d; delta2_expected(~isfinite(delta2_expected)) = 0;
        delta2_expected = reshape(delta2_expected, expected_shape);

        XXX_cdf_z_grid = reshape(result_cdf_z_grid, expected_shape);
        XXX_pdf_snp = reshape(XXX_pdf(snpi_in_chunk, :), expected_shape);

        power_denominator = power_denominator + options.calculate_z_cdf_weights(snps_in_chunk(snpi_in_chunk)) * sum(delta2_expected .* XXX_pdf_snp, 2);
        XXX_pdf_snp(abs(XXX_cdf_z_grid) < options.power_zthresh ) = 0;
        power_numerator = power_numerator + options.calculate_z_cdf_weights(snps_in_chunk(snpi_in_chunk)) * sum(delta2_expected .* XXX_pdf_snp, 2);
        plot(power_numerator ./ power_denominator)
        fprintf('.');
    end
end
