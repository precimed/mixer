function [cost, result] = BGMG_bivariate_cost(params, zmat, Hvec, Nmat, w_ld, ref_ld, options)
    % params    - struct with fields pi_vec, sig2_zero, sig2_beta, rho_zero, rho_beta
    %             pi_vec    --- row-vector, mixture weights, one per mixture component (excluding null component)
    %             sig2_zero --- col-vector, inflation from cryptic relatedness / sample stratification (sigma0^2), one per trait
    %             sig2_beta --- matrix,     discoverability, one per trait and mixture component (excluding null component)
    %             rho_zero  --- scalar,     spurious correlation of effect sizes due to sample overlap
    %             rho_beta  --- row-vector, correlation of effect sizes in each component
    %
    %             Example1: three-component mixture model, null + trait1 + trait2 (independent traits)
    %               struct('pi_vec', [0.1 0.2], 'sig2_zero', [1.05; 1.1], 'rho_zero', 0.2, ...
    %                                           'sig2_beta', [1e-3 0; 0 1e-4], 'rho_beta', 0.0)
    %
    %             Example1: two-component mixture model, null + pleiotropic component
    %               struct('pi_vec', [0.1 0.2], 'sig2_zero', [1.05; 1.1], 'rho_zero', 0.2, ...
    %                                           'sig2_beta', [1e-3 0; 0 1e-4], 'rho_beta', 0.0)
    %
    %             Example3:  saturated bivariate mixture model, null + trait1 + trait2 + pleio
    %               struct('pi_vec', [0.1 0.2 0.3], 'sig2_zero', [1.05; 1.1], 'rho_zero', 0.2, ...
    %                                               'sig2_beta', [1e-3 0 1e-3; 0 1e-4 1e-4], 'rho_beta', [0 0 0.5])
    %
    % mapparams - function that maps params into a vector (required to run fminsearch)
    % zmat      - matrix SNP x 2, z scores
    % Hvec      - vector SNP x 1, heterozigosity per SNP
    % Nmat      - number of subjects genotyped per SNP
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % ref_ld    - a struct with two fields
    %   .sum_r2 - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs
    %             ref_ld.sum_r2(j) = sum_i r^2_{ij}
    %   .chi_r4 - A fraction sum_r4 ./ sum_r2, where
    %             sum_r4(j) = sum_i r_{ij}^4
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'use_convolution'), options.use_convolution = false; end;  % experimental VERY SLOW option to calculate pdf via convolution
    if ~isfield(options, 'use_legacy_impl'), options.use_legacy_impl = false; end;  % legacy approximation implementation
    if ~isfield(options, 'use_poisson'), options.use_poisson = false; end;  % enable accurate pdf(z) estimation
    if ~isfield(options, 'poisson_sigma_grid_nodes') || ~isfinite(options.poisson_sigma_grid_nodes), options.poisson_sigma_grid_nodes = 5;  end;  % grid size of poisson projection
    if ~isfield(options, 'poisson_sigma_grid_scale'), options.poisson_sigma_grid_scale = 1;  end;
    if ~isfield(options, 'zmax'), options.zmax = +Inf; end;                 % ignore z scores larger than zmax

    % delta_hat_std_limit and delta_hat_std_step are used in posterior effect size
    % estimation. They express the grid to calculate posterior delta
    % distribution before taking mean or median statistic.
    if ~isfield(options, 'delta_hat_std_limit'), options.delta_hat_std_limit = 20; end;
    if ~isfield(options, 'delta_hat_std_step'),  options.delta_hat_std_step  = 0.5; end;
    if ~isfield(options, 'calculate_delta_hat'), options.calculate_delta_hat = false; end;
    if ~isfield(options, 'calculate_global_fdr'), options.calculate_global_fdr = false; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate

    % Validate input params
    if ~BGMG_util.validate_params(params); cost = nan; return; end;
    if isempty(ref_ld), error('ref_ld argument is required'); cost = nan; return; end;

    result = [];
    if ~options.use_poisson && (size(ref_ld.sum_r2, 2) ~= 1)
        ref_ld.sum_r2 = sum(ref_ld.sum_r2, 2);
        ref_ld.sum_r2_biased = sum(ref_ld.sum_r2_biased, 2);
        ref_ld.sum_r4_biased = sum(ref_ld.sum_r4_biased, 2);
    end

    ref_ld_sum_r2 = ref_ld.sum_r2;
    ref_ld_chi_r4 = ref_ld.sum_r4_biased ./ ref_ld.sum_r2_biased;
    ref_ld_chi_r4(ref_ld.sum_r4_biased == 0) = 0;

    defvec = isfinite(sum(zmat, 2) + Hvec + sum(Nmat, 2) + w_ld + sum(ref_ld_sum_r2, 2) + sum(ref_ld_chi_r4, 2)) & (Hvec > 0);
    defvec(any(abs(zmat) > options.zmax, 2)) = false;
    defvec(any(ref_ld_chi_r4 < 0, 2) | any(ref_ld_chi_r4 > 1, 2) | any(ref_ld_sum_r2 < 0, 2)) = false;

    zmat = zmat(defvec, :); Hvec = Hvec(defvec); Nmat = Nmat(defvec, :);
    w_ld = w_ld(defvec); ref_ld_sum_r2 = ref_ld_sum_r2(defvec, :); ref_ld_chi_r4 = ref_ld_chi_r4(defvec, :);

    % Weights as inverse TLD
    weights = 1 ./ w_ld;
    weights(w_ld < 1) = 1;

    num_traits = length(params.sig2_zero);  % number of traits in the anslysis
    num_mix    = length(params.pi_vec);  % number of components in the mixture
    num_snps   = length(Hvec);           % number of SNPs in the analysis
    num_ldr2bins = size(ref_ld_chi_r4, 2);  % number of bins in LDr2 histogram
    if num_traits ~= 2, error('Params must define two traits, unable to calculate bivariate cost.'); end;
    if num_mix > 3, error('Params must define more than mixtures, unable to calculate bivariate cost.'); end;
    if options.use_poisson && (num_mix ~= 3), error('Poisson is implemented specificly for 3 components, unable to calculate bivariate cost.'); end;

    sigma0_sqr = [params.sig2_zero(1), sqrt(prod(params.sig2_zero)) * params.rho_zero; ...
                                       sqrt(prod(params.sig2_zero)) * params.rho_zero, params.sig2_zero(2)];

    if ~options.use_legacy_impl && ~options.use_poisson
        % Approximation that preserves variance and kurtosis
        pi_vec     = repmat(params.pi_vec, [num_snps 1]);
        eta_factor = (pi_vec .* repmat(ref_ld_sum_r2, [1, num_mix]) + (1-pi_vec) .* repmat(ref_ld_chi_r4, [1, num_mix]));
        pi_vec     = (pi_vec .* repmat(ref_ld_sum_r2, [1, num_mix])) ./ eta_factor;

        % Safe-guard against weird cases
        pi_vec(~isfinite(pi_vec)) = 0;
        eta_factor(pi_vec == 0) = 0;
        
        % Multiply sig2_beta by eta_factor
        a11=zeros(num_snps, num_mix); a22=a11; a12=a11;
        for mixi = 1:num_mix
            a11(:, mixi) = eta_factor(:, mixi) .* Hvec .* Nmat(:, 1) * params.sig2_beta(1, mixi);
            a22(:, mixi) = eta_factor(:, mixi) .* Hvec .* Nmat(:, 2) * params.sig2_beta(2, mixi);
            a12(:, mixi) = sqrt(a11(:, mixi) .* a22(:, mixi)) * params.rho_beta(1, mixi);
        end
        
        if num_mix == 3, pdf = pdf_mix3(pi_vec, zmat, a11, a12, a22, sigma0_sqr);
        elseif num_mix == 2, pdf = pdf_mix2(pi_vec, zmat, a11, a12, a22, sigma0_sqr);
        elseif num_mix == 1, pdf = pdf_mix1(pi_vec, zmat, a11, a12, a22, sigma0_sqr);
        end
        % Calculate p.d.f. from 8-component mixture model
        % (analytical convolution of three 2-component gaussian mixtures)
        % Likelihood term, weighted by inverse TLD
        cost = sum(weights .* -log(pdf));

        if ~isfinite(cost), cost = NaN; end;
        if ~isreal(cost), cost = NaN; end;
    end
    
    if options.use_poisson
        Nmat_median = median(Nmat);

        % Poisson approximation + "fixed sigma grid"
        poisson_lambda = cell(num_mix, 1); poisson_sigma = cell(num_mix, 1);
        snp_contribution = 0;
        for mixi = 1:num_mix
            pi1u = params.pi_vec(mixi);
            poisson_lambda{mixi} = pi1u * ref_ld_sum_r2 ./ ((1-pi1u) * ref_ld_chi_r4);
            poisson_lambda{mixi}(~isfinite(poisson_lambda{mixi})) = 0;
            poisson_sigma{mixi}  = repmat(Hvec, [1 num_ldr2bins]) .* (1-pi1u) .* ref_ld_chi_r4;
            snp_contribution = snp_contribution + sum(poisson_lambda{mixi} .* poisson_sigma{mixi} * max(Nmat_median' .* params.sig2_beta(:, mixi)), 2);
        end

        snp_chunk_size  = floor(num_snps/20);  % TBD - what's reasonable number of SNPs chunks?
        [~, sorted_snps] = sort(snp_contribution);

        pdf = zeros(num_snps, 1);

        poisson_sigma_grid_nodes = (options.poisson_sigma_grid_scale.^2) * options.poisson_sigma_grid_nodes;

        % Make sure last chunk is as large as all previous
        snpi_vals = 1:snp_chunk_size:(num_snps-snp_chunk_size);
        snpj_vals = snpi_vals + (snp_chunk_size - 1); snpj_vals(end) = num_snps;

        has_speedup_info = isfield(options, 'speedup_info') && ~isempty(options.speedup_info);
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
            poisson_sigma_grid_delta = poisson_sigma_grid_limit / (poisson_sigma_grid_nodes - 1);
            poisson_sigma_grid       = 0:poisson_sigma_grid_delta:poisson_sigma_grid_limit;

            for mixi = 1:num_mix
            sigma_grid_pi{mixi}   = [];
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

                sigma_grid_pi{mixi} = convmatrix(sigma_grid_pi{mixi}, grid_coef_tmp);
            end
            end

            assert(params.sig2_beta(2, 1) == 0)
            assert(params.sig2_beta(1, 2) == 0)
            for k1=1:poisson_sigma_grid_nodes
            for k2=1:poisson_sigma_grid_nodes
            for k3=1:poisson_sigma_grid_nodes
                a11 = poisson_sigma_grid(k1) * params.sig2_beta(1, 1) * Nmat_median(1) + ...
                      poisson_sigma_grid(k3) * params.sig2_beta(1, 3) * Nmat_median(1) + ...
                      params.sig2_zero(1);
                a22 = poisson_sigma_grid(k2) * params.sig2_beta(2, 2) * Nmat_median(2) + ...
                      poisson_sigma_grid(k3) * params.sig2_beta(2, 3) * Nmat_median(2) + ...
                      params.sig2_zero(2);
                a12 = poisson_sigma_grid(k3) * sqrt(prod(params.sig2_beta(:, 3))) * params.rho_beta(:, 3) * sqrt(prod(Nmat_median)) + ...
                      sqrt(prod(params.sig2_zero)) * params.rho_zero;
                      
                pdf(snps_in_chunk) = pdf(snps_in_chunk) + ...
                    sigma_grid_pi{1}(:, k1) .* sigma_grid_pi{2}(:, k2) .* sigma_grid_pi{3}(:, k3) .* ...
                    fast_normpdf2(zmat(snps_in_chunk, 1), zmat(snps_in_chunk, 2), a11, a12, a22);
            end; end; end
        end

        % Likelihood term, weighted by inverse TLD
        weights = 1 ./ w_ld;
        weights(w_ld < 1) = 1;
        cost = sum(weights .* -log(pdf));
        if ~isfinite(cost), cost = NaN; end;
        if ~isreal(cost), cost = NaN; end;
    end

    if options.use_legacy_impl
        unique_non_zero = @(x)unique(x(x~=0));
        sig2_beta1 = unique_non_zero(params.sig2_beta(1, :));
        sig2_beta2 = unique_non_zero(params.sig2_beta(2, :));
        rho_beta   = unique_non_zero(params.rho_beta);

        % Covariance matrix [a11 a12; a12 a22]
        a11 = repmat(sig2_beta1, [num_snps 1]);
        a12 = repmat(sqrt(sig2_beta1*sig2_beta2) * rho_beta, [num_snps 1]);
        a22 = repmat(sig2_beta2, [num_snps 1]);
        if length(params.pi_vec) == 1
            pivec = repmat([0 0 params.pi_vec], [num_snps, 1]); pi0 = 1-sum(pivec, 2);
        else
            pivec = repmat(params.pi_vec, [num_snps, 1]); pi0 = 1-sum(pivec, 2);
        end

        sum_p = @(p1, p2)(1-(1-p1).*(1-p2));

        %new approximation
        pi1 = pivec(:, 1); pi2 = pivec(:, 2); pi3 = pivec(:, 3); pi0 = 1 - pi1 - pi2 - pi3;

        f1 = (pi1+pi3) .* ref_ld_sum_r2 + (pi0+pi2) .* ref_ld_chi_r4;
        f2 = (pi2+pi3) .* ref_ld_sum_r2 + (pi0+pi1) .* ref_ld_chi_r4;

        pi1_plus_pi3_eff = (pi1+pi3) .* ref_ld_sum_r2 ./ f1;
        pi2_plus_pi3_eff = (pi2+pi3) .* ref_ld_sum_r2 ./ f2;
        %pi3_eff = sum_p(pi3 .* ref_ld_r2 ./ (pi3 .* ref_ld_r2 + (pi0 + pi1 + pi2) .* r2eff), pi1_plus_pi3_eff .* pi2_plus_pi3_eff);
        pi3_eff = sum_p(ref_ld_sum_r2 .* ref_ld_chi_r4 .* (pi3 - pi1 .* pi2) ./ (f1 .* f2), pi1_plus_pi3_eff .* pi2_plus_pi3_eff);
        pi1_eff = max(pi1_plus_pi3_eff - pi3_eff, 0);
        pi2_eff = max(pi2_plus_pi3_eff - pi3_eff, 0);
        pivec = [pi1_eff pi2_eff pi3_eff]; pi0 = max(1-sum(pivec, 2), 0);
        a11 = a11 .* f1;
        a22 = a22 .* f2;
        %a12 = a12 .* pi3 .* ref_ld_r2 ./ (pi3_eff .* sqrt(a11 .* a22));
        a12 = sqrt(a11 .* a22) * rho_beta;

        % Covariance matrix [a11 a12; a12 a22]
        a11 = a11 .* Hvec .* Nmat(:, 1);
        a12 = a12 .* Hvec .* sqrt(Nmat(:, 1) .* Nmat(:, 2));
        a22 = a22 .* Hvec .* Nmat(:, 2);

        % Components
        pdf0 = pi0         .* fast_normpdf2(zmat(:, 1), zmat(:, 2),       sigma0_sqr(1,1),       sigma0_sqr(1,2),       sigma0_sqr(2,2));
        pdf1 = pivec(:, 1) .* fast_normpdf2(zmat(:, 1), zmat(:, 2), a11 + sigma0_sqr(1,1),       sigma0_sqr(1,2),       sigma0_sqr(2,2));
        pdf2 = pivec(:, 2) .* fast_normpdf2(zmat(:, 1), zmat(:, 2),       sigma0_sqr(1,1),       sigma0_sqr(1,2), a22 + sigma0_sqr(2,2));
        pdf3 = pivec(:, 3) .* fast_normpdf2(zmat(:, 1), zmat(:, 2), a11 + sigma0_sqr(1,1), a12 + sigma0_sqr(1,2), a22 + sigma0_sqr(2,2));
        pdf_legacy  = pdf0 + pdf1 + pdf2 + pdf3;
        %plot(pdf, pdf_legacy, '.')
        %plot(-log10(pdf), -log10(pdf_legacy), '.')

        cost_legacy = sum(weights .* -log(pdf_legacy));
        if ~isfinite(cost_legacy), cost_legacy = NaN; end;
        if ~isreal(cost_legacy), cost_legacy = NaN; end;
        if isnan(cost_legacy)
            fprintf('stop');
        end

        %fprintf('cost: %.5e vs %.5e, delta = %.5e\n', cost, cost_legacy, cost - cost_legacy);
        pdf = pdf_legacy;
        cost = cost_legacy;
    end

    if options.use_convolution
        distr   = @(p)(p ./ sum(p));
        dirac   = @(n)([zeros(1, (n-1)/2), 1, zeros(1,(n-1)/2)]);

        distr2   = @(p2)(p2 ./ sum(p2(:)));
        conv2_n  = @(p2, n)BGMG_util.op_power(p2, @(a,b)conv2(a, b, 'same'), n);
        diracX   = @(p)([zeros((length(p)-1)/2, length(p)); p; zeros((length(p)-1)/2, length(p))]);
        diracY   = @(p)(diracX(p)');
        dirac2   = @(n)(diracX(dirac(n)));
        mesh2X   = @(z)BGMG_util.colvec(meshgrid(z, z));
        mesh2Y   = @(z)BGMG_util.colvec(meshgrid(z, z)');
        mix2_c0  = @(z)(dirac2(length(z)));
        mix2_c1  = @(z, sig2)(diracX(distr(normpdf(z, 0, sig2))));
        mix2_c2  = @(z, sig2)(diracY(distr(normpdf(z, 0, sig2))));
        mix2_c3  = @(z, sig2, rho)(reshape(fast_normpdf2(mesh2X(z), mesh2Y(z), sig2(1), sqrt(prod(sig2)) * rho, sig2(2)), [length(z) length(z)]));
        mixture2 = @(z, pi_vec, sig2, rho)((1-sum(pi_vec))*mix2_c0(z) + pi_vec(1)*mix2_c1(z, sig2(1)) + ...
                                                                        pi_vec(2)*mix2_c2(z, sig2(2)) + ...
                                                                        pi_vec(3)*distr2(mix2_c3(z, sig2, rho)));
        % imagesc(-log10(conv2_n(mixture2(-6:0.1:6, [0.003 0.003 0.01], [1 1], 0.99), 200)), [0 5])

        pdf_by_conv = nan(size(pdf));

        ld_block = round(ref_ld_sum_r2./ref_ld_chi_r4);
        ld_sigm1 = sqrt(Nmat(:, 1).*Hvec.*ref_ld_chi_r4.*params.sig2_beta(1));
        ld_sigm2 = sqrt(Nmat(:, 2).*Hvec.*ref_ld_chi_r4.*params.sig2_beta(1));

        num_ld_block_bins = 10;
        num_ld_sigma_bins = 15;

        ld_sigm1_bins = [-Inf quantile(ld_sigm1, num_ld_sigma_bins) Inf];
        ld_sigm2_bins = [-Inf quantile(ld_sigm2, num_ld_sigma_bins) Inf];
        ld_block_bins = [-Inf unique(floor(logspace(0, log10(max(ld_block)+1), num_ld_block_bins)))];

        for ld_sigm1_bini = 2:length(ld_sigm1_bins)
        for ld_sigm2_bini = 2:length(ld_sigm2_bins)
            idx1 = (ld_sigm1 > ld_sigm1_bins(ld_sigm1_bini - 1)) & (ld_sigm1 <= ld_sigm1_bins(ld_sigm1_bini));
            idx2 = (ld_sigm2 > ld_sigm2_bins(ld_sigm2_bini - 1)) & (ld_sigm2 <= ld_sigm2_bins(ld_sigm2_bini));
            if sum(idx1 & idx2) == 0, continue; end;

        for ld_block_bini = 2:length(ld_block_bins)
            idx3 = (ld_block > ld_block_bins(ld_block_bini - 1)) & (ld_block <= ld_block_bins(ld_block_bini));
            idx = idx1 & idx2 & idx3;
            if (sum(idx) == 0), continue; end;

            fprintf('%i %i %i\n', ld_sigm1_bini, ld_sigm2_bini, ld_block_bini);

            mean_ld_block = round(mean(ld_block(idx)));
            mean_ld_sigm1 = mean(ld_sigm1(idx));
            mean_ld_sigm2 = mean(ld_sigm2(idx));

            z_max = max(100, round(max([zmat(idx, 1) ./ mean_ld_sigm1; zmat(idx, 2) ./ mean_ld_sigm2])) + 1);
            if z_max <= 100, step = 1;
            elseif z_max <= 200, step = 2;
            elseif z_max <= 500, step = 5;
            elseif z_max <= 1000, step = 10;
            else warning('too large z_max'); step = 10;
            end
            z = -z_max:step:z_max; z = z * min(mean_ld_sigm1, mean_ld_sigm2);

            table_to_pdf_factor = sum(sum(mix2_c3(z, [1 1], 0)));
            pdf_bini = table_to_pdf_factor * conv2(conv2_n(mixture2(z, params.pi_vec, [mean_ld_sigm1 mean_ld_sigm2], params.rho_beta(3)), mean_ld_block), distr2(mix2_c3(z, params.sig2_zero, params.rho_zero)), 'same');

            % Hack-hack to avoid concentrating the entire distribution at zero.
            pdf_bini((length(z)+1)/2, :) = pdf_bini((length(z)-1)/2, :);
            pdf_bini(:, (length(z)+1)/2) = pdf_bini(:, (length(z)-1)/2);
            pdf_bini((length(z)+1)/2, (length(z)+1)/2) = pdf_bini((length(z)-1)/2, (length(z)-1)/2);

            pdf_by_conv(idx) = interp2(z, z, pdf_bini, zmat(idx, 1), zmat(idx, 2));
            if any(isnan(pdf_by_conv(idx))),
                fprintf('error %i %i %i\n', ld_sigm1_bini, ld_sigm2_bini, ld_block_bini)
            end
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
    end

    if (nargout > 1)
        % probability density function
        result.pdf = pdf;
    end
    
    % TBD: all options below this line are broken - must be fixed.
    % Global FDR (False Discovery Rate)
    if false && (nargout > 1) && options.calculate_global_fdr
        fprintf('Estimate global false discovery rates...\n');
        snps=size(zmat, 1);
        z_grid = -options.delta_hat_std_limit:options.delta_hat_std_step:options.delta_hat_std_limit;
        
        result.condfdr_global = nan(size(zmat));
        result.conjfdr_global = nan(size(zmat, 1), 1);
        for snpi=1:snps
            if (mod(snpi, 1000) == 0),  fprintf('\tFinish %i SNPs out of %i\n', snpi, snps); end;
            for condi = 1:2
                if condi == 1, z1vec = z_grid; z2vec = repmat(zmat(snpi, 2), size(z_grid)); end;
                if condi == 2, z2vec = z_grid; z1vec = repmat(zmat(snpi, 1), size(z_grid)); end;

                % Conditional probability density functions for each component of the mixture
                pdf0 = pi0(snpi)      .* fast_normpdf2(z1vec, z2vec,             sigma0_sqr(1,1),             sigma0_sqr(1,2),             sigma0_sqr(2,2));
                pdf1 = pivec(snpi, 1) .* fast_normpdf2(z1vec, z2vec, a11(snpi) + sigma0_sqr(1,1),             sigma0_sqr(1,2),             sigma0_sqr(2,2));
                pdf2 = pivec(snpi, 2) .* fast_normpdf2(z1vec, z2vec,             sigma0_sqr(1,1),             sigma0_sqr(1,2), a22(snpi) + sigma0_sqr(2,2));
                pdf3 = pivec(snpi, 3) .* fast_normpdf2(z1vec, z2vec, a11(snpi) + sigma0_sqr(1,1), a12(snpi) + sigma0_sqr(1,2), a22(snpi) + sigma0_sqr(2,2));
                pdf  = pdf0 + pdf1 + pdf2 + pdf3;
                
                % Conditional cumulated distribution functions for each component of the mixture
                cdf0 = cumsum(pdf0) ./ sum(pdf); 
                cdf1 = cumsum(pdf1) ./ sum(pdf);
                cdf2 = cumsum(pdf2) ./ sum(pdf);
                cdf3 = cumsum(pdf3) ./ sum(pdf);
                cdf  = cdf0 + cdf1 + cdf2 + cdf3;

                % Calculate CDF from the tail to avoid values dangerously close to 1.0
                if interp1(z_grid, cdf, zmat(snpi, condi)) > 0.5
                    cdf0 = fliplr(cumsum(fliplr(pdf0))) ./ sum(pdf); 
                    cdf1 = fliplr(cumsum(fliplr(pdf1))) ./ sum(pdf);
                    cdf2 = fliplr(cumsum(fliplr(pdf2))) ./ sum(pdf);
                    cdf3 = fliplr(cumsum(fliplr(pdf3))) ./ sum(pdf);
                    cdf  = cdf0 + cdf1 + cdf2 + cdf3;
                end
                
                if condi == 1, FDR = (cdf0 + cdf2) ./ cdf; end;
                if condi == 2, FDR = (cdf0 + cdf1) ./ cdf; end;
                result.condfdr_global(snpi, condi) = interp1(z_grid, FDR, zmat(snpi, condi));
            end
        end
        result.conjfdr_global = max(result.condfdr_global, [], 2);
    end
    
    % Posterior effect size estimates
    if false && (nargout > 1) && options.calculate_delta_hat
        fprintf('Estimate posterior effect sizes...\n');
        snps=size(zmat, 1);
        result.delta_hat_mean = nan(size(zmat));
        result.delta_hat_var = nan(size(zmat));
        result.delta_hat_cov = nan(size(zmat, 1), 1);
        result.delta_hat_median = nan(size(zmat));
        for snpi=1:snps
            if (mod(snpi, 100) == 0),  fprintf('\tFinish %i SNPs out of %i\n', snpi, snps); end;
            delta1_grid = (-15:0.5:15);
            delta2_grid = (-15:0.5:15);
            [x, y] = meshgrid(delta1_grid, delta2_grid);
            middle_index = (length(delta1_grid)+1)/2;  % assume length(x) is odd
            if (mod(length(delta1_grid), 2) ~= 1), error('length(x) must be odd'); end;
            
            delta1vec = normpdf(delta1_grid, 0, sqrt(a11(snpi))); delta1vec=delta1vec/sum(delta1vec);
            delta2vec = normpdf(delta2_grid, 0, sqrt(a22(snpi))); delta2vec=delta2vec/sum(delta2vec);
            delta0 = zeros(size(x)); delta0(middle_index, middle_index) = pi0(snpi);
            delta1 = zeros(size(x)); delta1(:, middle_index) = pivec(snpi, 1) * delta1vec';
            delta2 = zeros(size(x)); delta2(middle_index, :) = pivec(snpi, 2) * delta2vec;
            delta3 = reshape(mvnpdf([x(:), y(:)], 0, [a11(snpi) a12(snpi); a12(snpi) a22(snpi)]), size(x)); delta3 = pivec(snpi, 3) * delta3 / sum(delta3(:));
            delta_prior = delta0 + delta1 + delta2 + delta3;  % imagesc(-log10(delta_prior))
            
            sigma0 = [params.sigma0(1).^2, params.rho0 .* prod(params.sigma0); params.rho0 .* prod(params.sigma0), params.sigma0(2).^2];
            like_func = reshape(mvnpdf(repmat(zmat(snpi, :), [numel(x), 1]), ...
                                [x(:), y(:)], ...
                                sigma0_sqr), size(x));
            delta_posterior = delta_prior .* like_func;
            delta_posterior = delta_posterior ./ sum(delta_posterior(:));  % normalize
            % imagesc(-log10(like_func))
            % imagesc(-log10(delta_posterior))

            result.delta_hat_mean(snpi, 1) = sum(delta_posterior(:) .* x(:));
            result.delta_hat_mean(snpi, 2) = sum(delta_posterior(:) .* y(:));
            result.delta_hat_var(snpi, 1)  = sum(delta_posterior(:) .* (x(:) - result.delta_hat_mean(snpi, 1)).^2);
            result.delta_hat_var(snpi, 2)  = sum(delta_posterior(:) .* (y(:) - result.delta_hat_mean(snpi, 2)).^2);
            result.delta_hat_cov(snpi, 1)      = sum(delta_posterior(:) .* (x(:) - result.delta_hat_mean(snpi, 1)) .* (y(:) - result.delta_hat_mean(snpi, 2)));
            result.delta_hat_median(snpi, 1) = delta1_grid(find(cumsum(sum(delta_posterior, 1))>0.5, 1, 'first'));  % calc median (works as a filter for noisy observations)
            result.delta_hat_median(snpi, 2) = delta2_grid(find(cumsum(sum(delta_posterior, 2))>0.5, 1, 'first'));  % calc median (works as a filter for noisy observations)
        end
    end
    
    % project pdf back into the original SNPs indexing (which may include undefined values)
    if exist('result', 'var'), result = restore_original_indexing(result, defvec); end

    if false && exist('result', 'var')
        result.params_per_snp.a11 = a11;
        result.params_per_snp.a12 = a12;
        result.params_per_snp.a22 = a22;
        result.params_per_snp.pi0 = pi0;
        result.params_per_snp.pivec = pivec;
    end
    
    if exist('result', 'var') && isfield(options, 'speedup_info'), result.speedup_info = options.speedup_info; end

    if options.verbose
        filt = @(x)unique(x(x~=0));
        fprintf('Bivariate : pi_vec=[%s], rho_beta=[%s], sig2_beta1=[%s], (eff. [%s]), sig2_beta2=[%s], (eff. [%s]), rho_zero=%.3f, sig2_zero=[%s], cost=%.3e\n', ...
            sprintf('%.3e ', params.pi_vec), ...
            sprintf('%.3f ', filt(params.rho_beta)), ...
            sprintf('%.2e ', filt(params.sig2_beta(1, :))), ...
            sprintf('%.2f ', filt(mean(a11))), ...
            sprintf('%.2e ', filt(params.sig2_beta(2, :))), ...
            sprintf('%.2f ', filt(mean(a22))), ...
            params.rho_zero, sprintf('%.3f ', params.sig2_zero), cost);

        if ~isnan(options.total_het)
            pleio_idx = all(params.sig2_beta ~= 0);
            fprintf('\th2=[%s], h2(pleio)=[%s]\n', ...
                sprintf('%.3f ', (params.sig2_beta*params.pi_vec')*options.total_het), ...
                sprintf('%.3f ', (params.sig2_beta(:, pleio_idx)*params.pi_vec(pleio_idx))*options.total_het));
        end
    end
end

function pdf = pdf_mix1(pi_vec, zmat, a11, a12, a22, sigma0_sqr)
        pdf = ...
            (1 - pi_vec(:, 1)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)                                    , ...
                    sigma0_sqr(1,2)                                    , ...
                    sigma0_sqr(2,2)                                    ) ...
            + ...
            (0 + pi_vec(:, 1)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1)                        , ...
                    sigma0_sqr(1,2) + a12(:, 1)                        , ...
                    sigma0_sqr(2,2) + a22(:, 1)                        );

end

function pdf = pdf_mix2(pi_vec, zmat, a11, a12, a22, sigma0_sqr)
        pdf = ...
            (1 - pi_vec(:, 1)) .* (1 - pi_vec(:, 2)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)                                    , ...
                    sigma0_sqr(1,2)                                    , ...
                    sigma0_sqr(2,2)                                    ) ...
            + ...
            (1 - pi_vec(:, 1)) .* (0 + pi_vec(:, 2)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)             + a11(:, 2)            , ...
                    sigma0_sqr(1,2)             + a12(:, 2)            , ...
                    sigma0_sqr(2,2)             + a22(:, 2)            ) ...
            + ...
            (0 + pi_vec(:, 1)) .* (1 - pi_vec(:, 2)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1)                        , ...
                    sigma0_sqr(1,2) + a12(:, 1)                        , ...
                    sigma0_sqr(2,2) + a22(:, 1)                        ) ...
            + ...
            (0 + pi_vec(:, 1)) .* (0 + pi_vec(:, 2)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1) + a11(:, 2)            , ...
                    sigma0_sqr(1,2) + a12(:, 1) + a12(:, 2)            , ...
                    sigma0_sqr(2,2) + a22(:, 1) + a22(:, 2)            );

end

function pdf = pdf_mix3(pi_vec, zmat, a11, a12, a22, sigma0_sqr)
        pdf = ...
            (1 - pi_vec(:, 1)) .* (1 - pi_vec(:, 2)) .* (1 - pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)                                    , ...
                    sigma0_sqr(1,2)                                    , ...
                    sigma0_sqr(2,2)                                    ) ...
            + ...
            (1 - pi_vec(:, 1)) .* (1 - pi_vec(:, 2)) .* (0 + pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)                         + a11(:, 3), ...
                    sigma0_sqr(1,2)                         + a12(:, 3), ...
                    sigma0_sqr(2,2)                         + a22(:, 3)) ...
            + ...
            (1 - pi_vec(:, 1)) .* (0 + pi_vec(:, 2)) .* (1 - pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)             + a11(:, 2)            , ...
                    sigma0_sqr(1,2)             + a12(:, 2)            , ...
                    sigma0_sqr(2,2)             + a22(:, 2)            ) ...
            + ...
            (1 - pi_vec(:, 1)) .* (0 + pi_vec(:, 2)) .* (0 + pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1)             + a11(:, 2) + a11(:, 3), ...
                    sigma0_sqr(1,2)             + a12(:, 2) + a12(:, 3), ...
                    sigma0_sqr(2,2)             + a22(:, 2) + a22(:, 3)) ...
            + ...
            (0 + pi_vec(:, 1)) .* (1 - pi_vec(:, 2)) .* (1 - pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1)                        , ...
                    sigma0_sqr(1,2) + a12(:, 1)                        , ...
                    sigma0_sqr(2,2) + a22(:, 1)                        ) ...
            + ...
            (0 + pi_vec(:, 1)) .* (1 - pi_vec(:, 2)) .* (0 + pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1)             + a11(:, 3), ...
                    sigma0_sqr(1,2) + a12(:, 1)             + a12(:, 3), ...
                    sigma0_sqr(2,2) + a22(:, 1)             + a22(:, 3)) ...
            + ...
            (0 + pi_vec(:, 1)) .* (0 + pi_vec(:, 2)) .* (1 - pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1) + a11(:, 2)            , ...
                    sigma0_sqr(1,2) + a12(:, 1) + a12(:, 2)            , ...
                    sigma0_sqr(2,2) + a22(:, 1) + a22(:, 2)            ) ...
            + ...
            (0 + pi_vec(:, 1)) .* (0 + pi_vec(:, 2)) .* (0 + pi_vec(:, 3)) .* ...
                fast_normpdf2(zmat(:, 1), zmat(:, 2), ...
                    sigma0_sqr(1,1) + a11(:, 1) + a11(:, 2) + a11(:, 3), ...
                    sigma0_sqr(1,2) + a12(:, 1) + a12(:, 2) + a12(:, 3), ...
                    sigma0_sqr(2,2) + a22(:, 1) + a22(:, 2) + a22(:, 3));
end

function pdf = fast_normpdf2(x1, x2, a11, a12, a22)
    % Calculation of log-likelihood and pdf, specific to bivariate normal
    % distribution with zero mean. It takes into account an explicit formula
    % for inverse 2x2 matrix, A = [a b; c d],  => A^-1 = [d -b; -c a] ./ det(A)
    dt = a11.*a22 - a12.*a12;  % det(A)
    log_exp = -0.5 * (a22.*x1.*x1 + a11.*x2.*x2 - 2.0*a12.*x1.*x2) ./ dt;
    log_dt  = -0.5 * log(dt);
    log_pi  = -1.0 * log(2*pi);
    pdf = exp(log_pi+log_dt+log_exp);
    if ~isreal(pdf)
        warning('complex numbers');
    end
    % loglike = sum(log_dt + log_exp) + length(dt) * log_pi;
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
