function [cost, result] = BGMG_univariate_cost2(params, zvec, Hvec, Nvec, w_ld, ref_ld, options)
    % INPUT:
    % params    - struct with fields pi_vec, sig2_zero, sig2_beta
    %             pi_vec --- mixture weights, one per mixture component (excluding null component)
    %             sig2_zero --- inflation from cryptic relatedness / sample stratification (sigma0^2)
    %             sig2_beta --- discoverability, one per mixture component (excluding null component)
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld    - a struct with fields sum_r2, sum_r2_biased, sum_r4_biased
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options
    %
    % OUTPUT
    % cost      - log-likelihood

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;    % enable or disable verbose logging
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate
    if ~isfield(options, 'zmax'), options.zmax = +Inf; end;           % ignore z scores above zmax

    if ~isfield(options, 'poisson_sigma_grid_nodes'), options.poisson_sigma_grid_nodes = 25;  end;  % grid size of poisson projection
    if ~isfield(options, 'poisson_sigma_grid_scale'), options.poisson_sigma_grid_scale = 1;  end;

    non_zero_mixi = (params.sig2_beta ~=0) & (params.pi_vec ~= 0);
    params.sig2_beta = params.sig2_beta(non_zero_mixi);
    params.pi_vec = params.pi_vec(non_zero_mixi);

    % Validate input params
    result=[];
    if ~BGMG_util.validate_params(params); cost = nan; return; end;

    ref_ld_sum_r2 = ref_ld.sum_r2;
    ref_ld_chi_r4 = ref_ld.sum_r4_biased ./ ref_ld.sum_r2_biased;
    ref_ld_chi_r4(ref_ld.sum_r4_biased == 0) = 0;

    defvec = isfinite(zvec + Hvec + Nvec + w_ld + sum(ref_ld_sum_r2, 2) + sum(ref_ld_chi_r4, 2)) & (Hvec > 0);
    defvec(any(abs(zvec) > options.zmax, 2)) = false;
    defvec(any(ref_ld_chi_r4 < 0, 2) | any(ref_ld_chi_r4 > 1, 2) | any(ref_ld_sum_r2 < 0, 2)) = false;

    zvec = zvec(defvec, :); Hvec = Hvec(defvec); Nvec = Nvec(defvec, :);
    w_ld = w_ld(defvec); ref_ld_sum_r2 = ref_ld_sum_r2(defvec, :); ref_ld_chi_r4 = ref_ld_chi_r4(defvec, :);
    
    num_traits = length(params.sig2_zero);  % number of traits in the anslysis
    num_mix    = length(params.pi_vec);  % number of components in the mixture
    num_snps   = length(Hvec);           % number of SNPs in the analysis
    num_ldr2bins = size(ref_ld_chi_r4, 2);  % number of bins in LDr2 histogram
    if num_traits ~= 1, error('Params define more than one trait, unable to calculate univariate cost.'); end;

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

    snp_chunk_size  = floor(num_snps/20);  % TBD - what's reasonable number of SNPs chunks?
    [~, sorted_snps] = sort(snp_contribution);

    pdf = zeros(num_snps, 1);

    poisson_sigma_grid_nodes = (options.poisson_sigma_grid_scale.^2) * options.poisson_sigma_grid_nodes;

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
        poisson_sigma_grid_delta = poisson_sigma_grid_limit / (poisson_sigma_grid_nodes - 1);
        poisson_sigma_grid       = 0:poisson_sigma_grid_delta:poisson_sigma_grid_limit;

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
    end

    % Likelihood term, weighted by inverse TLD
    weights = 1 ./ w_ld;
    weights(w_ld < 1) = 1;
    cost = sum(weights .* -log(pdf));
    if ~isfinite(cost), cost = NaN; end;
    if ~isreal(cost), cost = NaN; end;

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
