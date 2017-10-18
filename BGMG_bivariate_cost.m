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
    %   .sum_r4 - Similar to ref_ld.sum_r2, but contains a sum of r_{ij}^4
    %             ref_ld.sum_r4(j) = sum_i r^4_{ij}
    %             WARNING: the sum_r2 and sum_r4 must be consistent with each other
    % options   - options; see the below for a full list of configurable options

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'use_convolution'), options.use_convolution = false; end;  % experimental VERY SLOW option to calculate pdf via convolution

    % delta_hat_std_limit and delta_hat_std_step are used in posterior effect size
    % estimation. They express the grid to calculate posterior delta
    % distribution before taking mean or median statistic.
    if ~isfield(options, 'delta_hat_std_limit'), options.delta_hat_std_limit = 20; end;
    if ~isfield(options, 'delta_hat_std_step'),  options.delta_hat_std_step  = 0.5; end;
    if ~isfield(options, 'calculate_delta_hat'), options.calculate_delta_hat = true; end;
    if ~isfield(options, 'calculate_global_fdr'), options.calculate_global_fdr = true; end;
    if ~isfield(options, 'total_het'), options.total_het = nan; end;  % required for heritability estimate

    % Validate input params
    if ~BGMG_util.validate_params(params); cost = nan; return; end;
    if isempty(ref_ld), ref_ld_r2 = zero(size(w_ld)); ref_ld_r4 = zero(size(w_ld));
    else ref_ld_r2 = ref_ld.sum_r2; ref_ld_r4 = ref_ld.sum_r4; end;

    r2eff   = ref_ld_r4 ./ ref_ld_r2;

    defvec = isfinite(sum(zmat, 2) + Hvec + sum(Nmat, 2) + w_ld + ref_ld_r2 + ref_ld_r4) & (Hvec > 0);
    defvec((r2eff < 0) | (r2eff > 1) | (ref_ld_r4 < 0) | (ref_ld_r2 < 0)) = false;

    zmat = zmat(defvec, :); Hvec = Hvec(defvec); Nmat = Nmat(defvec, :);
    w_ld = w_ld(defvec); ref_ld_r2 = ref_ld_r2(defvec); ref_ld_r4 = ref_ld_r4(defvec); r2eff = r2eff(defvec);

    num_traits = length(params.sig2_zero);  % number of traits in the anslysis
    num_mix    = length(params.pi_vec);  % number of components in the mixture
    num_snps   = length(Hvec);           % number of SNPs in the analysis
    if num_traits ~= 2, error('Params must define two traits, unable to calculate bivariate cost.'); end;

    sigma0_sqr = [params.sig2_zero(1), sqrt(prod(params.sig2_zero)) * params.rho_zero; ...
                                   sqrt(prod(params.sig2_zero)) * params.rho_zero, params.sig2_zero(2)];

    % Approximation that preserves variance and kurtosis
    pi_vec     = repmat(params.pi_vec, [num_snps 1]);
    eta_factor = (pi_vec .* repmat(ref_ld_r2, [1, num_mix]) + (1-pi_vec) .* repmat(r2eff, [1, num_mix]));
    pi_vec     = (pi_vec .* repmat(ref_ld_r2, [1, num_mix])) ./ eta_factor;
        
    % Normalize pi vector, and compensate total variance by adjusting eta_factor
    [pi0, pi_vec_norm] = BGMG_util.normalize_pi_vec(pi_vec);
    eta_factor  = eta_factor .* (pi_vec ./ pi_vec_norm);
    eta_factor(pi_vec == 0) = 0;

    % Multiply sig2_beta by eta_factor
    a11=zeros(num_snps, num_mix); a22=a11; a12=a11;
    for mixi = 1:num_mix
        a11(:, mixi) = eta_factor(:, mixi) .* Hvec .* Nmat(:, 1) * params.sig2_beta(1, mixi);
        a22(:, mixi) = eta_factor(:, mixi) .* Hvec .* Nmat(:, 2) * params.sig2_beta(2, mixi);
        a12(:, mixi) = sqrt(a11(:, mixi) .* a22(:, mixi)) * params.rho_beta(1, mixi);
    end
    
    % Calculate p.d.f. from the mixture model
    pdf = pi0 .* fast_normpdf2(zmat(:, 1), zmat(:, 2), sigma0_sqr(1,1), sigma0_sqr(1,2), sigma0_sqr(2,2));
    for mixi = 1:num_mix
        pdf = pdf + pi_vec_norm(:, mixi) .* fast_normpdf2(zmat(:, 1), zmat(:, 2), a11(:, mixi) + sigma0_sqr(1,1), a12(:, mixi) + sigma0_sqr(1,2), a22(:, mixi) + sigma0_sqr(2,2));
    end

    % Likelihood term, weighted by inverse TLD
    weights = 1 ./ w_ld;
    weights(w_ld < 1) = 1;
    cost = sum(weights .* -log(pdf));
    if ~isfinite(cost), cost = NaN; end;
    if ~isreal(cost), cost = NaN; end;

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

        ld_block = round(ref_ld_r2./r2eff);
        ld_sigm1 = sqrt(Nmat(:, 1).*Hvec.*r2eff.*params.sig2_beta(1));
        ld_sigm2 = sqrt(Nmat(:, 2).*Hvec.*r2eff.*params.sig2_beta(1));

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

    % TBD: all options below this line are broken - must be fixed.
    if nargout > 1
        % probability density function
        result.pdf = pdf;
    
        % local conditional/conjunctive fdr
        result.condfdr_local = [pdf0 + pdf2, pdf0 + pdf1] ./ [pdf pdf];
        result.conjfdr_local  = (pdf0 + pdf1 + pdf2) ./ pdf;
    end
    
    % Global FDR (False Discovery Rate)
    if (nargout > 1) && options.calculate_global_fdr
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
    if (nargout > 1) && options.calculate_delta_hat
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

    if exist('result', 'var')
        result.params_per_snp.a11 = a11;
        result.params_per_snp.a12 = a12;
        result.params_per_snp.a22 = a22;
        result.params_per_snp.pi0 = pi0;
        result.params_per_snp.pivec = pivec;
    end

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

function pdf = fast_normpdf2(x1, x2, a11, a12, a22)
    % Calculation of log-likelihood and pdf, specific to bivariate normal
    % distribution with zero mean. It takes into account an explicit formula
    % for inverse 2x2 matrix, A = [a b; c d],  => A^-1 = [d -b; -c a] ./ det(A)
    dt = a11.*a22 - a12.*a12;  % det(A)
    log_exp = -0.5 * (a22.*x1.*x1 + a11.*x2.*x2 - 2.0*a12.*x1.*x2) ./ dt;
    log_dt  = -0.5 * log(dt);
    log_pi  = -1.0 * log(2*pi);
    pdf = exp(log_pi+log_dt+log_exp);
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
