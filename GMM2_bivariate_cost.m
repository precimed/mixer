function [cost, result] = GMM2_bivariate_cost(params, zmat, Hvec, Nmat, w_ld, ref_ld, mapparams, options)
    % params    - params struct with fields pivec, sigma_beta, rho_beta, sigma0, rho0
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
    if ~isfield(options, 'zmax'), options.zmax = 12; end;  % throw away z scores larger than this value
    if ~isfield(options, 'verbose'), options.verbose = false; end;  % enable or disable verbose logging
    if ~isfield(options, 'use_ref_ld'), options.use_ref_ld = true; end; % enable approximation that uses ref_ld

    % beta_hat_std_limit and beta_hat_std_step are used in posterior effect size
    % estimation. They express the grid to calculate posterior beta
    % distribution before taking mean or median statistic.
    % The values are expressed in 'normalized' form, e.g. in units of the
    % params.sigma_beta.
    if ~isfield(options, 'beta_hat_std_limit'), options.beta_hat_std_limit = 10; end;
    if ~isfield(options, 'beta_hat_std_step'),  options.beta_hat_std_step  = 0.1; end;
    if ~isfield(options, 'calculate_beta_hat'), options.calculate_beta_hat = true; end;
    if ~isfield(options, 'calculate_global_fdr'), options.calculate_global_fdr = true; end;

    if ~isstruct(params), params = mapparams(params); end;
    % params.sigma0     - vector with 2 elements in (0, +inf]
    % params.rho0       - scalar in [-1, 1]
    % params.sigma_beta - vector with 2 elements in (0, +inf]
    % params.rho_beta   - scalar in [-1, 1]
    % params.pivec      - vector with 3 elements in [0, 1]

    % Validate input params
    if any(params.sigma0 <= 0), warning('sigma0 can not be zero'); cost = nan; return; end;
    if any(params.sigma_beta) <= 0, warning('sigma_beta can not be zero'); cost = nan; return; end;
    if (params.rho0 < -1) || (params.rho0 > 1), warning('rho0 must be in [-1, 1]'); end;
    if (params.rho_beta < -1) || (params.rho_beta > 1), warning('rho_beta must be in [-1, 1]'); cost = nan; return; end;
    if any(params.pivec < 0), warning('pivec must be from 0 to 1'); cost = nan; return; end;
    if sum(params.pivec) > 1, warning('sum(pivec) can not exceed 1'); cost = nan; return; end;
    if isempty(ref_ld), ref_ld_r2 = zero(size(w_ld)); ref_ld_r4 = zero(size(w_ld));
    else ref_ld_r2 = ref_ld.sum_r2; ref_ld_r4 = ref_ld.sum_r4; end;

    r2eff   = ref_ld_r4 ./ ref_ld_r2;

    defvec = isfinite(sum(zmat, 2) + Hvec + sum(Nmat, 2) + w_ld + ref_ld_r2 + ref_ld_r4) & (Hvec > 0);
    defvec(any(abs(zmat) > options.zmax, 2)) = false;
    defvec((r2eff < 0) | (r2eff > 1) | (ref_ld_r4 < 0) | (ref_ld_r2 < 0)) = false;

    zmat = zmat(defvec, :); Hvec = Hvec(defvec); Nmat = Nmat(defvec, :);
    w_ld = w_ld(defvec); ref_ld_r2 = ref_ld_r2(defvec); ref_ld_r4 = ref_ld_r4(defvec); r2eff = r2eff(defvec);

    sigma0_sqr = [params.sigma0(1).^2, prod(params.sigma0) * params.rho0; prod(params.sigma0) * params.rho0, params.sigma0(2).^2];
    
    % Covariance matrix [a11 a12; a12 a22]
    a11 = repmat(params.sigma_beta(1).^2, size(Hvec));
    a12 = repmat(prod(params.sigma_beta) * params.rho_beta, size(Hvec)); 
    a22 = repmat(params.sigma_beta(2).^2, size(Hvec));

    pivec = repmat(params.pivec, size(Hvec)); pi0 = 1-sum(pivec, 2);

    if options.use_ref_ld,
        %previous approximation from GMM code
        %pi0 = pi0 .^ ref_ld_r2;
        %pivec = pivec .* repmat((1 - pi0) ./ sum(pivec, 2), [1 3]);
        
        sum_p = @(p1, p2)(1-(1-p1).*(1-p2));
        
        %new approximation
        pi1 = pivec(:, 1); pi2 = pivec(:, 2); pi3 = pivec(:, 3); pi0 = 1 - pi1 - pi2 - pi3;
        
        f1 = (pi1+pi3) .* ref_ld_r2 + (pi0+pi2) .* r2eff;
        f2 = (pi2+pi3) .* ref_ld_r2 + (pi0+pi1) .* r2eff;
        
        pi1_plus_pi3_eff = (pi1+pi3) .* ref_ld_r2 ./ f1;
        pi2_plus_pi3_eff = (pi2+pi3) .* ref_ld_r2 ./ f2;
        %pi3_eff = sum_p(pi3 .* ref_ld_r2 ./ (pi3 .* ref_ld_r2 + (pi0 + pi1 + pi2) .* r2eff), pi1_plus_pi3_eff .* pi2_plus_pi3_eff);
        pi3_eff = sum_p(ref_ld_r2 .* r2eff .* (pi3 - pi1 .* pi2) ./ (f1 .* f2), pi1_plus_pi3_eff .* pi2_plus_pi3_eff);
        pi1_eff = max(pi1_plus_pi3_eff - pi3_eff, 0);
        pi2_eff = max(pi2_plus_pi3_eff - pi3_eff, 0);
        pivec = [pi1_eff pi2_eff pi3_eff]; pi0 = max(1-sum(pivec, 2), 0);
        a11 = a11 .* f1;
        a22 = a22 .* f2;
        %a12 = a12 .* pi3 .* ref_ld_r2 ./ (pi3_eff .* sqrt(a11 .* a22));
        a12 = sqrt(a11 .* a22) * params.rho_beta;
    end

    % Covariance matrix [a11 a12; a12 a22]
    a11 = a11 .* Hvec .* Nmat(:, 1);
    a12 = a12 .* Hvec .* sqrt(Nmat(:, 1) .* Nmat(:, 2));
    a22 = a22 .* Hvec .* Nmat(:, 2);
    
    % Components
    pdf0 = pi0         .* fast_normpdf2(zmat(:, 1), zmat(:, 2),       sigma0_sqr(1,1),       sigma0_sqr(1,2),       sigma0_sqr(2,2));
    pdf1 = pivec(:, 1) .* fast_normpdf2(zmat(:, 1), zmat(:, 2), a11 + sigma0_sqr(1,1),       sigma0_sqr(1,2),       sigma0_sqr(2,2));
    pdf2 = pivec(:, 2) .* fast_normpdf2(zmat(:, 1), zmat(:, 2),       sigma0_sqr(1,1),       sigma0_sqr(1,2), a22 + sigma0_sqr(2,2));
    pdf3 = pivec(:, 3) .* fast_normpdf2(zmat(:, 1), zmat(:, 2), a11 + sigma0_sqr(1,1), a12 + sigma0_sqr(1,2), a22 + sigma0_sqr(2,2));
    pdf  = pdf0 + pdf1 + pdf2 + pdf3;
    
    % Likelihood term, weighted by inverse TLD
    weights = 1 ./ w_ld;
    weights(w_ld < 1) = 1;
    cost = sum(weights .* -log(pdf));
    if ~isfinite(cost), cost = NaN; end;
    if ~isreal(cost), cost = NaN; end;

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
        z_grid = -options.beta_hat_std_limit:options.beta_hat_std_step:options.beta_hat_std_limit;
        
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
    if (nargout > 1) && options.calculate_beta_hat
        fprintf('Estimate posterior effect sizes...\n');
        snps=size(zmat, 1);
        result.beta_hat_mean = nan(size(zmat));
        result.beta_hat_var = nan(size(zmat));
        result.beta_hat_cov = nan(size(zmat, 1), 1);
        result.beta_hat_median = nan(size(zmat));
        for snpi=1:snps
            if (mod(snpi, 1000) == 0),  fprintf('\tFinish %i SNPs out of %i\n', snpi, snps); end;
            beta1_grid = (-options.beta_hat_std_limit:options.beta_hat_std_step:options.beta_hat_std_limit) * params.sigma_beta(1);
            beta2_grid = (-options.beta_hat_std_limit:options.beta_hat_std_step:options.beta_hat_std_limit) * params.sigma_beta(2);
            [x, y] = meshgrid(beta1_grid, beta2_grid);
            middle_index = (length(beta1_grid)+1)/2;  % assume length(x) is odd
            if (mod(length(beta1_grid), 2) ~= 1), error('length(x) must be odd'); end;
            
            beta1vec = normpdf(beta1_grid, 0, sqrt(a11(snpi))); beta1vec=beta1vec/sum(beta1vec);
            beta2vec = normpdf(beta2_grid, 0, sqrt(a22(snpi))); beta2vec=beta2vec/sum(beta2vec);
            beta0 = zeros(size(x)); beta0(middle_index, middle_index) = pi0(snpi);
            beta1 = zeros(size(x)); beta1(:, middle_index) = pivec(snpi, 1) * beta1vec';
            beta2 = zeros(size(x)); beta2(middle_index, :) = pivec(snpi, 2) * beta2vec;
            beta3 = reshape(mvnpdf([x(:), y(:)], 0, [a11(snpi) a12(snpi); a12(snpi) a22(snpi)]), size(x)); beta3 = pivec(snpi, 3) * beta3 / sum(beta3(:));
            beta_prior = beta0 + beta1 + beta2 + beta3;  % imagesc(-log10(beta))
            
            like_func = reshape(mvnpdf(repmat(zmat(snpi, :), [numel(x), 1]), ...
                                [x(:) .* sqrt(Nmat(snpi, 1) * Hvec(snpi)), y(:) .* sqrt(Nmat(snpi, 2) * Hvec(snpi))], ...
                                params.sigma0), size(x));
            beta_posterior = beta_prior .* like_func;
            beta_posterior = beta_posterior ./ sum(beta_posterior(:));  % normalize
            % imagesc(-log10(beta_posterior))

            result.beta_hat_mean(snpi, 1) = sum(beta_posterior(:) .* x(:));  
            result.beta_hat_mean(snpi, 2) = sum(beta_posterior(:) .* y(:));  
            result.beta_hat_var(snpi, 1)  = sum(beta_posterior(:) .* (x(:) - result.beta_hat_mean(snpi, 1)).^2);
            result.beta_hat_var(snpi, 2)  = sum(beta_posterior(:) .* (y(:) - result.beta_hat_mean(snpi, 2)).^2);
            result.beta_hat_cov(snpi)      = sum(beta_posterior(:) .* (x(:) - result.beta_hat_mean(snpi, 1)) .* (y(:) - result.beta_hat_mean(snpi, 2)));
            result.beta_hat_median(snpi, 1) = beta1_grid(find(cumsum(sum(beta_posterior, 1))>0.5, 1, 'first'));  % calc median (works as a filter for noisy observations)
            result.beta_hat_median(snpi, 2) = beta2_grid(find(cumsum(sum(beta_posterior, 2))>0.5, 1, 'first'));  % calc median (works as a filter for noisy observations)
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
        fprintf('Bivariate : pi=[%s], rho_beta=%.3f, sigma_beta^2=[%s], (eff. [%s]), rho0=%.3f, sigma0^2=[%s], cost=%.3e\n', ...
            sprintf('%.3f ', params.pivec), params.rho_beta, sprintf('%.2e ', params.sigma_beta.^2), ...
            sprintf('%.2e ', [mean(a11),mean(a22)]), params.rho0, sprintf('%.3f ', params.sigma0.^2),cost);
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
