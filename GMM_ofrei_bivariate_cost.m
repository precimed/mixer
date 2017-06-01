function [cost, zprobvec] = GMM_ofrei_bivariate_cost(x, zvec, Hvec, Nvec, w_ld, ref_ld, mapparams, options)
    % x         - vectors of params
    % mapparams - function that maps x into a struct with fields pi1, sigma_beta, sigma0
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld    - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options

    params = mapparams(x);
    if params.sigma0 == 0, error('sigma0 can not be zero'); end;
    if params.sigma_beta == 0, error('sigma_beta can not be zero'); end;

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'x_limit'), options.x_limit = 20; end;             % tabulate probability distributions for [-xlim : xstep : xlim]
    if ~isfield(options, 'x_step'), options.x_step = 0.5; end;
    if ~isfield(options, 'ld_bins'), options.ld_bins = 30; end;             % number of LD bins to use

    defvec = isfinite(sum(zvec, 2) + Hvec + sum(Nvec, 2) + ref_ld + w_ld) & (Hvec > 0);
    defvec(any(abs(zvec) > 12, 2)) = false;
    %defvec(10000:end) = false;
    zvec = zvec(defvec, :); Hvec = Hvec(defvec); Nvec = Nvec(defvec, :); ref_ld = ref_ld(defvec); w_ld = w_ld(defvec);

    sigma0_sqr = [params.sigma0(1).^2, prod(params.sigma0) * params.rho0; prod(params.sigma0) * params.rho0, params.sigma0(2).^2];
    sbt_sqr    = [params.sigma_beta(1).^2, prod(params.sigma_beta) * params.rho_beta; prod(params.sigma_beta) * params.rho_beta, params.sigma_beta(2).^2];
    pivec = params.pivec;              % proportion of causal SNPs
    nsnp = size(zvec, 1);              % number of SNPs
    max_ref_ld = ceil(max(ref_ld));    % maximum TLD across SNPs

    conv2_n = @(x, n)op_power(x, @(a,b)conv2(a, b, 'same'), n);

    sbt_sqr = sbt_sqr .* [mean(Hvec .* Nvec(:, 1)), mean(Hvec .* sqrt(Nvec(:, 1) .* Nvec(:, 2))); mean(Hvec .* sqrt(Nvec(:, 1) .* Nvec(:, 2))), mean(Hvec .* Nvec(:, 2))]; % Hack-hack, use the mean here.

    % grid definition
    p = (-options.x_limit:options.x_step:options.x_limit);
    middle_index = (length(p)+1)/2;  % assume length(x) is odd
    if (mod(length(p), 2) ~= 1), error('length(x) must be odd'); end;
    [x,y] = meshgrid(p, p);
    mnvhist = @(sigma)reshape(mvnpdf([x(:), y(:)], 0, sigma), size(x));
    
    % prior beta distribution
    beta1vec = normpdf(p, 0, sqrt(sbt_sqr(1,1))); beta1vec=beta1vec/sum(beta1vec);
    beta2vec = normpdf(p, 0, sqrt(sbt_sqr(2,2))); beta2vec=beta2vec/sum(beta2vec);
    beta0 = zeros(size(x)); beta0(middle_index, middle_index) = 1 - sum(pivec);
    beta1 = zeros(size(x)); beta1(:, middle_index) = pivec(1) * beta1vec';
    beta2 = zeros(size(x)); beta2(middle_index, :) = pivec(2) * beta2vec;
    beta3 = mnvhist(sbt_sqr); beta3 = pivec(3) * beta3 / sum(beta3(:));
    beta = beta0 + beta1 + beta2 + beta3;  % imagesc(beta)
    
    % noise component
    noise = mnvhist(sigma0_sqr); noise = noise / sum(noise(:));    % imagesc(noise)

    % calculate probability density for $z = sum_{i=1}^ld beta_i + noise$
    % across exponentially-spaced set of ld values (see ld_table)
    ld_table = unique(ceil(logspace(0, log10(max_ref_ld + 1), options.ld_bins)));
    beta_conv = zeros(length(ld_table), size(x,1), size(x,2));
    for ld_index=1:length(ld_table)
        beta_conv(ld_index, :, :) = conv2(conv2_n(beta, ld_table(ld_index)), noise, 'same');
    end

    % perform lookup of z probabilities
    zprobvec = nan(nsnp, 1); 
    for ld_index=1:length(ld_table)
        index = ref_ld <= ld_table(ld_index);
        if (ld_index > 1), index = index & ref_ld > ld_table(ld_index - 1); end;
        zprobvec(index) = interp2(x, y, squeeze(beta_conv(ld_index, :, :)), zvec(index, 1), zvec(index, 2)) / ((p(2) - p(1)).^2);
    end
    if any(isnan(zprobvec)), error('failed to calculate p(z) for some SNPs'); end;

    % Likelihood term, weighted by inverse TLD
    weights = 1 ./ w_ld;
    cost = sum(weights .* -log(zprobvec));
    if ~isfinite(cost), cost = NaN; end;

    % project zprobvec back into the original SNPs indexing (which may include undefined values)
    tmp = zprobvec; zprobvec = nan(size(defvec)); zprobvec(defvec) = tmp;
    fprintf('BVT: pi=[%s], rho_beta=%.3f, sigma_beta^2=[%s], (eff. [%s]), rho0=%.3f, sigma0^2=[%s], cost=%.3e\n',sprintf('%.3f ', params.pivec), params.rho_beta, sprintf('%.2e ', params.sigma_beta.^2), sprintf('%.2e ', [sbt_sqr(1,1),sbt_sqr(2,2)]), params.rho0, sprintf('%.3f ', params.sigma0.^2),cost);
end


% TBD
%   - bin convolution calculations by sample size and heterozigosity (today just use mean)
%   - fix how we interpolate; use interp3 instead of interp2, because actually we lookup in 3D cube (x,y,LD)