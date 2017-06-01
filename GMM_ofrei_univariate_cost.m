function [cost, zprobvec] = GMM_ofrei_univariate_cost(x, zvec, Hvec, Nvec, ref_ld, w_ld, mapparams, options)
    % x         - vectors of params
    % mapparams - function that maps x into a struct with fields pi1, sigma_beta, sigma0
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld    - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options
    %
    % ====================
    % IMPLEMENTATION NOTES
    % ====================
    % There are many different approaches on how the cost can be
    % calculated. Some of the options that I evaluated include
    % - analytical approach (can be enabled by options.convolution = false)
    % - convolution approach, where p(z|ld) is calculated for a table of
    %   LD values. Currently I ignore dependence on heterozigosity and
    %   sample size.
    %   * active option: use table of LD values
    %   * commented out: all LD values from 1 to max(ref_ld)
    % - convolution approach, where I aim to calculate p(z|ld,het,N) 
    %   turned out to be really slow to calculate (seems o require a for
    %   loop across all SNPs),in addition it require complex logic to
    %   change the domain (values along X axis where I tabulate probability
    %   distributions).
    %   two options that I've evaluated:
    %   * z          ~ sqrt(NH) * sum r_i b_i + eps                (1)
    %   * z/sqrt(NH) ~            sum r_i b_i + eps/sqrt(NH)       (2)
    % =====================================================================
    
    params = mapparams(x);
    if params.sigma0 == 0, error('sigma0 can not be zero'); end;
    if params.sigma_beta == 0, error('sigma_beta can not be zero'); end;

    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'convolution'), options.convolution = true; end;   % a flag indicating whether to use convolution approach (default) or analytical approach
    if ~isfield(options, 'x_limit'), options.x_limit = 20; end;             % tabulate probability distributions for [-xlim : xstep : xlim]
    if ~isfield(options, 'x_step'), options.x_step = 0.01; end;
    if ~isfield(options, 'ld_bins'), options.ld_bins = 30; end;             % number of LD bins to use
    if ~isfield(options, 'mean_ld_r2'), options.mean_ld_r2 = 1.0; end;      % mean ld r2

    ref_ld = ref_ld ./ options.mean_ld_r2;
    defvec = isfinite(zvec + Hvec + Nvec + ref_ld + w_ld) & (Hvec > 0);
    defvec(abs(zvec) > 12) = false;
    %defvec(10000:end) = false;
    zvec = zvec(defvec); Hvec = Hvec(defvec); Nvec = Nvec(defvec); ref_ld = ref_ld(defvec); w_ld = w_ld(defvec);

    pivec = params.pivec;              % proportion of causal SNPs
    sbt = params.sigma_beta;           % 'discoverability', aka variance of the effect sizes
    sigma0 = params.sigma0;            % variance from non-heritable (noise) component
    nsnp = length(zvec);               % number of SNPs
    max_ref_ld = ceil(max(ref_ld));    % maximum TLD across SNPs

    if options.convolution % new convolution approach
        conv_n = @(x, n)op_power(x, @(a,b)conv(a, b, 'same'), n);

        sbt = sbt .* mean(sqrt(Nvec .* Hvec));  % Hack-hack, use the mean here.
        
        % prior distribution 
        x = (-options.x_limit:options.x_step:options.x_limit); % * sbt;   % <------ scale x by sbt.
        middle_index = (length(x)+1)/2;  % assume length(x) is odd
        if (mod(length(x), 2) ~= 1), error('length(x) must be odd'); end;
        beta = normpdf(x, 0, sqrt(options.mean_ld_r2) * sbt); beta = beta * pivec / sum(beta);
        beta(middle_index) = beta(middle_index) + (1 - pivec);

        % noise component
        noise = normpdf(x, 0, sigma0); noise = noise / sum(noise);
        
        % calculate probability density for $z = sum_{i=1}^ld beta_i + noise$
        % across exponentially-spaced set of ld values (see ld_table)
        ld_table = unique(ceil(logspace(0, log10(max_ref_ld + 1), options.ld_bins)));
        beta_conv = zeros(length(ld_table), length(x));
        for ld_index=1:length(ld_table)
            beta_conv(ld_index, :) = conv(conv_n(beta, ld_table(ld_index)), noise, 'same');
        end
        
        % perform lookup of z probabilities
        zprobvec = nan(size(zvec)); 
        for ld_index=1:length(ld_table)
            index = ref_ld <= ld_table(ld_index);
            if (ld_index > 1), index = index & ref_ld > ld_table(ld_index - 1); end;
            zprobvec(index) = interp1(x, beta_conv(ld_index, :), zvec(index)) / (x(2) - x(1));
        end
        if any(isnan(zprobvec)), error('failed to calculate p(z) for some SNPs'); end;
        
        if 0 % use all LD values (slow)
            % sum of LD beta distributions (for each LD up to the largest observed in the data)
            beta_conv = zeros(max_ref_ld, length(x));
            beta_conv(1, :) = conv(beta, noise, 'same');
            for ld=2:max_ref_ld
                beta_conv(ld, :) = conv(beta_conv(ld - 1, :), beta, 'same');
            end

            % perform lookup of z probabilities
            zprobvec = zeros(size(zvec));
            for ld=1:max_ref_ld
                index = (ld == max(1, ceil(ref_ld)));
                zprobvec(index) = interp1(x, beta_conv(ld, :), zvec(index)) / (x(2) - x(1));
            end
        end
        
        if 0 % interesting approach with scaling each Z, e.g. calc (2) instead of (1).
             %  z          ~ sqrt(NH) * sum r_i b_i + eps                (1)
             %  z/sqrt(NH) ~            sum r_i b_i + eps/sqrt(NH)       (2)
             % For this to work 'x' must be scaled
            zprobmat = beta_conv(max(1, ceil(ref_ld)), :);
            noisemat = normpdf(repmat(x, [nsnp, 1]), 0, repmat(sigma0 ./ sqrt(Nvec .* Hvec), [1 length(x)])); noisemat = noisemat ./ repmat(sum(noisemat, 2), [1 length(x)]);
            for i=1:nsnp, zprobmat(i, :) = conv(zprobmat(i, :), noisemat(i, :), 'same'); end;
            for i=1:nsnp, zprobvec(i   ) = interp1(x, zprobmat(i, :), zvec(i) ./ sqrt(Nvec(i) .* Hvec(i))); end;
            zprobvec = zprobvec * sbt ./ (x(2) - x(1));
            zprobvec(isnan(zprobvec)) = 0;
        end
        
        if 0 % interesting approach with re-scaling for each SNP. Correspond to (1) above.
            x_z = zvec; %x_z = -20:0.01:20;
            noise = normpdf(x_z, 0, sigma0); noise = noise / sum(noise);

            for i=1:nsnp
                ld_i = ceil(ref_ld(i));
                beta_i = beta_conv(ld_i, :);
                x_i = x * sqrt(Nvec(i) * Hvec(i));
                z_i = interp1(x_i, beta_i, x_z); z_i(isnan(z_i)) = 0; z_i = z_i / sum(z_i);
                z_i = conv(z_i, noise, 'same'); z_i = z_i / (x_z(2) - x_z(1));
                p = interp1(x_z, z_i, zvec(i));
                zprobvec = z_i; break;
            end
        end
        %image(-log10(beta_conv))
    else  % old 'analytical' version without convolution
        %if ~isfield(options, 'nbins'), options.nbins = 1; end;
        %params.pivec = params.pivec .* options.nbins;
        %params.sigma_beta = params.sigma_beta .* sqrt(options.nbins);
        %ref_ld = ref_ld ./ options.nbins;

        ts = min(50, max_ref_ld);          % table size; rows = total LD, columns = index K in summation (no more than 500, limit due to pi1 << 1)

        % Find pimat -- matrix of size "max_ref_ld X ts", containing binomial probabilities 
        pimat  = binopdf(repmat(0:ts, [max_ref_ld, 1]), repmat(1:max_ref_ld, [ts+1, 1])', pivec);
        if 0, image(-log10(pimat)); end;

        % Find zprobmat -- matrix of size "nsnp X ts", containing normal density
        zprobmat = normpdf(repmat(zvec, [1, ts+1]), 0, sqrt(sigma0 .^ 2 + repmat(Nvec .* Hvec, [1, ts+1]) .* repmat(0:ts, [nsnp, 1]) .* ((options.mean_ld_r2 * sbt).^2)));

        % Find zprobvec --- product of zprobmat and pimat
        zprobvec = sum(zprobmat .* pimat(max(1, ceil(ref_ld)), :), 2);
    end
    
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