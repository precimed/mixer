function [zvec, zprobvec_conv, zprobvec_amd] = GMM_ofrei_univariate_cost_point(x, Hvec, Nvec, ref_ld_hist, ref_ld_bins, mapparams, options)
    % x         - vectors of params
    % mapparams - function that maps x into a struct with fields pi1, sigma_beta, sigma0
    % zvec      - vector of z scores, one per SNP
    % Hvec      - heterozigosity per SNP
    % Nvec      - number of subjects genotyped per SNP
    % ref_ld_hist - LD score (total LD) per SNP, calculated across all genotypes or imputed SNPs,
    % ref_ld_bins - calculated as a histogram binned by levels of ld r2
    % w_ld      - LD score (total LD) per SNP, calculated across SNPs included in the template for the analysis
    % options   - options; see the below for a full list of configurable options

    params = mapparams(x);
    if params.sigma0 == 0, error('sigma0 can not be zero'); end;
    if params.sigma_beta == 0, error('sigma_beta can not be zero'); end;
    if length(Hvec) ~= 1 || length(Nvec) ~= 1 || size(ref_ld_hist, 1) ~= 1, error('Can not handle multiple SNPs in this code'); end;
    
    % List of all configurable options
    if ~exist('options', 'var'), options = struct(); end;
    if ~isfield(options, 'x_limit'), options.x_limit = 20; end;             % tabulate probability distributions for [-xlim : xstep : xlim]
    if ~isfield(options, 'x_step'), options.x_step = 0.01; end;
    if ~isfield(options, 'ld_bins'), options.ld_bins = 30; end;             % number of LD bins to use

    pivec = params.pivec;              % proportion of causal SNPs
    sbt = params.sigma_beta;           % 'discoverability', aka variance of the effect sizes
    sigma0 = params.sigma0;            % variance from non-heritable (noise) component

    conv_n = @(x, n)op_power(x, @(a,b)conv(a, b, 'same'), n);
    sbt = sbt * mean(sqrt(Nvec * Hvec));

    % prior distribution 
    x = (-options.x_limit:options.x_step:options.x_limit);
    middle_index = (length(x)+1)/2;  % assume length(x) is odd
    if (mod(length(x), 2) ~= 1), error('length(x) must be odd'); end;

    % noise component
    noise = normpdf(x, 0, sigma0); noise = noise / sum(noise);

    % convolutions according to LD histogram
    zvec = x; zprobvec_conv = noise;
    for tld_hist_bini = 1:size(ref_ld_hist, 2)
        mean_ld_r2 = mean(ref_ld_bins(tld_hist_bini, :));
        ref_ld = round(ref_ld_hist(:, tld_hist_bini) / mean_ld_r2);
        if ref_ld >= 1
            beta = normpdf(x, 0, sqrt(mean_ld_r2) * sbt); beta = beta * pivec / sum(beta);
            beta(middle_index) = beta(middle_index) + (1 - pivec);
            zprobvec_conv = conv(conv_n(beta, ref_ld), zprobvec_conv, 'same');
        end
    end
    
    % approximation by AMD
    Lvec = sum(ref_ld_hist); pi0effvec = (1-pivec).^Lvec;
    if pivec < 1
        pi1effvec = 1-pi0effvec;
        %pvec_tmp = pi0effvec.*normpdf(zvec,0,sig0)+pi1effvec.*normpdf(zvec,0,sqrt(sig0^2+Hvec*sigb^2));
        zprobvec_amd = pi0effvec.*normpdf(zvec,0,sigma0)+pi1effvec.*normpdf(zvec,0,sqrt(sigma0^2+sbt^2));
    else
        %pvec_tmp = normpdf(zvec,0,sqrt(sig0^2+Lvec.*Hvec*sigb^2));
        zprobvec_amd = normpdf(zvec,0,sqrt(sigma0^2+Lvec.*sbt^2));
    end

    fprintf('UVT: pi1 = %.3e, sigma_beta^2 = %.3e (eff. %.3e), sigma_0 = %.3e\n', pivec, params.sigma_beta.^2, sbt.^2, sigma0);
end
