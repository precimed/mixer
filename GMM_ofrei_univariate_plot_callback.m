function GMM_ofrei_univariate_plot_callback(options)
    if isfield(options, 'figure'), figure(options.figure); end;
    if isfield(options, 'axes'), axes(options.axes); end;
    subplot(2,2,1);

    % random subset of real LD histogram data.
     ref_ld_bins = [     0    0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000    0.9000;
                    0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000    0.9000    1.0000]';
     ref_ld_hist= [ ...
         0   13.6925   46.9704   21.2493   20.2187    1.6896    2.5692    2.3290   86.0306   44.6573;
         0    8.0521    3.9957    3.9949    1.8220    1.0022   16.6029    1.4494    2.5828   19.1663;
         0    1.6821    0.6669    0.3485         0    0.5856         0         0    1.6719    0.9598;
         0    8.7599   12.0511    4.1182    6.5419    3.5115    1.3193    3.1959         0    0.9706;
         0    1.8571    0.9014    0.3407    2.8109    6.2142    0.6986    3.8677   20.1679    1.8628;
         0   35.9966   28.7127   15.4653    1.2350    0.5100    1.2317         0    0.8113    9.9216;
         0   21.5186    3.0016    3.1965    0.9390    1.6970    9.7302    6.5146         0  112.9237;
         0   23.4158   22.4857    9.6881    4.9655    9.3562    7.2677    0.7011   27.7945   23.8091;
         0    4.4609    0.4738    1.0502    1.8982    2.5876         0         0         0    0.9199;
         0    3.6792    0.5480    0.6624    0.9595    6.8139    2.6086    0.7421    6.0515  287.2514;
         0    5.7819   11.4288         0    0.8915    0.5811    0.6819    1.5285    2.6350   26.3666;
         0   54.0642   40.9922   15.5923    2.8189    1.7074   17.7413   13.6445    8.0463   35.2606;
         0    3.8190    1.9604    0.3064    0.4468    1.5882    5.2965   21.4333   18.4806    2.8718;
         0    7.4225    8.6598   11.8033   15.4125    4.6817    2.5975    2.1383   43.7271   21.6123;
         0   17.7934    5.5324    9.1396    2.9075    5.4349    4.5940    9.7000    0.8014         0;
         0   20.3508   18.7842   22.8421   29.7822    3.3660   69.8273  167.5867    4.3151  235.8990;
         0   26.8614    7.2157    0.3276    1.8152    5.6400         0  118.1763    9.7634  189.5813;
         0   34.4090   30.6860   12.5587   14.0265   25.2630   25.2159    2.8914    0.8043   18.3686;
         0    8.0833    7.7403   18.8471   13.1084    5.4936    0.6762         0         0   23.5222;
         0   39.0505   25.1188    2.2550    0.9114   18.3475    3.7302    0.7504    1.6325    2.0000;
         0    9.5006    3.7541    3.0494    0.9238    1.6612    6.2244         0    1.6484    5.8065;
         0   13.4614    6.8903    6.2702         0    2.1391         0    5.3500    6.4546   46.6131;
         0    0.8345    3.1741    0.3769    0.4642    0.5045         0         0    0.8978    5.8751;
         0    6.5929    1.3540    1.9919         0         0    0.6620    0.7581    1.6560         0;
         0   12.0768    2.7138    1.8198    0.8972         0         0         0    0.8446    2.7609;
         0   16.9136    5.7712    1.6455    3.2308   11.6662    0.6804    0.7474    2.5459    4.8922;
         0    4.9505    3.3844    3.7933    5.6353   25.5442    3.8609    1.4929         0    7.8946;
         0   14.2447    2.7252    2.3738    4.3374   23.5081    1.2781         0   11.2067   78.0392;
         0    4.1394    2.0413         0    2.2031         0    2.6179         0         0   12.6291;
         0   13.4970   23.0710    3.0612   18.5326    3.0488    1.9270    0.7949   14.1193   18.8800;
         0    0.6225    0.2958    2.0907    4.9145   10.0530    0.6245    0.7907   12.5599    0.9438;
         0    3.6817    0.2910    0.3152         0         0         0         0         0    4.9361;
         0   41.6150   29.6188    9.4812   44.3233   25.3781    3.6133    1.5798    8.3071   14.5326;
         0    1.6568    2.4235    1.8768    3.2984    4.0995         0         0   15.6376   44.4551;
         0         0         0         0         0         0         0         0         0         0;
         0   23.0030    1.6078    1.8660    0.4332         0         0         0         0    7.8671;
         0   20.1490   11.0402   10.0918    1.9020    4.8620    1.3473    3.6202    0.8956   14.6881;
         0    4.6812    4.1642    6.9701    0.4376    1.6716   21.3510    5.9557    3.3386   22.6019;
         0         0    0.2780         0         0         0         0         0         0    5.9086;
         0   15.6773   13.2673    2.3190    1.3699    1.0922    1.3718         0         0   10.4867;
         0   27.2184    8.5346   18.2832   88.8391   37.7401   18.3397    0.7062         0    8.8545;
         0    6.8780   15.3547    4.7372    4.3856    6.6970    1.8733   20.0175    3.4146    4.9032;
         0   20.6507   27.0746   19.1184    7.8859   11.5453         0    2.2111    0.8971    0.9757;
         0   18.0836    2.3951   11.0134    1.3218         0    0.6346    1.5821    9.3207    6.6813;
         0    6.9015   19.9963    9.0308    9.6658         0    0.6855    2.2439    0.8572    0.9576;
         0   14.5537   17.6961    4.9067    1.7760    1.1030   24.1025    1.4801    3.4801   35.2977;
         0    1.2679    0.2169    0.7322    1.3777    2.0881   25.2647   24.4849   17.2439   94.4538;
         0    1.5933    2.9926    4.7706    0.8272    0.5965    3.2639    0.7939    3.4731    9.9711;
         0   22.1820    9.2705   13.7072    5.4083   10.9352    1.3551   16.2692   25.0309   38.3096;
         0    2.6307    2.0263    1.4247         0    1.5589         0         0    0.8735    7.7660];
         
    step=0.01;
    m_zvec  = (-10:step:10)'; 
    m_hvec = ones(size(m_zvec)) * 0.25;
    m_nvec = ones(size(m_zvec)) * options.N;
    m_ref_ld = ones(size(m_zvec)) * options.tld;
    m_w_ld = ones(size(m_zvec)) * options.tld;
    mapparams = @GMM_ofrei_univariate_mapparams;

    ref_ld_hist = ref_ld_hist(floor(options.mean_ld_r2 * size(ref_ld_hist, 1)), :);
    if options.tld > 0, ref_ld_hist = ref_ld_hist * options.tld / sum(ref_ld_hist); end;

    if 0
        fprintf('TLD=%i, mean ld r2=%.2f, N=%i\n', options.tld, options.mean_ld_r2, floor(options.N));
        [~, zprobvec] = GMM_ofrei_univariate_cost(mapparams(options.opts_struct), m_zvec, m_hvec, m_nvec, m_ref_ld, m_w_ld, mapparams, options);
        normprob = normpdf(m_zvec, 0, options.opts_struct.sigma0);

        if options.use_logscale, plot(m_zvec, log(zprobvec), m_zvec, log(normprob));
        else plot(m_zvec, zprobvec, m_zvec, normprob);
        end
        %legend('model', 'null', 'Location', 'NorthWest');
        xlabel('z');
        ylabel('p(z)');
        title(sprintf('TLD=%.1f, mean ld r2=%.2f, N=%.0f, pi1=%.2e, sbt=%.2e, sig0=%.2f', options.tld, options.mean_ld_r2, floor(options.N), options.opts_struct.pivec, options.opts_struct.sigma_beta, options.opts_struct.sigma0));
    else
        fprintf('TLD=%i, N=%i\n', sum(ref_ld_hist), floor(options.N));
        [m_zvec, zprobvec_conv, zprobvec_amd] = GMM_ofrei_univariate_cost_point(mapparams(options.opts_struct), mean(m_hvec), mean(m_nvec), ref_ld_hist, ref_ld_bins, mapparams, options);
        bar(mean(ref_ld_bins, 2), ref_ld_hist);
    	xlabel('r2');
        ylabel('TLD|r2');
        title(sprintf('TLD=%.1f, N=%.0f, pi1=%.2e, sbt=%.2e, sig0=%.2f', sum(ref_ld_hist), floor(options.N), options.opts_struct.pivec, options.opts_struct.sigma_beta, options.opts_struct.sigma0));
    end
    
    phi  = @(sig)normpdf(m_zvec, 0, sig);
    %pi=1e-3; sig0t=1.00; sig1t=1.3; pvec = (1-pi)*phi(sig0t)+pi*phi(sig1t);
    KL=@(p,q)sum(p .* (log(p) - log(q))); symKL=@(p,q)(KL(p,q)+KL(q,p));
    kurtosis=@(z,p)(sum(p .* z.^4) / sum(p .* z.^2)^2);
    variance=@(z,p)(sum(p .* z.^2));
    
    if 0
    cost = @(x)symKL((1-x(1))*phi(x(2))+x(1)*phi(x(3)), zprobvec_conv); % x = [pi1, sig0, sig1]
    x = fminsearch(cost, [0.1 1 2], struct('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8));
    if x(2)>x(3), tmp=x(2); x(2)=x(3); x(3)=tmp; x(1)=1-x(1); end; disp(x);
    zprobvec_fit=(1-x(1))*phi(x(2)) + x(1)*phi(x(3));
    elseif 0
        sig0 = options.opts_struct.sigma0;
    cost = @(x)symKL((1-x(1))*phi(sig0)+x(1)*phi(x(2)), zprobvec_conv); % x = [pi1, sig0, sig1]
    x = fminsearch(cost, [0.1 2], struct('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8));
    zprobvec_fit=(1-x(1))*phi(sig0) + x(1)*phi(x(2));
    x = [x(1) sig0 x(2)];
    else
        r2eff = sum(mean(ref_ld_bins, 2)' .* ref_ld_hist) ./ sum(ref_ld_hist);
        tld   = sum(ref_ld_hist);
        pi1   = options.opts_struct.pivec; pi0 = 1-pi1;
        x(1) = pi1 * tld / (pi1 * tld + pi0 * r2eff);
        x(2) = options.opts_struct.sigma0;
        x(3) = options.opts_struct.sigma_beta * sqrt(mean(m_hvec) * mean(m_nvec)) * sqrt(pi0 * r2eff + pi1 * tld);
        zprobvec_fit=(1-x(1))*phi(x(2)) + x(1)*phi(sqrt(x(2).^2 + x(3).^2));
    end

    %[sigma0, sigma_delta, pi1] = em(m_zvec,zprobvec_conv);
     %pf(m_zvec,zprobvec_conv);
    %costfunc = @(x)sum(...
   %     (log10(zprobvec_conv) - log10(x(1) * normpdf(m_zvec,0,x(2)) + (1-x(1)) * normpdf(m_zvec,0,x(3)))).^2 ...
   %     );
   %  x = fminsearch(costfunc, [0.01, 1, 2], struct('Display', 'on'))
     
    hs=subplot(2,2,2); cla(hs); % qq plot
    hold on
    qqlim=16;
    if exist('zprobvec', 'var'), plot(-log10(cumsum(zprobvec/sum(zprobvec))),-log10(normcdf(m_zvec,0,1)),'b'); end;
    if exist('zprobvec_conv', 'var'), plot(-log10(cumsum(zprobvec_conv/sum(zprobvec_conv))),-log10(normcdf(m_zvec,0,1)),'b'); end;
    if exist('zprobvec_amd', 'var'), plot(-log10(cumsum(zprobvec_amd/sum(zprobvec_amd))),  -log10(normcdf(m_zvec,0,1)),'r'); end;
    if exist('zprobvec_fit', 'var'), plot(-log10(cumsum(zprobvec_fit/sum(zprobvec_fit))),  -log10(normcdf(m_zvec,0,1)),'g'); end;
    plot([0 qqlim],[0 qqlim], 'k');
    xlim([0 10]); ylim([0 15]);
    legend('model-conv', 'model-amd', 'model-kl', 'null', 'Location', 'SouthEast');
    
    model_kl_str = sprintf('model-kl: kurtosis=%.3f, variance=%.3f',kurtosis(m_zvec, zprobvec_fit/sum(zprobvec_fit)),variance(m_zvec, zprobvec_fit/sum(zprobvec_fit)));
    model_conv_str = sprintf('model-conv: kurtosis=%.3f, variance=%.3f',kurtosis(m_zvec, zprobvec_conv/sum(zprobvec_conv)),variance(m_zvec, zprobvec_conv/sum(zprobvec_conv)));
    model_amd_str = sprintf('model-amd: kurtosis=%.3f, variance=%.3f',kurtosis(m_zvec, zprobvec_amd/sum(zprobvec_amd)),variance(m_zvec, zprobvec_amd/sum(zprobvec_amd)));
    model_kl_params = sprintf('model-kl: pi=%.4f, sigma0=%.3f, sigma_delta=%.3f',x(1),x(2),x(3));
    title(sprintf('%s\n%s\n%s\n%s', model_conv_str, model_amd_str, model_kl_str, model_kl_params));
    xlabel('Expected -log_1_0(P_v_a_l_u_e)');
    ylabel('Model -log_1_0(P_v_a_l_u_e)');
    
    
    
    %qq_x = norminv(logspace(-10,-0.1,100));
    %qq_x = -log10(normcdf(m_zvec,0,1));
    %qq_y = -log10(zprobvec);
    %qqlim=10; plot(qq_y, qq_x,[0 qqlim],[0 qqlim]); xlim([0 qqlim]); ylim([0 qqlim]);
end

function [sigma0, sigma_delta, pi1] = pf(z,p)
    to_pdf = @(x)(x / sum(x));
    
    sigma0 = sqrt(sum(p .* z.^2));
    p1 = to_pdf(abs(p - to_pdf(normpdf(z, 0, sigma0)))); 
    
    sigma1 = sqrt(sum(p1 .* z.^2));
    p2 = to_pdf(abs(p1 - to_pdf(normpdf(z, 0, sigma0))));
    
    sigma2 =  sqrt(sum(p2 .* z.^2));
end

function [sigma0, sigma_delta, pi1] = em(z,p)
    sigma0sqr = 1; sigma1sqr = 2; pi1 = 0.1;
    for i=1:1000
        p0 = (1-pi1) * normpdf(z, 0, sqrt(sigma0sqr));
        p1 = pi1 * normpdf(z, 0, sqrt(sigma1sqr));
        gamma = p1 ./ (p0 + p1);
        sigma0sqr = sum(p .* (1-gamma) .* z.^2) / sum(p .* (1-gamma));
        sigma1sqr = sum(p .* gamma .* z.^2) ./ sum(p .* gamma);
        pi1 = sum(p .* gamma)/sum(p);
    end
    
    sigma0 = sqrt(sigma0sqr);
    sigma_delta = sqrt(sigma1sqr - sigma0sqr);
end