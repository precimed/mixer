function GMM_ofrei_bivariate_plot_callback(options)
    %GMM_ofrei_univariate_plot_callback(options);
    %if isfield('options', 'axes2'), axes(options.axes2); end;
    subplot(2,2,2);
    
    xlimit = 10; step = 0.2;
    p  = (-xlimit:step:xlimit)';[z1,z2] = meshgrid(p, p);
    m_zvec = [z1(:), z2(:)];
    
    nsnp = size(m_zvec, 1);
    m_hvec = ones(nsnp, 1) * 0.25;
    m_nvec = ones(nsnp, 2) * options.N;
    m_ref_ld = ones(nsnp, 1) * options.tld;
    m_w_ld = ones(nsnp, 1) * options.tld;
    mapparams = @GMM_ofrei_bivariate_mapparams;

    options.x_step = 0.5;
    options.x_limit = 20;
    
    [~, zprobvec] = GMM_ofrei_bivariate_cost(mapparams(options.opts_struct2), m_zvec, m_hvec, m_nvec, m_ref_ld, m_w_ld, mapparams, options);
    zprobvec=reshape(zprobvec, size(z1));
    
    %if options.use_logscale,
    %    imagesc(log(zprobvec));colorbar;
    %else
        imagesc(zprobvec);colorbar;
    %end

    diag_elems  = diag(zprobvec); diag_elems = diag_elems / sum(diag_elems); 
    adiag_elems = diag(fliplr(zprobvec)); adiag_elems = adiag_elems / sum(adiag_elems);
    if options.use_logscale, diag_elems=log(diag_elems); adiag_elems = log(adiag_elems); end;
    subplot(2,2,1);
    xdiag = [-xlimit : 2 * xlimit / (length(diag_elems) - 1) : xlimit];
    plot(xdiag, diag_elems, xdiag, adiag_elems);
    if options.use_logscale, ylim([-20, 0]); xlim([-xlimit, xlimit]); end;
    

    %title(sprintf('TLD=%i, pi1=%.2e, sbt=%.2e, sig0=%.2f', tld, s.pivec, s.sigma_beta, s.sigma0));
end