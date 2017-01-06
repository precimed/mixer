function s = GMM_bivariate_show_struct(x_fit_struct)
    v2struct(x_fit_struct);

    s = sprintf('sig01=%.4f, sig02=%.4f, rho0=%.4f', sig01, sig02, rho0);
    if ~isempty(pivec), s = sprintf('%s, pivec=[%s]', s, sprintf('%.4f ', pivec)); end;
    if ~isempty(sig1vec), s = sprintf('%s, sig1vec=[%s]', s, sprintf('%.4f ', sig1vec)); end;
    if ~isempty(sig2vec), s = sprintf('%s, sig2vec=[%s]', s, sprintf('%.4f ', sig2vec)); end;
    if ~isempty(rhovec), s = sprintf('%s, rhovec=[%s]', s, sprintf('%.4f ', rhovec)); end;
end
