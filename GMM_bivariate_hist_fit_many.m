function gmm = GMM_bivariate_hist_fit_many(chrnumvec, posvec, annomat, mafvec, zmat, LDmat)

    % From http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/region.cgi?name=MHC&asm=GRCh37
    ivec_mhc = ((chrnumvec==6)&(posvec>=28477797)&(posvec<=33448354)); 

    TLDvec = annomat(:,end);
    Hvec = 2*mafvec.*(1-mafvec);
    genvar = Hvec;

    zmax = -norminv(5e-8/2); % Consider only non-GWAS significant SNPs
    %zmax = -norminv(1e-10/2);
    zedgevec = linspace(-zmax,zmax,30);

    nrep = 10;
    gvlist = linspace(0,0.5,6); ngvbins = length(gvlist)-1; 
    ldlist = [0 2 5 20 50 100 200]; nldbins = length(ldlist)-1;
    pruneflag = true;
    excludeMHCflag = false;
    params_fit = struct('edgevec_z',zedgevec,'gvlist',gvlist,'ldlist',ldlist,'pruneflag',pruneflag,'excludeMHCflag',excludeMHCflag);
    sumstats = cell(1,nrep);
    for repi = 1:nrep
      zvec1 = zmat(:,1); % Neff1 = Neffvec_fit(1);
      zvec2 = zmat(:,2); % Neff2 = Neffvec_fit(2);
      ld_gv_freq = NaN(length(ldlist)-1,length(gvlist)-1);
      hcmats = cell(nldbins,ngvbins);
      if excludeMHCflag, zvec1(ivec_mhc) = NaN; zvec2(ivec_mhc) = NaN; end
      if pruneflag, tmp = NaN(size(zvec1)); tmp(isfinite(zvec1+zvec2)) = rand(1,sum(isfinite(zvec1+zvec2))); tmp = GenStats_FastPrune(tmp,LDmat); zvec1(~isfinite(tmp)) = NaN; zvec2(~isfinite(tmp)) = NaN; end
      zmat_tmp = cat(2,zvec1,zvec2);
      defvec = isfinite(sum(zmat_tmp,2));
      for ldi = 1:nldbins
        for gvi = 1:ngvbins
          ivec = find(genvar>=gvlist(gvi)&genvar<=gvlist(gvi+1) & TLDvec>=ldlist(ldi)&TLDvec<ldlist(ldi+1) & defvec);
          gvfreq(gvi) = sum(ivec)/sum(defvec);
          hcmat = hist3(zmat_tmp(ivec,:),'Edges',{zedgevec,zedgevec}); hcmat = hcmat(1:end-1,1:end-1);
          hcmats{ldi,gvi} = hcmat;
        end
      end
      sumstats{repi} = struct('hcmats',{hcmats});
    end

    options_fit = struct();
    options_fit.MaxIter1 = 40;
    options_fit.MaxIter2 = 200;
    
    fprintf('Fit data to model (null distribution)\n');
    options_fit.mixtypevec = [];
    x0_struct = struct('sig01',1.000,'sig02',1.000,'rho0',0.000,'pivec',[],'sig1vec',[],'sig2vec',[],'rhovec',[]);
    gmm.null = fit_one(options_fit, x0_struct, sumstats, params_fit);
    
    % Extract parameters of the null distribution
    sig01 = gmm.null.x_fit_struct.sig01;
    sig02 = gmm.null.x_fit_struct.sig02;
    rho0  = gmm.null.x_fit_struct.rho0;
    sig1  = max(2, 1 + 3*(sig01 - 1));
    sig2  = max(2, 1 + 3*(sig02 - 1));
    rho   = 3 * rho0;
    pi   = 0.01;
    
    fprintf('Fit data to model (Each trait separately)\n');
    options_fit.mixtypevec = [1 2];
    x0_struct = struct('sig01',sig01,'sig02',sig02,'rho0',rho0, 'pivec', [pi pi], ...
                       'sig1vec', [sig1 0], 'sig2vec', [0 sig2],'rhovec', [0 0]);
    gmm.indep = fit_one(options_fit, x0_struct, sumstats, params_fit);
    
    fprintf('Fit data to model (pleiotropic component only)\n');
    options_fit = rmfield(options_fit, 'mixtypevec');
    x0_struct = struct('sig01',sig01,'sig02',sig02,'rho0',rho0, 'pivec', pi, ...
                       'sig1vec', sig1, 'sig2vec', sig2,'rhovec', rho);
    gmm.pleio = fit_one(options_fit, x0_struct, sumstats, params_fit);
    
    fprintf('Fit data to model (full model, including pleiotropic and independent components)\n');
    options_fit.mixtypevec = [1 2 3];
    x0_struct = struct('sig01',sig01,'sig02',sig02,'rho0',rho0, 'pivec', [pi pi pi], ...
                       'sig1vec', [sig1 0 sig1], 'sig2vec', [0 sig2 sig2],'rhovec', [0 0 rho]);
    gmm.full = fit_one(options_fit, x0_struct, sumstats, params_fit);

    fprintf('Fit data to model (independent components for the first trait, plust pleiotropic component)\n');
    options_fit.mixtypevec = [1 3];
    x0_struct = struct('sig01',sig01,'sig02',sig02,'rho0',rho0, 'pivec', [pi pi], ...
                       'sig1vec', [sig1 sig1], 'sig2vec', [0 sig2],'rhovec', [0 rho]);
    gmm.first = fit_one(options_fit, x0_struct, sumstats, params_fit);

    fprintf('Fit data to model (independent components for the second trait, plust pleiotropic component)\n');
    options_fit.mixtypevec = [2 3];
    x0_struct = struct('sig01',sig01,'sig02',sig02,'rho0',rho0, 'pivec', [pi pi], ...
                       'sig1vec', [0 sig1], 'sig2vec', [sig2 sig2],'rhovec', [0 rho]);
    gmm.second = fit_one(options_fit, x0_struct, sumstats, params_fit);
    
    fn = fieldnames(gmm);
    x = []; for i = 1:length(fn), x = [x gmm.(fn{i}).AIC]; end; 
    y = []; for i = 1:length(fn), y = [y gmm.(fn{i}).cost]; end; 
    for i = 1:length(fn)
        f = gmm.(fn{i});
        fprintf('%s cost=%.2f AIC=%.2f %s\n', fn{i}, f.cost-min(y), f.AIC-min(x), GMM_bivariate_show_struct(f.x_fit_struct));
    end
end

function result = fit_one(options_fit, x0_struct, sumstats, params_fit)
    clear GMM_bivariate_histcost
    options_fit.ngm = length(x0_struct.pivec);
    x0 = GMM_bivariate_mapparams(x0_struct,options_fit);
    resultStruct1 = GMM_bivariate_histfit(sumstats,params_fit,options_fit,x0); 
    result.cost = GMM_bivariate_histcost(resultStruct1.x_fit,sumstats,params_fit,options_fit);
    result.AIC = 2*length(resultStruct1.x_fit) + 2*result.cost;
    result.x_fit_struct= GMM_bivariate_mapparams(resultStruct1.x_fit,options_fit); 
end
