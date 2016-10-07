function [fitstruct] = GMM_bivariate_maketables(x,params,options,Neff1,Neff2)

ngm = options.ngm;

x_struct = GMM_bivariate_mapparams(x,options);

v2struct(x_struct);

sig0_d = 1e-6; % Probability mass at zero (should make code work for sig0_d = 0)

C0 = [sig01^2 sig01*sig02*rho0; sig01*sig02*rho0 sig02^2];
C0_d = [sig0_d^2 0; 0 sig0_d^2];
pi0 = 1-sum(pivec);

edgevec_z = params.edgevec_z; nzbins = length(edgevec_z)-1; meanvec_z = (edgevec_z(1:end-1)+edgevec_z(2:end))/2; 
gvlist = params.gvlist; ngvbins = length(gvlist)-1; ldlist = params.ldlist; nldbins = length(ldlist)-1;
[Z2 Z1] = meshgrid(edgevec_z,edgevec_z); zmat_edges = [Z1(:) Z2(:)]; [Z2 Z1] = meshgrid(meanvec_z,meanvec_z); zmat_means = [Z1(:) Z2(:)];
pmats_z = cell(nldbins,ngvbins); pmats_d = cell(nldbins,ngvbins); Cs_z = cell(nldbins,ngvbins); Cs_d = cell(nldbins,ngvbins); mu1mats = cell(nldbins,ngvbins); mu2mats = cell(nldbins,ngvbins);
for ldi = 1:nldbins
  for gvi = 1:ngvbins
    gv = (gvlist(gvi)+gvlist(gvi+1))/2; % Should probably take mean over actual gv-values
    ld = (ldlist(ldi)+ldlist(ldi+1))/2; % Should probably take mean over actual ld-values
    cmat_z = cell(1,ngm+1); pmat_z = cell(1,ngm+1);
    cmat_d = cell(1,ngm+1); pmat_d = cell(1,ngm+1);
    pmat_post = NaN(ngm+1,nzbins,nzbins);
    for gmi = 0:ngm, cmat_z{1+gmi} = NaN(nzbins+1); cmat_d{1+gmi} = NaN(nzbins+1); end
    pi0_eff = pi0.^ld; pi1_eff = 1-pi0_eff;
    pivec_eff = [pi0_eff pivec/sum(pivec)*pi1_eff];
    cmat_z{1}(:) = pi0_eff*BVNcdf(zmat_edges,[0 0],C0);
    cmat_d{1}(:) = pivec_eff(1)*BVNcdf(zmat_edges,[0 0],C0_d); Cs_z{1} = C0; Cs_d{1} = C0_d;
    for gmi = 1:ngm % Need to incorporate effects of Neff1 and Neff2
      if pi0>0 % Mixture?
        C = gv*[Neff1*sig1vec(gmi)^2 sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) Neff2*sig2vec(gmi)^2];
        cmat_z{1+gmi}(:) = pivec_eff(1+gmi)*BVNcdf(zmat_edges,[0 0],C0+C); 
        cmat_d{1+gmi}(:) = pivec_eff(1+gmi)*BVNcdf(zmat_edges,[0 0],C);
      else
        C = ld*gv*[Neff1*sig1vec(gmi)^2 sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) Neff2*sig2vec(gmi)^2];
        cmat_z{1+gmi}(:) = BVNcdf(zmat_edges,[0 0],C0+C); 
        cmat_d{1+gmi}(:) = BVNcdf(zmat_edges,[0 0],C);
      end
      Cs_z{1+gmi} = C0+C;
      Cs_d{1+gmi} = C;
    end
    pmat_z_sum = 0;
    for gmi = 0:ngm
      pmat_z{1+gmi} = max(eps,cm3_cdf2pdf(cmat_z{1+gmi}));
      pmat_z_sum = pmat_z_sum + pmat_z{1+gmi};
      pmat_d{1+gmi} = max(eps,cm3_cdf2pdf(cmat_d{1+gmi}));
    end
    pmat_d_post = 0; mu1mat = 0; mu2mat = 0;
    for gmi = 0:ngm
      pmat_post(1+gmi,:,:) = pmat_z{1+gmi}./pmat_z_sum;

if min(colvec(pmat_post(1+gmi,:,:)))<0, keyboard; end

      Cnew = inv(inv(C0)+inv(Cs_d{1+gmi}));
      pmat_d_tmp = NaN(nzbins,nzbins,nzbins,nzbins); mu1mat_tmp = NaN(nzbins,nzbins); mu2mat_tmp = NaN(nzbins,nzbins);
      for i1 = 1:nzbins
        for i2 = 1:nzbins
          tmp = reshape(mvnpdf(zmat_means,[0 0],Cs_d{1+gmi}),[nzbins nzbins]); pmat_d_tmp(i1,i2,:,:) = reshape(mvnpdf(zmat_means,[meanvec_z(i1) meanvec_z(i2)],C0),[nzbins nzbins]).*(tmp/sum(tmp(:)));
          munew = (Cnew*[meanvec_z(i1) meanvec_z(i2)]')'; % pmat2_d_tmp(i1,i2,:,:) = reshape(mvnpdf(zmat_means,munew,Cnew),[nzbins nzbins]);
          mu1mat_tmp(i1,i2) = munew(1); mu2mat_tmp(i1,i2) = munew(2);
        end
      end
      pmat_d_post = pmat_d_post + pivec_eff(1+gmi)*pmat_d_tmp; 
      mu1mat = mu1mat + squeeze(pmat_post(1+gmi,:,:)).*mu1mat_tmp; mu2mat = mu2mat + squeeze(pmat_post(1+gmi,:,:)).*mu2mat_tmp;
    end
    pmats_z{ldi,gvi} = pmat_z; pmats_z_sum{ldi,gvi} = pmat_z_sum; pmats_d{ldi,gvi} = pmat_d; pmats_post{ldi,gvi} = pmat_post; pmats_d_post{ldi,gvi} = pmat_d_post;
    mu1mats{ldi,gvi} = mu1mat; mu2mats{ldi,gvi} = mu2mat; 
  end
end

fitstruct = struct('Neff1',Neff1,'Neff2',Neff2,'gvlist',gvlist,'ldlist',ldlist,'edgevec_z',edgevec_z,'meanvec_z',meanvec_z,'pmats_z_sum',{pmats_z_sum},'pmats_z',{pmats_z},'mu1mats',{mu1mats},'mu2mats',{mu2mats},'pmats_post',{pmats_post},'pmats_d_post',{pmats_d_post});

