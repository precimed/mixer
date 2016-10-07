function [pvecs_post dvecs_post] = GMM_bivariate_lookup(fitstruct,params,options,zmat,TLDvec,Hvec)

ngm = options.ngm;
edgevec_z = params.edgevec_z; nzbins = length(edgevec_z)-1; meanvec_z = (edgevec_z(1:end-1)+edgevec_z(2:end))/2;
gvlist = fitstruct.gvlist; ngvbins = length(gvlist)-1; ldlist = fitstruct.ldlist; nldbins = length(ldlist)-1;

defvec = isfinite(sum(zmat,2));
ldivec = max(1,min(nldbins,floor(1+interp1(ldlist,[0:nldbins],TLDvec,'linear','extrap'))));
gvivec = max(1,min(ngvbins,floor(1+interp1(gvlist,[0:ngvbins],Hvec,'linear','extrap'))));
pvecs_post = NaN(size(zmat,1),ngm+1); dvecs_post = NaN(size(zmat));
for ldi = 1:nldbins
  for gvi = 1:length(gvlist)-1
    pmat_post = fitstruct.pmats_post{ldi,gvi};
    mu1mat = fitstruct.mu1mats{ldi,gvi};
    mu2mat = fitstruct.mu2mats{ldi,gvi};
    ivec = find((ldivec==ldi)&(gvivec==gvi)&defvec);
    indvec1 = max(1,min(nzbins,1+(zmat(ivec,1)-meanvec_z(1))/(meanvec_z(2)-meanvec_z(1))));
    indvec2 = max(1,min(nzbins,1+(zmat(ivec,2)-meanvec_z(1))/(meanvec_z(2)-meanvec_z(1))));
    dvecs_tmp = cat(2,interp2(mu1mat,indvec1,indvec2,'linear'),interp2(mu2mat,indvec1,indvec2,'linear'));
    pvecs_tmp = zeros(length(ivec),size(pvecs_post,2));
    for gmi = 0:ngm
      pvecs_tmp(:,1+gmi) = interp2(squeeze(pmat_post(1+gmi,:,:)),indvec1,indvec2,'linear'); 
    end 
    pvecs_post(ivec,:) = pvecs_tmp;
    dvecs_post(ivec,:) = dvecs_tmp;
  end
end

