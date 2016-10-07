function [cost fitstruct pvecs_post dvecs_post] = GMM_bivariate_histcost(x,sumstats,params,options,plotflag,zmat,TLDvec,Hvec)

if ~exist('plotflag','var'), plotflag = false; end

ngm = options.ngm;

x_struct = GMM_bivariate_mapparams(x,options);

v2struct(x_struct);

nrep = size(sumstats,2); ni = nrep;

sig0_d = 1e-6; % Probability mass at zero (should make code work for sig0_d = 0)

C0 = [sig01^2 sig01*sig02*rho0; sig01*sig02*rho0 sig02^2];
C0_d = [sig0_d^2 0; 0 sig0_d^2];
pi0 = 1-sum(pivec);

sig1vec = max(sig0_d,sig1vec); sig2vec = max(sig0_d,sig2vec);  % Avoid singular matrices -- should fix code to handle single-trait components

%Neff1 = 0.2; Neff2 = 0.2; % Half samples
Neff1 = 1; Neff2 = 1; % Full sample sizes

edgevec_z = params.edgevec_z; nzbins = length(edgevec_z)-1; meanvec_z = (edgevec_z(1:end-1)+edgevec_z(2:end))/2; 
gvlist = params.gvlist; ngvbins = length(gvlist)-1; ldlist = params.ldlist; nldbins = length(ldlist)-1;
cost = 0;
for i = 1:ni
  hcmats = sumstats{i}.hcmats;
  [Z2 Z1] = meshgrid(edgevec_z,edgevec_z); zmat_edges = [Z1(:) Z2(:)]; [Z2 Z1] = meshgrid(meanvec_z,meanvec_z); zmat_means = [Z1(:) Z2(:)];
  if nargout>1, pmats_z = cell(nldbins,ngvbins); pmats_d = cell(nldbins,ngvbins); Cs_z = cell(nldbins,ngvbins); Cs_d = cell(nldbins,ngvbins); mu1mats = cell(nldbins,ngvbins); mu2mats = cell(nldbins,ngvbins); end
  for ldi = 1:nldbins
    for gvi = 1:ngvbins
      hcmat = hcmats{ldi,gvi};
      gv = (gvlist(gvi)+gvlist(gvi+1))/2; % Should probably take mean over actual gv-values
      ld = (ldlist(ldi)+ldlist(ldi+1))/2; % Should probably take mean over actual ld-values
      cmat_z = cell(1,ngm+1); pmat_z = cell(1,ngm+1);
      cmat_d = cell(1,ngm+1); pmat_d = cell(1,ngm+1);
      pmat_post = NaN(ngm+1,nzbins,nzbins);
      for gmi = 0:ngm, cmat_z{1+gmi} = NaN(size(hcmat)+1); cmat_d{1+gmi} = NaN(size(hcmat)+1); end
      pi0_eff = pi0.^ld; pi1_eff = 1-pi0_eff;
      pivec_eff = [pi0_eff pivec/sum(pivec)*pi1_eff];
      cmat_z{1}(:) = pi0_eff*BVNcdf(zmat_edges,[0 0],C0);
      if nargout>1, cmat_d{1}(:) = pivec_eff(1)*BVNcdf(zmat_edges,[0 0],C0_d); Cs_z{1} = C0; Cs_d{1} = C0_d; end
      for gmi = 1:ngm % Need to incorporate effects of Neff1 and Neff2
        if pi0>0 % Mixture?
          C = gv*[Neff1*sig1vec(gmi)^2 sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) Neff2*sig2vec(gmi)^2];
          cmat_z{1+gmi}(:) = pivec_eff(1+gmi)*BVNcdf(zmat_edges,[0 0],C0+C); 
          if nargout>1, cmat_d{1+gmi}(:) = pivec_eff(1+gmi)*BVNcdf(zmat_edges,[0 0],C); end
        else
          C = ld*gv*[Neff1*sig1vec(gmi)^2 sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); sqrt(Neff1*Neff2)*rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) Neff2*sig2vec(gmi)^2];
          cmat_z{1+gmi}(:) = BVNcdf(zmat_edges,[0 0],C0+C); 
          if nargout>1, cmat_d{1+gmi}(:) = BVNcdf(zmat_edges,[0 0],C); end
        end
        if nargout>1
          Cs_z{1+gmi} = C0+C;
          Cs_d{1+gmi} = C;
        end
      end
      pmat_z_sum = 0;
      for gmi = 0:ngm
        pmat_z{1+gmi} = cm3_cdf2pdf(cmat_z{1+gmi});
        pmat_z_sum = pmat_z_sum + pmat_z{1+gmi};
        if nargout>1, pmat_d{1+gmi} = cm3_cdf2pdf(cmat_d{1+gmi}); end
      end
      cost = cost + -mnpdfln(rowvec(hcmat),rowvec(pmat_z_sum/sum(pmat_z_sum(:))))/ni; % Adjust sum of log-likelihood by number of independent prunings -- why is the cost not scaled by ni??
      if nargout>1 
        tic
        pmat_d_post = 0; mu1mat = 0; mu2mat = 0;
        for gmi = 0:ngm
          pmat_post(1+gmi,:,:) = pmat_z{1+gmi}./pmat_z_sum;
          Cnew = inv(inv(C0)+inv(Cs_d{1+gmi}));
          pmat_d_tmp = NaN(nzbins,nzbins,nzbins,nzbins); mu1mat_tmp = NaN(nzbins,nzbins); mu2mat_tmp = NaN(nzbins,nzbins);
          for i1 = 1:nzbins
            for i2 = 1:nzbins
              tmp = reshape(mvnpdf(zmat_means,[0 0],Cs_d{1+gmi}),[nzbins nzbins]); pmat_d_tmp(i1,i2,:,:) = reshape(mvnpdf(zmat_means,[meanvec_z(i1) meanvec_z(i2)],C0),[nzbins nzbins]).*(tmp/sum(tmp(:)));
              munew = (Cnew*[meanvec_z(i1) meanvec_z(i2)]')'; % pmat2_d_tmp(i1,i2,:,:) = reshape(mvnpdf(zmat_means,munew,Cnew),[nzbins nzbins]);
              mu1mat_tmp(i1,i2) = munew(1); mu2mat_tmp(i1,i2) = munew(2);
              if ismember(i1,[1 15 29]) & ismember(i2,[1 15 29]) & 0
                [meanvec_z(i2) meanvec_z(i1)]
                im = squeeze(pmat_d_tmp(i1,i2,:,:)); figure(1001); imagesc(log10(max(eps,im/max(im(:)))),[-5 0]); colormap(hot); axis xy; axis equal; axis tight; set(gca,'XTickLabel',{},'YTickLabel',{});
%                im2 = squeeze(pmat2_d_tmp(i1,i2,:,:)); figure(1002); imagesc(log10(max(eps,im2/max(im2(:)))),[-5 0]); colormap(hot); axis xy; axis equal; axis tight; set(gca,'XTickLabel',{},'YTickLabel',{});
                drawnow; pause;
              end
            end
          end
          pmat_d_post = pmat_d_post + pivec_eff(1+gmi)*pmat_d_tmp; 
          mu1mat = mu1mat + squeeze(pmat_post(1+gmi,:,:)).*mu1mat_tmp; mu2mat = mu2mat + squeeze(pmat_post(1+gmi,:,:)).*mu2mat_tmp;
        end
        pmats_z{ldi,gvi} = pmat_z; pmats_z_sum{ldi,gvi} = pmat_z_sum; pmats_d{ldi,gvi} = pmat_d; pmats_post{ldi,gvi} = pmat_post; pmats_d_post{ldi,gvi} = pmat_d_post;
        mu1mats{ldi,gvi} = mu1mat; mu2mats{ldi,gvi} = mu2mat; 
        toc
      end
    end
  end
end

if nargout==1, return; end

fitstruct = struct('Neff1',Neff1,'Neff2',Neff2,'gvlist',gvlist,'ldlist',ldlist,'edgevec_z',edgevec_z,'meanvec_z',meanvec_z,'pmats_z_sum',{pmats_z_sum},'pmats_z',{pmats_z},'mu1mats',{mu1mats},'mu2mats',{mu2mats},'pmats_post',{pmats_post},'pmats_d_post',{pmats_d_post});


% Plot observed & predicted histograms and means (should do this in separate function)
if plotflag
  muvec1_z = NaN(1,nzbins); muvec2_z = NaN(1,nzbins); muvec1_z_pred = NaN(1,nzbins); muvec2_z_pred = NaN(1,nzbins); muvec1_d_pred = NaN(1,nzbins); muvec2_d_pred = NaN(1,nzbins);
  figure(11); clf; title('Observed'); figure(12); clf; title('Predicted');
  figure(21); clf; figure(22); clf; % figure(31); clf; figure(32); clf;
  cnt = 0;
  for ldi = 1:nldbins
    for gvi = 1:length(gvlist)-1
      pmat_z = pmats_z_sum{ldi,gvi};
      pmat_d_post = pmats_d_post{ldi,gvi};
      hcmat = 0;
      for repi = 1:nrep
        hcmat = hcmat + sumstats{repi}.hcmats{ldi,gvi}/nrep;
      end
      cnt = cnt + 1;
      figure(11); subplot(nldbins,ngvbins,cnt); imagesc(log10(max(eps,hcmat/max(hcmat(:)))),[-5 0]); colormap(hot); axis xy; axis equal; axis tight; set(gca,'XTickLabel',{},'YTickLabel',{});
      figure(12); subplot(nldbins,ngvbins,cnt); imagesc(log10(max(eps,pmat_z/max(pmat_z(:)))),[-5 0]); colormap(hot); axis xy; axis equal; axis tight; set(gca,'XTickLabel',{},'YTickLabel',{});
      for j = 1:nzbins
        pvec1_z = rowvec(pmat_z(:,j)); pvec2_z = pmat_z(j,:); muvec1_z_pred(j) = sum(pvec1_z.*meanvec_z)/max(eps,sum(pvec1_z)); muvec2_z_pred(j) = sum(pvec2_z.*meanvec_z)/max(eps,sum(pvec2_z));
        pvec1_z = rowvec(hcmat(:,j)); pvec2_z = hcmat(j,:); muvec1_z(j) = condval(sum(pvec1_z)>100,sum(pvec1_z.*meanvec_z)/sum(pvec1_z),NaN); muvec2_z(j) = condval(sum(pvec2_z)>100,sum(pvec2_z.*meanvec_z)/sum(pvec2_z),NaN);
        pvec1_d = rowvec(sum(sum(pmat_d_post(:,j,:,:),4),1)); muvec1_d_pred(j) = sum(pvec1_d.*meanvec_z)/max(eps,sum(pvec1_d));
        pvec2_d = rowvec(sum(sum(pmat_d_post(j,:,:,:),3),2)); muvec2_d_pred(j) = sum(pvec2_d.*meanvec_z)/max(eps,sum(pvec2_d));
      end
      figure(21); subplot(nldbins,ngvbins,cnt); plot(meanvec_z,muvec1_z,meanvec_z,muvec1_z_pred,'LineWidth',3); ylim([-5 5]);
      figure(22); subplot(nldbins,ngvbins,cnt); plot(meanvec_z,muvec2_z,meanvec_z,muvec2_z_pred,'LineWidth',3); ylim([-5 5]);
%      figure(31); subplot(nldbins,ngvbins,cnt); plot(meanvec_z,muvec1_d_pred,meanvec_z,muvec1_z_pred,'LineWidth',3); ylim([-5 5]);
%      figure(32); subplot(nldbins,ngvbins,cnt); plot(meanvec_z,muvec2_d_pred,meanvec_z,muvec2_z_pred,'LineWidth',3); ylim([-5 5]);
      drawnow;
    end
  end
end


% Compute Bayesian p-values and conditional effect size expectancy values 

if exist('zmat','var')
  [pvecs_post dvecs_post] = GMM_bivariate_lookup(fitstruct,params,options,zmat,TLDvec,Hvec);
end

return

% ToDo
%   Separate out computation of lookup tables based on modified Neff1 and Neff2

% Make stratified z2z plots

pthresh = 0.05;
zvals_disc = linspace(-10.0,10.0,201);

dflist = [2 3 4 5]/10; ndf = length(dflist); nr = round(sqrt(ndf)); nc = ceil(ndf/nr);
pruneflag = false;
nrep = 100;
datastructmat = cell(ndf,nrep);
hist_logpvals = linspace(0,50,501);
figure(101); clf; figure(102); clf; figure(103); clf;
z2meanmat = NaN(ndf,nrep); lamgcmat = NaN(ndf,nrep); Neffmat = NaN(ndf,nrep);
time0 = now;
for dfi = 1:ndf
  dfrac = dflist(dfi);
  dmat = false(nrep,nsamp); rmat = false(nrep,nsamp);
  for repi = 1:nrep
    done = false;
    while ~done
      ivec = randperm(nsamp);
      mi = floor(dfrac*nsamp+rand);
      dmat(repi,:) = false; dmat(repi,ivec(1:mi)) = true; rmat(repi,:) = ~dmat(repi,:);
      done = abs((sum(Neffvec(dmat(repi,:)))/sum(Neffvec))-dfrac) < 0.10;
    end
  end
  zmat_repl_sum = 0; z2mat_repl_sum = 0; repratemat_sum = 0; cntmat_sum = 0;
  time1 = now;
  for repi = 1:nrep
    Neffvec_disc = rowvec(double(Neffvec(dmat(repi,:)))); Neffvec_repl = rowvec(double(Neffvec(rmat(repi,:))));
    Neff_disc = sum(Neffvec_disc); Neff_repl = sum(Neffvec_repl);
    zvec_disc = sum(zmat(:,dmat(repi,:)).*repmat(sqrt(Neffvec_disc),[size(zmat,1) 1]),2)/sqrt(sum(Neffvec_disc));
    zvec_repl = sum(zmat(:,rmat(repi,:)).*repmat(sqrt(Neffvec_repl),[size(zmat,1) 1]),2)/sqrt(sum(Neffvec_repl));
    z2vec_repl = zvec_repl.^2;
    pvec_disc = 2*normcdf(-abs(zvec_disc),0,1); logphistcount = hist(-log10(pvec_disc),hist_logpvals);
    pvec_repl = normcdf(-sign(zvec_disc).*zvec_repl,0,1); repvec = double(pvec_repl<pthresh);
    if pruneflag, tmp = NaN(size(zvec_disc)); tmp(isfinite(zvec_disc+zvec_repl)) = rand(1,sum(isfinite(zvec_disc+zvec_repl))); tmp = GenStats_FastPrune(tmp,LDmat); zvec_disc(~isfinite(tmp)) = NaN; prunevec = isfinite(zvec_disc); else prunevec = []; end
    pvec_disc = 2*normcdf(-abs(zvec_disc),0,1); logphistcount_pruned = hist(-log10(pvec_disc),hist_logpvals);
    for bini = 1:size(iimat,2)
      defvec = isfinite(zvec_disc+zvec_repl) & iimat(:,bini);
      [zvals_repl zvals_repl_cnt] = histmean([-zvec_disc(defvec); zvec_disc(defvec)],[-zvec_repl(defvec); zvec_repl(defvec)],zvals_disc);
      [z2vals_repl z2vals_repl_cnt] = histmean([-zvec_disc(defvec); zvec_disc(defvec)],[z2vec_repl(defvec); z2vec_repl(defvec)],zvals_disc);
      [repvals repvals_cnt] = histmean([-zvec_disc(defvec); zvec_disc(defvec)],[repvec(defvec); repvec(defvec)],zvals_disc);
      zmat_repl(:,bini) = zvals_repl;
      z2mat_repl(:,bini) = z2vals_repl;
      repratemat(:,bini) = repvals;
      cntmat(:,bini) = zvals_repl_cnt;
    end
    zmat_repl_sum = zmat_repl_sum + zmat_repl.*cntmat;
    z2mat_repl_sum = z2mat_repl_sum + z2mat_repl.*cntmat;
    repratemat_sum = repratemat_sum + repratemat.*cntmat;
    cntmat_sum = cntmat_sum + cntmat;
    datastructmat{dfi,repi} = struct('Neff_disc',Neff_disc,'Neff_repl',Neff_repl,'zmat_repl',zmat_repl,'cntmat',cntmat,'z2mat_repl',z2mat_repl,'repratemat',repratemat);
    figure(101); subplot(nr,nc,dfi); plot(zvals_disc,zmat_repl_sum(:,1:nbins)./cntmat_sum(:,1:nbins),zvals_disc,zmat_repl_sum(:,end)./cntmat_sum(:,end),'k-','LineWidt',2);
    h=xlabel('z_d_i_s_c'); set(h,'FontSize',22); h=ylabel('E(z_r_e_p_l | z_d_i_s_c)'); set(h,'FontSize',22); xlim([-6 6]); %ylim([-2 2]);
    figure(102); subplot(nr,nc,dfi); plot(zvals_disc,z2mat_repl_sum(:,1:nbins)./cntmat_sum(:,1:nbins),zvals_disc,z2mat_repl_sum(:,end)./cntmat_sum(:,end),'k-','LineWidt',2);
    h=xlabel('z_d_i_s_c'); set(h,'FontSize',22); h=ylabel('E(z_r_e_p_l^2 | z_d_i_s_c)'); set(h,'FontSize',22); xlim([-6 6]); %ylim([0 6]);
    figure(103); subplot(nr,nc,dfi); plot(zvals_disc,repratemat_sum(:,1:nbins)./cntmat_sum(:,1:nbins),zvals_disc,repratemat_sum(:,end)./cntmat_sum(:,end),'k-','LineWidt',2);
    h=xlabel('z_d_i_s_c'); set(h,'FontSize',22); h=ylabel('Expected Replication Rate'); set(h,'FontSize',22); ylim([0 1.1]); xlim([-6 6]);
    drawnow;
    fprintf(1,'rep %d of %d (now:%s done:%s)\n',repi,nrep,datestr(now,'HH:MM:SS'),datestr(time1+(now-time1)/repi*nrep,'HH:MM:SS'));
  end
  fprintf(1,'dfi %d of %d (now:%s done:%s)\n',dfi,ndf,datestr(now,'HH:MM:SS'),datestr(time0+(now-time0)/dfi*ndf,'HH:MM:SS'));
end

pi1 = 0.0512; sigd = 0.0248; sige = 0.0078; sig0 = 1.0078; % Best fit to z and z^2 replication SCZ
x0 = [logit(pi1,0) log(sigd) log(sige) log(sig0)]; % Full model
options_fminsearch = statset('MaxIter',10000,'MaxFunEvals',10000,'Display','iter','TolFun',1e-2,'TolTypeFun','abs','TolX',1e-3);

% Fit model (should also try fitting across bins, with roughness penalty or linear model for each param)
xfits = cell(1,size(iimat,2));
for bini = 1:size(iimat,2)
  costfun = @(x) GenStats_z2z_cost3(x,datastructmat(:),zvals_disc,genvar,[bini]);
  [xfit cost] = fminsearch(costfun,x0,options_fminsearch); x0 = xfit;
  xfits{bini} = xfit;
end


