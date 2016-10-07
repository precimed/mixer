function [cost fitstruct] = GMM_bivariate_cost(x,zmat,Hvec,Lvec,imat_prune,options)

x_struct = GMM_bivariate_mapparams(x,options);

v2struct(x_struct);

nsnp = size(zmat,1);
nprune = size(imat_prune,2);

Ledges = [0:100:500];
zmeans = linspace(-15,15,31); zstep = zmeans(2)-zmeans(1); zedges = [zmeans-zstep/2,zmeans(end)+zstep/2];
[tmp1 tmp2] = meshgrid(zmeans,zmeans); zmeangrid = [colvec(tmp1) colvec(tmp2)];

C0 = [sig01^2 sig01*sig02*rho0; sig01*sig02*rho0 sig02^2];
pi0 = 1-sum(pivec);

ngm = length(pivec);

cost = 0;
hcnt_sum = 0;
for prunei = 1:nprune
  tic
  pvec_tmp = pi0*mvnpdf(zmat(imat_prune(:,prunei),:),0,C0);
  for gmi = 1:ngm
    C = [sig1vec(gmi)^2 rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) sig2vec(gmi)^2];
    pvec_tmp = pvec_tmp + pivec(:,gmi).*mvnpdf(zmat(imat_prune(:,prunei),:),0,C0+C);
  end
  toc
  cost = cost + sum(-log(pvec_tmp));
  if nargout > 1
    tic
    valvec1 = cat(1,zmat(imat_prune(:,prunei),1),-zmat(imat_prune(:,prunei),1));
    valvec2 = cat(1,zmat(imat_prune(:,prunei),2),-zmat(imat_prune(:,prunei),2));
    [hcnt] = histc2d(valvec1,valvec2,zedges,zedges); hcnt = hcnt(1:end-1,1:end-1);
    hcnt_sum = hcnt_sum + hcnt;
    toc
  end
end

if nargout==1, return; end

ticks = find(mod(zmeans,5)==0); ticklabels = cellfun(@num2str,num2cell(zmeans(ticks)),'UniformOutput',false);
im = log10(hcnt_sum/max(hcnt_sum(:)));
figure(1); clf; imagesc(im,[-5 0]); colormap(hot); axis equal; axis xy; axis tight; set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels); drawnow;

im_sm = log10(max(0,real(smooth2d(max(0.5,hcnt_sum)/max(hcnt_sum(:)),3,3))));
figure(11); clf; imagesc(im_sm,[-5 0]); colormap(hot); axis equal; axis xy; axis tight; set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels); drawnow;
figure(12); contourf(im_sm); colormap(hot); axis equal; axis xy; axis tight; set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels); drawnow;

im_sm = real(smooth2d(log10(max(0.5,hcnt_sum)/max(hcnt_sum(:))),3,3));
figure(21); clf; imagesc(im_sm,[-5 0]); colormap(hot); axis equal; axis xy; axis tight; set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels); drawnow;
figure(22); contourf(im_sm); colormap(hot); axis equal; axis xy; axis tight; set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels); drawnow;

pdfmat_pred = GMM_bivariate_histpred(x_struct);

fitstruct = struct();

keyboard

return; 

% Not sure why visual fit to data gets worse with fitting

% Not sure why this doesn't work:

%pi0effvec = (pi0).^Lvec;
pi0effvec = repmat(pi0,size(Lvec)); % Fit z-scored directly (ignore L)
cost = 0;
hcnt_sum = 0;
for prunei = 1:nprune
  pi0effvec_tmp = pi0effvec(imat_prune(:,prunei));
  pieffmat = (1-pi0effvec_tmp)*pivec;
  Hmat_tmp = repmat(reshape(Hvec(imat_prune(:,prunei)),[1 1 length(pi0effvec_tmp)]),[2 2 1]); 
  C0mat = repmat(C0,[1 1 length(pi0effvec_tmp)]);
  pvec_tmp = pi0effvec_tmp.*mvnpdf(zmat(imat_prune(:,prunei),:),0,C0);
  for gmi = 1:ngm
    C = [sig1vec(gmi)^2 rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) sig2vec(gmi)^2];
    Cmat = Hmat_tmp.*repmat(C,[1 1 length(pi0effvec_tmp)]);
    pvec_tmp = pvec_tmp + pieffmat(:,gmi).*mvnpdf_amd(zmat(imat_prune(:,prunei),:),0,C0mat+Cmat);
  end
  if nargout > 1 
    valvec1 = cat(1,zmat(imat_prune(:,prunei),1),-zmat(imat_prune(:,prunei),1));
    valvec2 = cat(1,zmat(imat_prune(:,prunei),2),-zmat(imat_prune(:,prunei),2));
    [hcnt] = histc2d(valvec1,valvec2,zedges,zedges); hcnt = hcnt(1:end-1,1:end-1);
    hcnt_sum = hcnt_sum + hcnt;
  end
  cost = cost + sum(-log(pvec_tmp));
end


% ToDo
%   Turn back on modeling of effect of L
%   Compute predicted and observed bivariate histograms of z-scores, overall and by H-bin, and L-bin
%   Compute predicted and observed conditional expectancy of z-scores (from bivariate distribution)
%   Check w. Dominic / Chun re: up-to-date LDmat & TLD (L)
%   Extend to using RES, rather than TLD
%   Exclude MHC from fitting & pruning
%   Restrict range of TLD (L) in fitting
%   Allow for at least two causal SNPs in LD

%   Should try smoooth 2-d kernel / PDF estimator, allowing for scaling by H, and weighting by 1-pi0eff
%     First, based on raw z-scores
%     Then, based incorporating effects of H and L

%   Doesn't seem to fit data well

% Try manual fitting of generalized Gaussian mixture model (don't assume "Venn Diagram" structure):

x_struct2 = struct('sig01',x_struct.sig01,'sig02',x_struct.sig02,'rho0',x_struct.rho0,'pivec',[1e-2],'sig1vec',[1],'sig2vec',[1.0],'rhovec',[0.4]);
x_struct2 = struct('sig01',x_struct.sig01,'sig02',x_struct.sig02,'rho0',x_struct.rho0,'pivec',[1e-2],'sig1vec',[0.5],'sig2vec',[2.0],'rhovec',[0.9]);
x_struct2 = struct('sig01',x_struct.sig01,'sig02',x_struct.sig02,'rho0',x_struct.rho0,'pivec',[1e-2 1e-4],'sig1vec',[0.4 0.5],'sig2vec',[0.6 2.5],'rhovec',[0.4 0.95]);
pdfmat_pred = GMM_bivariate_histpred(x_struct2);

% ToDo
%   Implement data fit (start with fitting of z-scores (using mean H and L), then fit with predicted effects of L and H)
%   Also incorporate Neff1 and Neff2
%   Backup GMM_bivariate_fit, GMM_bivariate_cost, GMM_bivariate_mapparams; make new versions
%   Remove MHC from fit and pruning
%   Also compute & plot empirical and predited conditional expectancy values for z, delta or beta, replication rate)
%   Could implement as stepwise fit (adding Gaussian mixtures ony if fit improves significantly)


