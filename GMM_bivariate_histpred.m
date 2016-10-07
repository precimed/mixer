function pdfmat_pred = GMM_bivariate_histpred(x_struct)

nzbins = 31;
zmeans = linspace(-15,15,nzbins); zstep = zmeans(2)-zmeans(1); zedges = [zmeans-zstep/2,zmeans(end)+zstep/2];
[tmp1 tmp2] = meshgrid(zmeans,zmeans); zmeangrid = [colvec(tmp2) colvec(tmp1)];

v2struct(x_struct);

C0 = [sig01^2 sig01*sig02*rho0; sig01*sig02*rho0 sig02^2];
pi0 = 1-sum(pivec);

ngm = length(pivec);
%pdfmat_pred = zeros(length(zmeans));
pdfmat_pred = reshape(pi0*mvnpdf(zmeangrid,0,C0),[nzbins nzbins]);
for gmi = 1:ngm
  C = [sig1vec(gmi)^2 rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi); rhovec(gmi)*sig1vec(gmi)*sig2vec(gmi) sig2vec(gmi)^2];
  pdfmat_pred = pdfmat_pred + reshape(pivec(gmi)*mvnpdf(zmeangrid,0,C0+C),size(pdfmat_pred)); 
end

ticks = find(mod(zmeans,5)==0); ticklabels = cellfun(@num2str,num2cell(zmeans(ticks)),'UniformOutput',false);
figure(2); clf; imagesc(log10(pdfmat_pred/max(pdfmat_pred(:))),[-5 0]); colormap(hot); axis equal; axis xy; axis tight; set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels); drawnow;

