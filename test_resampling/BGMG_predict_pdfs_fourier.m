function [zvec Fpdfmat_z Fpdfmat_total Fpdfmat_indep Fpdfmat_pleio pdfmat_z pdfmat_total pdfmat_indep pdfmat_pleio] = BGMG_predict_pdfs_fourier(params,vals_r2,count_r2,delvals,Neff,H,dispflag)

% If Neff or H not specified, assume sigb is in units of del (already scaled)
if ~exist('Neff','var') | isempty(Neff), Neff = 1; end 
if ~exist('H','var') | isempty(H), H = 1; end
if ~exist('dispflag','var') | isempty(dispflag), dispflag = false; end

nsnp = sum(count_r2); nbins = length(vals_r2);

upsampfact = 1;
%upsampfact = 16; % This is needed to handle tiny sig_b * r^2 -- should find more computationally efficient solution
zvec = linspace(delvals(1),-delvals(1),upsampfact*length(delvals)+1); zvec = zvec(1:end-1);
zstep_orig = delvals(2)-delvals(1);
zstep = zvec(2)-zvec(1); Fs = 1/zstep; fvec = -1/2*Fs*zvec/zvec(1);

sigd1 = sqrt(Neff*H)*params.sigb_1;
sigd2 = sqrt(Neff*H)*params.sigb_2;

Fpdfvec1 = 1; Fpdfvec2 = 1;
for bini = 1:nbins
  if count_r2(bini)>0
    tmp1 = exp(-1/2*(2*pi*fvec*sigd1).^2*vals_r2(bini)); 
    tmp2 = exp(-1/2*(2*pi*fvec*sigd2).^2*vals_r2(bini)); 
%    tmp1 = (tmp1-tmp1(1))/(1-tmp1(1));
%    tmp2 = (tmp2-tmp2(1))/(1-tmp2(1));
    Fpdfvec1 = Fpdfvec1.*(params.pi_1*tmp1+(1-params.pi_1)).^count_r2(bini);
    Fpdfvec2 = Fpdfvec2.*(params.pi_2*tmp2+(1-params.pi_2)).^count_r2(bini);
    if 0
    figure(666); clf;
    subplot(2,2,1); plot(fvec,abs(tmp1));
    subplot(2,2,2); plot(fvec,abs(tmp2));
    subplot(2,2,3); plot(fvec,abs(Fpdfvec1));
    subplot(2,2,4); plot(fvec,abs(Fpdfvec2));
    figure(667); clf;
    subplot(2,2,1); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(tmp1)))))));
    subplot(2,2,2); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(tmp2)))))));
    subplot(2,2,3); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(Fpdfvec1)))))));
    subplot(2,2,4); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(Fpdfvec2)))))));
    end
  end
end

Fpdfmat_indep = colvec(Fpdfvec1)*rowvec(Fpdfvec2);

[zmat2 zmat1] = meshgrid(zvec);
[fmat2 fmat1] = meshgrid(fvec);
C = [params.sigb_1^2 params.sigb_1*params.sigb_2*params.rhob; params.sigb_1*params.sigb_2*params.rhob params.sigb_2^2];
tmp = 2*pi*[fmat1(:) fmat2(:)];
dist2list = sum((tmp*C).*tmp,2);
dist2mat = reshape(dist2list,size(zmat1));
Fpdfmat_pleio = 1;
for bini = 1:nbins
  if count_r2(bini)>0
    Fpdfmat_pleio = Fpdfmat_pleio.*(params.pi_3*exp(-1/2*dist2mat*vals_r2(bini))+(1-params.pi_3)).^count_r2(bini);
  end
end

C0 = [params.sig0_1^2 params.sig0_1*params.sig0_2*params.rho0; params.sig0_1*params.sig0_2*params.rho0 params.sig0_2^2];
tmp = 2*pi*[fmat1(:) fmat2(:)];
dist2list = sum((tmp*C0).*tmp,2);
dist2mat = reshape(dist2list,size(zmat1));
Fpdfmat_0 = exp(-1/2*dist2mat);

if dispflag % Show results without apodizing
  Fpdfmat_total = Fpdfmat_indep .* Fpdfmat_pleio;
  pdfmat_indep = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_indep))))); pdfmat_indep = pdfmat_indep/(sum(pdfmat_indep(:))*zstep^2)*(zstep/zstep_orig);
  pdfmat_pleio = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_pleio))))); pdfmat_pleio = pdfmat_pleio/(sum(pdfmat_pleio(:))*zstep^2)*(zstep/zstep_orig);
  pdfmat_total = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_total))))); pdfmat_total = pdfmat_total/(sum(pdfmat_total(:))*zstep^2)*(zstep/zstep_orig);

  figure(666); clf; crange = [-10 5];
  subplot(2,3,1); imagesc(Fpdfmat_indep); colormap(fire); axis equal tight xy;
  subplot(2,3,2); imagesc(Fpdfmat_pleio); colormap(fire); axis equal tight xy;
  subplot(2,3,3); imagesc(Fpdfmat_total); colormap(fire); axis equal tight xy;
  subplot(2,3,4); imagesc(log(pdfmat_indep),crange); colormap(fire); axis equal tight xy;
  subplot(2,3,5); imagesc(log(pdfmat_pleio),crange); colormap(fire); axis equal tight xy;
  subplot(2,3,6); imagesc(log(pdfmat_total),crange); colormap(fire); axis equal tight xy;
end

%apovec = rowvec(rowvec(hamming(length(zvec)+1)).^1); 
apovec = rowvec(chebwin(length(zvec)+1,160)); % Need sufficient sidelobe suppression to avoid "wrap"
apovec = apovec(1:end-1); apomat = colvec(apovec)*rowvec(apovec);
pi0_1 = (1-params.pi_1)^nsnp; pi0_2 = (1-params.pi_2)^nsnp;
Fpdfvec1 = pi0_1+apovec.*(Fpdfvec1-pi0_1); Fpdfvec2 = pi0_2+apovec.*(Fpdfvec2-pi0_2);
Fpdfmat_indep = colvec(Fpdfvec1)*rowvec(Fpdfvec2);
pi0_3 = (1-params.pi_3)^nsnp;
Fpdfmat_pleio = pi0_3+(apomat.*Fpdfmat_pleio-pi0_3);
Fpdfmat_total = Fpdfmat_indep .* Fpdfmat_pleio;
Fpdfmat_z = Fpdfmat_0 .* Fpdfmat_total;

if nargout <= 5, return; end % Compute inverse Fourier transform only if needed

pdfmat_indep = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_indep))))); pdfmat_indep = pdfmat_indep/(sum(pdfmat_indep(:))*zstep^2)*(zstep/zstep_orig);
pdfmat_pleio = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_pleio))))); pdfmat_pleio = pdfmat_pleio/(sum(pdfmat_pleio(:))*zstep^2)*(zstep/zstep_orig);
pdfmat_total = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_total))))); pdfmat_total = pdfmat_total/(sum(pdfmat_total(:))*zstep^2)*(zstep/zstep_orig);
pdfmat_z = max(0,real(ifftshift(ifft2(fftshift(Fpdfmat_z))))); pdfmat_z = pdfmat_z/(sum(pdfmat_z(:))*zstep^2)*(zstep/zstep_orig);

if dispflag % Show results with apodizing
  figure(667); clf; crange = [-10 5];
  subplot(2,3,1); imagesc(Fpdfmat_indep); colormap(fire); axis equal tight xy;
  subplot(2,3,2); imagesc(Fpdfmat_pleio); colormap(fire); axis equal tight xy;
  subplot(2,3,3); imagesc(Fpdfmat_total); colormap(fire); axis equal tight xy;
  subplot(2,3,4); imagesc(log(pdfmat_indep),crange); colormap(fire); axis equal tight xy;
  subplot(2,3,5); imagesc(log(pdfmat_pleio),crange); colormap(fire); axis equal tight xy;
  subplot(2,3,6); imagesc(log(pdfmat_total),crange); colormap(fire); axis equal tight xy;
end

return

% Extras 

Fpdfvec1 = 1; Fpdfvec2 = 1;
for bini = 1:nbins
  if count_r2(bini)>0
    tmp1 = exp(-1/2*(2*pi*fvec*sigd1).^2*vals_r2(bini)); tmp1 = (tmp1-tmp1(1))/(1-tmp1(1));
    tmp2 = exp(-1/2*(2*pi*fvec*sigd2).^2*vals_r2(bini)); tmp2 = (tmp2-tmp2(1))/(1-tmp2(1));
    Fpdfvec1 = Fpdfvec1.*(params.pi_1*tmp1+(1-params.pi_1)).^count_r2(bini);
    Fpdfvec2 = Fpdfvec2.*(params.pi_2*tmp2+(1-params.pi_2)).^count_r2(bini);
    figure(666); clf;
    subplot(2,2,1); plot(fvec,abs(tmp1));
    subplot(2,2,2); plot(fvec,abs(tmp2));
    subplot(2,2,3); plot(fvec,abs(Fpdfvec1));
    subplot(2,2,4); plot(fvec,abs(Fpdfvec2));
    figure(667); clf;
    subplot(2,2,1); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(tmp1)))))));
    subplot(2,2,2); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(tmp2)))))));
    subplot(2,2,3); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(Fpdfvec1)))))));
    subplot(2,2,4); plot(zvec,log(max(0,real(ifftshift(ifft(fftshift(Fpdfvec2)))))));
    pause
  end
end

apovec = rowvec(hamming(length(zvec)));
logFpdfvec1 = 0; logFpdfvec2 = 0;
for bini = 1:nbins
  if count_r2(bini)>0
    logFpdfvec1 = logFpdfvec1 + log(apovec.*(params.pi_1*exp(-1/2*(2*pi*fvec*sigd1).^2*vals_r2(bini))+(1-params.pi_1)))*count_r2(bini);
    logFpdfvec2 = logFpdfvec2 + log(apovec.*(params.pi_2*exp(-1/2*(2*pi*fvec*sigd2).^2*vals_r2(bini))+(1-params.pi_2)))*count_r2(bini);
  end
  Fpdfvec1 = exp(logFpdfvec1); Fpdfvec2 = exp(logFpdfvec2);
end

logFpdfmat_pleio = 0;
for bini = 1:nbins
  if count_r2(bini)>0
    logFpdfmat_pleio = logFpdfmat_pleio + log(params.pi_3*exp(-1/2*dist2mat*vals_r2(bini))+(1-params.pi_3))*count_r2(bini);
  end
end
Fpdfmat_pleio = exp(logFpdfmat_pleio);

