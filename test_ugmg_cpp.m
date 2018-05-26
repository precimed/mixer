%clear
loadlibrary('H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll', 'H:\GitHub\BGMG\src\bgmg_matlab.h');
%loadlibrary('H:\GitHub\BGMG\src\build_win\bin\Debug\bgmg.dll', 'H:\GitHub\BGMG\src\bgmg_matlab.h');

libfunctions('bgmg')
unloadlibrary('bgmg');

c1 = load('H:\work\hapgen_ldmat2_plink\bfile_merged_10K_ldmat_p01_SNPwind50k_chr1.ld.mat');
c1.index_A = c1.index_A + 1;
c1.index_B = c1.index_B + 1; % 0-based, comming from python
fprintf('%.i, %.i\n', min(c1.index_A), min(c1.index_B));

num_snp = max(c1.index_B);

dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\simu_h2=0.7_pi1u=1e-05_rep=1.trait1.mat');zvec = dat.zvec;
dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\simu_h2=0.1_pi1u=1e-05_rep=1.trait1.mat');zvec = dat.zvec;
dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\simu_h2=0.1_pi1u=0.0001_rep=1.trait1.mat');zvec = dat.zvec;
dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\simu_h2=0.7_pi1u=0.0001_rep=1.trait1.mat');zvec = dat.zvec;
dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\simu_h2=0.4_pi1u=0.001_rep=1.trait1.mat');zvec = dat.zvec;
dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\simu_h2=0.1_pi1u=0.01_rep=1.trait1.mat');zvec = dat.zvec;
%dat = load('H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'); zvec = dat.zvec(1:num_snp);
nvec = ones(size(zvec)) * 10000;

%ref = load('H:\Dropbox\shared\BGMG\1kG_phase3_EUR_11015883_reference_p01_SNPwind50k_per_allele_4bins.mat', 'w_ld', 'chrnumvec', 'posvec', 'mafvec');

ref = load('H:\Dropbox\shared\BGMG\HAPGEN_EUR_10K_chr1_gcta_reference_10Ksubj_50k_snps_biased_minr2bin2.mat');
weights = ref.w_ld; %(1:num_snp); 
hvec = ref.hvec; %2 * ref.mafvec(1:num_snp) .* ref.mafvec(1:num_snp);
 
def1= load('H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'); def1=def1.defvec; def1 = def1(1:num_snp);
def2= load('H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'); def2=def2.defvec; def2 = def2(1:num_snp);
defvec = def1 & def2 & isfinite(zvec) & (weights>=1); clear('def1', 'def2');
tag_indices = find(defvec);

trait=1;context=1;
kmax = 1;




m2c = @(x)(x-1); % convert matlab to c indices
check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
calllib('bgmg', 'bgmg_set_tag_indices', 0, length(defvec), sum(defvec), m2c(find(defvec)));  check();
calllib('bgmg', 'bgmg_set_zvec', 0, trait, sum(defvec), zvec(defvec));  check();
calllib('bgmg', 'bgmg_set_nvec', 0, trait, sum(defvec), nvec(defvec));  check();
calllib('bgmg', 'bgmg_set_weights', 0, sum(defvec), 1 ./ weights(defvec));  check();
calllib('bgmg', 'bgmg_set_option', 0,  'r2min', 0.01); check();
calllib('bgmg', 'bgmg_set_option', 0,  'kmax', kmax); check();
calllib('bgmg', 'bgmg_set_option', 0,  'max_causals', 10000); check();

calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, length(c1.r2), m2c(c1.index_A), m2c(c1.index_B), c1.r2);  check();
calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();
%calllib('bgmg', 'bgmg_set_hvec', 0, length(defvec), hvec);  check();

if 0
figure(2)
num_snp = length(zvec); 
pBuffer = libpointer('singlePtr', zeros(kmax*sum(defvec), 1, 'single'));
calllib('bgmg', 'bgmg_retrieve_tag_r2_sum', 0, 0, round(dat.causal_pi * num_snp), kmax*sum(defvec), pBuffer);  check(); 
tag_r2_sum = pBuffer.Value;
tag_r2_sum=reshape(tag_r2_sum, [kmax, sum(defvec)])';
clear pBuffer
end

clf;figure(1); hold on;
for repi=1:10
ax = gca;
ax.ColorOrderIndex = 1;
dat = load(['H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_normeff_chr1\' sprintf('simu_h2=0.7_pi1u=0.01_rep=%i.trait1.mat', repi)]);zvec = dat.zvec;
zvec(~isfinite(zvec)) = 100;
%sum(~isfinite(zvec(defvec)))
calllib('bgmg', 'bgmg_set_zvec', 0, trait, sum(defvec), zvec(defvec));  check();
% calculate model cdf
zgrid = single(0:0.1:15); 
pBuffer = libpointer('singlePtr', zeros(length(zgrid), 1, 'single'));
cost=calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, 1.0, dat.sigsq);  check(); 
calllib('bgmg', 'bgmg_calc_univariate_pdf', 0, dat.causal_pi, 1.0, dat.sigsq, length(zgrid), zgrid, pBuffer);  check(); 
pdf = pBuffer.Value'; clear pBuffer
pdf = pdf / sum(1./weights(defvec));
if (zgrid(1) == 0), zgrid = [-fliplr(zgrid(2:end)) zgrid];pdf = [fliplr(pdf(2:end)) pdf]; end
model_cdf = cumsum(pdf)  * (zgrid(2) - zgrid(1)) ;
X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
model_logpvec = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), zgrid(zgrid>=0))); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)

% calculate data cdf        
zvec_idx = zvec(defvec); weights_idx = 1./weights(defvec); 
weights_idx=weights_idx/sum(weights_idx);
[data_y, si] = sort(-log10(2*normcdf(-abs(zvec_idx))));
data_x=-log10(cumsum(weights_idx(si),1,'reverse'));
data_idx = ([data_y(2:end); +Inf] ~= data_y);
hv_logp = -log10(2*normcdf(-zgrid(zgrid >= 0)));
data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

% plot QQ plots
hData     = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
hModel    = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;

qq_options=[];
if ~isfield(qq_options, 'qqlimy'), qq_options.qqlimy = 20; end;
if ~isfield(qq_options, 'qqlimx'), qq_options.qqlimx = 7; end;
plot([0 qq_options.qqlimy],[0 qq_options.qqlimy], 'k--');
xlim([0 qq_options.qqlimx]); ylim([0 qq_options.qqlimy]);
end
            
%histogram(mean(tag_r2_sum))

% !!!!!!!!!!! TBD prepare test data. !!!!!!!!!!
% [done] validate this with synetic QQ plots.
% make this concurrent --- good idea

koef_vec = logspace(-1, 1, 31);
cost_pi = [];cost_sig2beta = []; cost_sig2zero = [];
for koef = koef_vec
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi * koef, 1.0, dat.sigsq);  check(); fprintf('%.3f\t', cost); cost_pi(end+1, 1) = cost;
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, 1.0 * koef, dat.sigsq);  check(); fprintf('%.3f\t', cost); cost_sig2zero(end+1, 1) = cost;
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, 1.0, dat.sigsq * koef);  check(); fprintf('%.3f\n', cost); cost_sig2beta(end+1, 1) = cost;
end
cost_pi(cost_pi > 1e99) = nan;
figure(1);plot(log10(koef_vec), cost_pi)
figure(2);plot(log10(koef_vec), cost_sig2beta)
figure(3);plot(log10(koef_vec), cost_sig2zero)

if 0
%for i=22
    fprintf('%i - loading...', i);
    cI = load(['H:\work\hapgen_ldmat2_plink\' sprintf('bfile_merged_10K_ldmat_p01_SNPwind50k_chr%i.ld.mat', i)]);
    fprintf(' transfer...\n');
    calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, length(cI.r2), cI.index_A, cI.index_B, cI.r2);  check();
end

fprintf('converting...\n');
calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();

calllib('bgmg', 'bgmg_find_snp_sample', 0); check();


% make sure to convert to zero-based indices

calllib('bgmg', 'bgmg_dispose', 0);

calllib('bgmg', 'bgmg_set_tag_indices', 0, 5, 3, [2 4 5]-1);  check();
calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, 3, [1 2 2]-1, [2 4 5]-1, [0.5 0.3 0.2]);  check();
calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();


DLL_PUBLIC int64_t bgmg_set_ld_r2_coo(int context_id, int length, int* snp_index, int* tag_index, float* r2);
  DLL_PUBLIC int64_t bgmg_set_ld_r2_csr(int context_id);


zvec = [1;2;3];
calllib('bgmg', 'bgmg_set_zvec', 0, length(defvec), sum(defvec), find(defvec));

sum_r2 = [1 2 3; 4 5 6];
sum_r4 = [11 12 13; 14 15 16];
calllib('bgmg', 'bgmg_set_ref_ld', 1, size(sum_r2, 1), size(sum_r2, 2), sum_r2', sum_r4');  % transpose to convert matlab's fortran layout for 2D matrices into standard C layout
fprintf('%s\n', calllib('bgmg', 'bgmg_get_last_error'));


return

p = 0.1;
n = 10;
i = 0:n;
t=(p .^ i);
sum(t) - (1 - p.^(n+1)) ./ (1-p)
sum(t .* i) - p ./ (1-p).^2
sum(t .* i .* i) - p*(p+1)/(1-p).^3






return

% fft
% sig2beta, sig2zero, pi
% rho0, rho_beta

% 1. ifft precision is limited to 1e-18 (while double-precision technically go up to 1e-323)
% 2. 


count_r2 = [0 0 00 0 0];
vals_r2  = [0.2 0.4 0.6 0.8 1.0];

nsnp = sum(count_r2); nbins = length(vals_r2);

delvals = linspace(-15,15,2^7+1); delstep = delvals(2)-delvals(1); 
tickvals = 5*[-10:10]; ticks = interp1(delvals,1:length(delvals),tickvals); tickvals = tickvals(isfinite(ticks)); ticks = ticks(isfinite(ticks));
if mod(length(delvals),2)~=0, delvals = delvals(1:end-1); end % Make sure there's an even number of zbins (for fft/fftshift)


if 0
   params_sig0=1.0;
   zvec=linspace(-15,15,2^9+2); zvec=zvec(1:end-1); 
   zstep = zvec(2)-zvec(1); Fs = 1/zstep; fvec = -1/2*Fs*zvec/zvec(1);
   Fpdfvec = exp(-1/2*(2*pi*fvec*params_sig0).^2); 
   pdfvec = max(0,real(ifftshift(ifft(fftshift(Fpdfvec)))))
   plot(zvec, -log10(pdfvec), zvec, -log10(zstep * normpdf(zvec, 0, params_sig0)));
end

upsampfact = 1;
%upsampfact = 16; % This is needed to handle tiny sig_b * r^2 -- should find more computationally efficient solution
zvec = linspace(delvals(1),-delvals(1),upsampfact*length(delvals)+1); zvec = zvec(1:end-1);
zstep_orig = delvals(2)-delvals(1);
zstep = zvec(2)-zvec(1); Fs = 1/zstep; fvec = -1/2*Fs*zvec/zvec(1);

params.pi_1 = 1e-13;
params.sigd1 = sqrt(1e-4);
params.sig0 = 1.0;

Fpdfvec1 = exp(-1/2*(2*pi*fvec*params.sig0).^2); 
if 0
for bini = 1:nbins
  if count_r2(bini)>0
    tmp1 = exp(-1/2*(2*pi*fvec*params.sigd1).^2*vals_r2(bini)); 
    Fpdfvec1 = Fpdfvec1.*(params.pi_1*tmp1+(1-params.pi_1)).^count_r2(bini);
  end
end
end

pdfmat_indep = max(0,real(ifftshift(ifft(fftshift(Fpdfvec1)))));
pdfmat_indep = pdfmat_indep/(sum(pdfmat_indep(:))*zstep)*(zstep/zstep_orig);

clf; hold on
plot(zvec, -log10(normpdf(zvec, 0, params.sig0)), 'r')
plot(zvec, -log10(pdfmat_indep))
plot(zvec, -log10(Fpdfvec1))

return
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
  subplot(2,3,1); imagesc(Fpdfmat_indep); colormap(winter); axis equal tight xy;
  subplot(2,3,2); imagesc(Fpdfmat_pleio); colormap(winter); axis equal tight xy;
  subplot(2,3,3); imagesc(Fpdfmat_total); colormap(winter); axis equal tight xy;
  subplot(2,3,4); imagesc(log(pdfmat_indep),crange); colormap(winter); axis equal tight xy;
  subplot(2,3,5); imagesc(log(pdfmat_pleio),crange); colormap(winter); axis equal tight xy;
  subplot(2,3,6); imagesc(log(pdfmat_total),crange); colormap(winter); axis equal tight xy;
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
  subplot(2,3,1); imagesc(Fpdfmat_indep); colormap(winter); axis equal tight xy;
  subplot(2,3,2); imagesc(Fpdfmat_pleio); colormap(winter); axis equal tight xy;
  subplot(2,3,3); imagesc(Fpdfmat_total); colormap(winter); axis equal tight xy;
  subplot(2,3,4); imagesc(log(pdfmat_indep),crange); colormap(winter); axis equal tight xy;
  subplot(2,3,5); imagesc(log(pdfmat_pleio),crange); colormap(winter); axis equal tight xy;
  subplot(2,3,6); imagesc(log(pdfmat_total),crange); colormap(winter); axis equal tight xy;
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
    logFpdfvec2     = logFpdfvec2 + log(apovec.*(params.pi_2*exp(-1/2*(2*pi*fvec*sigd2).^2*vals_r2(bini))+(1-params.pi_2)))*count_r2(bini);
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

