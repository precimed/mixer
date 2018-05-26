if ~libisloaded('bgmg'), loadlibrary('H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll', 'H:\GitHub\BGMG\src\bgmg_matlab.h'); end;

%loadlibrary('H:\GitHub\BGMG\src\build_win\bin\Debug\bgmg.dll', 'H:\GitHub\BGMG\src\bgmg_matlab.h');
%libfunctions('bgmg')
%unloadlibrary('bgmg');

ref = load('H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins.mat');
hits = sum(ref.pruneidxmat, 2); w_ld = size(ref.pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan; weights = 1./w_ld;
defvec = ref.defvec & (w_ld > 0); % & (ref.chrnumvec == 1);
tag_indices = find(defvec);

m2c = @(x)(x-1); % convert matlab to c indices
check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
context=0;kmax = 1;
calllib('bgmg', 'bgmg_set_tag_indices', 0, length(ref.defvec), length(tag_indices), m2c(tag_indices));  check();
calllib('bgmg', 'bgmg_set_option', 0,  'r2min', 0.01); check();
calllib('bgmg', 'bgmg_set_option', 0,  'kmax', kmax); check();
calllib('bgmg', 'bgmg_set_option', 0,  'max_causals', 200000); check();

for c=22:-1:1
    fprintf('chr %i...\n', c);
    c1 = load(['H:\work\hapgen_ldmat2_plink\', sprintf('bfile_merged_10K_ldmat_p01_SNPwind50k_chr%i.ld.mat', c)]);
    c1.index_A = c1.index_A + 1; 
    c1.index_B = c1.index_B + 1; % 0-based, comming from python
    calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, length(c1.r2), m2c(c1.index_A), m2c(c1.index_B), c1.r2);  check();
    clear('c1');
end

tic
calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();
hvec = ref.mafvec .* (1-ref.mafvec) * 2;
calllib('bgmg', 'bgmg_set_hvec', 0, length(hvec), hvec);  check();
toc





trait=1; 

%calllib('bgmg', 'bgmg_set_zvec', 0, trait, sum(defvec), zvec(defvec));  check();
%calllib('bgmg', 'bgmg_set_nvec', 0, trait, sum(defvec), nvec(defvec));  check();
calllib('bgmg', 'bgmg_set_weights', 0, sum(defvec), weights(defvec));  check();

return

dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_traits\simu_h2=0.1_pi1u=0.0001_rep=3.trait1.mat');zvec = dat.zvec;
nvec = ones(size(zvec)) * 100000;

%def1= load('H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'); def1=def1.defvec; def1 = def1(1:num_snp);
%def2= load('H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'); def2=def2.defvec; def2 = def2(1:num_snp);
%defvec = def1 & def2 & isfinite(zvec) & (weights>=1); clear('def1', 'def2');
%hvec = ref.hvec; %2 * ref.mafvec(1:num_snp) .* ref.mafvec(1:num_snp);


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

h2vec_str = {'0.1', '0.4', '0.7'};
pivec_str = {'1e-05', '0.0001', '0.001', '0.01'};

for h2_index = 1:length(h2vec_str)
for pi_index = 1:length(pivec_str)
for repi=1:1
   subplot(3, 4, (h2_index-1)*4 + pi_index);
   title(sprintf('h2=%s, pi=%s', h2vec_str{h2_index}, pivec_str{pi_index}));
%end;end;end;
ax = gca;
ax.ColorOrderIndex = 1;
dat= load(['H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_traits\', sprintf('simu_h2=%s_pi1u=%s_rep=%i.trait1.mat', h2vec_str{h2_index}, pivec_str{pi_index}, repi)]);
zvec = dat.zvec; nvec = ones(size(zvec)) * 100000;
zvec(~isfinite(zvec)) = 100;
%sum(~isfinite(zvec(defvec)))
calllib('bgmg', 'bgmg_set_zvec', 0, trait, sum(defvec), zvec(defvec));  check();
calllib('bgmg', 'bgmg_set_nvec', 0, trait, sum(defvec), nvec(defvec));  check();
% calculate model cdf
zgrid = single(0:0.1:15); 
if 1
pBuffer = libpointer('singlePtr', zeros(length(zgrid), 1, 'single'));
cost=calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, 1.0, dat.sigsq);  check(); 
calllib('bgmg', 'bgmg_calc_univariate_pdf', 0, dat.causal_pi, 1.0, dat.sigsq, length(zgrid), zgrid, pBuffer);  check(); 
pdf = pBuffer.Value'; clear pBuffer
pdf = pdf / sum(weights(defvec));
if (zgrid(1) == 0), zgrid = [-fliplr(zgrid(2:end)) zgrid];pdf = [fliplr(pdf(2:end)) pdf]; end
model_cdf = cumsum(pdf)  * (zgrid(2) - zgrid(1)) ;
X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
model_logpvec = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), zgrid(zgrid>=0))); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)
end
% calculate data cdf        
zvec_idx = zvec(defvec); weights_idx = weights(defvec); 
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
drawnow
end
end
end