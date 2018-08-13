bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll';
bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';

datafolder = 'H:\Dropbox\shared\BGMG\bgmgcpp';
plink_ld_bin = fullfile(datafolder, 'bfile_merged_ldmat_p01_SNPwind50k_chr1.ld.bin');
logfile = 'BGMG_cpp_example.bgmglib.log';
sumstat_file = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=1e-03_rep=5_tag1=customPolygenicOverlapAt0p375_tag2=evenPolygenicity';

load(fullfile(datafolder, sprintf('%s.trait1.mat', sumstat_file)));
params = load(fullfile(datafolder, sprintf('%s.params.mat', sumstat_file)));

load(fullfile(datafolder, 'mafvec.mat'))
load(fullfile(datafolder, 'posvec.mat'))
load(fullfile(datafolder, 'chrnumvec.mat'))
load(fullfile(datafolder, 'defvec_hapmap3.mat')); defvec = defvec & chrnumvec == 1 & isfinite(zvec);

BGMG_cpp.unload(); 
BGMG_cpp.load(bgmg_shared_library, bgmg_shared_library_header);
BGMG_cpp.init_log(logfile);

bgmglib = BGMG_cpp(); bgmglib.dispose()
bgmglib.defvec = defvec;
bgmglib.zvec1 = zvec(defvec);
bgmglib.nvec1 = 100000 * ones(sum(defvec), 1);

bgmglib.set_option('r2min', 0.01);
bgmglib.set_option('kmax', 1000);
bgmglib.set_option('max_causals', floor(0.01 * length(defvec)));
bgmglib.set_option('num_components', 1);  % use 3 for bivariate
bgmglib.set_option('cache_tag_r2sum', 1);
bgmglib.set_option('threads', 16);  % omp concurrency

bgmglib.mafvec = mafvec;
bgmglib.chrnumvec = chrnumvec;

bgmglib.set_ld_r2_coo_from_file(plink_ld_bin);
bgmglib.set_ld_r2_csr();

randprune_n = 16; randprune_r2 = 0.8;
bgmglib.set_weights_randprune(randprune_n, randprune_r2);
histogram(bgmglib.weights);  % histogram of weights induced by random pruning

bgmglib.set_option('diag', 0);  % print a bunch of useful into into log file

% plot univariate cost function (log likelihood) as a function of pivec
trait_index = 1;   % we've set zvec for the first trait (bgmglib.zvec1), hence trait_index=1.
pivec = sum(params.pi_vec_trait1) * (0.5:0.1:1.5);
cost = [];
for pival=pivec
    cost(end+1, 1) = bgmglib.calc_univariate_cost(trait_index, pival, 1.0, params.sig2_beta(1));
end
plot(pivec, cost, '.-')

% QQ plot
% Calculate data_logpvec
hv_z = linspace(0, min(max(abs(bgmglib.zvec1)), 38.0), 10000);
[data_y, si] = sort(-log10(2*normcdf(-abs(bgmglib.zvec1))));
data_weights = bgmglib.weights(si); data_weights = data_weights ./ sum(data_weights);
data_x=-log10(cumsum(data_weights,1,'reverse'));
data_idx = ([data_y(2:end); +Inf] ~= data_y);
hv_logp = -log10(2*normcdf(-hv_z));
data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

% Calculate model_logpvec
zgrid = single(0:0.25:15); 
pdf = bgmglib.calc_univariate_pdf(trait_index, sum(params.pi_vec_trait1), 1, params.sig2_beta(1), zgrid);
pdf = pdf / sum(bgmglib.weights);
if (zgrid(1) == 0), zgrid = [-fliplr(zgrid(2:end)) zgrid];pdf = [fliplr(pdf(2:end)) pdf]; end
model_cdf = cumsum(pdf)  * (zgrid(2) - zgrid(1)) ;
X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
model_logpvec = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), hv_z')); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)

hData  = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
hModel = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;
hNull  = plot(hv_logp,   hv_logp, 'k--', 'LineWidth', 1); hold on;

