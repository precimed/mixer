% This script runs Bivariate Causal Mixture for GWAS analysis. 

%{
% FULL LIST OF AVAILABLE PARAMETERS

% Full path to your libbgmg.so (linux), libbgmg.dylib (mac) or bgmg.dll (windows). See readme for how-to-build instructions.
bgmg_shared_library        = 'H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll';
bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';

% Input data
bim_file       = 'H:\GitHub\BGMG\LDSR\1000G_EUR_Phase3_plink\1000G.EUR.QC.@.bim';
frq_file       = 'H:\GitHub\BGMG\LDSR\1000G_EUR_Phase3_plink_freq\1000G.EUR.QC.@.frq';
plink_ld_bin   = 'H:\GitHub\BGMG\LDSR\1000G_EUR_Phase3_plink\1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'; 
chr_labels     = 1:22;

trait1_file    = 'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Data\PGC_SCZ_2014_EUR_qc_noMHC.sumstats.gz';
trait2_file    = 'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Data\PGC_BIP_2016_qc_noMHC.sumstats.gz';
trait1_params_file       = 'H:\GitHub\BGMG\LDSR\BGMG_results\PGC_SCZ_2014_EUR_qc_noMHC.ugmg.params.mat';
trait2_params_file       = 'H:\GitHub\BGMG\LDSR\BGMG_results\PGC_BIP_2016_qc_noMHC.ugmg.params.mat';
out_file                 = 'H:\GitHub\BGMG\LDSR\BGMG_results\PGC_SCZ_2014_EUR_qc_noMHC_vs_PGC_BIP_2016_qc_noMHC';

% Enable/disable features
DO_FIT_BGMG=true; 
STRATIFIED_QQ_PLOT=true; STRATIFIED_QQ_PLOT_DOWNSCALE=10; % enable/disable stratified QQ plot

% Optional parameters
exclude = '';                       % file containing SNP rs# to exclude from the anslysis
extract = '';                       % file containing SNP rs# to include in the anslysis
randprune_n=64; randprune_r2=0.1;   % random pruning options that define a weighting scheme on tag variants (avoid overcounting signal in large LD blocks)
kmax=5000;                          % number of sampling interation in pdf(z|params) model. Larger values => more accurate inference, but longer runtime, and larger memory usage
SEED=123;                           % seed for random number generator. Fix for reproducible results.
cache_tag_r2sum=1;                  % performance optimization. Set to 0 if you run out of RAM memory (but the model will run slower)
max_causal_fraction=0.03;           % upper threshold on polygenicity. This is required for technical reason - setting to 1.0 causes excesive memory usage.
r2min=0.05;                         % lower threshold for LD r2 values.
init_result_from_out_file='';       % path to .mat file with previous results. Use this together with DO_FIT_UGMG=false to make QQ plots on a larger set of variants, using previously fitted parameters.
CI_ALPHA=0.05;                      % enable confidence interval estimation
THREADS=-1;                         % specify how many threads to use (concurrency). "-1" means to use all available CPU power.
TolX = 1e-2; TolFun = 1e-2;         % fminserach tolerance (stop criteria)
z1max = nan; z2max = nan;           % enable right-censoring for z scores above certain threshold
qq_zgrid_lim = 38; qq_zgrid_step=0.25;  % control internal grid of z scores in stratified QQ plots

%}

if ~exist('out_file', 'var'), out_file = 'BGMG_result'; end;
if ~exist('bim_file', 'var'), error('bim_file is required'); end;
if ~exist('frq_file', 'var'), error('frq_file is required'); end;
if ~exist('trait1_file', 'var'), error('trait1_file is required'); end;
if ~exist('trait2_file', 'var'), error('trait2_file is required'); end;

% full path to bgmg shared library
if ~exist('bgmg_shared_library', 'var'), error('bgmg_shared_library is required'); end;
if ~exist('bgmg_shared_library_header', 'var'), [a,b,c]=fileparts(bgmg_shared_library); bgmg_shared_library_header = [fullfile(a, b), '.h']; clear('a', 'b','c'); end;

BGMG_cpp.unload(); 
BGMG_cpp.load(bgmg_shared_library, bgmg_shared_library_header);
BGMG_cpp.init_log([out_file, '.bgmglib.log']);
BGMG_cpp.log('out file: %s\n', out_file);

% plink_ld_mat is a file containing information about the LD structure. One
% file per chromosome. Each file should have index_A, index_B, r2
% variables. index_A and index_B corresponds to a zero-based index in the
% reference template. @ indicates the label of the chromosome.
if ~exist('plink_ld_bin', 'var'), error('plink_ld_bin is required'); end;
if ~exist('chr_labels', 'var'), chr_labels = 1:22; end;

if ~exist('randprune_n', 'var'), randprune_n = 64; end;
if ~exist('randprune_r2', 'var'), randprune_r2 = 0.1; end;
if ~exist('kmax', 'var'), kmax = 1000; end;
if ~exist('r2min', 'var'), r2min = 0.01; end;
if ~exist('max_causal_fraction', 'var'), max_causal_fraction = 0.03; end;
if ~exist('cache_tag_r2sum', 'var'), cache_tag_r2sum = 1; end;
if ~exist('SEED', 'var'), seed = nan; end;
if ~exist('z1max', 'var'), z1max = nan; end;  % 5.45 to remove GWS hits
if ~exist('z2max', 'var'), z2max = nan; end;
if ~exist('qq_zgrid_lim', 'var'), qq_zgrid_lim = 38; end;
if ~exist('qq_zgrid_step', 'var'), qq_zgrid_step = 0.25; end;

num_components = 3;  % bivariate

% The following three options control how to get bivariate params
% BGMG_cpp_run can take previous results file, and enhance it with
% additional features. Typical example would be to take a file produced
% after DO_FIT, and use it to produce QQ plots, etc.
if ~exist('DO_FIT_BGMG', 'var'), DO_FIT_BGMG = true; end;               % perform bivariate fit
if ~exist('init_from_params_file', 'var'), init_from_params_file = ''; end;

if isempty(init_from_params_file)
    if ~exist('trait1_params_file', 'var'), error('trait1_params_file is required'); end;
    if ~exist('trait2_params_file', 'var'), error('trait2_params_file is required'); end;

    trait1_params = load(trait1_params_file);
    if ~isfield(trait1_params, 'sig2_zero'), error('%s does not contain sig2_zero param', trait1_params_file); end;
    if ~isfield(trait1_params, 'sig2_beta'), error('%s does not contain sig2_beta param', trait1_params_file); end;
    if ~isfield(trait1_params, 'pi_vec'), error('%s does not contain pi_vec param', trait1_params_file); end;

    trait2_params = load(trait2_params_file);
    if ~isfield(trait2_params, 'sig2_zero'), error('%s does not contain sig2_zero param', trait2_params_file); end;
    if ~isfield(trait2_params, 'sig2_beta'), error('%s does not contain sig2_beta param', trait2_params_file); end;
    if ~isfield(trait2_params, 'pi_vec'), error('%s does not contain pi_vec param', trait2_params_file); end;
end

% Setup tolerance for fminsearch
if ~exist('TolX', 'var'), TolX = 1e-2; end;
if ~exist('TolFun', 'var'), TolFun = 1e-2; end;

if ~exist('FIT_FULL_MODEL', 'var'), FIT_FULL_MODEL = true; end;                % use full model (when false, use gaussian approximation)
if ~exist('STRATIFIED_QQ_PLOT', 'var'), STRATIFIED_QQ_PLOT = false; end;
if ~exist('STRATIFIED_QQ_PLOT_DOWNSCALE', 'var'), STRATIFIED_QQ_PLOT_DOWNSCALE = 1000; end;     % downscale #snps in stratified QQ plots (model prediction only)
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = 0.05; end;
if ~exist('THREADS', 'var'), THREADS = -1; end;
if ~exist('exclude', 'var'), exclude = ''; end;
if ~exist('extract', 'var'), extract = ''; end;

addpath('DERIVESTsuite');

bgmglib = BGMG_cpp();
bgmglib.dispose();
bgmglib.init(bim_file, frq_file, sprintf('%i ', chr_labels), trait1_file, trait2_file, exclude, extract);

bgmglib.set_option('r2min', r2min);
bgmglib.set_option('kmax', kmax);
bgmglib.set_option('max_causals', floor(max_causal_fraction * bgmglib.num_snp));
bgmglib.set_option('num_components', num_components);
bgmglib.set_option('cache_tag_r2sum', cache_tag_r2sum);
bgmglib.set_option('threads', THREADS);
if isfinite(SEED), bgmglib.set_option('seed', SEED); end;
if isfinite(z1max), bgmglib.set_option('z1max', z1max); end;
if isfinite(z2max), bgmglib.set_option('z2max', z2max); end;

for chr_index=1:length(chr_labels), 
    bgmglib.set_ld_r2_coo_from_file(strrep(plink_ld_bin,'@', sprintf('%i',chr_labels(chr_index)))); 
    bgmglib.set_ld_r2_csr(chr_labels(chr_index));
end;

bgmglib.set_weights_randprune(randprune_n, randprune_r2);

bgmglib.set_option('diag', 0);
% Preparation is done, BGMG library is fully setup. Now we can use it to
% calculate model QQ plots and univariate or bivariate cost function.

options = [];
options.total_het = sum(2*bgmglib.mafvec.*(1-bgmglib.mafvec));
options.verbose = true;
options.ci_alpha = CI_ALPHA;
options.title = TITLE;
options.fit_full_model = FIT_FULL_MODEL;
options.trait1_nval = median(bgmglib.nvec1);
options.trait2_nval = median(bgmglib.nvec2);

disp(options)

if ~isempty(init_from_params_file)
    BGMG_cpp.log('Loading params from %s...\n', init_from_params_file);
    trait12_params = load(init_from_params_file);
    BGMG_cpp.log('Univariate and bivariate params load from initial file');
    trait12_params = BGMG_cpp_fit_bivariate_fast_constrained(trait12_params)
    clear('tmp');
else
    % if there is no initial approximation, setup it from fast model
    trait12_params = BGMG_cpp_fit_bivariate_fast(trait1_params, trait2_params);
end

% Fit bivariate or univariate model to the data
bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);

result = [];
if DO_FIT_BGMG
    result.bivariate = BGMG_cpp_fit_bivariate(trait12_params, options);
    trait12_params = result.bivariate.params;
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);  % restore fast_cost flag as BGMG_cpp_fit_univariate may change it (CI intervals are based on fast cost function)
end

result.trait1_file = trait1_file;
result.trait2_file = trait2_file;
result.bim_file = bim_file;
result.frq_file = frq_file;
result.options = options;
result.params = trait12_params;

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.params.mat'], '-struct', 'trait12_params');
BGMG_cpp.log('Params saved to %s.params.mat\n', out_file);

bgmglib.set_option('diag', 0);

try
if STRATIFIED_QQ_PLOT
    options.downscale = STRATIFIED_QQ_PLOT_DOWNSCALE;
	if isfinite(qq_zgrid_lim), options.qq_zgrid_lim = qq_zgrid_lim; end;
	if isfinite(qq_zgrid_step), options.qq_zgrid_step = qq_zgrid_step; end;

    [figures, plot_data] = BGMG_cpp_stratified_qq_plot(trait12_params, options);
    result.bivariate.stratified_qq_plot_fit_data.trait1 = plot_data(1, :);
    result.bivariate.stratified_qq_plot_fit_data.trait2 = plot_data(2, :);
    %print(figures.tot{1}, sprintf('%s.trait1.stratqq.pdf', out_file), '-dpdf')
    %print(figures.tot{2}, sprintf('%s.trait2.stratqq.pdf', out_file), '-dpdf')
end
catch err
    BGMG_cpp.log_error(err)
    result.bivariate.stratified_qq_plot_fit_data_error = err;
end

bgmglib.set_option('diag', 0);

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.mat'], 'result');
BGMG_cpp.log('Results saved to %s.mat\n', out_file);

if exist('jsonencode')
    str=jsonencode(result);
    fileID = fopen([out_file '.json'] , 'w');
    fprintf(fileID, str);
    fclose(fileID);
    BGMG_cpp.log('Results saved to %s.json\n', out_file);
else
    warning('jsonencode does not exist, cannot convert resulting mat files to json')
end
