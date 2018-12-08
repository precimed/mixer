% This script runs Univariate Causal Mixture for GWAS analysis. 
% To start analysis from command line run:
%
%   matlab -nodisplay -nosplash -nodesktop -r "bim_file='chr@.bim'; frq_file='chr@.frq'; trait1_file='PGC_SCZ_2014_noMHC.sumstats.gz'; ...; BGMG_run; exit;"
%
% See below for the FULL LIST OF AVAILABLE PARAMETERS.
%
% Run flow:
% 1. Load and initialize BGMG library (bgmglib, native c++ plugin for Matlab)
% 2. Load reference pannel (SNP/CHR/BP/A1/A2) from .bim file(s), plink format
% 3. Load allele frequencies for reference pannel (FRQ) from .frq file(s) plink format
% 4. Load trait data (zvec, nvec) from .sumstats.gz LDSR-formatted file
%    (can be output of https://github.com/bulik/ldsc/munge_sumstats.py)
% 5. Load LD structure from binary files
%    (binary files should be produced from plink ***.ld.gz files using bgmg-cli converter)
% 6. Setup weights based on random pruning.
% 7. Setup params ('initial approximation'). Two options available
%       - load params from previous run
%       - fit params using gaussian approximation ('fit from scratch')
% 8. [optional] fit UGMG parameters using full model for GWAS z scores
% 9. [optional] produce QQ plots
% 10. [optional] produce partitioned QQ plots (bins of MAF and LD score)
% 11. [optional] produce power curves (proportion of heritability explained as a function of GWAS sample size)
% 12. Save results to <out_file>.[mat, pdf, log]

% NB! you may uncomment the following lines (e.i. remove '%{' and '%}' symbold) if you wish to set parameters here in .m file, and not from command line.

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
out_file       = 'H:\GitHub\BGMG\LDSR\BGMG_results\PGC_SCZ_2014_EUR_qc_noMHC.ugmg';

trait1_file    = 'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Data\PGC_BIP_2016_qc_noMHC.sumstats.gz';
out_file       = 'H:\GitHub\BGMG\LDSR\BGMG_results\PGC_BIP_2016_qc_noMHC.ugmg';

% Enable/disable features
DO_FIT_UGMG=true; 
QQ_PLOT=true; QQ_PLOT_DOWNSCALE = 10;         % enable/disable QQ plots
QQ_PLOT_BINS=true; QQ_PLOT_BINS_DOWNSCALE=10; % enable/disable partitioned QQ plot (maf/ldscore bins)
POWER_PLOT=true; POWER_PLOT_DOWNSCALE=10;     % enable/disable power plots

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
z1max = nan;                        % enable right-censoring for z scores above certain threshold

%}

if ~exist('out_file', 'var'), out_file = 'UGMG_result'; end;
if ~exist('bim_file', 'var'), error('bim_file is required'); end;
if ~exist('frq_file', 'var'), error('frq_file is required'); end;
if ~exist('trait1_file', 'var'), error('trait1_file is required'); end;

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
num_components = 1;  % univariate

% The following three options control how to get univariate & bivariate params
% BGMG_cpp_run can take previous results file, and enhance it with
% additional features. Typical example would be to take a file produced
% after DO_FIT, and use it to produce QQ plots, etc.
% Allows to fir BGMG with already fitted univaraite estimates
if ~exist('DO_FIT_UGMG', 'var'), DO_FIT_UGMG = true; end;               % perform univariate fit
if ~exist('init_from_params_file', 'var'), init_from_params_file = ''; end;

% Setup tolerance for fminsearch
if ~exist('TolX', 'var'), TolX = 1e-2; end;
if ~exist('TolFun', 'var'), TolFun = 1e-2; end;

if ~exist('FIT_FULL_MODEL', 'var'), FIT_FULL_MODEL = true; end;                % use full model (when false, use gaussian approximation)
if ~exist('QQ_PLOT', 'var'), QQ_PLOT = false; end;   % make QQ plots
if ~exist('QQ_PLOT_DOWNSCALE', 'var'), QQ_PLOT_DOWNSCALE = 100; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('QQ_PLOT_BINS', 'var'), QQ_PLOT_BINS = false; end;   % make QQ plots
if ~exist('QQ_PLOT_BINS_DOWNSCALE', 'var'), QQ_PLOT_BINS_DOWNSCALE = 10; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('POWER_PLOT', 'var'), POWER_PLOT = false; end;  % make power plots with fitted parameters
if ~exist('POWER_PLOT_DOWNSCALE', 'var'), POWER_PLOT_DOWNSCALE = 10; end;  % make power plots with fitted parameters
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = 0.05; end;
if ~exist('THREADS', 'var'), THREADS = -1; end;
if ~exist('exclude', 'var'), exclude = ''; end;
if ~exist('extract', 'var'), extract = ''; end;

addpath('DERIVESTsuite');

bgmglib = BGMG_cpp();
bgmglib.dispose();
bgmglib.init(bim_file, frq_file, sprintf('%i ', chr_labels), trait1_file, '', exclude, extract);

bgmglib.set_option('r2min', r2min);
bgmglib.set_option('kmax', kmax);
bgmglib.set_option('max_causals', floor(max_causal_fraction * bgmglib.num_snp));
bgmglib.set_option('num_components', num_components);
bgmglib.set_option('cache_tag_r2sum', cache_tag_r2sum);
bgmglib.set_option('threads', THREADS);
if isfinite(SEED), bgmglib.set_option('seed', SEED); end;
if isfinite(z1max), bgmglib.set_option('z1max', z1max); end;

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

disp(options)

% Load params from previous runs.
if ~isempty(init_from_params_file)
    BGMG_cpp.log('Loading params from %s...\n', init_from_params_file);
    trait1_params = load(init_from_params_file);
    trait1_params = BGMG_cpp_fit_univariate_fast_constrained(1, trait1_params);
    BGMG_cpp.log('Univariate params load from initial file');
    clear('tmp');
else
    % If there is no initial approximation, setup it from fast model
    trait1_params = BGMG_cpp_fit_univariate_fast(1);
end

% Fit bivariate or univariate model to the data
bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);

result = [];
if DO_FIT_UGMG
    result.univariate{1} = BGMG_cpp_fit_univariate(1, trait1_params, options);
    trait1_params = result.univariate{1}.params;
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);  % restore fast_cost flag as BGMG_cpp_fit_univariate may change it (CI intervals are based on fast cost function)
end

result.trait1_file = trait1_file;
result.bim_file = bim_file;
result.frq_file = frq_file;
result.options = options;
result.params = trait1_params;

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.params.mat'], '-struct', 'trait1_params');
BGMG_cpp.log('Params saved to %s.params.mat\n', out_file);

bgmglib.set_option('diag', 0);

% Produce power plots
try
if POWER_PLOT
    options.downscale = POWER_PLOT_DOWNSCALE;trait_index=1;
    figures.tot = figure;
    plot_data = BGMG_cpp_power_plot(trait1_params, trait_index, options);
    result.univariate{trait_index}.power_plot_data = plot_data;
    print(figures.tot, sprintf('%s.power.pdf', out_file), '-dpdf')
end
catch err
    BGMG_cpp.log_error(err)
    result.univariate{trait_index}.power_plot_data_error = err;
end

% Produce QQ plots
try
if QQ_PLOT
    options.downscale = QQ_PLOT_DOWNSCALE;trait_index=1;
    figures.tot = figure;
    plot_data = BGMG_cpp_qq_plot(trait1_params, trait_index, options);

    % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
    result.univariate{trait_index}.qq_plot_data = plot_data;
    print(figures.tot, sprintf('%s.qq.pdf', out_file), '-dpdf')
end
catch err
    BGMG_cpp.log_error(err)
    result.univariate{trait_index}.qq_plot_data_error = err;
end

try
if QQ_PLOT_BINS
    options.downscale = QQ_PLOT_BINS_DOWNSCALE;trait_index=1;
    mafvec = bgmglib.mafvec(bgmglib.defvec);
    ldscore = bgmglib.ld_tag_r2_sum;
    maf_bins = [-inf quantile(mafvec, 2) inf];
    tld_bins = [-inf quantile(ldscore, 2) inf];
    total_fig = figure;set(gcf, 'Position', get(0, 'Screensize'));
    options.full_annotation = false;
    for i=1:3
        for j=1:3
            subplot(3,3,(i-1)*3+j);
            options.title = sprintf('$$ maf \\in [%.3f,%.3f) $$ \n $$ L \\in [%.3f,%.3f) $$', maf_bins(i), maf_bins(i+1), tld_bins(j), tld_bins(j+1));
            options.mask = (mafvec > maf_bins(i)) & (mafvec <= maf_bins(i+1)) & (ldscore > tld_bins(j)) & (ldscore <= tld_bins(j+1));
            plot_data = BGMG_cpp_qq_plot(trait1_params, trait_index, options);
            result.univariate{trait_index}.qq_plot_bins_data{i,j} = plot_data;
        end
    end
    set(total_fig,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(total_fig, sprintf('%s.qq.bins.pdf', out_file), '-dpdf')
end
catch err
    BGMG_cpp.log_error(err)
    result.univariate{trait_index}.qq_plot_bins_data_error = err;
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
