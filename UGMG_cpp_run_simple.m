% This script can be run univariate mixture analysis from command line like this:
%
%   matlab -nodisplay -nosplash -nodesktop -r "trait1_file='PGC_SCZ_2014.mat'; reference_file='1kG_phase3_EUR_11015883_reference_holland.mat'; BGMG_run; exit;"
%
% You may use this script with one trait (univariate analysis).
% The results are saved into a text file as well as into .mat file.
% Currently you may need to modify this file to pass additional
% parameters, later all parameters will be exposed as a text config file.
%
% Run flow:
% 1. Load and initialize BGMG library (bgmglib)
% 2. Load reference file (posvec, chrnumvec, mafvec)
% 3. Load trait data (zvec, nvec)
% 6. Load LD structure from binary files directly into bgmglib. Setup weights based on random pruning.
% 7. Setup params ('initial approximation'). Two options available
%       - load params from previous run
%       - fit params using gaussian approximation ('fit from scratch')
% 8. [optional] fit UGMG model
% 9. [optional] produce QQ plots
% 10. Save results to <out_file>.[mat, pdf, log]

if 0
% example of the input
bgmg_shared_library        = 'H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll';
bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';
plink_ld_bin   = 'H:\work\hapgen_ldmat2_plink\1000Genome_ldmat_p05_SNPwind50k_chr@.ld.bin'; chr_labels = 1:22;
bim_file       = 'H:\work\1000Genome\chr@.bim';
frq_file       = 'H:\work\1000Genome\chr@.frq';
trait1_file    = 'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Data\PGC_SCZ_2014_EUR_qc_noMHC.sumstats.gz';
out_file       = 'PGC_SCZ_2014_EUR_qc_noMHC.ugmg';
max_causal_fraction=0.01;
cache_tag_r2sum=true;
kmax=1000;
SEED=123;
r2min=0.05;
DO_FIT_UGMG=true; QQ_PLOT=true; QQ_PLOT_DOWNSCALE = 10;
end

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
num_components = 1;  % univariate

% The following three options control how to get univariate & bivariate params
% BGMG_cpp_run can take previous results file, and enhance it with
% additional features. Typical example would be to take a file produced
% after DO_FIT, and use it to produce QQ plots, etc.
% Allows to fir BGMG with already fitted univaraite estimates
if ~exist('DO_FIT_UGMG', 'var'), DO_FIT_UGMG = true; end;               % perform univariate fit
if ~exist('init_result_from_out_file', 'var'), init_result_from_out_file = ''; end;

% Setup tolerance for fminsearch
if ~exist('TolX', 'var'), TolX = 1e-2; end;
if ~exist('TolFun', 'var'), TolFun = 1e-2; end;

if ~exist('FIT_FULL_MODEL', 'var'), FIT_FULL_MODEL = true; end;                % use full model (when false, use gaussian approximation)
if ~exist('QQ_PLOT', 'var'), QQ_PLOT = false; end;   % make QQ plots
if ~exist('QQ_PLOT_DOWNSCALE', 'var'), QQ_PLOT_DOWNSCALE = 100; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('QQ_PLOT_BINS', 'var'), QQ_PLOT_BINS = false; end;   % make QQ plots
if ~exist('QQ_PLOT_BINS_DOWNSCALE', 'var'), QQ_PLOT_BINS_DOWNSCALE = 10; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('UGMG_LOGLIKE_PLOT', 'var'), UGMG_LOGLIKE_PLOT = false; end;
if ~exist('POWER_PLOT', 'var'), POWER_PLOT = false; end;  % make power plots with fitted parameters
if ~exist('POWER_PLOT_DOWNSCALE', 'var'), POWER_PLOT_DOWNSCALE = 10; end;  % make power plots with fitted parameters
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;
if ~exist('THREADS', 'var'), THREADS = -1; end;

addpath('DERIVESTsuite');
addpath('PolyfitnTools');

bgmglib = BGMG_cpp();
bgmglib.dispose();
bgmglib.init(bim_file, frq_file, sprintf('%i ', chr_labels), trait1_file, '');

bgmglib.set_option('r2min', r2min);
bgmglib.set_option('kmax', kmax);
bgmglib.set_option('max_causals', floor(max_causal_fraction * bgmglib.num_snp));
bgmglib.set_option('num_components', num_components);
bgmglib.set_option('cache_tag_r2sum', cache_tag_r2sum);
bgmglib.set_option('threads', THREADS);
if isfinite(SEED), bgmglib.set_option('seed', SEED); end;

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

params = [];

% Load params from previous runs.
if ~isempty(init_result_from_out_file)
    BGMG_cpp.log('Loading init file from %s...\n', init_result_from_out_file);
    tmp = load(init_result_from_out_file);
    params.univariate{1} = tmp.result.univariate{1}.params;
    BGMG_cpp.log('Univariate params load from initial file');
    clear('tmp');
end

% If there is no initial approximation, setup it from fast model
if ~isfield(params, 'univariate'),
    params.univariate{1} = BGMG_cpp_fit_univariate_fast(1);
end

% Fit bivariate or univariate model to the data
bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);

result = [];
if DO_FIT_UGMG
    result.univariate{1} = BGMG_cpp_fit_univariate(1, params.univariate{1}, options);
    params.univariate{1} = result.univariate{1}.params;
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);  % restore fast_cost flag as BGMG_cpp_fit_univariate may change it (CI intervals are based on fast cost function)
end

result.trait1_file = trait1_file;
result.bim_file = bim_file;
result.frq_file = frq_file;
result.options = options;
result.params = params;

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.preliminary.mat'], 'result');
BGMG_cpp.log('Results saved to %s.preliminary.mat\n', out_file);

bgmglib.set_option('diag', 0);

% Produce power plots
if POWER_PLOT
    options.downscale = POWER_PLOT_DOWNSCALE;trait_index=1;
    figures.tot = figure;
    plot_data = BGMG_cpp_power_plot(params.univariate{trait_index}, trait_index, options);
    result.univariate{trait_index}.power_plot_data = plot_data;
    print(figures.tot, sprintf('%s.power.pdf', out_file), '-dpdf')
end

% Produce QQ plots
if QQ_PLOT
    options.downscale = QQ_PLOT_DOWNSCALE;trait_index=1;
    figures.tot = figure;
    plot_data = BGMG_cpp_qq_plot(params.univariate{trait_index}, trait_index, options);

    % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
    result.univariate{trait_index}.qq_plot_data = plot_data;
    print(figures.tot, sprintf('%s.qq.pdf', out_file), '-dpdf')
end

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
            plot_data = BGMG_cpp_qq_plot(params.univariate{trait_index}, trait_index, options);
            result.univariate{trait_index}.qq_plot_bins_data{i,j} = plot_data;
        end
    end
    set(total_fig,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(total_fig, sprintf('%s.qq.bins.pdf', out_file), '-dpdf')
end

bgmglib.set_option('diag', 0);

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.mat'], 'result');
BGMG_cpp.log('Results saved to %s.mat\n', out_file);
