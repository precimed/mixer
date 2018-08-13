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
% 3. Load trait data (zvec, nvec) and user-provided defvecs
% 4. [optional] (for debugging)     Restrict analysis to chr1
% 5. [optional] (hardprune feature) Restrict set of SNPs based on LD structure
% 6. Load LD structure from binary files directly into bgmglib. Setup weights based on random pruning.
% 7. Setup params ('initial approximation'). Three options available
%       - load true params from simu file (sig2zero is fitted in this case)
%       - load params from previous run
%       - fit params using gaussian approximation ('fit from scratch')
% 8. [optional] fit UGMG model
% 9. [optional] produce QQ plots
% 10. Save results to <out_file>.[mat, pdf, log]

if 0
bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll';
bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';
plink_ld_bin = 'H:\work\hapgen_ldmat2_plink\bfile_merged_ldmat_p01_SNPwind50k_chr@.ld.bin'; chr_labels = 1:22;
reference_file = 'H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins_wld.mat';
filename = 'simu_h2=0.4_pi1u=0.001_rep=1';
trait1_file = ['H:\GitHub\BGMG\' filename '.trait1.mat']; trait1_nvec=100000;

kmax=1000; cache_tag_r2sum=true; r2min=0.05; SEED=123;
out_folder = 'results_2018_08_11.kmax=1000';

% QQ plots with true params
out_file = fullfile(out_folder, [filename '.true']); 
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'};
simu_params_file = ['H:\GitHub\BGMG\' filename '.params.mat']; init_result_from_out_file='';
DO_FIT_UGMG=false; QQ_PLOT=true; QQ_PLOT_DOWNSCALE = 10; QQ_PLOT_BINS=1; QQ_PLOT_BINS_DOWNSCALE = 1; UGMG_cpp_run;

% Fit UGMG parameters
out_file = fullfile(out_folder, [filename '.fit']);
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3_hardprune_p1.mat'};
simu_params_file = ''; init_result_from_out_file = '';
DO_FIT_UGMG=true; QQ_PLOT=true; QQ_PLOT_DOWNSCALE = 10; UGMG_cpp_run;

% QQ plots with fitted params
out_file = fullfile(out_folder, [filename '.fit.test']);
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'};
simu_params_file = ''; init_result_from_out_file = fullfile(out_folder, [filename '.fit.mat']);
DO_FIT_UGMG=false; QQ_PLOT=true; QQ_PLOT_DOWNSCALE = 10; UGMG_cpp_run

% Test hardprune feature
out_file = fullfile(out_folder, [filename '.hardprune']);
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'};
simu_params_file = ''; init_result_from_out_file = '';
hardprune_r2 = 0.1; hardprune_plink_ld_bin = plink_ld_bin;
DO_FIT_UGMG=false; QQ_PLOT=false; QQ_PLOT_DOWNSCALE = 10; UGMG_cpp_run
end


if ~exist('out_file', 'var'), out_file = 'UGMG_result'; end;
if exist('DO_FIT', 'var'), error('DO_FIT is deprecated; use DO_FIT_UGMG'); end;

% full path to bgmg shared library
if ~exist('bgmg_shared_library', 'var'), error('bgmg_shared_library is required'); end;
if ~exist('bgmg_shared_library_header', 'var'), [a,b,c]=fileparts(bgmg_shared_library); bgmg_shared_library_header = [fullfile(a, b), '.h']; clear('a', 'b','c'); end;

BGMG_cpp.unload(); 
BGMG_cpp.load(bgmg_shared_library, bgmg_shared_library_header);
BGMG_cpp.init_log([out_file, '.bgmglib.log']);
BGMG_cpp.log('out file: %s\n', out_file);

% defvec_file determine the set of tag S\NPs. It must have a binary variable "defvec" of the same length as defined by reference file ( 1 = tag snp, 0 = exclude from the analysis ). 
% When multiple defvec_files are defined we take an overlap (e.i. consider only tag SNPs defined across all defvec files).
if ~exist('defvec_files', 'var'), defvec_files = {}; end;

% plink_ld_mat is a file containing information about the LD structure. One
% file per chromosome. Each file should have index_A, index_B, r2
% variables. index_A and index_B corresponds to a zero-based index in the
% reference template. @ indicates the label of the chromosome.
if ~exist('plink_ld_bin', 'var'), error('plink_ld_bin is required'); end;
if ~exist('chr_labels', 'var'), chr_labels = 1:22; end;

if ~exist('trait1_file', 'var'), error('trait1_file is required'); end;
if ~exist('trait1_nvec', 'var'), trait1_nvec = nan; end;

if ~exist('randprune_n', 'var'), randprune_n = 10; end;
if ~exist('randprune_r2', 'var'), randprune_r2 = 0.8; end;
if ~exist('hardprune_r2', 'var'), hardprune_r2 = nan; end; % 1 iterations at this threshold to reduce set of SNPs (as in defvec)
if ~exist('hardprune_plink_ld_bin', 'var'), hardprune_plink_ld_bin = ''; end;

if ~exist('kmax', 'var'), kmax = 1000; end;
if ~exist('r2min', 'var'), r2min = 0.01; end;
if ~exist('max_causal_fraction', 'var'), max_causal_fraction = 0.03; end;
if ~exist('cache_tag_r2sum', 'var'), cache_tag_r2sum = 1; end;
if ~exist('SEED', 'var'), seed = nan; end;

% The following three options control how to get univariate & bivariate params
% BGMG_cpp_run can take previous results file, and enhance it with
% additional features. Typical example would be to take a file produced
% after DO_FIT, and use it to produce QQ plots, etc.
% Allows to fir BGMG with already fitted univaraite estimates
if ~exist('DO_FIT_UGMG', 'var'), DO_FIT_UGMG = true; end;               % perform univariate fit
if ~exist('simu_params_file', 'var'), simu_params_file = ''; end                  % load true parameters from 'simu' (simulation tool)
if ~exist('init_result_from_out_file', 'var'), init_result_from_out_file = ''; end;
if ~isempty(simu_params_file) && ~isempty(init_result_from_out_file), error('simu_params_file can not be used with init_result_from_out_file'); end;

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
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;
if ~exist('THREADS', 'var'), THREADS = nan; end;

if POWER_PLOT, error('not yet implemented in c++ version'); end;

% reference file containing mafvec, chrnumvec and posvec for all SNPs to consider in this analysis. 
if ~exist('reference_file', 'var'), error('reference_file is required'); end;
BGMG_cpp.log('Loading reference from %s... ', reference_file); ref = load(reference_file, 'mafvec', 'posvec', 'chrnumvec'); fprintf('OK.\n');
clear('defvec'); defvec_tmp = true(length(ref.mafvec), 1);
 
for i=1:length(defvec_files)
    BGMG_cpp.log('Loading %s... ', defvec_files{i}); cur_defvec = load(defvec_files{i}); defvec_tmp = defvec_tmp & cur_defvec.defvec; fprintf('OK.\n');
    BGMG_cpp.log('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));
end; clear('i');

BGMG_cpp.log('Loading %s...', trait1_file); trait1_data = load(trait1_file); fprintf('OK.\n'); num_components = 1;
if (length(trait1_data.zvec) ~= length(ref.chrnumvec)), error('trait1_file is incompatible with the reference'); end;
if isfinite(trait1_nvec), trait1_data.nvec = ones(size(trait1_data.zvec)) * trait1_nvec; end;
if ~isfield(trait1_data, 'nvec'), error('nvec is not available in trait1_file, and trait1_nvec parameter is not set'); end;
cur_defvec.defvec = isfinite(trait1_data.zvec + trait1_data.nvec); defvec_tmp = defvec_tmp & cur_defvec.defvec;
BGMG_cpp.log('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));

if (length(chr_labels) == 1) && (chr_labels(1) == 1) && ( all(ref.chrnumvec == 1) || (find(ref.chrnumvec == 1, 1, 'last' ) < find(ref.chrnumvec ~= 1, 1 )))
    chrlabel = chr_labels(1);
    BGMG_cpp.log('Reduce reference to a signle chromosome (chr%i).\n', chrlabel);
    trait1_data.zvec = trait1_data.zvec(ref.chrnumvec == chrlabel);
    trait1_data.nvec = trait1_data.nvec(ref.chrnumvec == chrlabel);
    defvec_tmp = defvec_tmp(ref.chrnumvec == chrlabel);
    ref.mafvec = ref.mafvec(ref.chrnumvec == chrlabel);
    ref.posvec = ref.posvec(ref.chrnumvec == chrlabel);
    ref.chrnumvec = ref.chrnumvec(ref.chrnumvec == chrlabel);  % must be the last statement as we replace ref.chrnumvec
    clear('chrlabel');
end
  
addpath('DERIVESTsuite');
addpath('PolyfitnTools');

if isfinite(hardprune_r2)
    if ~exist('hardprune_plink_ld_bin', 'var'), error('hardprune_plink_ld_bin is required'); end;
    defvec_tmp = BGMG_util.find_hardprune_indices(defvec_tmp, hardprune_r2, ref.mafvec, hardprune_plink_ld_bin, chr_labels);
end

% finalize defvec, from here it must not change.
defvec = defvec_tmp; clear('defvec_tmp', 'cur_defvec');
tag_indices = find(defvec);

BGMG_cpp.log('%i tag SNPs will go into fit and/or qq plots, etc\n', length(tag_indices));

bgmglib = BGMG_cpp();
bgmglib.dispose()
bgmglib.defvec = defvec;
bgmglib.chrnumvec = ref.chrnumvec;

bgmglib.set_option('r2min', r2min);
bgmglib.set_option('kmax', kmax);
bgmglib.set_option('max_causals', floor(max_causal_fraction * length(defvec)));
bgmglib.set_option('num_components', num_components);
bgmglib.set_option('cache_tag_r2sum', cache_tag_r2sum);
bgmglib.set_option('threads', THREADS);
if isfinite(SEED), bgmglib.set_option('seed', SEED); end;

bgmglib.mafvec = ref.mafvec;

for chr_index=1:length(chr_labels), 
    bgmglib.set_ld_r2_coo_from_file(strrep(plink_ld_bin,'@', sprintf('%i',chr_labels(chr_index)))); 
    bgmglib.set_ld_r2_csr(chr_labels(chr_index));
end;

bgmglib.set_weights_randprune(randprune_n, randprune_r2);

bgmglib.zvec1 = trait1_data.zvec(defvec);
bgmglib.nvec1 = trait1_data.nvec(defvec);

bgmglib.set_option('diag', 0);

% Preparation is done, BGMG library is fully setup. Now we can use it to
% calculate model QQ plots and univariate or bivariate cost function.

options = [];
options.total_het = sum(2*ref.mafvec.*(1-ref.mafvec));
options.verbose = true;
options.ci_alpha = CI_ALPHA;
options.title = TITLE;
options.fit_full_model = FIT_FULL_MODEL;

disp(options)

params = [];
% Load true params for simulated data. Do a quick fit to initialize sig2_zero 
if ~isempty(simu_params_file),
    tmp_params = load(simu_params_file);    trait_index = 1;
    params.univariate{trait_index}.pi_vec = tmp_params.causal_pi;
    params.univariate{trait_index}.sig2_beta = tmp_params.sigsq;

    bgmglib.set_option('fast_cost', 1); % ~FIT_FULL_MODEL);   % <- shortcut, initialize sig2zero based on fast model
    fitfunc = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), struct('Display', 'on', 'TolX', TolX, 'TolFun', TolFun)));
    fit_sig2_zero = fitfunc(struct('sig2_zero', 1), @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', params.univariate{trait_index}.pi_vec, 'sig2_beta', params.univariate{trait_index}.sig2_beta)));
    params.univariate{trait_index}.sig2_zero = fit_sig2_zero.sig2_zero;
    true_params = params;  % save, just in case

    BGMG_cpp.log('Params loaded from intput file (synthetic data).\n');
end;

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
end

result.trait1_file = trait1_file;
result.reference_file = reference_file;
result.options = options;
result.params = params;
if ~isempty(simu_params_file), result.true_params = true_params; end;

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.preliminary.mat'], 'result');
BGMG_cpp.log('Results saved to %s.preliminary.mat\n', out_file);

bgmglib.set_option('diag', 0);

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
