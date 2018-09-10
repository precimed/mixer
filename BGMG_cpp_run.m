% This script can be run univariate or bivariate mixture analysis from command line like this:
%
%   matlab -nodisplay -nosplash -nodesktop -r "trait1_file='PGC_SCZ_2014.mat'; trait2_file='LIPIDS_TG_2013.mat'; reference_file='1kG_phase3_EUR_11015883_reference_holland.mat'; BGMG_run; exit;"
%
% You may use this script with one trait (univariate analysis) or with two
% traits (bivariate analysis).
% The results are saved into a text file as well as into .mat file.
% Currently you may need to modify this file to pass additional
% parameters, later all parameters will be exposed as a text config file.
%
% For information about various parameters look into BGMG_fit.m.
%
% Some typical parameters (for my own machine)

% variables=who; for i=1:length(variables), eval(variables{i}); end

if 0
bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll';
bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';
plink_ld_bin = 'H:\work\hapgen_ldmat2_plink\bfile_merged_ldmat_p01_SNPwind50k_chr@.ld.bin'; chr_labels = 1; % 1:22;
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3_hardprune_p1.mat'};
%defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'};
%filename = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=1e-03_rep=5_tag1=customPolygenicOverlapAt0p375_tag2=evenPolygenicity';
filename = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=10_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity';
trait1_file = ['H:\work\simu_9pi_params\' filename '.trait1.mat']; trait1_nvec=100000;
trait2_file = ['H:\work\simu_9pi_params\' filename '.trait2.mat']; trait2_nvec=100000;
simu_params_file = ['H:\work\simu_9pi_params\' filename '.params.mat'];

kmax=1000;

reference_file = 'H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins_wld.mat';
DO_FIT_UGMG=true; DO_FIT_BGMG=true;
FIT_FULL_MODEL=true;
QQ_PLOT=false;STRATIFIED_QQ_PLOT=false;BGMG_LOGLIKE_PLOT=false;
cache_tag_r2sum=false;
MAF_THRESH=0.01;
r2min=0.05;
CI_ALPHA=0.05;
out_file = ['results_2018_09_06/' filename];
end

if ~exist('out_file', 'var'), out_file = 'BGMG_result'; end;
if exist('DO_FIT', 'var'), error('DO_FIT is deprecated; use DO_FIT_UGMG or DO_FIT_BGMG'); end;

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
if ~exist('MAF_THRESH', 'var'), MAF_THRESH = nan; end;
if ~exist('EXCLUDE_MHC', 'var'), EXCLUDE_MHC = true; end;

% plink_ld_mat is a file containing information about the LD structure. One
% file per chromosome. Each file should have index_A, index_B, r2
% variables. index_A and index_B corresponds to a zero-based index in the
% reference template. @ indicates the label of the chromosome.
if ~exist('plink_ld_bin', 'var'), error('plink_ld_bin is required'); end;
if ~exist('chr_labels', 'var'), chr_labels = 1:22; end;

if ~exist('trait1_file', 'var'), error('trait1_file is required'); end;
if ~exist('trait2_file', 'var'), trait2_file = ''; end;
if ~exist('trait1_nvec', 'var'), trait1_nvec = nan; end;
if ~exist('trait2_nvec', 'var'), trait2_nvec = nan; end;

if ~exist('randprune_n', 'var'), randprune_n = 16; end;
if ~exist('randprune_r2', 'var'), randprune_r2 = 0.8; end;
if ~exist('hardprune_r2', 'var'), hardprune_r2 = nan; end; % 1 iterations at this threshold to reduce set of SNPs (as in defvec)
if ~exist('hardprune_plink_ld_bin', 'var'), hardprune_plink_ld_bin = ''; end;

if ~exist('kmax', 'var'), kmax = 1000; end;
if ~exist('r2min', 'var'), r2min = 0.01; end;
if ~exist('max_causal_fraction', 'var'), max_causal_fraction = 0.008; end;
if ~exist('cache_tag_r2sum', 'var'), cache_tag_r2sum = 1; end;

% The following three options control how to get univariate & bivariate params
% BGMG_cpp_run can take previous results file, and enhance it with
% additional features. Typical example would be to take a file produced
% after DO_FIT, and use it to produce QQ plots, etc.
% Allows to fir BGMG with already fitted univaraite estimates
if ~exist('DO_FIT_UGMG', 'var'), DO_FIT_UGMG = true; end;               % perform univariate fit
if ~exist('DO_FIT_BGMG', 'var'), DO_FIT_BGMG = true; end;               % perform bivariate fit
if ~exist('simu_params_file', 'var'), simu_params_file = ''; end                  % load true parameters from 'simu' (simulation tool)
if ~exist('init_result_from_out_file', 'var'), init_result_from_out_file = ''; end;
if ~exist('init_trait1_from_out_file', 'var'), init_trait1_from_out_file = ''; end;
if ~exist('init_trait2_from_out_file', 'var'), init_trait2_from_out_file = ''; end;

% Setup tolerance for fminsearch
if ~exist('TolX', 'var'), TolX = 1e-2; end;
if ~exist('TolFun', 'var'), TolFun = 1e-2; end;

if ~exist('FIT_FULL_MODEL', 'var'), FIT_FULL_MODEL = true; end;                % use full model (when false, use gaussian approximation)
if ~exist('QQ_PLOT', 'var'), QQ_PLOT = false; end;   % make QQ plots
if ~exist('QQ_PLOT_DOWNSCALE', 'var'), QQ_PLOT_DOWNSCALE = 100; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('STRATIFIED_QQ_PLOT', 'var'), STRATIFIED_QQ_PLOT = false; end;
if ~exist('STRATIFIED_QQ_PLOT_DOWNSCALE', 'var'), STRATIFIED_QQ_PLOT_DOWNSCALE = 1000; end;     % downscale #snps in stratified QQ plots (model prediction only)
if ~exist('BGMG_LOGLIKE_PLOT', 'var'), BGMG_LOGLIKE_PLOT = false; end;
if ~exist('UGMG_LOGLIKE_PLOT', 'var'), UGMG_LOGLIKE_PLOT = false; end;
if ~exist('POWER_PLOT', 'var'), POWER_PLOT = false; end;  % make power plots with fitted parameters
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;
if ~exist('THREADS', 'var'), THREADS = -1; end;

% reference file containing mafvec, chrnumvec and posvec for all SNPs to consider in this analysis. 
if ~exist('reference_file', 'var'), error('reference_file is required'); end;
BGMG_cpp.log('Loading reference from %s... ', reference_file); ref = load(reference_file, 'mafvec', 'posvec', 'chrnumvec'); fprintf('OK.\n');
clear('defvec'); defvec_tmp = true(length(ref.mafvec), 1);
 
for i=1:length(defvec_files)
    BGMG_cpp.log('Loading %s... ', defvec_files{i}); cur_defvec = load(defvec_files{i}); defvec_tmp = defvec_tmp & cur_defvec.defvec; fprintf('OK.\n');
    BGMG_cpp.log('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));
end; clear('i');
if isfinite(MAF_THRESH), defvec_tmp = defvec_tmp & (ref.mafvec >= MAF_THRESH); BGMG_cpp.log('Exclude %i variants due to mafvec (%i variants remain)\n', sum(ref.mafvec < MAF_THRESH), sum(defvec_tmp)); end
if EXCLUDE_MHC, defvec_mhc = ~((ref.chrnumvec == 6) & (ref.posvec > 25e6) & (ref.posvec < 35e6)); defvec_tmp = defvec_tmp & defvec_mhc; BGMG_cpp.log('Exclude %i variants in MHC region (%i variants remain)\n', sum(~defvec_mhc), sum(defvec_tmp)); end

BGMG_cpp.log('Loading %s...', trait1_file); trait1_data = load(trait1_file); fprintf('OK.\n'); num_components = 1;
if (length(trait1_data.zvec) ~= length(ref.chrnumvec)), error('trait1_file is incompatible with the reference'); end;
if isfinite(trait1_nvec), trait1_data.nvec = ones(size(trait1_data.zvec)) * trait1_nvec; end;
if isfield(trait1_data, 'NCASE') && isfield(trait1_data, 'NCONTROL'), trait1_data.nvec = 4./(1./trait1_data.NCASE + 1./trait1_data.NCONTROL); end;
if ~isfield(trait1_data, 'nvec'), error('nvec is not available in trait1_file, and trait1_nvec parameter is not set'); end;
cur_defvec.defvec = isfinite(trait1_data.zvec + trait1_data.nvec); defvec_tmp = defvec_tmp & cur_defvec.defvec;
BGMG_cpp.log('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));

if ~isempty(trait2_file),
    BGMG_cpp.log('Loading %s...', trait2_file); trait2_data = load(trait2_file); fprintf('OK.\n'); num_components = 3;
    if (length(trait2_data.zvec) ~= length(ref.chrnumvec)), error('trait2_file is incompatible with the reference'); end;
    if isfinite(trait2_nvec), trait2_data.nvec = ones(size(trait2_data.zvec)) * trait2_nvec; end;
    if isfield(trait2_data, 'NCASE') && isfield(trait2_data, 'NCONTROL'), trait2_data.nvec = 4./(1./trait2_data.NCASE + 1./trait2_data.NCONTROL); end;
    if ~isfield(trait2_data, 'nvec'), error('nvec is not available in trait2_file, and trait2_nvec parameter is not set'); end;
    cur_defvec.defvec = isfinite(trait2_data.zvec + trait2_data.nvec); defvec_tmp = defvec_tmp & cur_defvec.defvec;
    BGMG_cpp.log('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));
end

if (length(chr_labels) == 1) && (chr_labels(1) == 1) && ( all(ref.chrnumvec == 1) || (find(ref.chrnumvec == 1, 1, 'last' ) < find(ref.chrnumvec ~= 1, 1 )))
    chrlabel = chr_labels(1);
    BGMG_cpp.log('Reduce reference to a signle chromosome (chr%i).\n', chrlabel);
    trait1_data.zvec = trait1_data.zvec(ref.chrnumvec == chrlabel);
    trait1_data.nvec = trait1_data.nvec(ref.chrnumvec == chrlabel);
    if ~isempty(trait2_file),
        trait2_data.zvec = trait2_data.zvec(ref.chrnumvec == chrlabel);
        trait2_data.nvec = trait2_data.nvec(ref.chrnumvec == chrlabel);
    end
    defvec_tmp = defvec_tmp(ref.chrnumvec == chrlabel);
    ref.mafvec = ref.mafvec(ref.chrnumvec == chrlabel);
    ref.posvec = ref.posvec(ref.chrnumvec == chrlabel);
    ref.chrnumvec = ref.chrnumvec(ref.chrnumvec == chrlabel);  % must be the last statement as we replace ref.chrnumvec
    clear('chrlabel');
end
  
addpath('DERIVESTsuite');
addpath('PolyfitnTools');

if isfinite(hardprune_r2)
    if ~exist('hardprune_plink_ld_mat', 'var'), error('randprune_r2_plink_ld_mat is required'); end;
    defvec_tmp = BGMG_util.find_hardprune_indices(defvec_tmp, hardprune_r2, ref.mafvec, ref.chrnumvec, hardprune_plink_ld_bin, chr_labels);
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

bgmglib.mafvec = ref.mafvec;

for chr_index=1:length(chr_labels),
    bgmglib.set_ld_r2_coo_from_file(strrep(plink_ld_bin,'@', sprintf('%i',chr_labels(chr_index)))); 
    bgmglib.set_ld_r2_csr(chr_labels(chr_index));
end;

bgmglib.set_weights_randprune(randprune_n, randprune_r2);

bgmglib.zvec1 = trait1_data.zvec(defvec);
bgmglib.nvec1 = trait1_data.nvec(defvec);

if ~isempty(trait2_file),
    bgmglib.zvec2 = trait2_data.zvec(defvec);
    bgmglib.nvec2 = trait2_data.nvec(defvec);
end

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
if ~isempty(simu_params_file),
    tmp_params = load(simu_params_file);
    tmp_params.sig2_zero = [1;1];
    % Make sure that first component is specific to trait1, second
    % component - specific to trait2, and third is pleiotropic component
    idx = tmp_params.pi_vec_trait1 ~= 0 & tmp_params.pi_vec_trait2 == 0; assert(sum(idx) <= 1);
    tmp_params.pi_vec_tmp(1) = sum(tmp_params.pi_vec_trait1(idx));
    tmp_params.rho_beta_tmp(1) = sum(tmp_params.rho_beta(idx));
    
    idx = tmp_params.pi_vec_trait2 ~= 0 & tmp_params.pi_vec_trait1 == 0; assert(sum(idx) <= 1);
    tmp_params.pi_vec_tmp(2) = sum(tmp_params.pi_vec_trait2(idx));
    tmp_params.rho_beta_tmp(2) = sum(tmp_params.rho_beta(idx));
    
    idx = tmp_params.pi_vec_trait2 ~= 0 & tmp_params.pi_vec_trait1 ~= 0; assert(sum(idx) <= 1);
    tmp_params.pi_vec_tmp(3) = sum(tmp_params.pi_vec_trait2(idx));
    tmp_params.rho_beta_tmp(3) = sum(tmp_params.rho_beta(idx));
    
    tmp_params.pi_vec = tmp_params.pi_vec_tmp;
    tmp_params.rho_beta = tmp_params.rho_beta_tmp;
    tmp_params = rmfield(tmp_params, {'pi_vec_trait2', 'pi_vec_trait1', 'pi_vec_tmp', 'rho_beta_tmp'});
    
    if all(size(tmp_params.sig2_beta) == [1 2])
        tmp_params.sig2_beta = [tmp_params.sig2_beta(1) 0 tmp_params.sig2_beta(1); 0 tmp_params.sig2_beta(2) tmp_params.sig2_beta(2)];
    end
    
    if (tmp_params.pi_vec(end)<1e-15), tmp_params.pi_vec(end) = prod(tmp_params.pi_vec(1:end-1)); end;

    true_params = [];
    true_params.univariate{1} = struct('sig2_zero', tmp_params.sig2_zero(1), 'pi_vec', sum(tmp_params.pi_vec([1 3])), 'sig2_beta', tmp_params.sig2_beta(1,1));
    true_params.univariate{2} = struct('sig2_zero', tmp_params.sig2_zero(2), 'pi_vec', sum(tmp_params.pi_vec([2 3])), 'sig2_beta', tmp_params.sig2_beta(2,2));
    true_params.bivariate = tmp_params;
    
    % Do a quick fit to initialize sig2_zero and rho_zero
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);
     for trait_index = 1:2
         fitfunc = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), struct('Display', 'on', 'TolX', TolX, 'TolFun', TolFun)));
         fit_sig2_zero  = fitfunc(struct('sig2_zero', 1), @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', true_params.univariate{trait_index}.pi_vec, 'sig2_beta', true_params.univariate{trait_index}.sig2_beta)));
         true_params.univariate{trait_index}.sig2_zero = fit_sig2_zero.sig2_zero;
     end
     
     zmattmp = [bgmglib.zvec1, bgmglib.zvec2];
     corrmat = corr(zmattmp(all(~isnan(zmattmp), 2), :));
     fitfunc = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), struct('Display', 'on', 'TolX', TolX, 'TolFun', TolFun)));
     fit_rho_zero = fitfunc(...
         struct('rho_zero', corrmat(1,2)), ...
         @(x)BGMG_util.BGMG_mapparams3(x, struct(...
            'sig2_zero', [true_params.univariate{1}.sig2_zero, true_params.univariate{2}.sig2_zero], ... 
            'sig2_beta', true_params.bivariate.sig2_beta, ...
            'rho_beta', true_params.bivariate.rho_beta, ...
            'pi_vec', true_params.bivariate.pi_vec)));
    true_params.bivariate.rho_zero = fit_rho_zero.rho_zero;
    true_params.bivariate.sig2_zero = [true_params.univariate{1}.sig2_zero, true_params.univariate{2}.sig2_zero];

    BGMG_cpp.log('Params loaded from intput file (synthetic data)\n');
    params = true_params;
end;

if ~isempty(init_trait1_from_out_file) || ~isempty(init_trait2_from_out_file)
    if isempty(trait2_file), error('init_traitN_from_out_file is specific for BGMG fit'); end;
    if isempty(init_trait1_from_out_file) || isempty(init_trait2_from_out_file), error('init_traitN_from_out_file should be specified for both traits (N=1 and N=2)'); end;
    if ~isempty(init_result_from_out_file), error('init_result_from_out_file is incompatible with init_traitN_from_out_file'); end;

    BGMG_cpp.log('Loading existing params for trait1 from %s...\n', init_trait1_from_out_file);
    tmp = load(init_trait1_from_out_file); params.univariate{1} = tmp.univariate{1}.params;
    if length(tmp.univariate) ~= 1 || isfield(tmp, 'bivariate'), error('init_trait1_from_out_file does not appear to be a result of univariate analysis'); end;
    BGMG_cpp.log('Loading existing params for trait2 from %s...\n', init_trait2_from_out_file);
    tmp = load(init_trait2_from_out_file); params.univariate{2} = tmp.univariate{1}.params;
    if length(tmp.univariate) ~= 1 || isfield(tmp, 'bivariate'), error('init_trait2_from_out_file does not appear to be a result of univariate analysis'); end;
    clear('tmp');
    BGMG_cpp.log('Univariate params load from initial file');
end

if ~isempty(init_result_from_out_file)
    BGMG_cpp.log('Loading init file from %s...\n', init_result_from_out_file);
    tmp = load(init_result_from_out_file);
    params.univariate{1} = tmp.result.univariate{1}.params;
    params.univariate{2} = tmp.result.univariate{2}.params;
    params.bivariate = tmp.result.bivariate.params;
    BGMG_cpp.log('Univariate and bivariate params load from initial file');
    clear('tmp');
end

% if there is no initial approximation, setup it from fast model
if ~isfield(params, 'univariate')
    if ~isempty(trait1_file), params.univariate{1} = BGMG_cpp_fit_univariate_fast(1); end
    if ~isempty(trait2_file), params.univariate{2} = BGMG_cpp_fit_univariate_fast(2); end
end
if ~isfield(params, 'bivariate') && ~isempty(trait2_file)
    params.bivariate = BGMG_cpp_fit_bivariate_fast(params.univariate{1}, params.univariate{2});
end

% Fit bivariate or univariate model to the data
bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);

result = [];
if DO_FIT_UGMG
    result.univariate{1} = BGMG_cpp_fit_univariate(1, params.univariate{1}, options);
    params.univariate{1} = result.univariate{1}.params;
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);  % restore fast_cost flag as BGMG_cpp_fit_univariate may change it (CI intervals are based on fast cost function)
    if ~isempty(trait2_file),
        result.univariate{2} = BGMG_cpp_fit_univariate(2, params.univariate{2}, options);
        params.univariate{2} = result.univariate{2}.params;
        bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);  % restore fast_cost flag as BGMG_cpp_fit_univariate may change it (CI intervals are based on fast cost function)
        
        % Update initial approximation for bivariate fit (sig2_zero, sig2_beta)
        params.bivariate.sig2_zero = [params.univariate{1}.sig2_zero; params.univariate{2}.sig2_zero];
        params.bivariate.sig2_beta = [params.univariate{1}.sig2_beta 0 params.univariate{1}.sig2_beta; ...
                                      0 params.univariate{2}.sig2_beta params.univariate{2}.sig2_beta];
                                  
        % Also, update pi_vec. Preserve univariate estimates, and set
        % pi12frac to what it was in params.bivariate (e.i. either what we
        % load, or what we fit with BGMG_cpp_fit_bivariate_fast
        tmp.pi_vec = params.bivariate.pi_vec;
        tmp.pi_frac = tmp.pi_vec(3) / min(sum(tmp.pi_vec([1, 3])), sum(tmp.pi_vec([2, 3])));
        tmp.pi_vec = [params.univariate{1}.pi_vec, params.univariate{2}.pi_vec];
        tmp.pi12 = min(tmp.pi_vec) * tmp.pi_frac;
        params.bivariate.pi_vec = [tmp.pi_vec(1) - tmp.pi12, tmp.pi_vec(2) - tmp.pi12, tmp.pi12];
        clear('tmp');
    end
end

if DO_FIT_BGMG
    result.bivariate = BGMG_cpp_fit_bivariate(params.bivariate, options);
    params.bivariate = result.bivariate.params;
end

result.trait1_file = trait1_file;
result.trait2_file = trait2_file;
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
    for trait_index = 1:(1 + ~isempty(trait2_file))
        options.downscale = QQ_PLOT_DOWNSCALE;
        figures.tot = figure;
        plot_data = BGMG_cpp_qq_plot(params.univariate{trait_index}, trait_index, options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_true_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.pdf', out_file, trait_index), '-dpdf')
    end
end

if STRATIFIED_QQ_PLOT && ~isempty(trait2_file)
    options.downscale = STRATIFIED_QQ_PLOT_DOWNSCALE;
    [figures, plot_data] = BGMG_cpp_stratified_qq_plot(params.bivariate, options);
    result.bivariate.stratified_qq_plot_fit_data = plot_data;
    print(figures.tot, sprintf('%s.stratqq.pdf', out_file), '-dpdf')
end

if UGMG_LOGLIKE_PLOT
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);
    for trait_index = 1:(1 + ~isempty(trait2_file))
        bgmglib.clear_loglike_cache();
        x0 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(result.params.univariate{trait_index});
        x_range = abs(x0)*0.0 + 1; x_low = x0 - x_range; x_up  = x0 + x_range;
        arg_index=3; %1:length(x0)
        xgrid = linspace(x_low(arg_index), x_up(arg_index), 20)';
        for xi = 1:length(xgrid)
            x=x0; x(arg_index)=xgrid(xi);
            BGMG_util.UGMG_fminsearch_cost(BGMG_util.UGMG_mapparams1_decorrelated_parametrization(x), trait_index);
        end
        ut = bgmglib.extract_univariate_loglike_trajectory()
        ut.cost = ut.cost-min(ut.cost);
        figure(10+trait_index);
        
        def = isfinite(xgrid+ut.cost);
        x = xgrid(def);
        y = ut.cost(def);
        [curve2, gof2] = fit(x, y, 'poly2');
        [curve3, gof3] = fit(x, y, 'poly3');
        fprintf('%.3f vs %.3f\n', gof2.rsquare, gof3.rsquare)

        plot(x, [(curve2(x)-y).^2, (curve3(x)-y).^2], '.')
        x0opt = fminsearch(@(x)curve3(x), x0(arg_index));

        clf; 
        subplot(1,2,1); plot(ut.pivec, [ut.cost, curve2(xgrid), curve3(xgrid)], '.-');
        subplot(1,2,2); plot(xgrid, [ut.cost, curve2(xgrid), curve3(xgrid)], '.-'); hold on; plot(x0opt, curve3(x0opt), '*k');
    end
end

if BGMG_LOGLIKE_PLOT && ~isempty(trait2_file)
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);
    bgmglib.clear_loglike_cache();

    x0 = BGMG_util.BGMG_mapparams3_decorrelated_parametrization(result.params.bivariate);
    x_range = 2; x_low = x0 - x_range; x_up  = x0 + x_range;
    arg_index=9;
    xgrid = linspace(x_low(arg_index), x_up(arg_index), 20)';
    for xi = 1:length(xgrid)
        x=x0; x(arg_index)=xgrid(xi);
        cost = BGMG_util.BGMG_fminsearch_cost(BGMG_util.BGMG_mapparams3_decorrelated_parametrization(x));
    end
    bt = bgmglib.extract_bivariate_loglike_trajectory();
    
    def = isfinite(xgrid+bt.cost);
    x = xgrid(def);
    y = bt.cost(def);
    curve3 = fit(x, y, 'poly3');
    x0opt = fminsearch(@(x)curve3(x), x0(arg_index));
        
    x0(arg_index) = x0opt;
    %result.params = options.mapparams(x0);
        
    figure(20);
    subplot(1,2,1); plot(bt.pivec(:, 3) ./ min(sum(bt.pivec(:, [1 3]), 2), sum(bt.pivec(:, [2 3]), 2)), [bt.cost, curve3(xgrid)])
    subplot(1,2,2); plot(xgrid, [bt.cost, curve3(xgrid)])
    
    %result.bivariate.loglike_plot_data_fit = plots_data;
    %result.bivariate.loglike_plot_trajectory_fit = bgmglib.extract_bivariate_loglike_trajectory();
    %print(figures.tot, sprintf('%s.loglike.pdf', out_file), '-dpdf')
end

bgmglib.set_option('diag', 0);

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.mat'], 'result');
BGMG_cpp.log('Results saved to %s.mat\n', out_file);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper code below 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    % Helper code to save all results to a text file
    % SIMU_BGMG2_2018_06_17
    dirs=dir('H:\work\SIMU_BGMG_9pifrac_20000\*.bgmg.mat');
    fileID = fopen(['H:\work\SIMU_BGMG_9pifrac_20000\SIMU_BGMG_9pifrac_20000.csv'], 'w');
    has_header = false;
    for i=1:length(dirs)
        try
        x = load(fullfile('H:\work\SIMU_BGMG_9pifrac_20000', dirs(i).name));
        if ~isfield(x.result, 'options'), x.result.options = []; end;
        if ~isfield(x.result.bivariate, 'params'), continue; end;
        catch
            fprintf('error: %s\n', dirs(i).name);
            continue
        end
        [header, data] = BGMG_util.result2str_point_estimates(x.result, x.result.options);
        if ~has_header
            fprintf(fileID, 'trait1_file\ttrait2_file\tbgmg_mat_file\treference_file\t%s\n', header);
            has_header=1;
        end
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\n', x.result.trait1_file, x.result.trait2_file, dirs(i).name, x.result.reference_file, data);
    end

    if fileID ~= 2, fclose(fileID); end;
end

if 0
    pBuffer = libpointer('singlePtr', zeros(sum(defvec), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', 0, sum(defvec), pBuffer);  check(); weights_bgmg = pBuffer.Value;
    calllib('bgmg', 'bgmg_retrieve_ld_tag_r2_sum', 0, sum(defvec), pBuffer);  check(); bgmg_sum_r2 = pBuffer.Value;
    calllib('bgmg', 'bgmg_retrieve_ld_tag_r4_sum', 0, sum(defvec), pBuffer);  check(); bgmg_sum_r4 = pBuffer.Value;
    save('test_bgmg_sum_r2_sum_r4.mat', 'defvec', 'bgmg_sum_r2', 'bgmg_sum_r4');
    ref_hapgen_100k = load('H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins.mat');
    sum_r2 = sum(ref_hapgen_100k.sum_r2_biased(defvec, :), 2);
    sum_r4 = sum(ref_hapgen_100k.sum_r4_biased(defvec, :), 2);
    clear pBuffer
end

if 0
    
pi1u = '1e-02';
files=    dir(['H:\work\SIMU_BGMG2_loglike_2018_06_09\*pi1u=', pi1u, '*rep=2*.bgmg.mat'])
for i=1:length(files), files(i).dir = 'H:\work\SIMU_BGMG2_loglike_2018_06_09\'; end

files2=    dir(['H:\work\SIMU_BGMG2_loglike_2018_06_09_kmax1000\to_plot\*pi1u=',pi1u,'*'])
for i=1:length(files2), files2(i).dir = 'H:\work\SIMU_BGMG2_loglike_2018_06_09_kmax1000\to_plot\'; end

files = [files; files2];

f = figure;
legends = {};
for i=1:length(files)
    fprintf('%s\n', files(i).name)
    try
        x=  load([files(i).dir, files(i).name]);
        %if ~isempty(findstr(files(i).name, 'simu_h2=0.4_rg=0.0_pi1u=1e-04_pi2u=1e-04_pi12=3e-05_rep=1_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity.cpp.bgmg.mat')) continue; end;
        %x.result.bivariate.params;
        x.result.bivariate.loglike_plot_fit_data
    catch
        continue
    end
    
    if isfield(x.result.bivariate, 'params')
        params=x.result.bivariate.params;
        if any(params.sig2_zero>1.2), warning('non-converged'); continue; end;
    end
    plots_data = x.result.bivariate.loglike_plot_fit_data;
    
        
    fn = files(i).name;
    type = find([~isempty(findstr(fn, 'random')) ~isempty(findstr(fn, 'partial')) ~isempty(findstr(fn, 'complete'))]);
    if ~isempty(findstr(fn, 'kmax1000')) plot_symbol='-', else plot_symbol='.-'; end;
    
    
    fnames = fieldnames(plots_data);
    for fni=1:length(fnames)
        subplot(3,3,fni);hold on; 
         ax = gca;ax.ColorOrderIndex = type;
   
        plot_name=fnames{fni};
        y = plots_data.(plot_name).y; y(y>1e99)=nan;
        plot(plots_data.(plot_name).x, y-min(y), [plot_symbol]);
        plot_name(plot_name=='_')=' ';title(plot_name);
    end
    leg = files(i).name; leg(leg=='_')=' '; legends{end+1, 1} = leg;
    subplot(3,3,8:9);  ax = gca;ax.ColorOrderIndex = type; hold on; plot(0,0,plot_symbol);
end
    subplot(3,3,8:9); legend(legends);
    set(f,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(f, sprintf('loglike_pi1u=%s_joint',pi1u), '-dpdf')
end

% TBD: re-test confidence intervals estimation


if 0
    % find centromere location
    % 1:121535434- 124535434
        sum((ref.chrnumvec==1) & (ref.posvec> 121535434) & (ref.posvec < 124535434))
        sum((ref.chrnumvec==2) & (ref.posvec> 92326171) & (ref.posvec < 95326171))

    
    figure(1);clf;hold on;
    for chri=1:22
        p = ref.posvec((ref.chrnumvec==chri));  % defvec == hm3 as defined by LDSR
        diff = p(2:end) - p(1:end-1);
        if sum(diff>5e5)==0
            continue;
        end
        subplot(6,4,chri); plot(diff); title(sprintf('chr%i', chri));
        for i=find(diff>1e6)'
            fprintf('chr=%i %i - %i (%.1f MB)\n', chri, p(i), p(i+1), single(diff(i))/1e6 )
        end
    end
end


if 0
plink_ld_mat = 'H:\work\hapgen_ldmat2_plink\bfile_merged_10K_ldmat_p01_SNPwind50k_chr@.ld.mat'; chr_labels = 1:22;
for chr_index=length(chr_labels):-1:1
    plink_ld_mat_chr = strrep(plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
    fprintf('Loading %s...', plink_ld_mat_chr); tmp = load(plink_ld_mat_chr); fprintf('OK.\n');
    tmp.index_A = tmp.index_A(tmp.r2 >= 0.1);
    tmp.index_B = tmp.index_B(tmp.r2 >= 0.1);
    tmp.r2      = tmp.r2     (tmp.r2 >= 0.1);
    plink_ld_mat_chr_90 = strrep('H:\work\hapgen_ldmat2_plink\bfile_merged_10K_ldmat_p10_SNPwind50k_chr@.ld.mat', '@', sprintf('%i',chr_labels(chr_index)));
    save(plink_ld_mat_chr_90, '-struct', 'tmp', '-v7');
end
end

if 0
    % quick and dirty fminsearch
    p0 = result.univariate{1}.params;
    x0 = [p0.pi_vec, p0.sig2_zero, p0.sig2_beta];
  
    fminsearch_options = struct('Display', 'iter');
    mapparams = @(x)struct('pi_vec', x(1), 'sig2_zero', x(2), 'sig2_beta', x(3));
    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', 1); check();
    UGMG_fminsearch_cost = @(ov)calllib('bgmg', 'bgmg_calc_univariate_cost_with_deriv', 0, trait_index, ov.pi_vec, ov.sig2_zero, ov.sig2_beta, 3, pBuffer);
    UGMG_fminsearch_cost2 = @(ov)calllib('bgmg', 'bgmg_calc_univariate_cost', 0, trait_index, ov.pi_vec, ov.sig2_zero, ov.sig2_beta);
    fit = @(x0)fminsearch(@(x)UGMG_fminsearch_cost(mapparams(x)), x0, fminsearch_options);
    
    tic;UGMG_fminsearch_cost(p0);toc
    pBuffer.Value
    tic;UGMG_fminsearch_cost2(p0);toc
    
    fit(x0)

    options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','off');
    fminunc(@(x)BGMG_util.UGMG_fminsearch_cost_with_gradient(mapparams(x)), x0, options)
    
    fminsearch_options = struct('Display', 'off'); 
    fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x)), x0, fminsearch_options);
    
    clear pBuffer

end

if 0
    options.params0 = struct('sig2_zero', 1, 'pi_vec', sum(trait1_data.causal_pi), 'sig2_beta', trait1_data.sigsq);
    result.univariate{1} = BGMG_cpp_fit_univariate(trait1_data.zvec, trait1_data.nvec, options);
end

if 0
    params = struct('sig2_zero', 1, 'pi_vec', sum(trait1_data.causal_pi), 'sig2_beta', trait1_data.sigsq);
    BGMG_cpp_fit_kl_univariate(trait1_data.zvec(defvec), trait1_data.nvec(defvec), struct('params0', params));
    
    % Retrieve weights from c++ library
    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
    pBuffer = libpointer('singlePtr', zeros(sum(defvec), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', 0, sum(defvec), pBuffer);  check(); 
    weights_bgmg = pBuffer.Value;
    clear pBuffer
    
    data_weights = weights_bgmg; data_weights = data_weights/sum(data_weights);
    model_weights = weights_bgmg;
    
    % fit KL divergence between model & data pdf
    % Calculate data_logpvec
    zvec = trait1_data.zvec(defvec);
    hv_z = linspace(0, min(max(abs(zvec)), 38.0), 10000);
    [data_y, si] = sort(-log10(2*normcdf(-abs(zvec))));
    data_x=-log10(cumsum(data_weights(si),1,'reverse'));
    data_idx = ([data_y(2:end); +Inf] ~= data_y);
    hv_logp = -log10(2*normcdf(-hv_z));
    data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

    % Calculate model_logpvec
    
    zgrid = single(0:0.05:15); 
    pBuffer = libpointer('singlePtr', zeros(length(zgrid), 1, 'single'));
    calllib('bgmg', 'bgmg_calc_univariate_pdf', 0, trait_index, params.pi_vec, params.sig2_zero, params.sig2_beta, length(zgrid), zgrid, pBuffer);  check(); 
    pdf = pBuffer.Value'; clear pBuffer
    pdf = pdf / sum(model_weights);
    if (zgrid(1) == 0), zgrid = [-fliplr(zgrid(2:end)) zgrid];pdf = [fliplr(pdf(2:end)) pdf]; end
    model_cdf = cumsum(pdf)  * (zgrid(2) - zgrid(1)) ;
    X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
    model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
    model_logpvec = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), hv_z')); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)

    clf; hold on
    hData  = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
    hModel = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;
 
end


if 0
    % re-save LD info in C++ format
    plink_ld_mat='1000Genome_ldmat_p05_SNPwind50k_chr@.ld.mat';chr_labels=1:22;
    for chr_index=1:length(chr_labels)
        plink_ld_mat_chr = strrep(plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
        fprintf('Loading %s...', plink_ld_mat_chr); tmp = load(plink_ld_mat_chr); fprintf('OK.\n');
        fileID = fopen([plink_ld_mat_chr(1:end-4) '.bin'],'w');
        fwrite(fileID,uint64(length(tmp.index_A)),'uint64');
        fwrite(fileID,int32(tmp.index_A),'int32');
        fwrite(fileID,int32(tmp.index_B),'int32');
        fwrite(fileID,single(tmp.r2),'single');
        fclose(fileID);
    end
end

if 0
    % new plots for loglike trajectory - a as a function of the same params
    % as we use in fminsearch
    t = result.bivariate.loglike_plot_trajectory_fit;
    t = result.bivariate.loglike_fit_trajectory;
    p0 = result.bivariate.params;
    func_map_params = @(x)BGMG_util.BGMG_mapparams3(x, struct('sig2_zero', p0.sig2_zero, 'sig2_beta', p0.sig2_beta(:, end), 'rho_zero', p0.rho_zero));
    
    x0 = func_map_params(p0)
    x = []
    for i=1:length(t.cost)
        x = [x; func_map_params(struct('pi_vec', t.pivec(i, :), 'rho_beta', [0 0 t.rho_beta(i)]))];
    end
    
    idx = t.cost < quantile(t.cost, 0.5);
    
    figure(1); clf;
    [~, x0id] = min(sum((x  - repmat(x0, [size(x, 1), 1])).^2, 2));
    for param_index=1:4
        xu = unique(x(idx, param_index))
        cu = nan(size(xu));
        for i=1:length(xu)
            cu(i) = min(t.cost(x(:, param_index) == xu(i)));
        end
        subplot(2,2,param_index);  hold on;
        %plot(xu, cu, '.');
        plot(x0(param_index), t.cost(x0id), '*');
        plot(x(idx, param_index), t.cost(idx), '.');
    end
    
    
    polymodel = polyfitn(x(idx, :),t.cost(idx), 2)
    cost2 = polyvaln(polymodel, x);
    plot(t.cost(idx), cost2(idx), '.')
    optfunc = @(xarg)polyvaln(polymodel, xarg);
    fminsearch(optfunc, x0, struct('Display', 'iter'))

    
    
end

if 0 
    %datafolder = 'H:\\GitHub\\BGMG\\tes2017_07_10';
    bgmglib.verbose=false;
    datafolder = '/home/oleksandr/precimed/BGMG/test_cost_func';
    outtag_vec= {'fitFromTrue'; 'fitFromScratch'}
    for repi=1:10
        trait1_data = load(fullfile(datafolder, sprintf('simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=%i_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity.trait1.mat', repi)));
        trait2_data = load(fullfile(datafolder, sprintf('simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=%i_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity.trait2.mat', repi)));
        bgmglib.zvec1=trait1_data.zvec(bgmglib.defvec);
        bgmglib.zvec2=trait2_data.zvec(bgmglib.defvec);

    for outtagi=1:2;
        outtag=outtag_vec{outtagi};
        try
            load(fullfile(datafolder, sprintf('simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=%i_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity_outtag=%s.hardprune.bgmg.mat', repi,outtag)))
        catch
            fprintf('\n');
            continue
        end
    ov = result.params.bivariate;
    cost_recalculated = bgmglib.calc_bivariate_cost(ov.pi_vec, ov.sig2_beta(:, 3), ov.rho_beta(3), ov.sig2_zero, ov.rho_zero);

    %pi_vec = result.bivariate.params.pi_vec;
    %fprintf('%s%i: pi_vec=%.3e, %.3e, %.3e\n', strpad(outtag, 15, 'post', ' '), repi, pi_vec(1), pi_vec(2), pi_vec(3) );
    [~, id] = min(result.bivariate.loglike_fit_trajectory.cost);
    pi_vec2=result.bivariate.loglike_fit_trajectory.pivec(id, :);
    strpad=@(x)[x repmat(' ', [1, 15-length(x)])];
    fprintf('%s%i: pi_vec=%.3e, %.3e, %.3e, cost=%.8e, cost_recalculated=%.8e, diff=%.4f \n', strpad(outtag), repi, pi_vec2(1), pi_vec2(2), pi_vec2(3), result.bivariate.loglike_fit_trajectory.cost(id), cost_recalculated, cost_recalculated-result.bivariate.loglike_fit_trajectory.cost(id) );

    plotfile = sprintf('%s_%s_v2_simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=%i_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity.trait1.mat', out_file, outtag, repi)
    for trait_index = 1:2
        options.downscale = 10;
        [figures, plot_data] = BGMG_cpp_qq_plot(result.univariate{trait_index}.params, trait_index, options);
        print(figures.tot, sprintf('%s.trait%i.qq.pdf', plotfile, trait_index), '-dpdf')
        save(sprintf('%s.trait%i.qq.mat', plotfile, trait_index), 'plot_data', 'result');
    end

    options.downscale = 1000;
    [figures, plot_data] = BGMG_cpp_stratified_qq_plot(result.params.bivariate, options);
    print(figures.tot, sprintf('%s.stratqq.pdf', plotfile), '-dpdf')
    save(sprintf('%s.stratqq.mat', plotfile), 'plot_data', 'result')
    end
    end
end

if 0 
    % summarize results of fit-from-scratch runs
    folders = {'SIMU_BGMG_9pifrac_withSampleOverlap50k', 'SIMU_BGMG_9pifrac_withoutSampleOverlap'};
    for folder_index=1:length(folders)
    for j=1:4
    figure(folder_index);
    %load(fullfile('H:\work\SIMU_BGMG_9pifrac_run3', sprintf('simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=%i_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity_outtag=fitFromScratch.kmax5000.hardprune.run3.bgmg.mat', j)))
    load(fullfile('H:\work', folders{folder_index}, sprintf('simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=%i_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity_outtag=fitFromScratch.kmax5000.hardprune.run4.bgmg.mat', j)))
    for i=1:2
    r=result.univariate{i};f=r.loglike_adj_trajectory;
    subplot(4,3,i+3*(j-1)); plot(f.pivec, f.cost, '.-', r.params.pi_vec, min(f.cost), '*')
     xlabel(sprintf('pi%iu', i)); 
     if i==1, ylabel(sprintf('iter %i', j)); end;
    end
    r=result.bivariate;f=r.loglike_adj_trajectory;
    pifrac = f.pivec(:, 3) ./ (min(f.pivec(:, 1:2)')' + f.pivec(:, 3));
    pifrac0 = r.params.pi_vec(3)/(min(r.params.pi_vec([1:2])) + r.params.pi_vec(3));    
    subplot(4,3,3+3*(j-1)); plot(pifrac, f.cost, '.-', pifrac0, min(f.cost), '*');
     xlabel('pifrac = pi12/min(pi1u, pi2u)');
     %xlim([0 0.25]);
    end
    f=figure(folder_index);
    set(f,'PaperOrientation','portrait','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(f, fullfile('H:\Dropbox\analysis\2018_07_07_BGMG_fit_procedure',sprintf('loglike_pi1u=3e-3_bgmg_fromScratch_%s.pdf', folders{folder_index})), '-dpdf')
    end
end