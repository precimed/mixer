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
plink_ld_bin = 'H:\work\hapgen_ldmat2_plink\bfile_merged_ldmat_p01_SNPwind50k_chr@.ld.bin'; chr_labels = 1; %1:22;
%randprune_r2_plink_ld_mat = ''; randprune_r2_defvec_threshold = nan;
%randprune_r2_plink_ld_mat = 'H:\work\hapgen_ldmat2_plink\bfile_merged_ldmat_p01_SNPwind50k_chr@.ld.mat';
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'};
%defvec_files = {};
%defvec_files = {'H:\Dropbox\shared\BGMG\defvecs\defvec_11m_infinium-omniexpress-24-v1-3-a1.mat', ...
%                'H:\Dropbox\shared\BGMG\defvecs\defvec_highld_regions.mat', ... 
%                'H:\Dropbox\shared\BGMG\defvecs\centromere_plus_3cm_locations.mat' };

filename = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=1e-03_rep=5_tag1=customPolygenicOverlapAt0p375_tag2=evenPolygenicity'
%filename = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=0e+00_rep=10_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity'
trait1_file = ['H:\work\simu_9pi_params\' filename '.trait1.mat']; trait1_nvec=100000;
trait2_file = ['H:\work\simu_9pi_params\' filename '.trait2.mat']; trait2_nvec=100000;
params_file = ['H:\work\simu_9pi_params\' filename '.params.mat'];

kmax=1000;

%trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
%trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait2_nvec=100000;
%trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=3e-03_rep=10_tag1=completePolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
%trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=3e-03_rep=10_tag1=completePolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait1_nvec=100000;
%trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=8e-04_rep=10_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
%trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=8e-04_rep=10_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait1_nvec=100000;

reference_file = 'H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins_wld.mat';
DO_FIT=false;FIT_FULL_MODEL=false;STRATIFIED_QQ_PLOT_FIT=false;QQ_PLOT_TRUE=true;LOGLIKE_PLOT_TRUE=false;QQ_PLOT_FIT=false;cache_tag_r2sum=true;
MAF_THRESH=0.01;
out_file = 'tmptesting6';
init_result_from_out_file = ''; %'tmptesting5';
%out_file = 'BGMG_random_overlap_chr1_pi1u=3e-03';
%out_file = 'BGMG_full_overlap_chr1_pi1u=3e-03';

bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\Debug\bgmg.dll';
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'}
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_1kG_phase3_EUR.mat' };
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_1kG_phase3_EUR.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'}
trait1_file = 'H:\Dropbox\shared\BGMG\p_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc\simu_h2=0.40_pi1u=1e-3_rep=1.ugmg.mat'
trait1_file = 'H:\GitHub\BGMG\GIANT_HEIGHT_2014_lift.mat';
trait1_file = 'H:\GitHub\BGMG\PGC_SCZ_2014.mat';
reference_file = 'H:\Dropbox\shared\BGMG\1kG_phase3_EUR_11015883_reference_p01_SNPwind50k_per_allele_4bins.mat'
DO_FIT = false; QQ_PLOT_TRUE = false;
trait1_file  = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat';
trait2_file  = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat';
end

addpath('DERIVESTsuite');

if ~exist('out_file', 'var'), out_file = 'BGMG_result'; end;

% BGMG_cpp_run can take previous results file, and enhance it with
% additional features. Typical example would be to take a file produced
% after DO_FIT, and use it to produce QQ plots, etc.
if ~exist('init_result_from_out_file', 'var'), init_result_from_out_file = ''; end;

% Allows to fir BGMG with already fitted univaraite estimates
if ~exist('init_trait1_from_out_file', 'var'), init_trait1_from_out_file = ''; end;
if ~exist('init_trait2_from_out_file', 'var'), init_trait2_from_out_file = ''; end;

% full path to bgmg shared library
if ~exist('bgmg_shared_library', 'var'), error('bgmg_shared_library is required'); end;
if ~exist('bgmg_shared_library_header', 'var'), [a,b,c]=fileparts(bgmg_shared_library); bgmg_shared_library_header = [fullfile(a, b), '.h']; clear('a', 'b','c'); end;

BGMG_cpp.unload(); 
BGMG_cpp.load(bgmg_shared_library, bgmg_shared_library_header);
BGMG_cpp.init_log([out_file, '.bgmglib.log']);

% reference file containing mafvec, chrnumvec and posvec for all SNPs to consider in this analysis. 
if ~exist('reference_file', 'var'), error('reference_file is required'); end;
%fprintf('Loading reference from %s... ', reference_file); ref = load(reference_file, 'mafvec', 'chrnumvec', 'posvec'); fprintf('OK.\n');
fprintf('Loading reference from %s... ', reference_file); ref = load(reference_file, 'mafvec', 'posvec', 'chrnumvec'); fprintf('OK.\n');
clear('defvec'); defvec_tmp = true(length(ref.mafvec), 1);
    
% defvec_file determine the set of tag SNPs. It must have a binary variable "defvec" of the same length as defined by reference file ( 1 = tag snp, 0 = exclude from the analysis ). 
% When multiple defvec_files are defined we take an overlap (e.i. consider only tag SNPs defined across all defvec files).
if ~exist('defvec_files', 'var'), defvec_files = {}; end;
if ~exist('MAF_THRESH', 'var'), MAF_THRESH = nan; end;
if ~exist('EXCLUDE_MHC', 'var'), EXCLUDE_MHC = true; end;

for i=1:length(defvec_files)
    fprintf('Loading %s... ', defvec_files{i}); cur_defvec = load(defvec_files{i}); defvec_tmp = defvec_tmp & cur_defvec.defvec; fprintf('OK.\n');
    fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));
end; clear('i');
if isfinite(MAF_THRESH), defvec_tmp = defvec_tmp & (ref.mafvec >= MAF_THRESH); fprintf('Exclude %i variants due to mafvec (%i variants remain)\n', sum(ref.mafvec < MAF_THRESH), sum(defvec_tmp)); end
if EXCLUDE_MHC, defvec_mhc = ~((ref.chrnumvec == 6) & (ref.posvec > 25e6) & (ref.posvec < 35e6)); defvec_tmp = defvec_tmp & defvec_mhc; fprintf('Exclude %i variants in MHC region (%i variants remain)\n', sum(~defvec_mhc), sum(defvec_tmp)); end

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
if ~exist('params_file', 'var'), params_file = ''; end

if ~isempty(params_file),
    params = load(params_file);
    params.sig2_zero = [1;1];
    % Make sure that first component is specific to trait1, second
    % component - specific to trait2, and third is pleiotropic component
    idx = params.pi_vec_trait1 ~= 0 & params.pi_vec_trait2 == 0; assert(sum(idx) <= 1);
    params.pi_vec_tmp(1) = sum(params.pi_vec_trait1(idx));
    params.rho_beta_tmp(1) = sum(params.rho_beta(idx));
    
    idx = params.pi_vec_trait2 ~= 0 & params.pi_vec_trait1 == 0; assert(sum(idx) <= 1);
    params.pi_vec_tmp(2) = sum(params.pi_vec_trait2(idx));
    params.rho_beta_tmp(2) = sum(params.rho_beta(idx));
    
    idx = params.pi_vec_trait2 ~= 0 & params.pi_vec_trait1 ~= 0; assert(sum(idx) <= 1);
    params.pi_vec_tmp(3) = sum(params.pi_vec_trait2(idx));
    params.rho_beta_tmp(3) = sum(params.rho_beta(idx));
    
    params.pi_vec = params.pi_vec_tmp;
    params.rho_beta = params.rho_beta_tmp;
    params = rmfield(params, {'pi_vec_trait2', 'pi_vec_trait1', 'pi_vec_tmp', 'rho_beta_tmp'});
    
    if all(size(params.sig2_beta) == [1 2])
        params.sig2_beta = [params.sig2_beta(1) 0 params.sig2_beta(1); 0 params.sig2_beta(2) params.sig2_beta(2)];
    end

    true_params = [];
    true_params.univariate{1} = struct('sig2_zero', params.sig2_zero(1), 'pi_vec', sum(params.pi_vec([1 3])), 'sig2_beta', params.sig2_beta(1,1));
    true_params.univariate{2} = struct('sig2_zero', params.sig2_zero(2), 'pi_vec', sum(params.pi_vec([2 3])), 'sig2_beta', params.sig2_beta(2,2));
    true_params.bivariate = params;
end

fprintf('Loading %s...', trait1_file); trait1_data = load(trait1_file); fprintf('OK.\n'); num_components = 1;
if (length(trait1_data.zvec) ~= length(ref.chrnumvec)), error('trait1_file is incompatible with the reference'); end;
if isfinite(trait1_nvec), trait1_data.nvec = ones(size(trait1_data.zvec)) * trait1_nvec; end;
if ~isfield(trait1_data, 'nvec'), error('nvec is not available in trait1_file, and trait1_nvec parameter is not set'); end;
cur_defvec.defvec = isfinite(trait1_data.zvec + trait1_data.nvec); defvec_tmp = defvec_tmp & cur_defvec.defvec;
fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));

if ~isempty(trait2_file),
    fprintf('Loading %s...', trait2_file); trait2_data = load(trait2_file); fprintf('OK.\n'); num_components = 3;
    if (length(trait2_data.zvec) ~= length(ref.chrnumvec)), error('trait2_file is incompatible with the reference'); end;
    if isfinite(trait2_nvec), trait2_data.nvec = ones(size(trait2_data.zvec)) * trait2_nvec; end;
    if ~isfield(trait2_data, 'nvec'), error('nvec is not available in trait2_file, and trait2_nvec parameter is not set'); end;
    cur_defvec.defvec = isfinite(trait2_data.zvec + trait2_data.nvec); defvec_tmp = defvec_tmp & cur_defvec.defvec;
    fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec_tmp));
end

if ~exist('randprune_n', 'var'), randprune_n = 16; end;
if ~exist('randprune_r2', 'var'), randprune_r2 = 0.8; end;
if ~exist('randprune_r2_weight_threshold', 'var'), randprune_r2_weight_threshold = nan; end; % SNPs with weight below this threshold are excluded from tag snps
if ~exist('randprune_r2_defvec_threshold', 'var'), randprune_r2_defvec_threshold = nan; end; % 1 iterations at this threshold to reduce set of SNPs (as in defvec)
if ~exist('randprune_r2_plink_ld_mat', 'var'), randprune_r2_plink_ld_mat = ''; end;

if ~exist('kmax', 'var'), kmax = 1000; end;
if ~exist('r2min', 'var'), r2min = 0.01; end;
if ~exist('max_causal_fraction', 'var'), max_causal_fraction = 0.02; end;
if ~exist('cache_tag_r2sum', 'var'), cache_tag_r2sum = 1; end;

if ~exist('DO_FIT', 'var'), DO_FIT = true; end;                % perform fitting
if ~exist('FIT_INIT_TRUE_PARAMS', 'var'), FIT_INIT_TRUE_PARAMS = true; end;    % on synthetic data, init fit with true parameters
if ~exist('FIT_FULL_MODEL', 'var'), FIT_FULL_MODEL = true; end;                % use full model (when false, use gaussian approximation)
if ~exist('FIT_WITH_CONSTRAINS', 'var'), FIT_WITH_CONSTRAINS = true; end;      % fit bivariate model with univariate constrains
if ~exist('QQ_PLOT_TRUE', 'var'), QQ_PLOT_TRUE = false; end;   % make QQ plots with true parameters
if ~exist('QQ_PLOT_FIT', 'var'), QQ_PLOT_FIT = false; end;     % make QQ plots with fitted parameters
if ~exist('QQ_PLOT_DOWNSCALE', 'var'), QQ_PLOT_DOWNSCALE = 1; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('STRATIFIED_QQ_PLOT_FIT', 'var'), STRATIFIED_QQ_PLOT_FIT = false; end;
if ~exist('STRATIFIED_QQ_PLOT_TRUE', 'var'), STRATIFIED_QQ_PLOT_TRUE = false; end;
if ~exist('STRATIFIED_QQ_PLOT_DOWNSCALE', 'var'), STRATIFIED_QQ_PLOT_DOWNSCALE = 10; end;     % downscale #snps in stratified QQ plots (model prediction only)
if ~exist('LOGLIKE_PLOT_FIT', 'var'), LOGLIKE_PLOT_FIT = false; end;
if ~exist('LOGLIKE_PLOT_TRUE', 'var'), LOGLIKE_PLOT_TRUE = false; end;
if ~exist('POWER_PLOT_FIT', 'var'), POWER_PLOT_FIT = false; end;  % make power plots with fitted parameters
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;
if ~exist('THREADS', 'var'), THREADS = nan; end;

if POWER_PLOT_FIT, error('not yet implemented in c++ version'); end;

result = [];
if ~isempty(init_trait1_from_out_file) || ~isempty(init_trait2_from_out_file)
    if isempty(trait2_file), error('init_traitN_from_out_file is specific for BGMG fit'); end;
    if isempty(init_trait1_from_out_file) || isempty(init_trait2_from_out_file), error('init_traitN_from_out_file should be specified for both traits (N=1 and N=2)'); end;
    if ~isempty(init_result_from_out_file), error('init_result_from_out_file is incompatible with init_traitN_from_out_file'); end;

    fprintf('Loading existing params for trait1 from %s...\n', init_trait1_from_out_file);
    tmp = load(init_trait1_from_out_file); result.univariate{1} = tmp.univariate{1};
    if length(tmp.univariate) ~= 1 || isfield(tmp, 'bivariate'), error('init_trait1_from_out_file does not appear to be a result of univariate analysis'); end;
    fprintf('Loading existing params for trait2 from %s...\n', init_trait2_from_out_file);
    tmp = load(init_trait2_from_out_file); result.univariate{2} = tmp.univariate{1};
    if length(tmp.univariate) ~= 1 || isfield(tmp, 'bivariate'), error('init_trait2_from_out_file does not appear to be a result of univariate analysis'); end;
    clear('tmp');
end

if ~isempty(init_result_from_out_file)
    if DO_FIT, error('init_result_from_out_file is incompatible with DO_FIT'); end;
    fprintf('Loading init file from %s...\n', init_result_from_out_file);
    tmp = load(init_result_from_out_file); result = tmp.result;
    clear('tmp');
end


if (length(chr_labels) == 1) && (chr_labels(1) == 1) && ( all(ref.chrnumvec == 1) || (find(ref.chrnumvec == 1, 1, 'last' ) < find(ref.chrnumvec ~= 1, 1 )))
    chrlabel = chr_labels(1);
    fprintf('Reduce reference to a signle chromosome (chr%i).\n', chrlabel);
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
  
if isfinite(randprune_r2_defvec_threshold)
    % Use hard threshold to exlude sinonimous SNPs from fit. Just one
    % iteration of random pruning with very high r2 threshold. Non-selected
    % SNPs are excluded.
    if ~exist('randprune_r2_plink_ld_mat', 'var') error('randprune_r2_plink_ld_mat is required'); end;
    fprintf('Excluding variants based on random pruning at %.3f threshold...\n', randprune_r2_defvec_threshold);
    tag_indices_tmp = find(defvec_tmp);
    bgmglib=BGMG_cpp(1);
    bgmglib.dispose();
    bgmglib.defvec = defvec_tmp;
    bgmglib.hvec = ref.mafvec .* (1-ref.mafvec) * 2;

    for chr_index=1:length(chr_labels)
        plink_ld_mat_chr = strrep(randprune_r2_plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
        tmp = load(plink_ld_mat_chr); tmp.index_A = tmp.index_A + 1; tmp.index_B = tmp.index_B + 1; % 0-based, comming from python
        bgmglib.set_ld_r2_coo_from_matlab_indices(tmp.index_A, tmp.index_B, tmp.r2);
        printf('OK.\n');
    end
    bgmglib.set_ld_r2_csr();
    bgmglib.set_weights_randprune(1, randprune_r2_defvec_threshold);
    weights_bgmg = bgmglib.weights;
    defvec_tmp(tag_indices_tmp(weights_bgmg==0)) = false;
    fprintf('Exclude %i variants after random pruning at %.3f threshold (%i variants rmain)\n', sum(weights_bgmg == 0), randprune_r2_defvec_threshold, sum(defvec_tmp));
end

if isfinite(randprune_r2_weight_threshold)
    % Use soft threshold on random pruning weight to exlude low-weight
    % SNPs. Perform usual pruning. Then sort SNPs according to their weight, and
    % exclude those that contribute to the bottom 10% (randprune_r2_weight_threshold)
    if ~exist('randprune_r2_plink_ld_mat', 'var') error('randprune_r2_plink_ld_mat is required'); end;
    fprintf('Excluding variants with too low weight, based on random pruning at %.3f threshold...\n', randprune_r2);
    tag_indices_tmp = find(defvec_tmp);
    bgmglib=BGMG_cpp(1);
    bgmglib.dispose();
    bgmglib.defvec = defvec_tmp;
    bgmglib.hvec = ref.mafvec .* (1-ref.mafvec) * 2;

    for chr_index=1:length(chr_labels)
        plink_ld_mat_chr = strrep(randprune_r2_plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
        tmp = load(plink_ld_mat_chr); tmp.index_A = tmp.index_A + 1; tmp.index_B = tmp.index_B + 1; % 0-based, comming from python
        bgmglib.set_ld_r2_coo_from_matlab_indices(tmp.index_A, tmp.index_B, tmp.r2);
        printf('OK.\n');
    end
    bgmglib.set_ld_r2_csr();
    bgmglib.set_weights_randprune(randprune_n, randprune_r2);
    weights_bgmg = bgmglib.weights;
    bgmglib.dispose();
    
    weights_bgmg_sorted = sort(weights_bgmg);
    weights_bgmg_cumsum = cumsum(weights_bgmg_sorted/sum(weights_bgmg_sorted));
    weights_bgmg_thresh = weights_bgmg_sorted(find(weights_bgmg_cumsum > randprune_r2_weight_threshold, 1, 'first'));
    
    defvec_tmp(tag_indices_tmp(weights_bgmg < weights_bgmg_thresh)) = false;
    fprintf('Exclude %i variants (%i variants remain)\n', sum(weights_bgmg < weights_bgmg_thresh), sum(defvec_tmp));
end

% finalize defvec, from here it must not change.
defvec = defvec_tmp; clear('defvec_tmp', 'cur_defvec');
tag_indices = find(defvec);

fprintf('%i tag SNPs will go into fit and/or qq plots, etc\n', length(tag_indices));

bgmglib = BGMG_cpp();
bgmglib.dispose()
bgmglib.defvec = defvec;
bgmglib.zvec1 = trait1_data.zvec(defvec);
bgmglib.zvec2 = trait2_data.zvec(defvec);

bgmglib.set_option('r2min', r2min);
bgmglib.set_option('kmax', kmax);
bgmglib.set_option('max_causals', floor(max_causal_fraction * length(defvec)));
bgmglib.set_option('num_components', num_components);
bgmglib.set_option('cache_tag_r2sum', cache_tag_r2sum);
bgmglib.set_option('threads', THREADS);

bgmglib.hvec = ref.mafvec .* (1-ref.mafvec) * 2;

for chr_index=1:length(chr_labels)
    plink_ld_bin_chr = strrep(plink_ld_bin,'@', sprintf('%i',chr_labels(chr_index)));
    fprintf('Loading %s...', plink_ld_bin_chr); bgmglib.set_ld_r2_coo_from_file(plink_ld_bin_chr); fprintf('OK.\n'); 
    clear('plink_ld_bin_chr');
end; clear('chr_index');
bgmglib.set_ld_r2_csr();

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
options.fit_with_constrains = FIT_WITH_CONSTRAINS;

if exist('true_params', 'var'), 
    % Do a quick fit to initialize sig2_zero and rho_zero
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);
     for trait_index = 1:2
         fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.UGMG_fminsearch_cost(mapparams(x), trait_index), mapparams(x0), struct('Display', 'on')));
         fit_sig2_zero  = fit(struct('sig2_zero', 1), @(x)BGMG_util.UGMG_mapparams1(x, struct('pi_vec', true_params.univariate{trait_index}.pi_vec, 'sig2_beta', true_params.univariate{trait_index}.sig2_beta)));
         true_params.univariate{trait_index}.sig2_zero = fit_sig2_zero.sig2_zero;
     end
     
     fit = @(x0, mapparams)mapparams(fminsearch(@(x)BGMG_util.BGMG_fminsearch_cost(mapparams(x)), mapparams(x0), struct('Display', 'on')));
     fit_rho_zero = fit(...
         struct('rho_zero', 0.1), ...
         @(x)BGMG_util.BGMG_mapparams3(x, struct(...
            'sig2_zero', [true_params.univariate{1}.sig2_zero, true_params.univariate{2}.sig2_zero], ... 
            'sig2_beta', true_params.bivariate.sig2_beta, ...
            'rho_beta', true_params.bivariate.rho_beta, ...
            'pi_vec', true_params.bivariate.pi_vec)));
    true_params.bivariate.rho_zero = fit_rho_zero.rho_zero;

    options.true_params = true_params;
end;

disp(options)

% Fit bivariate or univariate model to the data
if DO_FIT
    if isempty(trait2_file),
        result.univariate{1} = BGMG_cpp_fit_univariate(1, options);
    else
        if ~isfield(result, 'univariate')
            result.univariate{1} = BGMG_cpp_fit_univariate(1, options);
            result.univariate{2} = BGMG_cpp_fit_univariate(2, options);
        end
        result.bivariate = BGMG_cpp_fit_bivariate(...
            result.univariate{1}.params, ...
            result.univariate{2}.params, ...
            options);
    end

    result.trait1_file = trait1_file;
    result.trait2_file = trait2_file;
    result.reference_file = reference_file;
    result.options = options;
    
    % Save the result in .mat file
    % (this overrides previously saved file)
    save([out_file '.preliminary.mat'], 'result');
    fprintf('Results saved to %s.mat\n', out_file);
end

bgmglib.set_option('diag', 0);

% Produce QQ plots with true params (only works for synthetic data, of course)
if QQ_PLOT_TRUE
    for trait_index = 1:(1 + ~isempty(trait2_file))
        options.downscale = QQ_PLOT_DOWNSCALE;
        [figures, plot_data] = BGMG_cpp_qq_plot(true_params.univariate{trait_index}, trait_index, options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_true_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.true.pdf', out_file, trait_index), '-dpdf')
    end
end

% Produce QQ plots with fitted params
if QQ_PLOT_FIT
    for trait_index = 1:length(result.univariate)
        options.downscale = QQ_PLOT_DOWNSCALE;
        [figures, plot_data] = BGMG_cpp_qq_plot(result.univariate{trait_index}.params, trait_index, options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_fit_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.fit.pdf', out_file, trait_index), '-dpdf')
    end
end

if STRATIFIED_QQ_PLOT_TRUE && ~isempty(params_file) && ~isempty(trait2_file)
    options.downscale = STRATIFIED_QQ_PLOT_DOWNSCALE;
    [figures, plot_data] = BGMG_cpp_stratified_qq_plot(true_params.bivariate, options);
    result.bivariate.stratified_qq_plot_fit_data = plot_data;
    print(figures.tot, sprintf('%s.stratqq.fit.pdf', out_file), '-dpdf')
end

if STRATIFIED_QQ_PLOT_FIT && (DO_FIT || ~isempty(init_result_from_out_file)) && ~isempty(trait2_file)
    options.downscale = STRATIFIED_QQ_PLOT_DOWNSCALE;
    [figures, plot_data] = BGMG_cpp_stratified_qq_plot(result.bivariate.params, options);
    result.bivariate.stratified_qq_plot_fit_data = plot_data;
    print(figures.tot, sprintf('%s.stratqq.fit.pdf', out_file), '-dpdf')
end

if LOGLIKE_PLOT_FIT && (DO_FIT || ~isempty(init_result_from_out_file)) && ~isempty(trait2_file)
    f = figure; hold on;
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);
    bgmglib.clear_loglike_cache();
    [figures, plots_data] = BGMG_cpp_loglike_plot(result.bivariate.params);
    result.bivariate.loglike_plot_data_fit = plots_data;
    result.bivariate.loglike_plot_trajectory_fit = bgmglib.extract_bivariate_loglike_trajectory();
    print(figures.tot, sprintf('%s.loglike.fit.pdf', out_file), '-dpdf')
end

if LOGLIKE_PLOT_TRUE && ~isempty(trait2_file) && ~isempty(params_file)
    bgmglib.set_option('fast_cost', ~FIT_FULL_MODEL);
    bgmglib.clear_loglike_cache();
    [figures, plots_data] = BGMG_cpp_loglike_plot(true_params.bivariate);
    result.bivariate.loglike_plot_data_true = plots_data;
    result.bivariate.loglike_plot_trajectory_true = bgmglib.extract_bivariate_loglike_trajectory();
    print(figures.tot, sprintf('%s.loglike.true.pdf', out_file), '-dpdf')
end

bgmglib.set_option('diag', 0);

if ~exist('result', 'var')
    error('No options selected; please enable DO_FIT or QQ_PLOT_TRUE');
end

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.mat'], 'result');
fprintf('Results saved to %s.mat\n', out_file);

if 0
    % Helper code to save all results to a text file
    % SIMU_BGMG2_2018_06_17
    dirs=dir('H:\work\SIMU_BGMG_9pifrac_2018_06_21\*.bgmg.mat');
    fileID = fopen(['H:\work\SIMU_BGMG_9pifrac_2018_06_21.csv'], 'w');
    has_header = false;
    for i=1:length(dirs)
        try
        x = load(fullfile('H:\work\SIMU_BGMG_9pifrac_2018_06_21', dirs(i).name));
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