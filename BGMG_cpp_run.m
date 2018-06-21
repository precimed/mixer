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
if 0
bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll';
bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';
plink_ld_mat = 'H:\work\hapgen_ldmat2_plink\bfile_merged_10K_ldmat_p01_SNPwind50k_chr@.ld.mat'; chr_labels = 1;
randprune_r2_plink_ld_mat = ''; randprune_r2_defvec_threshold = nan;
randprune_r2_plink_ld_mat = 'H:\work\hapgen_ldmat2_plink\bfile_merged_10K_ldmat_p10_SNPwind50k_chr@.ld.mat';
%defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'};
defvec_files = {};
defvec_files = {'H:\Dropbox\shared\BGMG\defvecs\defvec_11m_infinium-omniexpress-24-v1-3-a1.mat', ...
                'H:\Dropbox\shared\BGMG\defvecs\defvec_highld_regions.mat', ... 
                'H:\Dropbox\shared\BGMG\defvecs\centromere_plus_3cm_locations.mat' };
trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait2_nvec=100000;
trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=3e-04_rep=10_tag1=completePolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=3e-04_rep=10_tag1=completePolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait2_nvec=100000;


%trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
%trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait2_nvec=100000;
%trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=3e-03_rep=10_tag1=completePolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
%trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=3e-03_rep=10_tag1=completePolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait1_nvec=100000;
%trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=8e-04_rep=10_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
%trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=8e-04_rep=10_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait1_nvec=100000;

reference_file = 'H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins_wld.mat';
DO_FIT=false;FIT_FULL_MODEL=true;STRATIFIED_QQ_PLOT_FIT=true;QQ_PLOT_TRUE=true;LOGLIKE_PLOT_TRUE=true;QQ_PLOT_FIT=true;cache_tag_r2sum=true;
MAF_THRESH=0.01;
out_file = 'tmptesting5_ext';
init_result_from_out_file = 'tmptesting5';
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

% full path to bgmg shared library
if ~exist('bgmg_shared_library', 'var'), error('bgmg_shared_library is required'); end;
if ~exist('bgmg_shared_library_header', 'var'), [a,b,c]=fileparts(bgmg_shared_library); bgmg_shared_library_header = [fullfile(a, b), '.h']; clear('a', 'b','c'); end;
if libisloaded('bgmg'), unloadlibrary('bgmg'); end;
if ~libisloaded('bgmg'), fprintf('Loading bgmg library: %s, %s... ', bgmg_shared_library, bgmg_shared_library_header); loadlibrary(bgmg_shared_library, bgmg_shared_library_header, 'alias', 'bgmg');  fprintf('OK.\n'); end;
calllib('bgmg', 'bgmg_init_log', [out_file, '.bgmglib.log']);

% reference file containing mafvec, chrnumvec and posvec for all SNPs to consider in this analysis. 
if ~exist('reference_file', 'var'), error('reference_file is required'); end;
%fprintf('Loading reference from %s... ', reference_file); ref = load(reference_file, 'mafvec', 'chrnumvec', 'posvec'); fprintf('OK.\n');
fprintf('Loading reference from %s... ', reference_file); ref = load(reference_file); fprintf('OK.\n');
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
if ~exist('plink_ld_mat', 'var'), error('plink_ld_mat is required'); end;
if ~exist('chr_labels', 'var'), chr_labels = 1:22; end;

if ~exist('trait1_file', 'var'), error('trait1_file is required'); end;
if ~exist('trait2_file', 'var'), trait2_file = ''; end;
if ~exist('trait1_nvec', 'var'), trait1_nvec = nan; end;
if ~exist('trait2_nvec', 'var'), trait2_nvec = nan; end;

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

if ~exist('randprune_n', 'var'), randprune_n = 100; end;
if ~exist('randprune_r2', 'var'), randprune_r2 = 0.8; end;
if ~exist('randprune_r2_weight_threshold', 'var'), randprune_r2_weight_threshold = nan; end; % SNPs with weight below this threshold are excluded from tag snps
if ~exist('randprune_r2_defvec_threshold', 'var'), randprune_r2_defvec_threshold = nan; end; % 1 iterations at this threshold to reduce set of SNPs (as in defvec)
if ~exist('randprune_r2_plink_ld_mat', 'var'), randprune_r2_plink_ld_mat = ''; end;

if ~exist('kmax', 'var'), kmax = 1000; end;
if ~exist('r2min', 'var'), r2min = 0.01; end;
if ~exist('max_causal_fraction', 'var'), max_causal_fraction = 0.02; end;
if ~exist('cache_tag_r2sum', 'var'), cache_tag_r2sum = 1; end;

if ~exist('DO_FIT', 'var'), DO_FIT = true; end;                % perform fitting
if ~exist('FIT_FULL_MODEL', 'var'), FIT_FULL_MODEL = true; end;                % use full model (when false, use gaussian approximation)
if ~exist('FIT_WITH_CONSTRAINS', 'var'), FIT_WITH_CONSTRAINS = true; end;      % fit bivariate model with univariate constrains
if ~exist('QQ_PLOT_TRUE', 'var'), QQ_PLOT_TRUE = false; end;   % make QQ plots with true parameters
if ~exist('QQ_PLOT_FIT', 'var'), QQ_PLOT_FIT = false; end;     % make QQ plots with fitted parameters
if ~exist('QQ_PLOT_DOWNSCALE', 'var'), QQ_PLOT_DOWNSCALE = 1; end;     % downscale #snps in QQ plots (model prediction only)
if ~exist('STRATIFIED_QQ_PLOT_FIT', 'var'), STRATIFIED_QQ_PLOT_FIT = false; end;
if ~exist('STRATIFIED_QQ_PLOT_DOWNSCALE', 'var'), STRATIFIED_QQ_PLOT_DOWNSCALE = 10; end;     % downscale #snps in stratified QQ plots (model prediction only)
if ~exist('LOGLIKE_PLOT_FIT', 'var'), LOGLIKE_PLOT_FIT = false; end;
if ~exist('LOGLIKE_PLOT_TRUE', 'var'), LOGLIKE_PLOT_TRUE = false; end;
if ~exist('POWER_PLOT_FIT', 'var'), POWER_PLOT_FIT = false; end;  % make power plots with fitted parameters
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;

if POWER_PLOT_FIT, error('not yet implemented in c++ version'); end;

result = [];
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
  
m2c = @(x)(x-1); % convert matlab to c indices
check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
check_for_context = @(context)fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', context));

if isfinite(randprune_r2_defvec_threshold)
    % Use hard threshold to exlude sinonimous SNPs from fit. Just one
    % iteration of random pruning with very high r2 threshold. Non-selected
    % SNPs are excluded.
    if ~exist('randprune_r2_plink_ld_mat', 'var') error('randprune_r2_plink_ld_mat is required'); end;
    fprintf('Excluding variants based on random pruning at %.3f threshold...\n', randprune_r2_defvec_threshold);
    context = 1; tag_indices_tmp = find(defvec_tmp);
    calllib('bgmg', 'bgmg_set_tag_indices', context, length(defvec_tmp), length(tag_indices_tmp), m2c(tag_indices_tmp));  check_for_context(context);

    for chr_index=1:length(chr_labels)
        plink_ld_mat_chr = strrep(randprune_r2_plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
        tmp = load(plink_ld_mat_chr); tmp.index_A = tmp.index_A + 1; tmp.index_B = tmp.index_B + 1; % 0-based, comming from python
        calllib('bgmg', 'bgmg_set_ld_r2_coo', context, length(tmp.r2), m2c(tmp.index_A), m2c(tmp.index_B), tmp.r2); fprintf('OK.\n'); check(); 
    end
    calllib('bgmg', 'bgmg_set_ld_r2_csr', context);  check_for_context(context);
    calllib('bgmg', 'bgmg_set_weights_randprune', context, 1, randprune_r2_defvec_threshold);  check_for_context(context);
    
    pBuffer = libpointer('singlePtr', zeros(sum(defvec_tmp), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', context, sum(defvec_tmp), pBuffer);  check_for_context(context); weights_bgmg = pBuffer.Value;
    clear pBuffer
    calllib('bgmg', 'bgmg_dispose', context);  check_for_context(context);

    defvec_tmp(tag_indices_tmp(weights_bgmg==0)) = false;
    fprintf('Exclude %i variants after random pruning at %.3f threshold (%i variants remain)\n', sum(weights_bgmg == 0), randprune_r2_defvec_threshold, sum(defvec_tmp));
end

if isfinite(randprune_r2_weight_threshold)
    % Use soft threshold on random pruning weight to exlude low-weight
    % SNPs. Perform usual pruning. Then sort SNPs according to their weight, and
    % exclude those that contribute to the bottom 10% (randprune_r2_weight_threshold)
    if ~exist('randprune_r2_plink_ld_mat', 'var') error('randprune_r2_plink_ld_mat is required'); end;
    fprintf('Excluding variants with too low weight, based on random pruning at %.3f threshold...\n', randprune_r2);
    context = 1; tag_indices_tmp = find(defvec_tmp);
    calllib('bgmg', 'bgmg_set_tag_indices', context, length(defvec_tmp), length(tag_indices_tmp), m2c(tag_indices_tmp));  check_for_context(context);

    for chr_index=1:length(chr_labels)
        plink_ld_mat_chr = strrep(randprune_r2_plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
        tmp = load(plink_ld_mat_chr); tmp.index_A = tmp.index_A + 1; tmp.index_B = tmp.index_B + 1; % 0-based, comming from python
        calllib('bgmg', 'bgmg_set_ld_r2_coo', context, length(tmp.r2), m2c(tmp.index_A), m2c(tmp.index_B), tmp.r2); fprintf('OK.\n'); check(); 
    end
    calllib('bgmg', 'bgmg_set_ld_r2_csr', context);  check_for_context(context);
    calllib('bgmg', 'bgmg_set_weights_randprune', context, randprune_n, randprune_r2);  check_for_context(context);
    
    pBuffer = libpointer('singlePtr', zeros(sum(defvec_tmp), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', context, sum(defvec_tmp), pBuffer);  check_for_context(context); weights_bgmg = pBuffer.Value;
    clear pBuffer
    calllib('bgmg', 'bgmg_dispose', context);  check_for_context(context);
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

calllib('bgmg', 'bgmg_set_tag_indices', 0, length(defvec), length(tag_indices), m2c(tag_indices));  check();
calllib('bgmg', 'bgmg_set_option', 0,  'r2min', r2min); check();
calllib('bgmg', 'bgmg_set_option', 0,  'kmax', kmax); check();
calllib('bgmg', 'bgmg_set_option', 0,  'max_causals', floor(max_causal_fraction * length(defvec))); check();  
calllib('bgmg', 'bgmg_set_option', 0,  'num_components', num_components); check();
calllib('bgmg', 'bgmg_set_option', 0,  'cache_tag_r2sum', cache_tag_r2sum); check();

hvec = ref.mafvec .* (1-ref.mafvec) * 2; 
calllib('bgmg', 'bgmg_set_hvec', 0, length(hvec), hvec);  check(); clear('hvec');

for chr_index=1:length(chr_labels)
    plink_ld_mat_chr = strrep(plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
    fprintf('Loading %s...', plink_ld_mat_chr); tmp = load(plink_ld_mat_chr); fprintf('OK.\n');
    tmp.index_A = tmp.index_A + 1; 
    tmp.index_B = tmp.index_B + 1; % 0-based, comming from python
    fprintf('Importing %s to bgmg library...', plink_ld_mat_chr); 
    
    % chunk all LD r2 elements into pieces, each no more than 1e6.
    chunks_up = length(tmp.r2); chunks_step = 10000000;
    chunks = 1:chunks_step:chunks_up; 
    if (length(chunks) >= 3), chunks_start = chunks(1:end-1); chunks_end = [chunks(2:(end-1))-1 chunks_up];
    else chunks_start = 1; chunks_end=chunks_up; end;
    if (sum(chunks_end - chunks_start + 1) ~= chunks_up), error('error in chunking logig'); end;
    
    for chunk_index = 1:length(chunks_start)
        a = chunks_start(chunk_index); b = chunks_end(chunk_index);
        calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, b-a+1, m2c(tmp.index_A(a:b)), m2c(tmp.index_B(a:b)), tmp.r2(a:b)); fprintf('OK.\n'); check(); 
    end
    
    clear('tmp', 'plink_ld_mat_chr', 'a', 'b', 'chunk_index', 'chunk_start', 'chunk_end', 'chunks_step', 'chunk_up');
end; clear('chr_index');

calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();
calllib('bgmg', 'bgmg_set_weights_randprune', 0, randprune_n, randprune_r2);  check();

calllib('bgmg', 'bgmg_set_zvec', 0, 1, sum(defvec), trait1_data.zvec(defvec));
calllib('bgmg', 'bgmg_set_nvec', 0, 1, sum(defvec), trait1_data.nvec(defvec));

if ~isempty(trait2_file),
    calllib('bgmg', 'bgmg_set_zvec', 0, 2, sum(defvec), trait2_data.zvec(defvec));
    calllib('bgmg', 'bgmg_set_nvec', 0, 2, sum(defvec), trait2_data.nvec(defvec));
end

calllib('bgmg', 'bgmg_set_option', 0,  'diag', 0); check();
% Preparation is done, BGMG library is fully setup. Now we can use it to
% calculate model QQ plots and univariate or bivariate cost function.

options = [];
options.total_het = sum(2*ref.mafvec.*(1-ref.mafvec));
options.verbose = true;
options.ci_alpha = CI_ALPHA;
options.title = TITLE;
options.fit_full_model = FIT_FULL_MODEL;
options.fit_with_constrains = FIT_WITH_CONSTRAINS;

% Save true parameters (if available)
if isfield(trait1_data, 'causal_pi'), options.causal_pi_T1 = trait1_data.causal_pi; end;
if isfield(trait1_data, 'sigsq'), options.sigsq_T1 = trait1_data.sigsq; end;
if ~isempty(trait2_file),
    if isfield(trait2_data, 'causal_pi'), options.causal_pi_T2 = trait2_data.causal_pi; end;
    if isfield(trait2_data, 'sigsq'), options.sigsq_T2 = trait2_data.sigsq; end;
end

disp(options)

% Fit bivariate or univariate model to the data
if DO_FIT
    if ~isempty(trait2_file),
        zmat = [trait1_data.zvec, trait2_data.zvec];
        nmat = [trait1_data.nvec, trait2_data.nvec];
        result2 = BGMG_cpp_fit(zmat, nmat, options);
    else
        zvec = [trait1_data.zvec];
        nvec = [trait1_data.nvec];
        result2 = BGMG_cpp_fit(zvec, nvec, options);
    end
    
    fnames = fieldnames(result2);
    for findex = 1:length(fnames)
        result.(fnames{findex}) = result2.(fnames{findex});
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

calllib('bgmg', 'bgmg_set_option', 0,  'diag', 0); check();

% Produce QQ plots with true params (only works for synthetic data, of course)
if QQ_PLOT_TRUE
    for trait_index = 1:(1 + ~isempty(trait2_file))
        if trait_index==1
            qq_params = struct('sig2_zero', 1, 'pi_vec', sum(trait1_data.causal_pi), 'sig2_beta', trait1_data.sigsq);
            qq_data = trait1_data;
        else
            qq_params = struct('sig2_zero', 1, 'pi_vec', sum(trait2_data.causal_pi), 'sig2_beta', trait2_data.sigsq);
            qq_data = trait2_data;
        end
        options.downscale = QQ_PLOT_DOWNSCALE;
        [figures, plot_data] = BGMG_cpp_qq_plot(qq_params, qq_data.zvec(defvec), options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_true_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.true.pdf', out_file, trait_index), '-dpdf')
    end
end

% Produce QQ plots with fitted params
if QQ_PLOT_FIT
    for trait_index = 1:length(result.univariate)
        if trait_index==1
            qq_params = result.univariate{trait_index}.params;
            qq_data = trait1_data;
        else
            qq_params = result.univariate{trait_index}.params;
            qq_data = trait2_data;
        end
        options.downscale = QQ_PLOT_DOWNSCALE;
        [figures, plot_data] = BGMG_cpp_qq_plot(qq_params, qq_data.zvec(defvec), options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_fit_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.fit.pdf', out_file, trait_index), '-dpdf')
    end
end

if STRATIFIED_QQ_PLOT_FIT && (DO_FIT || ~isempty(init_result_from_out_file)) && ~isempty(trait2_file)
    %ov=[]; ov.pi_vec = [0, 0, trait1_data.causal_pi];     ov.sig2_beta = [trait1_data.sigsq, trait2_data.sigsq];     ov.sig2_zero = [1, 1];    ov.rho_zero = 0; ov.rho_beta = 0;
    zmat = [trait1_data.zvec, trait2_data.zvec]; 
    options.downscale = STRATIFIED_QQ_PLOT_DOWNSCALE;
    [figures, plot_data] = BGMG_cpp_stratified_qq_plot(result.bivariate.params, zmat(defvec, :), options);
    result.bivariate.stratified_qq_plot_fit_data = plot_data;
    print(figures.tot, sprintf('%s.stratqq.fit.pdf', out_file), '-dpdf')
end

if LOGLIKE_PLOT_FIT && (DO_FIT || ~isempty(init_result_from_out_file)) && ~isempty(trait2_file)
    f = figure; hold on;
    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', ~FIT_FULL_MODEL); check();
    [figures, plots_data] = BGMG_cpp_loglike_plot(result.bivariate.params);
    result.bivariate.loglike_plot_fit_data = plots_data;
    print(figures.tot, sprintf('%s.loglike.fit.pdf', out_file), '-dpdf')
    %subplot(3,3,7); legend('random overlap', 'full overlap');
end

if LOGLIKE_PLOT_TRUE && ~isempty(trait2_file)
    if ~DO_FIT && isempty(init_result_from_out_file)
        fprintf('Fit model with fit_full_model=false to find sig2zero and rho_zero params\n');
        options.fit_full_model = false;
        result_sig0_rho0 = BGMG_cpp_fit([trait1_data.zvec, trait2_data.zvec], [trait1_data.nvec, trait2_data.nvec], options);
        options.fit_full_model = FIT_FULL_MODEL;
    else
        result_sig0_rho0 = result;
    end

    true_params = [];
    pi1u = trait1_data.causal_pi;
    pi2u = trait2_data.causal_pi;
    if ~isempty(strfind(trait1_file, 'random')) % <-------- TBD; this needs to be fixed, we need to save true pi_vec somewhere.
        pi12 = sum(pi1u)*sum(pi2u);
    else
        pi12 = sum(pi1u(pi1u>0 & pi2u>0));
    end
    pi1u=sum(pi1u); pi2u=sum(pi2u);
    true_params.pi_vec = [pi1u - pi12, pi2u - pi12, pi12];
    true_params.sig2_beta = [trait1_data.sigsq,0,trait1_data.sigsq; 0, trait1_data.sigsq, trait2_data.sigsq];
    true_params.rho_beta = [0.0 0.0 0.0];  % <-------- TBD; this needs to be fixed, we need to save true rho_beta somewhere.
    true_params.sig2_zero = result_sig0_rho0.bivariate.params.sig2_zero;
    true_params.rho_zero = result_sig0_rho0.bivariate.params.rho_zero;
    
    calllib('bgmg', 'bgmg_set_option', 0, 'fast_cost', ~FIT_FULL_MODEL); check();
    [figures, plots_data] = BGMG_cpp_loglike_plot(true_params);
    result.bivariate.loglike_plot_fit_data = plots_data;
    print(figures.tot, sprintf('%s.loglike.true.pdf', out_file), '-dpdf')
    %subplot(3,3,7); legend('random overlap', 'full overlap');
end

calllib('bgmg', 'bgmg_set_option', 0,  'diag', 0); check();

if ~exist('result', 'var')
    error('No options selected; please enable DO_FIT or QQ_PLOT_TRUE');
end

% Save the result in .mat file
% (this overrides previously saved file)
save([out_file '.mat'], 'result');
fprintf('Results saved to %s.mat\n', out_file);

if 0
    % Helper code to save all results to a text file
    dirs=dir('H:\work\SIMU_BGMG2_2018_06_17\*.bgmg.mat');
    fileID = fopen(['H:\work\SIMU_BGMG2_2018_06_18.csv'], 'w');
    has_header = false;
    for i=1:length(dirs)
        try
        x = load(fullfile('H:\work\SIMU_BGMG2_2018_06_17', dirs(i).name));
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