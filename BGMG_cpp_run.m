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
plink_ld_mat = 'H:\work\hapgen_ldmat2_plink\bfile_merged_10K_ldmat_p01_SNPwind50k_chr@.ld.mat'; chr_labels = 1:22;
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'};
trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; trait1_nvec=100000;
trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; trait2_nvec=100000;
reference_file = 'H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins_wld.mat';
DO_FIT=true;QQ_PLOT_TRUE=true;QQ_PLOT_FIT=true;

bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\Debug\bgmg.dll';
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'}
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_1kG_phase3_EUR.mat' };
defvec_files = {'H:\Dropbox\shared\BGMG\defvec_1kG_phase3_EUR.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'}
trait1_file = 'H:\Dropbox\shared\BGMG\p_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc\simu_h2=0.40_pi1u=1e-3_rep=1.ugmg.mat'
trait1_file = 'H:\GitHub\BGMG\GIANT_HEIGHT_2014_lift.mat';
trait1_file = 'H:\GitHub\BGMG\PGC_SCZ_2014.mat';
reference_file = 'H:\Dropbox\shared\BGMG\1kG_phase3_EUR_11015883_reference_p01_SNPwind50k_per_allele_4bins.mat'
DO_FIT = false; QQ_PLOT_TRUE = true;
init_file = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.bgmg.mat';
trait1_file  = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat';
trait2_file  = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat';
end

addpath('DERIVESTsuite');

if ~exist('out_file', 'var'), out_file = 'BGMG_result'; end;

% full path to bgmg shared library
if ~exist('bgmg_shared_library', 'var'), error('bgmg_shared_library is required'); end;
if ~exist('bgmg_shared_library_header', 'var'), [a,b,c]=fileparts(bgmg_shared_library); bgmg_shared_library_header = [fullfile(a, b), '.h']; clear('a', 'b','c'); end;
if libisloaded('bgmg'), unloadlibrary('bgmg'); end;
if ~libisloaded('bgmg'), fprintf('Loading bgmg library: %s, %s... ', bgmg_shared_library, bgmg_shared_library_header); loadlibrary(bgmg_shared_library, bgmg_shared_library_header, 'alias', 'bgmg');  fprintf('OK.\n'); end;
calllib('bgmg', 'bgmg_init_log', [out_file, '.bgmglib.log']);

% reference file containing mafvec, chrnumvec and posvec for all SNPs to consider in this analysis. 
if ~exist('reference_file', 'var'), error('reference_file is required'); end;
fprintf('Loading reference from %s... ', reference_file); ref = load(reference_file, 'mafvec', 'chrnumvec', 'posvec'); fprintf('OK.\n');
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
if isfinite(MAF_THRESH), defvec_tmp = defvec_tmp & (mafvec >= MAF_THRESH); fprintf('Exclude %i variants due to mafvec (%i variants remain)\n', sum(mafvec < MAF_THRESH), sum(defvec_tmp)); end
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
if ~exist('kmax', 'var'), kmax = 1000; end;
if ~exist('r2min', 'var'), r2min = 0.01; end;
if ~exist('max_causal_fraction', 'var'), max_causal_fraction = 0.02; end;

if ~exist('DO_FIT', 'var'), DO_FIT = true; end;                % perform fitting
if ~exist('QQ_PLOT_TRUE', 'var'), QQ_PLOT_TRUE = false; end;   % make QQ plots with true parameters
if ~exist('QQ_PLOT_FIT', 'var'), QQ_PLOT_FIT = false; end;     % make QQ plots with fitted parameters
if ~exist('POWER_PLOT_FIT', 'var'), POWER_PLOT_FIT = false; end;  % make power plots with fitted parameters
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;

if POWER_PLOT_FIT, error('not yet implemented in c++ version'); end;

if (length(chr_labels) == 1) && (chr_labels(1) == 1) && (find(ref.chrnumvec == 1, 1, 'last' ) < find(ref.chrnumvec ~= 1, 1 ))
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
    
% finalize defvec, from here it must not change.
defvec = defvec_tmp; clear('defvec_tmp', 'cur_defvec');
tag_indices = find(defvec);

m2c = @(x)(x-1); % convert matlab to c indices
check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
calllib('bgmg', 'bgmg_set_tag_indices', 0, length(defvec), length(tag_indices), m2c(tag_indices));  check();
calllib('bgmg', 'bgmg_set_option', 0,  'r2min', r2min); check();
calllib('bgmg', 'bgmg_set_option', 0,  'kmax', kmax); check();
calllib('bgmg', 'bgmg_set_option', 0,  'max_causals', floor(max_causal_fraction * length(defvec))); check();  
calllib('bgmg', 'bgmg_set_option', 0,  'num_components', num_components); check();

for chr_index=1:length(chr_labels)
    plink_ld_mat_chr = strrep(plink_ld_mat,'@', sprintf('%i',chr_labels(chr_index)));
    fprintf('Loading %s...', plink_ld_mat_chr); tmp = load(plink_ld_mat_chr); fprintf('OK.\n');
    tmp.index_A = tmp.index_A + 1; 
    tmp.index_B = tmp.index_B + 1; % 0-based, comming from python
    fprintf('Importing %s to bgmg library...', plink_ld_mat_chr); 
    calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, length(tmp.r2), m2c(tmp.index_A), m2c(tmp.index_B), tmp.r2); fprintf('OK.\n'); check(); 
    clear('tmp', 'plink_ld_mat_chr');
end; clear('chr_index');

calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();
calllib('bgmg', 'bgmg_set_weights_randprune', 0, randprune_n, randprune_r2);  check();
hvec = ref.mafvec .* (1-ref.mafvec) * 2; calllib('bgmg', 'bgmg_set_hvec', 0, length(hvec), hvec);  check(); clear('hvec');

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
        result = BGMG_cpp_fit(zmat, nmat, options);
    else
        zvec = [trait1_data.zvec];
        nvec = [trait1_data.nvec];
        result = BGMG_cpp_fit(zvec, nvec, options);
    end

    result.trait1_file = trait1_file;
    result.trait2_file = trait2_file;
    result.reference_file = reference_file;
    result.options = options;
end

calllib('bgmg', 'bgmg_set_option', 0,  'diag', 0); check();

% Produce QQ plots with true params (only works for synthetic data, of course)
if QQ_PLOT_TRUE
    for trait_index = 1:(1 + ~isempty(trait2_file))
        if trait_index==1
            qq_params = struct('sig2_zero', 1, 'pi_vec', trait1_data.causal_pi, 'sig2_beta', trait1_data.sigsq);
            qq_data = trait1_data;
        else
            qq_params = struct('sig2_zero', 1, 'pi_vec', trait2_data.causal_pi, 'sig2_beta', trait2_data.sigsq);
            qq_data = trait2_data;
        end
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
        [figures, plot_data] = BGMG_cpp_qq_plot(qq_params, qq_data.zvec(defvec), options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_fit_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.fit.pdf', out_file, trait_index), '-dpdf')
    end
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
    dirs=dir('H:\work\simu_2018_03_12.bgmg.mat.tar\*.bgmg.mat');
    fileID = fopen(['bgmg_results.csv'], 'w');
    for i=1:length(dirs)
        x = load(fullfile('H:\work\simu_2018_03_12.bgmg.mat.tar', dirs(i).name));
        [header, data] = BGMG_util.result2str_point_estimates(x.result, x.result.options);
        if i==1,
            fprintf(fileID, 'trait1_file\ttrait2_file\treference_file\t%s\n', header);
        end
        fprintf(fileID, '%s\t%s\t%s\t%s\n', x.result.trait1_file, x.result.trait2_file, x.result.reference_file, data);
    end

    fclose(fileID);
end

% TBD: re-test confidence intervals estimation
