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

addpath('DERIVESTsuite');

if ~exist('reference_file', 'var'), error('reference_file is required'); end;
if ~exist('trait1_file', 'var'), error('trait1_file is required'); end;
if ~exist('trait2_file', 'var'), trait2_file = ''; end;
if ~exist('defvec_files', 'var'), defvec_files = {}; end;
if ~exist('LDmat_file', 'var'), LDmat_file = ''; end;
if ~exist('LDmat_file_variable', 'var'), LDmat_file_variable = 'LDr2_p8sparse'; end;
if ~exist('nprune', 'var'), nprune = 10; end;
if ~exist('trait1_nvec', 'var'), trait1_nvec = 100000; end;
if ~exist('trait2_nvec', 'var'), trait2_nvec = 100000; end;
if ~exist('out_file', 'var'), out_file = 'BGMG_result'; end;
if ~exist('init_file', 'var'), init_file = ''; end;

if ~exist('USE_POISSON', 'var'), USE_POISSON = true; end;
if ~exist('USE_POISSON_BGMG', 'var'), USE_POISSON_BGMG = false; end;

if ~exist('DO_FIT', 'var'), DO_FIT = true; end;                % perform fitting
if ~exist('QQ_PLOT_TRUE', 'var'), QQ_PLOT_TRUE = false; end;   % make QQ plots with true parameters
if ~exist('QQ_PLOT_FIT', 'var'), QQ_PLOT_FIT = false; end;     % make QQ plots with fitted parameters
if ~exist('POWER_PLOT_FIT', 'var'), POWER_PLOT_FIT = false; end;  % make power plots with fitted parameters
if ~exist('BGMG_RELAX_ALL', 'var'), BGMG_RELAX_ALL = false; end;
if ~exist('TITLE', 'var'), TITLE = 'title'; end;
if ~exist('CI_ALPHA', 'var'), CI_ALPHA = nan; end;
if ~exist('BIVARIATE_PENALTY_FACTOR', 'var'), BIVARIATE_PENALTY_FACTOR = 10; end;
if ~exist('SCAD_PENALTY', 'var'), SCAD_PENALTY = nan; end;
if ~exist('MAF_THRESH', 'var'), MAF_THRESH = nan; end;
if ~exist('FORCE_W_LD', 'var'), FORCE_W_LD = false; end;

if ~exist('plot_HL_bins', 'var'), plot_HL_bins = false; end;

% defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'}
% defvec_files = {'H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'}
% defvec_files = {'H:\Dropbox\shared\BGMG\defvec_1kG_phase3_EUR.mat' };
% defvec_files = {'H:\Dropbox\shared\BGMG\defvec_1kG_phase3_EUR.mat', 'H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'}
% LDmat_file = 'H:\NORSTORE\oleksanf\11015833\hapgen\LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat';
% LDmat_file = 'H:\Dropbox\shared\BGMG\ldmat_p8_BPwind10M_n503.mat'; LDmat_file_variable='LDmat';
% trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'
% trait1_file = 'H:\Dropbox\shared\BGMG\p_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc\simu_h2=0.40_pi1u=1e-3_rep=1.ugmg.mat'
% trait1_file = 'H:\GitHub\BGMG\GIANT_HEIGHT_2014_lift.mat';
% trait1_file = 'H:\GitHub\BGMG\PGC_SCZ_2014.mat';
% reference_file = 'H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins_wld.mat'
% reference_file = 'H:\Dropbox\shared\BGMG\1kG_phase3_EUR_11015883_reference_p01_SNPwind50k_per_allele_4bins.mat'
% DO_FIT = false; QQ_PLOT_TRUE = true;
% init_file = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.bgmg.mat';
% trait1_file  = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat';
% trait2_file  = 'H:\work\SIMU_BGMG_random_pi12\test\simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat';

if ~isempty(init_file) && isempty(trait2_file), error('init_file makes sence only for bivariate'); end;
if ~isempty(init_file) && isfinite(CI_ALPHA), error('init_file is incompatible with CI_ALPHA'); end;
if ~isempty(LDmat_file) && FORCE_W_LD, error('LDmat_file must not be used together with FORCE_W_LD'); end;
if POWER_PLOT_FIT && ~isempty(trait2_file), error('power plots are only supported for univariate traits'); end;
if POWER_PLOT_FIT && ~DO_FIT, error('power plots are only supported with DO_FIT=true'); end;

clear('pruneidxmat_or_w_ld', 'pruneidxmat', 'w_ld', 'hvec');
fprintf('Loading reference from %s...\n', reference_file);
load(reference_file); num_snps = length(mafvec);
if exist('pruneidxmat', 'var') && isempty(LDmat_file) && ~FORCE_W_LD
    fprintf('w_ld was changed to random pruning-based weighting\n');
    hits = sum(pruneidxmat, 2); w_ld = size(pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan;
end

if ~exist('hvec', 'var'),
    hvec = 2 .* mafvec .* (1-mafvec);
    warning('hvec does not exist in the reference; hvec is set based on mafvec.');
end

if exist('w_ld', 'var') && (size(w_ld, 2) > 1), w_ld = sum(w_ld, 2); end;
if ~exist('w_ld', 'var') && isempty(LDmat_file), error('Binary pruning matrix is required for reference without w_ld'); end;
if ~isempty(LDmat_file), w_ld = zeros(num_snps, 1); end;  % ignore w_ld because we will re-generate it based on random pruning
defvec = true(num_snps, 1);
fprintf('%i SNPs in the reference\n', num_snps);
ref_ld  = struct('sum_r2', sum_r2, 'sum_r4_biased', sum_r4_biased, 'sum_r2_biased', sum_r2_biased);
cur_defvec = struct('defvec', isfinite(sum(sum_r2, 2) + sum(sum_r4_biased, 2) + sum(sum_r2_biased, 2) + hvec + mafvec + w_ld));
defvec = defvec & cur_defvec.defvec;
fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec));

for i=1:length(defvec_files)
    fprintf('Loading %s...\n', defvec_files{i});
    cur_defvec = load(defvec_files{i});
    defvec = defvec & cur_defvec.defvec;
    fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec));
end

if isfinite(MAF_THRESH), defvec = defvec & (mafvec >= MAF_THRESH);
    fprintf('Exclude %i variants due to mafvec (%i variants remain)\n', sum(mafvec < MAF_THRESH), sum(defvec));
end

fprintf('Loading %s...\n', trait1_file);
trait1_data = load(trait1_file);
if ~isfield(trait1_data, 'nvec'), trait1_data.nvec = ones(size(trait1_data.zvec)) * trait1_nvec; end;
cur_defvec.defvec = isfinite(trait1_data.zvec + trait1_data.nvec);
defvec = defvec & cur_defvec.defvec;
fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec));

if ~isempty(trait2_file),
    fprintf('Loading %s...\n', trait2_file);
    trait2_data = load(trait2_file);
    if ~isfield(trait2_data, 'nvec'), trait2_data.nvec = ones(size(trait2_data.zvec)) * trait2_nvec; end;
    cur_defvec.defvec = isfinite(trait2_data.zvec + trait2_data.nvec);
    defvec = defvec & cur_defvec.defvec;
    fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec));
end

if ~isempty(LDmat_file)
    fprintf('Loading %s...\n', LDmat_file);
    LDmat = load(LDmat_file);
    if isfield(LDmat, 'mafvec')
        cur_defvec.defvec = isfinite(LDmat.mafvec);
        defvec = defvec & cur_defvec.defvec;
        fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec));
    end
    
    % perform random pruning
    fprintf('Perform %i iterations of random pruning ', nprune);
    pruneidxmat = false(size(defvec,1), nprune);
    for prune_repi=1:nprune,
        tmp=rand(size(chrnumvec,1),1);
        tmp(~defvec) = NaN;
        pruneidxmat(:,prune_repi) = isfinite(FastPrune(tmp, LDmat.(LDmat_file_variable)));
        fprintf('.');
    end;
    fprintf('done.\n');

    fprintf('Effective number of SNPs on each iteration of random pruning:\n');
    sum(pruneidxmat)

    fprintf('w_ld was changed to random pruning-based weighting\n');
    hits = sum(pruneidxmat, 2); w_ld = size(pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan;
    
    cur_defvec.defvec = isfinite(w_ld);
    defvec = defvec & cur_defvec.defvec;
    fprintf('Exclude %i variants (%i variants remain)\n', sum(~cur_defvec.defvec), sum(defvec));
end

if exist('pruneidxmat', 'var'),
    pruneidxmat_or_w_ld = pruneidxmat;
else
    pruneidxmat_or_w_ld = w_ld;
end;

trait1_data.zvec(~defvec) = nan;
if ~isempty(trait2_file), trait2_data.zvec(~defvec) = nan; end;

options = [];
options.total_het = total_het;
options.verbose = true;
options.ci_alpha = CI_ALPHA;
options.use_poisson = USE_POISSON;
options.use_poisson_bgmg = USE_POISSON_BGMG;
options.title = TITLE;
options.relax_all = BGMG_RELAX_ALL;
options.bivariate_penalty_factor = BIVARIATE_PENALTY_FACTOR;
options.scad_penalty = SCAD_PENALTY;

if ~isempty(init_file)
    fprintf('Loading %s...\n', init_file);
    init_point = load(init_file);
    options.params1 = init_point.result.univariate{1}.params;
    options.params2 = init_point.result.univariate{2}.params;
end

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
        result = BGMG_fit([trait1_data.zvec, trait2_data.zvec], hvec, [trait1_data.nvec, trait2_data.nvec], w_ld, ref_ld, options);
    else
        result = BGMG_fit(trait1_data.zvec, hvec, trait1_data.nvec, w_ld, ref_ld, options);
    end

    result.trait1_file = trait1_file;
    result.trait2_file = trait2_file;
    result.reference_file = reference_file;
    result.options = options;
end

% Produce power plots (Proportion of narrow-sense chip heritability, S(N),
% captured by genome-wide significant SNPs as a function of sample size)
if POWER_PLOT_FIT && DO_FIT
    trait_index=1;
    power_plot_data = UGMG_power_plot(result.univariate{trait_index}.params, trait1_data.zvec, hvec, trait1_data.nvec, ref_ld, options);
    result.univariate{trait_index}.power_plot_fit_data = power_plot_data;
    power_plot_figure = figure; hold on;
    plot(log10(power_plot_data.nvalues), power_plot_data.prop_h2, '-');
    ax = gca; ax.ColorOrderIndex = 1;
    plot(log10(power_plot_data.nvalue_current), power_plot_data.prop_h2_current, 'o');
    set(gca, 'xtick', [3, 4, 5, 6, 7, 8]);
    xlabel('log_1_0(N)');
    ylabel('S(N)');
    legend(sprintf('N=%.0f, S(N)=%.1f%%', power_plot_data.nvalue_current, 100*power_plot_data.prop_h2_current), 'Location', 'NorthWest')
    trait1_name = strsplit(trait1_file, {'\\', '/'}); trait1_name=trait1_name{end};trait1_name(trait1_name=='_')='-';
    title(trait1_name);
    print(power_plot_figure, sprintf('%s.trait%i.power.fit.pdf', out_file, trait_index), '-dpdf')
end

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
        options.plot_HL_bins = plot_HL_bins;
        options.mafvec = mafvec;
        [figures, plot_data] = UGMG_qq_plot(qq_params, qq_data.zvec, hvec, qq_data.nvec, pruneidxmat_or_w_ld, ref_ld, options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_true_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.true.pdf', out_file, trait_index), '-dpdf')

        if (plot_HL_bins)
            set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
            print(figures.bin, sprintf('%s.trait%i.qq-bin.true.pdf', out_file, trait_index),'-dpdf')
        end
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
        options.plot_HL_bins = plot_HL_bins;
        options.mafvec = mafvec;
        [figures, plot_data] = UGMG_qq_plot(qq_params, qq_data.zvec, hvec, qq_data.nvec, pruneidxmat_or_w_ld, ref_ld, options);

        % To reproduce the same curve: plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
        result.univariate{trait_index}.qq_plot_fit_data = plot_data;
        print(figures.tot, sprintf('%s.trait%i.qq.fit.pdf', out_file, trait_index), '-dpdf')

        if (plot_HL_bins)
            set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
            print(figures.bin, sprintf('%s.trait%i.qq-bin.fit.pdf', out_file, trait_index),'-dpdf')
        end
    end
end

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
 
if 0
    pp={};legend_text ={};
    pp{end+1,1} =load('GIANT_HEIGHT_2014_lift.mat.mat'); legend_text{end+1, 1} = 'HEIGHT';
    pp{end+1,1}=load('PGC_SCZ_2014.result.mat');legend_text{end+1, 1} = 'SCZ';
    pp{end+1,1}=load('LIPIDS_HDL_2013.mat.mat');legend_text{end+1, 1} = 'HDL';
    figure; hold on;
    for i=1:length(pp), pp{i}=pp{i}.result.univariate{1}.power_plot_fit_data; end;
    for i=1:length(pp), plot(log10(pp{i}.nvalues), pp{i}.prop_h2, '-'); end
    ax = gca; ax.ColorOrderIndex = 1;
    for i=1:length(pp), plot(log10(pp{i}.nvalue_current), pp{i}.prop_h2_current, 'o'); end
    xlabel('log_1_0(N)'); ylabel('S(N)');
    set(gca, 'xtick', [3, 4, 5, 6, 7, 8]);
    for i=1:length(pp), legend_text{i}=sprintf('%s: N=%.0f, S(N)=%.1f%%',legend_text{i}, pp{i}.nvalue_current, 100*pp{i}.prop_h2_current); end
    legend(legend_text, 'Location', 'NorthWest')
end
