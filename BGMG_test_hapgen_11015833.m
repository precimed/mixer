% You are expected to change active directry to the folder that contains this script
addpath('DERIVESTsuite');

if ~exist('reference_path', 'var'), reference_path = 'H:\work\HAPGEN_EUR_100K_11015883_reference_holland'; end;
if ~exist('USE_HAPMAP', 'var'), USE_HAPMAP = false; end;
if ~exist('USE_POISSON', 'var'), USE_POISSON = true; end;
if ~exist('USE_SAVED_PRUNING_INDICES', 'var'), USE_SAVED_PRUNING_INDICES = false; end;
if ~exist('r2_aggregated_bins', 'var'), r2_aggregated_bins = 4; end;
if ~exist('trait1_file', 'var'), trait1_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait1.mat'; end;
if ~exist('trait2_file', 'var'), trait2_file = 'H:\work\SIMU_HAPGEN_EUR_100K_11015883_traits\simu_h2=0.7_rg=0.0_pi1u=3e-04_pi2u=3e-04_pi12=9e-08_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity.trait2.mat'; end;
if ~exist('sample_size', 'var'), sample_size = 100000; end;
 
if USE_SAVED_PRUNING_INDICES && USE_HAPMAP
    error('USE_SAVED_PRUNING_INDICES is incompatible with USE_HAPMAP')
end

if USE_HAPMAP && ~exist('hapmap', 'var')
    load(fullfile(reference_path, 'snpIDs_11p3m.mat'))                                          % snpIDs
    hapmap = load(fullfile(reference_path, 'infomat_hapmap3.mat'));      % A1vec, A2vec, chrnumvec, posvec, snpidlist
    [is_in_11m, index_to_11m] = ismember(hapmap.snpidlist, snpIDs);
    mask_in_11m = false(length(snpIDs), 1); mask_in_11m(index_to_11m(index_to_11m ~= 0)) = true;                     % 1184120 SNPs in the mask
    defvec = mask_in_11m; save('H:\Dropbox\shared\BGMG\defvec_hapmap3.mat', 'defvec');
end

if ~USE_SAVED_PRUNING_INDICES && ~exist('LDr2_p8sparse', 'var') 
    load(fullfile(reference_path, 'LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat'))  % LDr2_p8sparse
end

if ~exist('LDr2hist', 'var')
    load(fullfile(reference_path, 'LDr2hist_zontr2_p5_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat'))     % LDr2hist
end

load(fullfile(reference_path, 'mafvec.mat'))                                                    % mafvec
load(fullfile(reference_path, 'chrpos_11015833.mat'))                                           % chrnumvec, posvec
load(fullfile(reference_path, 'TLDr2_zontr2_p1_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat'));       % tldvec

num_snps = length(chrnumvec);  % 11015833

total_het = nansum(2 .* mafvec .* (1-mafvec));
Hvec      = 2 .* mafvec .* (1 - mafvec);

ldr2binsNum     = size(LDr2hist,2);  %  100
ldr2binsEdges   = linspace( 0,1,ldr2binsNum+1);
ldr2binsCenters = nan(1,ldr2binsNum); for i=1:ldr2binsNum, ldr2binsCenters(i) = (ldr2binsEdges(i+1)+ldr2binsEdges(i))/2; end

minr2 = 0.05; minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
LDr2 = LDr2hist .* repmat(ldr2binsCenters, [num_snps, 1]);
LDr4 = LDr2hist .* repmat(ldr2binsCenters.^2, [num_snps, 1]);
numSNPsInLDr2_gt_r2min_vec = sum(LDr2hist(:,minr2bin:end),2);

r2_aggregated_bin_size = ldr2binsNum / r2_aggregated_bins;
assert(mod(ldr2binsNum,r2_aggregated_bins) == 0);
assert(r2_aggregated_bin_size > minr2bin);
sum_r2=zeros(num_snps, r2_aggregated_bins );
sum_r4=zeros(num_snps, r2_aggregated_bins );
for r2_bini=1:r2_aggregated_bins
    r2_bin_index = max(minr2bin, 1+(r2_bini-1)*r2_aggregated_bin_size) : (r2_bini*r2_aggregated_bin_size);
    sum_r2(:, r2_bini) = sum(LDr2(:, r2_bin_index), 2);
    sum_r4(:, r2_bini) = sum(LDr4(:, r2_bin_index), 2);
end
shape_param = sum_r4 ./ sum_r2;
ref_ld  = struct('sum_r2', sum_r2, 'sum_r4_biased', sum_r4, 'sum_r2_biased', sum_r2);

TLD_MAX=600;
LD_BLOCK_SIZE_MAX = 2000;
MAF_THRESH = 0.01;
mhcvec = chrnumvec==6 & posvec >= 25e6 & posvec <= 35e6;

trait1_data = load(trait1_file);
trait2_data = load(trait2_file);

defvec = true(num_snps, 1); fprintf('%i SNPs in the template\n', sum(defvec));
defvec = defvec & isfinite(Hvec) & isfinite(tldvec) & isfinite(trait1_data.zvec + trait2_data.zvec);  fprintf('%i SNPs left after filtering missing values (maf, tld, zvec, etc)\n', sum(defvec));
defvec = defvec & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX); fprintf('%i SNPs left after filtering large LD blocks (<= %.2f)\n', sum(defvec), LD_BLOCK_SIZE_MAX);
defvec = defvec & (tldvec <= TLD_MAX); fprintf('%i SNPs left after filtering high LD SNPs (<= %.2f)\n', sum(defvec), TLD_MAX);
defvec = defvec & ~mhcvec; fprintf('%i SNPs left after filtering MHC\n', sum(defvec));
defvec = defvec & (mafvec >= MAF_THRESH); fprintf('%i SNPs left after filtering low MAF SNPs (>= %.3f)\n', sum(defvec), MAF_THRESH);
% save('H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat', 'defvec')
if USE_HAPMAP
  defvec = defvec & mask_in_11m; fprintf('%i SNPs left after excluding all non-HapMap3 SNPs\n', sum(defvec));
end

if ~USE_SAVED_PRUNING_INDICES 
    % Perform random pruning at LDr2 0.8 threshold....
    nprune = 10;
    fprintf('Perform %i iterations of random pruning ', nprune);
    pruneidxmat = false(size(defvec,1), nprune);
    for prune_repi=1:nprune,
        tmp=rand(size(chrnumvec,1),1);
        tmp(~defvec) = NaN;
        pruneidxmat(:,prune_repi) = isfinite(FastPrune(tmp, LDr2_p8sparse));
        fprintf('.');
    end;
    fprintf('done.\n');
    
    % This is how pre-saved pruning indices were created:
    % save('/space/syn03/1/data/oleksandr/SynGen2/11015833/pruneidxmat_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat', 'pruneidxmat', 'defvec');
else
    saved_pruning_indices = load(fullfile(reference_path, 'pruneidxmat_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat'));
    if any(defvec ~= saved_pruning_indices.defvec), error('pre-generated random pruning indices are incompatible with current defvec'); end
    pruneidxmat = saved_pruning_indices.pruneidxmat;
end

fprintf('Effective number of SNPs on each iteration of random pruning:\n');
sum(pruneidxmat)
hits = sum(pruneidxmat, 2); w_ld = size(pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan;

trait1_data.zvec(~defvec) = nan;
trait2_data.zvec(~defvec) = nan;

options = [];
options.total_het = total_het;  % Total heterozigosity across all SNPs in 1kG phase3
options.verbose = true;
options.ci_alpha = nan;
options.use_poisson = USE_POISSON;
options.title = 'title';

sum_r2 = ref_ld.sum_r2;
sum_r4_biased = ref_ld.sum_r4_biased;
sum_r2_biased = ref_ld.sum_r2_biased;
save('H:\work\HAPGEN_EUR_100K_11015883_reference_holland_cleaned\reference.mat', 'sum_r2', 'sum_r2_biased', 'sum_r4_biased', 'w_ld', 'mafvec', 'total_het', 'chrnumvec', 'posvec');

result = BGMG_fit([trait1_data.zvec, trait2_data.zvec], Hvec, ones(num_snps, 2) * sample_size, w_ld, ref_ld, options);
result.trait1 = trait1;
result.trait2 = trait2;
result.pruneidxmat = pruneidxmat;
result.defvec = defvec;
result.options = options;

% Save the result in .mat file
fname = sprintf('BGMG_run_%s-%s', trait1_name, trait2_name);
if exist('out_file', 'var'), fname = out_file; end;
save([fname '.mat'], 'result');

close all




%save(sprintf('/home/oleksandr/space/SynGen2/11015833/BGMG_results/task_pi=%s_h2=%s_rep=%i_ldr2bins=%i.mat', pivec_str{pi_index}, h2vec_str{h2_index}, rep_index, task.ref_ld_bins), 'task');

task.options.plot_HL_bins = false;
%task.options.poisson_sigma_grid_nodes = 25;
%params.sig2_zero = task.params.sig2_zero;
%params.pi_vec = [task.params.pi_vec/2, task.params.pi_vec/2];
%params.sig2_beta = [2* task.params.sig2_beta, task.params.sig2_beta];
[figures, plot_data] = UGMG_qq_plot(task.params, task.zvec, task.Hvec, task.nvec, task.pruneidxmat, task.ref_ld, task.options);
%[figures, plot_data] = DATA_qq_plot(             task.zvec, task.Hvec,            task.pruneidxmat, task.ref_ld, task.options);
% To reproduce the same curve:
% plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
save(sprintf('plot_simu_data_%s_%i_r2bin%i.mat', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)), 'plot_data');
print(figures.tot, sprintf('%s_%i_r2bin%i.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')
%set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
%print(figures.bin, sprintf('%s_%i_r2bin%i_HL.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')

if 0
    % Perform fitting of the parameters
    hits = sum(task.pruneidxmat, 2); w_ld = size(task.pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan;
    task.options.use_poisson = false;
    result = BGMG_fit(task.zvec, task.Hvec, task.nvec, w_ld, task.ref_ld, task.options);
    [figures, plot_data] = UGMG_qq_plot(result.univariate{1}.params, task.zvec, task.Hvec, task.nvec, task.pruneidxmat, task.ref_ld, task.options);
%    plot_data_2018_02_19 = plot_data;
    save(sprintf('plot_data_%s_%i_r2bin%i_fitted_gaussian.mat', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)), 'plot_data');

    task.options.use_poisson = true;
    result = BGMG_fit(task.zvec, task.Hvec, task.nvec, w_ld, task.ref_ld, task.options);
    [figures, plot_data] = UGMG_qq_plot(result.univariate{1}.params, task.zvec, task.Hvec, task.nvec, task.pruneidxmat, task.ref_ld, task.options);
    save(sprintf('plot_data_%s_%i_r2bin%i_fitted_poisson.mat', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)), 'plot_data');

    %print(figures.tot, sprintf('%s_%i_r2bin%i_fitted.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')
    %set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    %print(figures.bin, sprintf('%s_%i_r2bin%i_fitted_HL.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')
end


% Load and display previously generated results
if 0

for pi_index = 1:4
for h2_index = 1:3

    pi = pivec_str{pi_index};
    h2 = h2vec_str{h2_index};
    r2bins = 4;

    RECALCULATE_MODEL = false;
    if RECALCULATE_MODEL  % re-calculate model prediction from tasks
    load(sprintf('C:/Storage/UGMG_tasks/task_pi=%s_h2=%s_rep=X_ldr2bins=%i.mat', pi, h2, r2bins));
    task.options.poisson_sigma_grid_limit = nan;
    task.options.poisson_sigma_grid_nodes = 10;
    task.options.plot_HL_bins = false;
    [figures, my_plot_data] = UGMG_qq_plot(task.params, task.zvec, task.Hvec, task.nvec, task.pruneidxmat, task.ref_ld, task.options);
    end

    figure(666); hold on;
    myAxis{pi_index, h2_index} = subplot(3,4, pi_index + (h2_index-1) * 4); hold on;

    loaded = 0;
    files = dir(sprintf('C:/Storage/UGMG_results_11M/plot_data_HAPGEN pi=%s h2=%s rep=*_9*_r2bin%i.mat', pi, h2, r2bins));
    for file = files'
        try
            load(['C:/Storage/UGMG_results_11M/' file.name]);
            loaded = loaded + 1;
            set(gca,'ColorOrderIndex',1)
            if RECALCULATE_MODEL
                plot(plot_data.data_logpvec, plot_data.hv_logp, my_plot_data.model_logpvec, my_plot_data.hv_logp);
            else
                plot(plot_data.data_logpvec, plot_data.hv_logp,    plot_data.model_logpvec,    plot_data.hv_logp);
            end
        catch
        end
    end
    if loaded == 0, error('No files found'); end;
    title(sprintf('pi=%s h2=%s r2bins=%i #rep=%i', pi, h2, r2bins, loaded));
    qq_options.qqlimy = 20;
    qq_options.qqlimx = 7;
    plot([0 qq_options.qqlimy],[0 qq_options.qqlimy], 'k--');
    xlim([0 qq_options.qqlimx]); ylim([0 qq_options.qqlimy]);
    drawnow
end
end

for pi_index = 1:4
for h2_index = 1:3
set(myAxis{pi_index,h2_index}, 'Position', myAxis{pi_index,h2_index}.Position - [0.005 * (pi_index - 1)  0.00 0 0])
end
end

set(666,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(666, 'test.pdf', '-dpdf')

end

if 0
   qq_options.qqlimy=20;
   close all; f=figure; 
    for pi_index=1:4
        for h2_index=1:3
            subplot(3,4,(pi_index)+4*(h2_index-1));hold on;
            is_file = exist(sprintf('SynGen11mBeta1_4Alex_15-Feb-2018_H_10_L_10_B_1/resultsSim_pi1i_%i_h2i_%i_betai_1.mat', pi_index, h2_index), 'file');
            if is_file, load(sprintf('SynGen11mBeta1_4Alex_15-Feb-2018_H_10_L_10_B_1/resultsSim_pi1i_%i_h2i_%i_betai_1.mat', pi_index, h2_index)); end;
            load(sprintf('2018_02_19_plot_hapgen/plot_data_HAPGEN pi=%s h2=%s rep=1_9510654_r2bin4.mat', pivec_str{pi_index}, h2vec_str{h2_index}));
            %figure;
            %plot([resstructEstF.data_logq'], resstructEstF.data_logp);
            %plot(plot_data.data_logpvec, plot_data.hv_logp);

            if is_file, plot(resstructEstF.data_logq', resstructEstF.data_logp, '-k', resstructEstF.modelFit_logq, resstructEstF.data_logp, ':b', resstructEstFtrue.modelTru_logq, resstructEstF.data_logp, '-b'); end;
%            plot([resstructEstF.data_logq'], resstructEstF.data_logp);
 %           plot([resstructEstFtrue.modelTru_logq], resstructEstF.data_logp);
            plot(plot_data.data_logpvec, plot_data.hv_logp, '-k', plot_data.model_logpvec, plot_data.hv_logp, '-r')

  %          plot(plot_data_2018_02_19.data_logpvec, plot_data_2018_02_19.hv_logp, plot_data_2018_02_19.model_logpvec, plot_data_2018_02_19.hv_logp)
            ylim([0 qq_options.qqlimy]); xlim([0 7])
            title(sprintf('pi=%s h2=%s',pivec_str{pi_index}, h2vec_str{h2_index}))
            plot([0 qq_options.qqlimy],[0 qq_options.qqlimy], 'k--');
            if (pi_index==1 && h2_index==1)
            legend('DH-data', 'DH-fit', 'DH-true', 'OF-data', 'OF-true', 'null')
            end
        end
    end
    set(f,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(f, 'comparison.pdf', '-dpdf')

end
