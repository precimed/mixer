addpath('DERIVESTsuite');
if ~exist('LDr2_p8sparse', 'var') 
load('/space/syn03/1/data/oleksandr/hapgen/LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat')  % LDr2_p8sparse
load('/space/syn03/1/data/oleksandr/hapgen/LDr2hist_zontr2_p5_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat')     % LDr2hist
load('/space/syn03/1/data/oleksandr/hapgen/setUpParameterResults_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat')  % parameterResults
load('/space/syn03/1/data/oleksandr/hapgen/p_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat')                      % pvecs (cell 10x4x3),  h2list, pi1list
load('/space/syn03/1/data/oleksandr/SynGen2/11015833/mafvec.mat')                                          % mafvec
load('/space/syn03/1/data/oleksandr/SynGen2/11015833/chrpos_11015833.mat')                                 % chrnumvec, posvec
load('/space/syn03/1/data/oleksandr/hapgen/snpIDs_11p3m.mat')                                              % snpIDs
load('/space/syn03/1/data/oleksandr/hapgen/TLDr2_zontr2_p1_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat');       % tldvec
hapmap = load('/space/syn03/1/data/GWAS/SUMSTAT/LDSR/MATLAB_Annot/infomat.mat');                           % A1vec, A2vec, chrnumvec, posvec, snpidlist 

num_snps = length(chrnumvec);  % 11015833

[is_in_11m, index_to_11m] = ismember(hapmap.snpidlist, snpIDs);
mask_in_11m = false(num_snps, 1); mask_in_11m(index_to_11m(index_to_11m ~= 0)) = true;  % 1184120 SNPs in the mask

total_het = nansum(2 .* mafvec .* (1-mafvec));
Hvec      = 2 .* mafvec .* (1 - mafvec);

ldr2binsNum     = size(LDr2hist,2);  %  100
ldr2binsEdges   = linspace( 0,1,ldr2binsNum+1);
ldr2binsCenters = nan(1,ldr2binsNum); for i=1:ldr2binsNum, ldr2binsCenters(i) = (ldr2binsEdges(i+1)+ldr2binsEdges(i))/2; end

minr2 = 0.05; minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
LDr2 = LDr2hist .* repmat(ldr2binsCenters, [num_snps, 1]);
LDr4 = LDr2hist .* repmat(ldr2binsCenters.^2, [num_snps, 1]);
numSNPsInLDr2_gt_r2min_vec = sum(LDr2hist(:,minr2bin:end),2);
end

r2_aggregated_bins = 4;
if ~exist('ref_ld') || (size(ref_ld.sum_r2, 2) ~= r2_aggregated_bins)
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
end

TLD_MAX=600;
LD_BLOCK_SIZE_MAX = 2000;
MAF_THRESH = 0.01;
mhcvec = chrnumvec==6 & posvec >= 25e6 & posvec <= 35e6;
USE_HAPMAP = true;
z = @(logpvec, zvec) -norminv(10.^-logpvec/2).*sign(zvec);

for rep_index = 1:10
for h2_index = [2 1 3]
for pi_index = 1:4
try

task = [];
task.zvec = z(-log10(pvecs{rep_index,pi_index,h2_index}), randn(num_snps, 1));
task.gwasParams = parameterResults{rep_index,pi_index,h2_index};
task.nvec   = 100000 * ones(size(task.zvec));
task.ref_ld = ref_ld;
task.Hvec   = Hvec;

defvec = true(num_snps, 1); fprintf('%i SNPs in the template\n', sum(defvec));
defvec = defvec & isfinite(Hvec) & isfinite(tldvec) & isfinite(task.zvec);  fprintf('%i SNPs left after filtering missing values (maf, tld, zvec, etc)\n', sum(defvec));
defvec = defvec & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX); fprintf('%i SNPs left after filtering large LD blocks (<= %.2f)\n', sum(defvec), LD_BLOCK_SIZE_MAX);
defvec = defvec & (tldvec <= TLD_MAX); fprintf('%i SNPs left after filtering high LD SNPs (<= %.2f)\n', sum(defvec), TLD_MAX);
defvec = defvec & ~mhcvec; fprintf('%i SNPs left after filtering MHC\n', sum(defvec));
defvec = defvec & (mafvec >= MAF_THRESH); fprintf('%i SNPs left after filtering low MAF SNPs (>= %.3f)\n', sum(defvec), MAF_THRESH);
if USE_HAPMAP
  defvec = defvec & mask_in_11m; fprintf('%i SNPs left after excluding all non-HapMap3 SNPs\n', sum(defvec));
end

%task.defvec = defvec;
task.defvec = true(sum(defvec), 1);
task.zvec = task.zvec(defvec);
task.nvec = task.nvec(defvec);
task.Hvec = task.Hvec(defvec);
task.ref_ld = struct('sum_r2', ref_ld.sum_r2(defvec, :), 'sum_r4_biased', ref_ld.sum_r4_biased(defvec, :), 'sum_r2_biased', ref_ld.sum_r2_biased(defvec, :));
task.ref_ld_bins = size(task.ref_ld.sum_r2, 2);

% Perform random pruning at LDr2 0.8 threshold....
nprune = 10;
fprintf('Perform %i iterations of random pruning ', nprune);
task.pruneidxmat = false(size(defvec,1), nprune);
for prune_repi=1:nprune,
    tmp=rand(size(chrnumvec,1),1);
    tmp(~defvec) = NaN;
    task.pruneidxmat(:,prune_repi) = isfinite(FastPrune(tmp, LDr2_p8sparse));
    fprintf('.');
end;
task.pruneidxmat = task.pruneidxmat(defvec, :);
fprintf('done.\nEffective number of SNPs on each iteration of random pruning:\n');
sum(task.pruneidxmat)

pivec_str = {'1e-5', '1e-4', '1e-3', '1e-2'};
h2vec_str = {'0.10', '0.40', '0.70'};

task.options = [];
task.options.total_het = total_het;  % Total heterozigosity across all SNPs in 1kG phase3
task.options.verbose = true;
task.options.ci_alpha = nan;
task.options.use_poisson = 1;
task.options.title = sprintf('HAPGEN pi=%s h2=%s rep=%i', pivec_str{pi_index}, h2vec_str{h2_index}, rep_index);

task.params = struct('sig2_zero', 1, 'pi_vec', task.gwasParams.pi1True, 'sig2_beta', task.gwasParams.sig2betaEst);
disp(task.options.title)
disp(task.options)
disp(task)
close all
save(sprintf('/home/oleksandr/space/SynGen2/11015833/BGMG_results/task_pi=%s_h2=%s_rep=%i_ldr2bins=%i.mat', pivec_str{pi_index}, h2vec_str{h2_index}, rep_index, task.ref_ld_bins), 'task');

[figures, plot_data] = UGMG_qq_plot(task.params, task.zvec, task.Hvec, task.nvec, task.pruneidxmat, task.ref_ld, task.options);
% To reproduce the same curve:
% plot(plot_data.data_logpvec, plot_data.hv_logp, plot_data.model_logpvec, plot_data.hv_logp)
save(sprintf('plot_data_%s_%i_r2bin%i.mat', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)), 'plot_data');
print(figures.tot, sprintf('%s_%i_r2bin%i.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_r2bin%i_HL.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')

if 0
    % Perform fitting of the parameters
    hits = sum(task.pruneidxmat, 2); w_ld = 1./hits; w_ld(hits==0) = nan;
    result = BGMG_fit(task.zvec, task.Hvec, task.nvec, w_ld, task.ref_ld, task.options);
    figures = UGMG_qq_plot(result.univariate{1}.params, task.zvec, task.Hvec, task.nvec, task.pruneidxmat, task.ref_ld, task.options);
    print(figures.tot, sprintf('%s_%i_r2bin%i_fitted.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(figures.bin, sprintf('%s_%i_r2bin%i_fitted_HL.pdf', task.options.title, sum(isfinite(task.zvec)), size(task.ref_ld.sum_r2, 2)),'-dpdf')
end

catch
    disp('error');
end

end;end;end

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
    subplot(3,4, pi_index + (h2_index-1) * 4); hold on;

    loaded = 0;
    files = dir(sprintf('C:/Storage/UGMG_results/plot_data/plot_data_HAPGEN pi=%s h2=%s rep=*_*_r2bin%i.mat', pi, h2, r2bins));
    for file = files'
        try
            load(['C:/Storage/UGMG_results/plot_data/' file.name]);
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

set(666,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(666, 'test.pdf', '-dpdf')

end


