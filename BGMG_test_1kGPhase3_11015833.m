addpath('DERIVESTsuite');
USE_HAPMAP = false;
if USE_HAPMAP && ~exist('hapmap', 'var')
    load('/space/syn03/1/data/oleksandr/hapgen/snpIDs_11p3m.mat')                                              % snpIDs
    hapmap = load('/space/syn03/1/data/GWAS/SUMSTAT/LDSR/MATLAB_Annot/infomat.mat');                           % A1vec, A2vec, chrnumvec, posvec, snpidlist
    [is_in_11m, index_to_11m] = ismember(hapmap.snpidlist, snpIDs);
    mask_in_11m = false(num_snps, 1); mask_in_11m(index_to_11m(index_to_11m ~= 0)) = true;                     % 1184120 SNPs in the mask
end

if ~exist('LDr2_p8sparse', 'var')
load('/space/syn03/1/data/oleksandr/SynGen2/11015833/1000Genome_ldmat/ldmat_p8_BPwind10M_n503.mat'); LDr2_p8sparse=LDmat;        % LDr2_p8sparse
load('/space/md10/8/data/holland/genetics/LDhistAndTLD_1kGPhase3_hg19/LDr2hist_zontr2_p5_1kGHG1pc.mat')    % LDr2hist
load('/space/md10/8/data/holland/genetics/LDhistAndTLD_1kGPhase3_hg19/TLDr2_zontr2_p5_1kGHG1pc');          % tldvec
load('/space/md10/8/data/holland/genetics/LDhistAndTLD_1kGPhase3_hg19/mafvec_1kGPIII14p9m_n_HG1pc');       % mafvec
load('/space/syn03/1/data/oleksandr/SynGen2/11015833/chrpos_11015833.mat')                                 % chrnumvec, posvec

mafvec = mafvec(:);

GWAS_STUDIES = [];
GWAS_FOLDER = '/space/syn03/1/data/oleksandr/SynGen2/11015833/SUMSTAT/';
GWAS_NAMES = {'PGC_SCZ_2014', 'PGC_BIP_2016_qc'};
for gwas_index = 1:length(GWAS_NAMES)
	fname = fullfile(GWAS_FOLDER, GWAS_NAMES{gwas_index});
	fprintf('Loading %s...\n', fname);
	GWAS_STUDIES.(GWAS_NAMES{gwas_index}) = load(fname);
end

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

gwas_indices = [1 2];

for gwas_index = 1:size(gwas_indices, 1)
try
gwas_name1 = GWAS_NAMES{gwas_indices(gwas_index, 1)};
gwas_name2 = GWAS_NAMES{gwas_indices(gwas_index, 2)};
gwas1 = GWAS_STUDIES.(gwas_name1);
gwas2 = GWAS_STUDIES.(gwas_name2);

task = [];
task.zmat = [gwas1.zvec gwas2.zvec]; 
task.nmat   = [gwas1.nvec gwas2.nvec];
task.ref_ld = ref_ld;
task.Hvec   = Hvec;

defvec = true(num_snps, 1); fprintf('%i SNPs in the template\n', sum(defvec));
defvec = defvec & isfinite(Hvec) & isfinite(tldvec) & isfinite(sum(task.nmat + task.zmat, 2));  fprintf('%i SNPs left after filtering missing values (maf, tld, zmat, nmat, etc)\n', sum(defvec));
defvec = defvec & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX); fprintf('%i SNPs left after filtering large LD blocks (<= %.2f)\n', sum(defvec), LD_BLOCK_SIZE_MAX);
defvec = defvec & (tldvec <= TLD_MAX); fprintf('%i SNPs left after filtering high LD SNPs (<= %.2f)\n', sum(defvec), TLD_MAX);
defvec = defvec & ~mhcvec; fprintf('%i SNPs left after filtering MHC\n', sum(defvec));
defvec = defvec & (mafvec >= MAF_THRESH); fprintf('%i SNPs left after filtering low MAF SNPs (>= %.3f)\n', sum(defvec), MAF_THRESH);
if USE_HAPMAP
  defvec = defvec & mask_in_11m; fprintf('%i SNPs left after excluding all non-HapMap3 SNPs\n', sum(defvec));
end

%task.defvec = defvec;
task.defvec = true(sum(defvec), 1);
task.zmat = task.zmat(defvec, :);
task.nmat = task.nmat(defvec, :);
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

task.options = [];
task.options.total_het = total_het;  % Total heterozigosity
task.options.verbose = true;
task.options.ci_alpha = nan;
task.options.use_poisson = true;
task.options.title = sprintf('%s - %s', gwas_name1, gwas_name2);

disp(task.options.title)
disp(task.options)
disp(task)
close all
% save(sprintf('/home/oleksandr/space/SynGen2/11015833/BGMG_GWAS_results/%s_%i_r2bin%i.mat', gwas_name, sum(isfinite(task.zmat)), size(task.ref_ld.sum_r2, 2)), 'task');

% Perform fitting of the parameters
hits = sum(task.pruneidxmat, 2); w_ld = size(task.pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan;
result = BGMG_fit(task.zmat, task.Hvec, task.nmat, w_ld, task.ref_ld, task.options);

% Show bivariate results, and univariate results in context of poisson and gaussian mixtures
disp(result.bivariate.params)
sum(result.bivariate.params.pi_vec([1, 3]))
sum(result.bivariate.params.pi_vec([2, 3]))

disp(result.univariate{1}.params)
disp(result.univariate{2}.params)

for plot_trait_index = 1:2
[figures, plot_data] = UGMG_qq_plot(result.univariate{plot_trait_index}.params, task.zmat(:, plot_trait_index), task.Hvec, task.nmat(:, plot_trait_index), task.pruneidxmat, task.ref_ld, task.options);
save(sprintf('BGMG_plot_data_%s_%i_r2bin%i.mat',  GWAS_NAMES{gwas_indices(gwas_index, plot_trait_index)} , sum(isfinite(task.zmat(:, plot_trait_index))), size(task.ref_ld.sum_r2, 2)), 'plot_data');
print(figures.tot, sprintf('BGMG_%s_%i_r2bin%i_fitted.pdf',GWAS_NAMES{gwas_indices(gwas_index, plot_trait_index)} , sum(isfinite(task.zmat(:, plot_trait_index))), size(task.ref_ld.sum_r2, 2)),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('BGMG_%s_%i_r2bin%i_fitted_HL.pdf', GWAS_NAMES{gwas_indices(gwas_index, plot_trait_index)}, sum(isfinite(task.zmat(:, plot_trait_index))), size(task.ref_ld.sum_r2, 2)),'-dpdf')
end

catch e
    fprintf(1, 'The identifier was:\n%s',e.identifier);
    fprintf(1, 'There was an error! The message was:\n%s',e.message);
    disp('error');
end

end
