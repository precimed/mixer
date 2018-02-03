% This script can be run univariate or bivariate mixture analysis from command line like this:
%
%   matlab -nodisplay -nosplash -nodesktop -r "trait1='PGC_SCZ_2014'; trait2='LIPIDS_TG_2013'; BGMG_run; exit;"
%
% You may use this script with one trait (univariate analysis) or with two
% traits (bivariate analysis).
% The results are saved into a text file as well as into .mat file.
% Currently you may need to modify this file to pass additional
% parameters, later all parameters will be exposed as a text config file.
%
% You may also set trait1 to a file containing Nx2 zvec and nvec,
% this will trigger bivariate analysis even though trait2 is not set.
% This is handfull, for example, in simulations when two GWAS results
% are saved in one file.
%
% For information about various parameters look into BGMG_fit.m.
%try

if 0
    trait1 = 'PGC_BIP_2016_qc.mat';
    trait2 = 'PGC_SCZ_2014.mat';
    reference_data='reference_data_11M.mat';
    data_path='H:/GitHub/BGMG/reference_data_11M/SUMSTAT';
    USE_HAPMAP=true;
end

fprintf('trait1: %s\n', trait1);
if exist('trait2', 'var'), fprintf('trait2: %s\n', trait2); end;

[~, trait1_name, ~] = fileparts(trait1);
if exist('trait2', 'var'), [~, trait2_name, ~] = fileparts(trait2); end;

assert(exist('data_path', 'var') ~= 0);
assert(exist('reference_data', 'var') ~= 0);
addpath('DERIVESTsuite');

if ~exist('LDr2_p8sparse', 'var')
    fprintf('Loading reference data from %s...\n', reference_data);
    load(reference_data);
end

num_snps = length(chrnumvec);
TLD_MAX=600;
LD_BLOCK_SIZE_MAX = 2000;
MAF_THRESH = 0.01;
Hvec = 2*mafvec .* (1-mafvec);  
mhcvec = (chrnumvec==6) & (posvec > 25e6) & (posvec < 35e6);

data1  = load(fullfile(data_path, trait1));
if exist('data1_neff', 'var'), data1.nvec = ones(size(data1.zvec)) * data1_neff; end;
if exist('trait2', 'var'),
    data2 = load(fullfile(data_path, trait2)); 
    if exist('data2_neff', 'var'), data2.nvec = ones(size(data2.zvec)) * data2_neff; end;
end

defvec = true(num_snps, 1); fprintf('%i SNPs in the template\n', sum(defvec));
defvec = defvec & isfinite(data1.zvec + data1.nvec);  fprintf('%i SNPs left after filtering missing zvec and nvec values, trait %s\n', sum(defvec), trait1);
if exist('trait2', 'var'), defvec = defvec & isfinite(data2.zvec + data2.nvec);  fprintf('%i SNPs left after filtering missing zvec and nvec values, trait %s\n', sum(defvec), trait2); end;
defvec = defvec & isfinite(Hvec) & isfinite(tldvec);  fprintf('%i SNPs left after filtering missing values (maf, tld, etc)\n', sum(defvec));
defvec = defvec & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX); fprintf('%i SNPs left after filtering large LD blocks (<= %.2f)\n', sum(defvec), LD_BLOCK_SIZE_MAX);
defvec = defvec & (tldvec <= TLD_MAX); fprintf('%i SNPs left after filtering high LD SNPs (<= %.2f)\n', sum(defvec), TLD_MAX);
defvec = defvec & ~mhcvec; fprintf('%i SNPs left after filtering MHC\n', sum(defvec));
defvec = defvec & (mafvec >= MAF_THRESH); fprintf('%i SNPs left after filtering low MAF SNPs (>= %.3f)\n', sum(defvec), MAF_THRESH);
if exist('USE_HAPMAP', 'var') && USE_HAPMAP
  defvec = defvec & hapmap3_mask; fprintf('%i SNPs left after excluding all non-HapMap3 SNPs\n', sum(defvec));
end

nprune = 10;
if ~exist('last_prune_options', 'var') || ...
        last_prune_options.nprune ~= nprune || ...
        any(last_prune_options.defvec ~= defvec)
    fprintf('Perform %i iterations of random pruning ', nprune);
    pruneidxmat = false(size(defvec,1), nprune);
    for prune_repi=1:nprune,
        tmp=rand(size(chrnumvec,1),1);
        tmp(~defvec) = NaN;
        pruneidxmat(:,prune_repi) = isfinite(FastPrune(tmp, LDr2_p8sparse));
        fprintf('.');
    end;
    fprintf('done.\nEffective number of SNPs on each iteration of random pruning:\n');
    sum(pruneidxmat)
    last_prune_options.nprune = nprune;
    last_prune_options.defvec = defvec; 
else
    fprintf('Reuse pruning indices from previous run.\n');
end
hits = sum(pruneidxmat, 2); w_ld = size(pruneidxmat, 2) ./ hits; w_ld(hits==0) = nan;

options = [];
options.total_het = 2*sum(mafvec .* (1-mafvec));
options.verbose = true;
options.ci_alpha = nan;
options.use_poisson = false;
options.title = sprintf('%s - %s', trait1_name, trait2_name);

data1.zvec(~defvec) = nan;
if exist('trait2', 'var'),
    data2.zvec(~defvec) = nan;
    result = BGMG_fit([data1.zvec, data2.zvec], Hvec, [data1.nvec, data2.nvec], w_ld, ref_ld, options);
    result.trait1 = trait1;
    result.trait2 = trait2;
    result.pruneidxmat = pruneidxmat;
    result.defvec = defvec;
    result.options = options;

    % Save the result in .mat file
    fname = sprintf('BGMG_run_%s-%s', trait1_name, trait2_name);
    save([fname '.mat'], 'result');
else
    result = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options);
    result.trait1 = trait1;
    
    if 0
    % Infinitesimal model (pi1 constrained to 1)
    options_inf = options; options_inf.fit_infinitesimal=true;
    result_inf = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options_inf);

    % 3-component mixture (null + two causal components)
    options_mix2 = options; options_mix2.fit_two_causal_components=true;
    result_mix2 = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options_mix2);
    end
    
    % Save the result in .mat file
    if size(data1.zvec, 2) == 1
        fname = sprintf('UGMG_run_%s', trait1_name);
    else
        fname = sprintf('BGMG_run_%s', trait1_name);
    end
    save([fname '.mat'], 'result');
end

fileID = fopen([fname '.txt'], 'w');
BGMG_util.result2str(1,      result)
BGMG_util.result2str(fileID, result)
fclose(fileID);

if ~exist('trait2', 'var')
    options.title = trait1; options.title(options.title=='_')= '-';
    options.plot_HL_bins = false;
    [figures] = UGMG_qq_plot(result.univariate{1}.params, data1.zvec, Hvec, data1.nvec, pruneidxmat, ref_ld, options);
    print(figures.tot, sprintf('%s', fname),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(figures.bin, sprintf('%s_HL', fname),'-dpdf')

    if 0
    options.title = sprintf('%s (infinitesimal model)', trait1); options.title(options.title=='_')= '-';
    [figures] = UGMG_qq_plot(result_inf.univariate{1}.params, data1.zvec, Hvec, data1.nvec, pruneidxmat, ref_ld, options);
    print(figures.tot, sprintf('%s_inf', fname),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3');
    print(figures.bin, sprintf('%s_inf_HL', fname),'-dpdf')

    fileID = fopen([fname '.inf.txt'], 'w');
    BGMG_util.result2str(1,      result_inf)
    BGMG_util.result2str(fileID, result_inf)
    fclose(fileID);
    
    options.title = sprintf('%s (null+small+large)', trait1); options.title(options.title=='_')= '-';
    [figures] = UGMG_qq_plot(result_mix2.univariate{1}.params, data1.zvec, Hvec, data1.nvec, pruneidxmat, ref_ld, options);
    print(figures.tot, sprintf('%s_mix2', fname),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3');
    print(figures.bin, sprintf('%s_mix2_HL', fname),'-dpdf')
    
    fileID = fopen([fname '.mix2.txt'], 'w');
    BGMG_util.result2str(1,      result_mix2)
    BGMG_util.result2str(fileID, result_mix2)
    fclose(fileID);
    end
end

if 0
    rbp = result.bivariate.params;
    n1 = nanmedian(data1.nvec);
    n2 = nanmedian(data2.nvec);
    
    a=0; b=1;
    params_plus.pi_vec = rbp.pi_vec;
    params_plus.sig2_zero = a^2 * rbp.sig2_zero(1) + 2 * rbp.rho_zero * sqrt(a*b*prod(rbp.sig2_zero)) + b^2 * rbp.sig2_zero(2);
    params_plus.sig2_beta = [a^2 * n1*rbp.sig2_beta(1, 1), ...
        b^2 * n2*rbp.sig2_beta(2, 2)];  
     x= (a^2 * n1*rbp.sig2_beta(1, 3)) + 2 *a*b* rbp.rho_beta(1, 3) * ...
        sqrt(n1*n2*prod(rbp.sig2_beta(:, 3))) +(b^2)* n2*rbp.sig2_beta(2, 3);
     params_plus.sig2_beta(1,3) = x;

    params_minus.pi_vec = rbp.pi_vec;
    params_minus.sig2_zero = rbp.sig2_zero(1) - 2 * rbp.rho_zero * sqrt(prod(rbp.sig2_zero)) + rbp.sig2_zero(2);
    params_minus.sig2_beta = [n1*rbp.sig2_beta(1, 1), n2*rbp.sig2_beta(2, 2), n1*rbp.sig2_beta(1, 3) - 2 * rbp.rho_beta(1, 3) * sqrt(n1*n2*prod(rbp.sig2_beta(:, 3))) + n2*rbp.sig2_beta(2, 3)];

    options.poisson_sigma_grid_scale=1;
    [figures, plot_data] = UGMG_qq_plot(params_plus, a*data1.zvec+b*data2.zvec, Hvec, ones(size(Hvec)), pruneidxmat, ref_ld, options);

end