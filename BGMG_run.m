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
fprintf('trait1: %s\n', trait1);
if exist('trait2', 'var'), fprintf('trait2: %s\n', trait2); end;

[~, trait1_name, ~] = fileparts(trait1);
if exist('trait2', 'var'), [~, trait2_name, ~] = fileparts(trait2); end;

if ~exist('data_path', 'var')
t = 'C:\Users\Oleksandr\Documents\GitHub\BGMG\reference_data';         if exist(t, 'dir'), data_path = t; end;
t = '/work/users/oleksanf/NORSTORE/GWAS_SUMSTAT/LDSR/MATLAB_Data';     if exist(t, 'dir'), data_path = t; end;
t = '/space/syn03/1/data/GWAS_SUMSTAT/LDSR/MATLAB_Data';               if exist(t, 'dir'), data_path = t; end;
end
if ~exist('data_path', 'var'), error('Unable to locate data folder'); end;

addpath('DERIVESTsuite');
addpath('wprctile');

if ~exist('reference_data', 'var'), reference_data = 'reference_data'; end;
fprintf('Loading reference data from %s...\n', reference_data);
load(fullfile(reference_data, 'mafvec.mat'),'mafvec');
load(fullfile(reference_data, 'chrnumvec.mat'),'chrnumvec');
load(fullfile(reference_data, 'posvec.mat'),'posvec');
ref_ld2 = load(fullfile(reference_data, 'ref_l2.mat'));
biased_ref_ld4 = load(fullfile(reference_data, 'biased_ref_l4.mat'));
biased_ref_ld2 = load(fullfile(reference_data, 'biased_ref_l2.mat'));
w_ld2 = load(fullfile(reference_data, 'w_ld.mat'));

r2bins = 20;
LDr2_binned = nan(length(chrnumvec), r2bins);
for r2bini=1:r2bins
    tmp = load([reference_data sprintf('_LDr2_hist/base-r2bin%i.1m.l2.mat', r2bini)]);
    LDr2_binned(:, r2bini) = tmp.annomat;
end

ref_ld = struct('sum_r2', ref_ld2.annomat, 'chi_r4', biased_ref_ld4.annomat ./ biased_ref_ld2.annomat, 'LDr2_binned', LDr2_binned);
Hvec = 2*mafvec .* (1-mafvec);  w_ld  = w_ld2.annomat;

% Remember to exclude MHC, also for summary stats called ending with
% no_MHC.mat. no_MHC suffix imply exclusion of 26-34, as done by LD score
% regresison. In some cases this seems to be not enough. 25-35 is the 
% standard from Ricopilli pipeline.
mhc = (chrnumvec==6) & (posvec > 25e6) & (posvec < 35e6);

% mafvec_all = load('mafvec_all_snps_1000G_Phase3_frq.mat')
% options.total_het = 2*sum(mafvec_all.mafvec .* (1-mafvec_all.mafvec))
options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs
options.verbose = true;
options.use_convolution = false;

data1  = load(fullfile(data_path, trait1)); data1.zvec(mhc, :) = nan; if isfield(data1, 'infovec'), data1.zvec(data1.infovec < 0.9, :) = nan; end; defvec = isfinite(data1.zvec); 
if exist('data1_neff', 'var'), data1.nvec = ones(size(data1.zvec)) * data1_neff; end;
if exist('trait2', 'var'),
    data2 = load(fullfile(data_path, trait2)); data2.zvec(mhc, :) = nan;  if isfield(data2, 'infovec'), data2.zvec(data2.infovec < 0.9, :) = nan; end; defvec = isfinite(data1.zvec + data2.zvec); 
    if exist('data2_neff', 'var'), data2.nvec = ones(size(data2.zvec)) * data2_neff; end;
end

weights = 1./w_ld;
weights(w_ld < 1) = 1;
weights(~isfinite(weights)) = 0;

if 1
    % Do binning on (N*H), TLD and chi2.
    bin_num = [10, 10, 100];
    bin_factor = {data1.nvec .* Hvec, ref_ld.sum_r2, ref_ld.chi_r4};

    %bins_idx = false(length(chrnumvec), prod(bin_num));
    bin_edges{1} = [-Inf, wquantile(bin_factor{1}, bin_num(1)-1, weights), +Inf];
    
    bined.ref_ld_chi_r4 = nan(prod(bin_num), 2);
    bined.ref_ld_sum_r2 = nan(prod(bin_num), 2);
    bined.Hvec = nan(prod(bin_num), 2);
    bined.Nvec = nan(prod(bin_num), 2);
    bined.bin_size = nan(prod(bin_num), 1);
    
    
    for bini1 = 1:bin_num(1)
        fprintf('.');
        idx1 = (bin_factor{1} >= bin_edges{1}(bini1)) & (bin_factor{1} < bin_edges{1}(bini1+1));
        bin_edges{2} = [-Inf, wquantile(bin_factor{2}(idx1), bin_num(2)-1, weights(idx1)), +Inf];
        for bini2 = 1:bin_num(2)
            idx2 = (bin_factor{2} >= bin_edges{2}(bini2)) & (bin_factor{2} < bin_edges{2}(bini2+1));
            bin_edges{3} = [-Inf, wquantile(bin_factor{3}(idx1&idx2), bin_num(3)-1, weights(idx1&idx2)), +Inf];
            for bini3 = 1:bin_num(3)
                idx3 = (bin_factor{3} >= bin_edges{3}(bini3)) & (bin_factor{3} < bin_edges{3}(bini3+1));
                %bins_idx(:, sub2ind(bin_num, bini1, bini2, bini3)) = idx1 & idx2 & idx3;
                bin_idx = idx1 & idx2 & idx3;
                
                i = sub2ind(bin_num, bini1, bini2, bini3);
                bined.bin_size(i) = sum(bin_idx);                
                % TBD: fill binned container with data.
                
                
                ref_ld.chi_r4(bin_idx) = mean(ref_ld.chi_r4(bin_idx));
                ref_ld.sum_r2(bin_idx) = mean(ref_ld.sum_r2(bin_idx));
                data1.nvec(bin_idx) = mean(data1.nvec(bin_idx));
                Hvec(bin_idx) = mean(Hvec(bin_idx));
            end
        end
    end
    fprintf('\n');
%    bins_idx = bins_idx(:, sum(bins_idx) ~= 0);  % make sure all bins are non-empty
    %assert(all(sum(bins_idx, 2)==1))             % make sure bins cover all SNPs
    %histogram(sum(bins))                        % inspect distribution of bin sizes

    if 0
    for bini = 1:size(bins_idx, 2)
        bin_idx = bins_idx(:, bini);
        ref_ld.chi_r4(bin_idx) = mean(ref_ld.chi_r4(bin_idx));
        ref_ld.sum_r2(bin_idx) = mean(ref_ld.sum_r2(bin_idx));
        data1.nvec(bin_idx) = mean(data1.nvec(bin_idx));
        Hvec(bin_idx) = mean(Hvec(bin_idx));
    end
    end
end

if exist('trait2', 'var'),
    result = BGMG_fit([data1.zvec, data2.zvec], Hvec, [data1.nvec, data2.nvec], w_ld, ref_ld, options);
    result.trait1 = trait1;
    result.trait2 = trait2;

    % Save the result in .mat file
    fname = sprintf('BGMG_run_%s-%s', trait1_name, trait2_name);
    save([fname '.mat'], 'result');
else
    result = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options);
    result.trait1 = trait1;
    result.univariate{1}.cdf = [];
    if 0
    % Infinitesimal model (pi1 constrained to 1)
    options_inf = options; options_inf.fit_infinitesimal=true;
    result_inf = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options_inf);
    result_inf.univariate{1}.cdf = [];

    % 3-component mixture (null + two causal components)
    options_mix2 = options; options_mix2.fit_two_causal_components=true;
    result_mix2 = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options_mix2);
    result_mix2.univariate{1}.cdf = [];
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
    [result.univariate{1}.cdf, figures] = UGMG_qq_plot(result.univariate{1}.params, data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options, result.univariate{1}.cdf);
    print(figures.tot, sprintf('%s', fname),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(figures.bin, sprintf('%s_HL', fname),'-dpdf')
if 0
    options.title = sprintf('%s (infinitesimal model)', trait1); options.title(options.title=='_')= '-';
    [result_inf.univariate{1}.cdf, figures] = UGMG_qq_plot(result_inf.univariate{1}.params, data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options, result_inf.univariate{1}.cdf);
    print(figures.tot, sprintf('%s_inf', fname),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3');
    print(figures.bin, sprintf('%s_inf_HL', fname),'-dpdf')

    fileID = fopen([fname '.inf.txt'], 'w');
    BGMG_util.result2str(1,      result_inf)
    BGMG_util.result2str(fileID, result_inf)
    fclose(fileID);
    
    options.title = sprintf('%s (null+small+large)', trait1); options.title(options.title=='_')= '-';
    [result_mix2.univariate{1}.cdf, figures] = UGMG_qq_plot(result_mix2.univariate{1}.params, data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options, result_mix2.univariate{1}.cdf);
    print(figures.tot, sprintf('%s_mix2', fname),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3');
    print(figures.bin, sprintf('%s_mix2_HL', fname),'-dpdf')
    
    fileID = fopen([fname '.mix2.txt'], 'w');
    BGMG_util.result2str(1,      result_mix2)
    BGMG_util.result2str(fileID, result_mix2)
    fclose(fileID);
end
end
