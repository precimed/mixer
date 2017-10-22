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
% For information about various parameters look into BGMG_fit.m.

fprintf('trait1: %s\n', trait1);
if exist('trait2', 'var'), fprintf('trait2: %s\n', trait2); end;

if ~exist('data_path', 'var')
t = 'H:\NORSTORE\MMIL\SUMSTAT\LDSR\MATLAB_Data';                       if exist(t, 'dir'), data_path = t; end;
t = '/work/users/oleksanf/NORSTORE/GWAS_SUMSTAT/LDSR/MATLAB_Data';     if exist(t, 'dir'), data_path = t; end;
t = '/space/syn03/1/data/GWAS_SUMSTAT/LDSR/MATLAB_Data';               if exist(t, 'dir'), data_path = t; end;
end
if ~exist('data_path', 'var'), error('Unable to locate data folder'); end;

addpath('DERIVESTsuite');

load('reference_data/mafvec.mat','mafvec');
load('reference_data/chrnumvec.mat','chrnumvec');
load('reference_data/posvec.mat','posvec');
ref_ld2 = load('reference_data/ref_l2.mat');
biased_ref_ld4 = load('reference_data/biased_ref_l4.mat');
biased_ref_ld2 = load('reference_data/biased_ref_l2.mat');
w_ld2 = load('reference_data/w_ld.mat');

ref_ld = struct('sum_r2', biased_ref_ld2.annomat, 'chi_r4', biased_ref_ld4.annomat ./ biased_ref_ld2.annomat);
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

data1  = load(fullfile(data_path, trait1)); data1.zvec(mhc) = nan;
if exist('trait2', 'var'),
    data2 = load(fullfile(data_path, trait2)); data2.zvec(mhc) = nan;
    result = BGMG_fit([data1.zvec, data2.zvec], Hvec, [data1.nvec, data2.nvec], w_ld, ref_ld, options);
    result.trait1 = trait1;
    result.trait2 = trait2;

    % Save the result in .mat file
    fname = sprintf('BGMG_run_%s_vs_%s', trait1, trait2);
    save([fname '.mat'], 'result');
else
    result = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options);
    result.trait1 = trait1;

    % Save the result in .mat file
    fname = sprintf('UGMG_run_%s', trait1);
    save([fname '.mat'], 'result');
end

fileID = fopen([fname '.txt'], 'w');
BGMG_util.result2str(1,      result)
BGMG_util.result2str(fileID, result)
fclose(fileID);
