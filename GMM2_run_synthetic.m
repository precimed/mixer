% This script can be run univariate or bivariate mixture analysis from command line like this:
%
%   matlab -nodisplay -nosplash -nodesktop -r "file='gwaslogp_h2=0.50_pi=1e-03_pi12=1e-04_rho12=0.50_n12=100000_tp3cf96be5_403d_4411_9b4e_317dcaefcab9.mat'; GMM2_run_synthetic; exit;"
%
% You may use this script with one trait (univariate analysis) or with two
% traits (bivariate analysis).
% The results are saved into a text file as well as into .mat file.
% Currently you may need to modify this file to pass additional
% parameters, later all parameters will be exposed as a text config file.
%
% For information about various parameters look into GMM2_fit.m.

fprintf('file: %s\n', file);
load(file);  % -> frame

outdir = 'figures'; if ~exist(outdir,'dir'), mkdir(outdir); end

addpath('DERIVESTsuite');

% Load data
t = 'H:\NORSTORE\GWAS_SUMSTAT\LDSR';                       if exist(t, 'dir'), ldsr_path = t; end;
t = '/work/users/oleksanf/NORSTORE/GWAS_SUMSTAT/LDSR';     if exist(t, 'dir'), ldsr_path = t; end;
t = '/space/syn03/1/data/GWAS_SUMSTAT/LDSR';               if exist(t, 'dir'), ldsr_path = t; end;
if ~exist('ldsr_path', 'var'), error('Unable to locate data folder'); end;

load(fullfile(ldsr_path, 'MATLAB_Annot/ldmat_1m_p8.mat'),'mafvec');
load(fullfile(ldsr_path, 'MATLAB_Annot/infomat.mat'));
ref_ld2 = load(fullfile(ldsr_path, 'MATLAB_Annot/1000G_EUR_Phase3_ref_ld2.mat'));
ref_ld4 = load(fullfile(ldsr_path, 'MATLAB_Annot/1000G_EUR_Phase3_ref_ld4.mat'));
w_ld2 = load(fullfile(ldsr_path, 'MATLAB_Annot/1000G_EUR_Phase3_w_ld2.mat'));

ref_ld = struct();
ref_ld.sum_r2 = ref_ld2.annomat;
ref_ld.sum_r4 = ref_ld4.annomat;

z = @(logpvec, zvec) -norminv(10.^-logpvec/2).*sign(zvec);
Hvec = 2*mafvec .* (1-mafvec);  w_ld  = w_ld2.annomat;
mhc = (chrnumvec==6) & (posvec > 25e6) & (posvec < 35e6);
frame.gwaszscore(mhc, :) = nan;

TotalHET = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs

% Run the model
options = struct('verbose', false, 'plot_costlines', false, 'plot_costmaps', false, 'total_het', TotalHET);

result = GMM2_fit(frame.gwaszscore, Hvec, frame.gwassize, w_ld, ref_ld, options);
result.trait1 = 'trait1';
result.trait2 = 'trait2';
result.file = file;

% Save the result in .mat file
fname = sprintf('GMM2_run_%s', file);
save([fname '.mat'], 'result');

fileID = fopen([fname '.txt'], 'w');
% Describe results in a plain text file
for i=1:length(result.univariate)
    if i==1, fprintf(fileID,'Univariate analysis for %s:\n', result.trait1); else  fprintf(fileID,'Univariate analysis for %s:\n', result.trait2); end;
    u = result.univariate{i}.params.uncertainty;
    fn = fieldnames(u);
    for j=1:length(fn)
        p = u.(fn{j});
        fprintf(fileID,'\t%s=%.3g, se=%.3g, ci@95%%=[%.3g,%.3g], wald=%.3f, H0="%s", Pvalue=%.3g\n', fn{j}, p.value, p.se, p.ci(1), p.ci(2), p.wald, p.H0, p.pvalue);
    end
end
if isfield(result, 'bivariate')
    fprintf(fileID,'Bivariate analysis for %s vs %s:\n', result.trait1, result.trait2);
    u = result.bivariate.params.uncertainty;
    fn = fieldnames(u);
    for j=1:length(fn)
        p = u.(fn{j});
        fprintf(fileID,'\t%s=%.3g, se=%.3g, ci@95%%=[%.3g,%.3g], wald=%.3f, H0="%s", Pvalue=%.3g\n', fn{j}, p.value, p.se, p.ci(1), p.ci(2), p.wald, p.H0, p.pvalue);
    end
end
fclose(fileID);
