addpath('DERIVESTsuite');
addpath('H:\NORSTORE\oleksanf\holland_genetics\GPSIM');  % for GenStats_meta
addpath('H:\NORSTORE\oleksanf\holland_matlab\qq');

%workspace_data_file = 'data_2017_11_15.mat';
workspace_data_file = 'data_2017_12_06.mat';
if 0
    % Test consistency with Dominic's model
    info1m = load('H:\NORSTORE\MMIL\SUMSTAT\LDSR\MATLAB_Annot\infomat.mat');
    info9m = load('H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Annot\infomat.mat', 'chrnumvec', 'posvec', 'snpidlist');
    chrnumvec = info9m.chrnumvec; posvec = info9m.posvec;

    [is_in_9m, index_to_9m] = ismember(info1m.snpidlist, info9m.snpidlist);
    mask_in_9m = false(size(info9m.chrnumvec)); mask_in_9m(index_to_9m(index_to_9m ~= 0)) = true;

    % Load reference data from Dominic, align to 9M template
    load('H:\NORSTORE\oleksanf\holland_genetics\LDhistAndTLD_1kGPhase3_hg19\LDr2hist_zontr2_p5_1kGHG1pc.mat')
    load('H:\NORSTORE\oleksanf\holland_genetics\LDhistAndTLD_1kGPhase3_hg19\mafvec_1kGPIII14p9m_n_HG1pc.mat', 'mafvec')
    load('H:\NORSTORE\oleksanf\holland_genetics\LDhistAndTLD_1kGPhase3_hg19\TLDr2_zontr2_p5_1kGHG1pc.mat')
    load('H:\NORSTORE\oleksanf\holland_genetics\LDhistAndTLD_1kGPhase3_hg19\info9m_n_1kGHG1pc_RemappingIndices.mat');

    defvecForAlex1m_new = load('defvecForAlex1m_new.mat');
    defvecForAlexBIP6m = load('defvecForAlexBIP6m.mat');
    defvecForAlexScz5m = load('defvecForAlexScz5m.mat');

    nsnps9m = length(info9m.chrnumvec);
    mafvec_9m   = nan(nsnps9m,1);
    tldvec_9m   = nan(nsnps9m,1);
    LDr2hist_9m = zeros(nsnps9m,size(LDr2hist,2)); %   % FIX UP AND UNCOMMENT!
    LDr2hist_9m(:,end) = 1;                        % Initialize to ones: "missing" SNPs will have r^2 with self equal to 1.  % FIX UP AND UNCOMMENT!

    tldvec_9m(i9m_9m_n_1kGHG1pc)     = tldvec(i1kGHG1pc_9m_n_1kGHG1pc);
    mafvec_9m(i9m_9m_n_1kGHG1pc)     = mafvec(i1kGHG1pc_9m_n_1kGHG1pc);
    LDr2hist_9m(i9m_9m_n_1kGHG1pc,:) = LDr2hist(i1kGHG1pc_9m_n_1kGHG1pc,:);  % FIX UP AND UNCOMMENT!
    Hvec_9m = 2*mafvec_9m.*(1-mafvec_9m);
    clear('tldvec', 'mafvec', 'LDr2hist');

    ldr2binsNum     = size(LDr2hist_9m,2);  %  100
    ldr2binsEdges   = linspace( 0,1,ldr2binsNum+1);
    ldr2binsCenters = nan(1,ldr2binsNum); for i=1:ldr2binsNum, ldr2binsCenters(i) = (ldr2binsEdges(i+1)+ldr2binsEdges(i))/2; end

    %minr2 = 0.05; % If r^2<0.05, ignore --> noise. USUALLY HAVE USED THIS ONE.
    minr2 = 0.1;  % If r^2<0.05, ignore --> noise.
    minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    numSNPsInLDr2_gt_r2min_vec = sum(LDr2hist_9m(:,minr2bin:end),2);
    
    % calculate chi2 for two r2 bins ("low" and "high" r2)
    minr2 = 0.05; minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    r2_bin_index = false(2, ldr2binsNum);
    r2_bin_index(1, minr2bin:(ldr2binsNum/2)) = true;
    r2_bin_index(2, ((ldr2binsNum/2) + 1):end) = true;
    
    LDr2tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters, [size(LDr2hist_9m, 1), 1]);
    LDr4tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters.^2, [size(LDr2hist_9m, 1), 1]);
    
    alex_sum_r2=zeros(nsnps9m, 2);
    alex_sum_r4=zeros(nsnps9m, 2);
    for r2_bini=1:2
        alex_sum_r2(:, r2_bini) = sum(LDr2tot_9m(:, r2_bin_index(r2_bini, :)), 2);
        alex_sum_r4(:, r2_bini) = sum(LDr4tot_9m(:, r2_bin_index(r2_bini, :)), 2);
    end
    
    fullFileName = 'H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Data\PGC2_SCZ52_multisites.mat';    
    scz_meta = load(fullFileName,'logpmat','logormat','N_A_vec','N_U_vec');
    scz_meta.Neffvec = 4./(1./scz_meta.N_A_vec+1./scz_meta.N_U_vec);   % CHANGE THE 2 to 4.    % Neffvec is the 1x52 vector of sample sizes over the 52 sub-studies.
    scz_meta.zmat = -norminv(10.^-scz_meta.logpmat/2).*sign(scz_meta.logormat);
    [scz.zvec_meta, scz.logpvec_meta, scz.Neff] = GenStats_meta(scz_meta.zmat,scz_meta.Neffvec);   % Neff = sum(Neffvec);  where Neffvec is the vector of substudy sample sizes.
    scz.Neff = repmat(scz.Neff, size(scz.zvec_meta));
    clear('scz_meta');

    fullFileName = 'H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Data\PGC2_BIP32_multisites_9m.mat';
    bip_meta = load(fullFileName,'logpmat','betamat','Nvec_A','Nvec_U');
    bip_meta.Neffvec = 4./(1./bip_meta.Nvec_A+1./bip_meta.Nvec_U);  %  sum(Neffvec) == 46582
    bip_meta.zmat = -norminv(10.^-bip_meta.logpmat/2).*sign(bip_meta.betamat);
    [bip.zvec_meta, bip.logpvec_meta, bip.Neff] = GenStats_meta(bip_meta.zmat,bip_meta.Neffvec);   % Neff = sum(Neffvec);  where Neffvec is the vector of substudy sample sizes.
    bip.Neff = repmat(bip.Neff, size(bip.zvec_meta));
    clear('bip_meta');

    scz_noMETA = load('H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Data\PGC2_SCZ52.mat');
    scz.zvec = -norminv(10.^-scz_noMETA.logpvec_pgc2_scz/2).*sign(scz_noMETA.zvec_pgc2_scz);
    scz.zvec(~isfinite(scz.zvec_meta)) = nan;

    bip_noMETA = load('H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Data\PGC_BIP32_9m.mat');
    bip.zvec = -norminv(10.^-bip_noMETA.logpvec_pgc_bip32/2).*sign(bip_noMETA.zvec_pgc_bip32);
    bip.zvec(~isfinite(bip.zvec_meta)) = nan;
    
    bmi    = load('H:\NORSTORE\MMIL\SUMSTAT\MAT_9M\GIANT_BMI_2015_EUR_lift.mat');
    height = load('H:\NORSTORE\MMIL\SUMSTAT\MAT_9M\GIANT_HEIGHT_2014_lift.mat');
    scz3   = load('H:\NORSTORE\MMIL\SUMSTAT\MAT_9M\PGC_SCZ_0917.mat');

    % Align Alex's reference data
    if 0
        ref_ld2 = load(fullfile('reference_data', 'ref_l2.mat'));
        biased_ref_ld4 = load(fullfile('reference_data', 'biased_ref_l4.mat'));
        biased_ref_ld2 = load(fullfile('reference_data', 'biased_ref_l2.mat'));
        w_ld2 = load(fullfile('reference_data', 'w_ld.mat'));
        mafvec_1m = load(fullfile('reference_data', 'mafvec.mat'),'mafvec');
        scz_1m_withMHC = load('H:\work\bgmg-annot-1000G-Phase3-EUR\1m_with_MHC\PGC_SCZ_2014.mat');
        bip_1m_withMHC = load('H:\work\bgmg-annot-1000G-Phase3-EUR\1m_with_MHC\PGC_BIP_2016_qc.mat');

        index_to_9m_nnz = index_to_9m(index_to_9m ~= 0);
        alex.ref_ld2 = nan(size(chrnumvec)); alex.ref_ld2(index_to_9m_nnz) = ref_ld2.annomat(is_in_9m);
        alex.biased_ref_ld4 = nan(size(chrnumvec)); alex.biased_ref_ld4(index_to_9m_nnz) = biased_ref_ld4.annomat(is_in_9m);
        alex.biased_ref_ld2 = nan(size(chrnumvec)); alex.biased_ref_ld2(index_to_9m_nnz) = biased_ref_ld2.annomat(is_in_9m);
        alex.w_ld2 = nan(size(chrnumvec)); alex.w_ld2(index_to_9m_nnz) = w_ld2.annomat(is_in_9m);
        alex.mafvec = nan(size(chrnumvec)); alex.mafvec(index_to_9m_nnz) = mafvec_1m.mafvec(is_in_9m);
        alex.bip_zvec = nan(size(chrnumvec)); alex.bip_zvec(index_to_9m_nnz) = bip_1m_withMHC.zvec(is_in_9m);
        alex.bip_nvec = nan(size(chrnumvec)); alex.bip_nvec(index_to_9m_nnz) = bip_1m_withMHC.nvec(is_in_9m);
        alex.scz_zvec = nan(size(chrnumvec)); alex.scz_zvec(index_to_9m_nnz) = scz_1m_withMHC.zvec(is_in_9m);
        alex.scz_nvec = nan(size(chrnumvec)); alex.scz_nvec(index_to_9m_nnz) = scz_1m_withMHC.nvec(is_in_9m);

        corrcoef(alex.mafvec, mafvec_9m, 'rows', 'complete')   % => 0.9999
        corrcoef(alex.ref_ld2, tldvec_9m, 'rows', 'complete')   % => 0.9868
        corrcoef(alex.bip_zvec, bip.zvec, 'rows', 'complete')   % => 0.9791
        corrcoef(alex.scz_zvec, scz.zvec, 'rows', 'complete')   % => 0.9750
        floor([nanmean(alex.bip_nvec) nanmean(bip.Neff)])       %   24683       46582
        floor([nanmean(alex.scz_nvec) nanmean(scz.Neff)])       %   40373       76326
    end
    
    TLD_MAX=600;
    LD_BLOCK_SIZE_MAX = 2000;
    MAF_THRESH = 0.01;
    %mhcvec = chrnumvec==6 & posvec >= 25e6 & posvec <= 35e6;

    % Save defvec
    defvec0 = isfinite(Hvec_9m + tldvec_9m) & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX) & (tldvec_9m <= TLD_MAX) & (mafvec_9m >= MAF_THRESH);
    scz.defvec_1M = defvec0 & mask_in_9m & isfinite(scz.zvec); fprintf('sum(scz.defvec_1M) = %i\n', sum(scz.defvec_1M));
    bip.defvec_1M = defvec0 & mask_in_9m & isfinite(bip.zvec); fprintf('sum(bip.defvec_1M) = %i\n', sum(bip.defvec_1M));
    scz.defvec_5M = defvec0 &              isfinite(scz.zvec); fprintf('sum(scz.defvec_5M) = %i\n', sum(scz.defvec_5M));
    bip.defvec_5M = defvec0 &              isfinite(bip.zvec); fprintf('sum(bip.defvec_5M) = %i\n', sum(bip.defvec_5M));
    
    bmi.defvec_1M    = defvec0 & mask_in_9m & isfinite(bmi.zvec);    fprintf('sum(bmi.defvec_1M) = %i\n',    sum(bmi.defvec_1M));
    height.defvec_1M = defvec0 & mask_in_9m & isfinite(height.zvec); fprintf('sum(height.defvec_1M) = %i\n', sum(height.defvec_1M));    
    scz3.defvec_1M = defvec0 & mask_in_9m & isfinite(scz3.zvec); fprintf('sum(height.defvec_1M) = %i\n', sum(scz3.defvec_1M));    

    if any(scz.defvec_1M ~= defvecForAlex1m_new.defvec), warning('scz.defvec_1M mismatch'); end;
    if any(scz.defvec_5M ~= defvecForAlexScz5m.defvec), warning('scz.defvec_5M mismatch'); end;
    if any(bip.defvec_5M ~= defvecForAlexBIP6m.defvec), warning('bip.defvec_5M mismatch'); end;
    
    % Generate random pruning matrices
    load('H:\NORSTORE\oleksanf\holland_genetics\LD_1kGPhase3_hg19_compatibleWith_info9m_hg18\LD_r2_gt_p8_chrs_1_22_from_1kGPhase3_hg19_1kGHG1pc_compatibleWith_info9m_hg18.mat');
    LDmat = LDr2_p8sparse;   clear LDr2_p8sparse;   % 9279485x9279485
    nprune = 10;
    scz.pruneidxmat_1M = false(size(scz.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~scz.defvec_1M) = NaN; scz.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    bip.pruneidxmat_1M = false(size(bip.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~bip.defvec_1M) = NaN; bip.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');

    scz.pruneidxmat_5M = false(size(scz.defvec_5M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~scz.defvec_5M) = NaN; scz.pruneidxmat_5M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    bip.pruneidxmat_5M = false(size(bip.defvec_5M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~bip.defvec_5M) = NaN; bip.pruneidxmat_5M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');

    bmi.pruneidxmat_1M = false(size(bmi.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~bmi.defvec_1M) = NaN; bmi.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    height.pruneidxmat_1M = false(size(height.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~height.defvec_1M) = NaN; height.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    scz3.pruneidxmat_1M = false(size(scz3.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~scz3.defvec_1M) = NaN; scz3.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');

    %save(data_file, 'scz', 'bip', 'alex',                     'chrnumvec', 'posvec', 'mask_in_9m', 'index_to_9m', 'mafvec_9m', 'Hvec_9m', 'tldvec_9m', 'numSNPsInLDr2_gt_r2min_vec');
    save(workspace_data_file, 'scz', 'bip', 'bmi', 'height', 'scz3',  'alex_sum_r2', 'alex_sum_r4', 'chrnumvec', 'posvec', 'mask_in_9m', 'index_to_9m', 'mafvec_9m', 'Hvec_9m', 'tldvec_9m', 'numSNPsInLDr2_gt_r2min_vec', 'nprune');
    keyboard
    clear
    return
else
    load(workspace_data_file);
    %ref_ld  = struct('sum_r2', tldvec_9m, 'sum_r2_biased', tldvec_9m, 'sum_r4_biased', alex_shape_param .* tldvec_9m);
    ref_ld  = struct('sum_r2', alex_sum_r2, 'sum_r2_biased', alex_sum_r2, 'sum_r4_biased', alex_sum_r4);
end




%task=scz; disp_trait1_name='SCZ'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_1M;clear('task');
%task=scz; disp_trait1_name='SCZ'; zvec=task.zvec; zvec(~task.defvec_5M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_5M;clear('task');
task=bip; disp_trait1_name='BIP'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_1M;clear('task');
%task=bip; disp_trait1_name='BIP'; zvec=task.zvec; zvec(~task.defvec_5M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_5M;clear('task');

%task=bmi; disp_trait1_name='BMI'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.nvec; pruneidxmat=task.pruneidxmat_1M;clear('task');
%task=height; disp_trait1_name='HEIGHT'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.nvec; pruneidxmat=task.pruneidxmat_1M;clear('task');
%task=scz3; disp_trait1_name='SCZsep17'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.nvec; pruneidxmat=task.pruneidxmat_1M;clear('task');

if 1
    options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs in 1kG phase3
    options.verbose = true;
    options.ci_alpha = nan;

    % Standard 2-component mixture (null+causal)
    hits = sum(pruneidxmat, 2); w_ld = 1./hits; w_ld(hits==0) = nan;
    result = BGMG_fit(zvec, Hvec_9m, nvec, w_ld, ref_ld, options);
    result_cdf = [];
    
    if 0
    % Infinitesimal model (only causal component, without null)
    options_inf = options; options_inf.fit_infinitesimal=true;
    result_inf = BGMG_fit(zvec, Hvec_9m, nvec, w_ld, ref_ld, options_inf);
    result_inf_cdf = [];

    % 3-component mixture (null + two causal components)
    options_mix2 = options; options_mix2.fit_two_causal_components=true;
    result_mix2 = BGMG_fit(zvec, Hvec_9m, nvec, w_ld, ref_ld, options_mix2);
    result_mix2_cdf = [];
    end
end
close all
options.title = disp_trait1_name;
options.use_poisson = true;
[result_cdf, figures] = UGMG_qq_plot(result.univariate{1}.params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_cdf);
print(figures.tot, sprintf('%s_poisson_with_r2bins_%i', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_poisson_with_r2bins_%i_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')

if 0
options.title = sprintf('%s (infinitesimal model)', disp_trait1_name);
[result_inf_cdf, figures] = UGMG_qq_plot(result_inf.univariate{1}.params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_inf_cdf);
print(figures.tot, sprintf('%s_%i_inf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_inf_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
 
options.title = sprintf('%s (null+small+large)', disp_trait1_name);
[result_mix2_cdf, figures] = UGMG_qq_plot(result_mix2.univariate{1}.params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_mix2_cdf);
print(figures.tot, sprintf('%s_%i_mix2', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_mix2_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
end
if 0
    params = struct('sig2_zero', 1.239, 'pi_vec', 0.004727, 'sig2_beta', 0.000041); result_dominic_cdf = []; % 1M SCZ
    params = struct('sig2_zero', 1.204, 'pi_vec', 0.004795, 'sig2_beta', 0.000040); result_dominic_cdf = []; % 5M SCZ
    options.title = sprintf('%s', disp_trait1_name);
    [result_dominic_cdf, figures] = UGMG_qq_plot(params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_dominic_cdf);
    print(figures.tot, sprintf('%s_%i_Dominic_params', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
    set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
    print(figures.bin, sprintf('%s_%i_Dominic_params_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
end
