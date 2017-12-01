addpath('DERIVESTsuite');
addpath('H:\NORSTORE\oleksanf\holland_genetics\GPSIM');  % for GenStats_meta
addpath('H:\NORSTORE\oleksanf\holland_matlab\qq');

if 0
    % Align 11M LD matrix to 9M template
    load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat')
    [iLD11, jLD11] = find(LDr2_p8sparse);
    idx = ismember(iLD11, i1kGHG1pc_9m_n_1kGHG1pc) & ismember(jLD11, i1kGHG1pc_9m_n_1kGHG1pc);
    iLD11 = iLD11(idx); jLD11 = jLD11(idx);
    [~, iLD9_idx] = ismember(iLD11, i1kGHG1pc_9m_n_1kGHG1pc);
    [~, jLD9_idx] = ismember(jLD11, i1kGHG1pc_9m_n_1kGHG1pc);
    iLD9 = i9m_9m_n_1kGHG1pc(iLD9_idx);
    jLD9 = i9m_9m_n_1kGHG1pc(jLD9_idx);
    LDmat = sparse(double(iLD9),double(jLD9),true,double(nsnps9m),double(nsnps9m));
    LDmat = LDmat | speye(double(nsnp));
    LDmat = LDmat | (LDmat - LDmat');
    LDr2_p8sparse = LDmat;
    save('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc_converted_to_9M.mat', 'LDr2_p8sparse', '-v7.3');
    
    % validation - looks OK
    hgLDmat = load('H:\NORSTORE\oleksanf\holland_genetics\LD_1kGPhase3_hg19_compatibleWith_info9m_hg18\LD_r2_gt_p8_chrs_1_22_from_1kGPhase3_hg19_1kGHG1pc_compatibleWith_info9m_hg18.mat');
    idx = 8000000:8003000;
    clf;hold on;
    spy(LDr2_p8sparse(idx, idx), 'r');
    spy(hgLDmat.LDr2_p8sparse(idx, idx));
    
    clf;hold on;
    spy(hgLDmat.LDr2_p8sparse(mhcvec, mhcvec));
    spy(LDr2_p8sparse(mhcvec, mhcvec), 'r');
end

workspace_data_file = 'data_2017_12_01_syn.mat';
if 0
    % Test consistency with Dominic's model
    info1m = load('H:\NORSTORE\MMIL\SUMSTAT\LDSR\MATLAB_Annot\infomat.mat');
    info9m = load('H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Annot\infomat.mat', 'chrnumvec', 'posvec', 'snpidlist');
    chrnumvec = info9m.chrnumvec; posvec = info9m.posvec;

    [is_in_9m, index_to_9m] = ismember(info1m.snpidlist, info9m.snpidlist);
    mask_in_9m = false(size(info9m.chrnumvec)); mask_in_9m(index_to_9m(index_to_9m ~= 0)) = true;

    % Load reference data from Dominic, align to 9M template
    load('H:\NORSTORE\oleksanf\holland_genetics\SynGenNov2017\simDataForAlex.mat')
    load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\LDr2hist_zontr2_p5_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat')
    load('H:\NORSTORE\oleksanf\holland_genetics\LDhistAndTLD_1kGPhase3_hg19\info9m_n_1kGHG1pc_RemappingIndices.mat');
    load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc_converted_to_9M.mat');

    total_het = nansum(2 .* mafvec .* (1-mafvec));
    causals_het = nansum(2 .* mafvec(betaTrueVec~=0) .* (1-mafvec(betaTrueVec~=0)));
    
    nsnps9m = length(info9m.chrnumvec);
    mafvec_9m   = nan(nsnps9m,1);
    tldvec_9m   = nan(nsnps9m,1);
    LDr2hist_9m = zeros(nsnps9m,size(LDr2hist,2)); %   % FIX UP AND UNCOMMENT!
    LDr2hist_9m(:,end) = 1;                        % Initialize to ones: "missing" SNPs will have r^2 with self equal to 1.  % FIX UP AND UNCOMMENT!

    tldvec_9m(i9m_9m_n_1kGHG1pc)     = tldvec(i1kGHG1pc_9m_n_1kGHG1pc);
    mafvec_9m(i9m_9m_n_1kGHG1pc)     = mafvec(i1kGHG1pc_9m_n_1kGHG1pc);
    LDr2hist_9m(i9m_9m_n_1kGHG1pc,:) = LDr2hist(i1kGHG1pc_9m_n_1kGHG1pc,:);  % FIX UP AND UNCOMMENT!
    Hvec_9m = 2*mafvec_9m.*(1-mafvec_9m);
 
    ldr2binsNum     = size(LDr2hist_9m,2);  %  100
    ldr2binsEdges   = linspace( 0,1,ldr2binsNum+1);
    ldr2binsCenters = nan(1,ldr2binsNum); for i=1:ldr2binsNum, ldr2binsCenters(i) = (ldr2binsEdges(i+1)+ldr2binsEdges(i))/2; end

    %minr2 = 0.05; % If r^2<0.05, ignore --> noise. USUALLY HAVE USED THIS ONE.
    minr2 = 0.1;  % If r^2<0.05, ignore --> noise.
    minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    numSNPsInLDr2_gt_r2min_vec = sum(LDr2hist_9m(:,minr2bin:end),2);
    
    minr2 = 0.05; minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    LDr2tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters, [size(LDr2hist_9m, 1), 1]);
    LDr4tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters.^2, [size(LDr2hist_9m, 1), 1]);
    
    % calculate chi2 for two r2 bins ("low" and "high" r2)
    r2_bin_index = false(2, ldr2binsNum);
    r2_bin_index(1, minr2bin:(ldr2binsNum/2)) = true;
    r2_bin_index(2, ((ldr2binsNum/2) + 1):end) = true;
    
    alex_sum_r2=zeros(nsnps9m, 2);
    alex_sum_r4=zeros(nsnps9m, 2);
    for r2_bini=1:2
        alex_sum_r2(:, r2_bini) = sum(LDr2tot_9m(:, r2_bin_index(r2_bini, :)), 2);
        alex_sum_r4(:, r2_bini) = sum(LDr4tot_9m(:, r2_bin_index(r2_bini, :)), 2);
    end
    alex_shape_param = alex_sum_r4 ./ alex_sum_r2;

    TLD_MAX=600;
    LD_BLOCK_SIZE_MAX = 2000;
    MAF_THRESH = 0.01;
    mhcvec = chrnumvec==6 & posvec >= 25e6 & posvec <= 35e6;

    task.zvec   = nan(nsnps9m,1);
    task.zvec(i9m_9m_n_1kGHG1pc) = zvec(i1kGHG1pc_9m_n_1kGHG1pc);
    task.nvec   = 100000 * ones(size(task.zvec));
    
    % Save defvec
    defvec0 = isfinite(Hvec_9m + tldvec_9m) & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX) & (tldvec_9m <= TLD_MAX) & (mafvec_9m >= MAF_THRESH);
    task.defvec = defvec0 & mask_in_9m & ~mhcvec & isfinite(task.zvec); fprintf('sum(defvec_1M) = %i\n', sum(task.defvec));

    nprune = 10;
    task.pruneidxmat = false(size(task.defvec,1), nprune);
    LDmat = LDr2_p8sparse;
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~task.defvec) = NaN; task.pruneidxmat(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    sum(task.pruneidxmat)
    
    save(workspace_data_file, 'task', 'alex_shape_param', 'alex_sum_r2', 'alex_sum_r4', 'gwasParams',  'total_het', 'causals_het', 'chrnumvec', 'posvec', 'mask_in_9m', 'index_to_9m', 'mafvec_9m', 'Hvec_9m', 'tldvec_9m', 'numSNPsInLDr2_gt_r2min_vec', 'nprune');
    keyboard
    clear
    return
else
    load(workspace_data_file);
    ref_ld  = struct('sum_r2', alex_sum_r2, 'chi_r4', alex_shape_param);
end


disp_trait1_name='synthetic'; zvec=task.zvec; zvec(~task.defvec) = nan; nvec=task.nvec; pruneidxmat=task.pruneidxmat;

options.total_het = total_het;  % Total heterozigosity across all SNPs in 1kG phase3
options.verbose = true;
options.ci_alpha = nan;

if 0
    % Standard 2-component mixture (null+causal)
    hits = sum(pruneidxmat, 2); w_ld = 1./hits; w_ld(hits==0) = nan;
    result = BGMG_fit(zvec, Hvec_9m, nvec, w_ld, ref_ld, options);
    result_cdf = [];
    true_sig2_zero = result.univariate{1}.params.sig2_zero;
end

close all
if 0
options.title = disp_trait1_name;
[result_cdf, figures] = UGMG_qq_plot(result.univariate{1}.params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_cdf);
print(figures.tot, sprintf('%s_%i', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
end

params = struct('sig2_zero', 1, 'pi_vec', gwasParams.pi1True, 'sig2_beta', gwasParams.sig2betaEst); result_true_params_cdf = [];
options.title = sprintf('%s', disp_trait1_name);
[result_true_params_cdf, figures] = UGMG_qq_plot(params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_true_params_cdf);
print(figures.tot, sprintf('%s_%i_poisson_with_r2bins_true_params', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_poisson_with_r2bins_true_params_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')

if 0
true_sig2_zero = 1.088;
params = struct('sig2_zero', true_sig2_zero, 'pi_vec', gwasParams.pi1True, 'sig2_beta', gwasParams.sig2betaEst); result_true_model_sig2zero_params_cdf = [];
options.title = sprintf('%s', disp_trait1_name);
[result_true_model_sig2zero_params_cdf, figures] = UGMG_qq_plot(params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_true_model_sig2zero_params_cdf);
print(figures.tot, sprintf('%s_%i_true_params_model_sig2zero', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_true_params_model_sig2zero_HL', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
end