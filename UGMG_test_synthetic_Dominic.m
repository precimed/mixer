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

USE_DATA_DEC_17 = 1;

%workspace_data_file = 'data_2017_12_05a_syn.mat';
%workspace_data_file = 'data_2017_12_15_synA.mat';
workspace_data_file = 'data_2018_01_22.mat';
if 0
    % Test consistency with Dominic's model
    info1m = load('H:\NORSTORE\MMIL\SUMSTAT\LDSR\MATLAB_Annot\infomat.mat');
    info9m = load('H:\NORSTORE\MMIL\new_9m_SNPs_processed_summary\GWAS_Annot\infomat.mat', 'chrnumvec', 'posvec', 'snpidlist');
    chrnumvec = info9m.chrnumvec; posvec = info9m.posvec;
    dataDec17 = load('H:\NORSTORE\oleksanf\holland_genetics\SynGen_14-Dec-2017\sim_pi1_1em3_h2_p4_pruned_nrep_1_minr2_p05.mat');
    [is_in_9m, index_to_9m] = ismember(info1m.snpidlist, info9m.snpidlist);
    mask_in_9m = false(size(info9m.chrnumvec)); mask_in_9m(index_to_9m(index_to_9m ~= 0)) = true;

    % Load reference data from Dominic, align to 9M template
    simDataForAlex=load('H:\NORSTORE\oleksanf\holland_genetics\SynGenNov2017\simDataForAlex.mat');
    if 0, load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\LDr2hist_zontr2_p5_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat'); end;
    LDr2=load('H:\NORSTORE\oleksanf\11015833\LDr2_biased_10Ksubj_1cm_100bins.mat');
    LDr4=load('H:\NORSTORE\oleksanf\11015833\LDr4_biased_10Ksubj_1cm_100bins.mat');
    
    load('H:\NORSTORE\oleksanf\holland_genetics\LDhistAndTLD_1kGPhase3_hg19\info9m_n_1kGHG1pc_RemappingIndices.mat');
    load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\LD_r2_gt_p8_chrs_1_22_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc_converted_to_9M.mat');
    mafvec_plink = load('H:\NORSTORE\oleksanf\11015833\mafvec.mat');
    if ~USE_DATA_DEC_17 && ~exist('p_HG80p3m', 'var'), p_HG80p3m=load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2\p_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat'); end;
    parameterResults=load('H:\NORSTORE\oleksanf\holland_GWAS_data_SynGen2/setUpParameterResults_HG80p3m_HG80p3m_n_1kGPIII14p9m_1pc.mat');

    %total_het = nansum(2 .* simDataForAlex.mafvec .* (1-simDataForAlex.mafvec));
    total_het = nansum(2 .* mafvec_plink.mafvec .* (1-mafvec_plink.mafvec));
    
    nsnps9m = length(info9m.chrnumvec);
    mafvec_9m   = nan(nsnps9m,1);
    tldvec_9m   = nan(nsnps9m,1);
    if 0
    LDr2hist_9m = zeros(nsnps9m,size(LDr2hist,2)); %   % FIX UP AND UNCOMMENT!
    LDr2hist_9m(:,end) = 1;                        % Initialize to ones: "missing" SNPs will have r^2 with self equal to 1.  % FIX UP AND UNCOMMENT!
    end
    LDr2_9m = zeros(nsnps9m,size(LDr2.LDr2,2));
    LDr4_9m = zeros(nsnps9m,size(LDr4.LDr4, 2));

    if 0
    tldvec_9m(i9m_9m_n_1kGHG1pc)     = simDataForAlex.tldvec(i1kGHG1pc_9m_n_1kGHG1pc);
    mafvec_9m(i9m_9m_n_1kGHG1pc)     = simDataForAlex.mafvec(i1kGHG1pc_9m_n_1kGHG1pc);
    LDr2hist_9m(i9m_9m_n_1kGHG1pc,:) = LDr2hist(i1kGHG1pc_9m_n_1kGHG1pc,:);  % FIX UP AND UNCOMMENT!
    end
    mafvec_9m(i9m_9m_n_1kGHG1pc)     = mafvec_plink.mafvec(i1kGHG1pc_9m_n_1kGHG1pc);
    LDr2_9m(i9m_9m_n_1kGHG1pc,:) = LDr2.LDr2(i1kGHG1pc_9m_n_1kGHG1pc, :);
    LDr4_9m(i9m_9m_n_1kGHG1pc,:) = LDr4.LDr4(i1kGHG1pc_9m_n_1kGHG1pc, :);
    tldvec_9m  = sum(LDr2_9m, 2);   
    Hvec_9m = 2*mafvec_9m.*(1-mafvec_9m);
 
    ldr2binsNum     = size(LDr2_9m,2);  %  100
    ldr2binsEdges   = linspace( 0,1,ldr2binsNum+1);
    ldr2binsCenters = nan(1,ldr2binsNum); for i=1:ldr2binsNum, ldr2binsCenters(i) = (ldr2binsEdges(i+1)+ldr2binsEdges(i))/2; end

    LDr2hist_9m = (LDr2_9m .* LDr2_9m ./ LDr4_9m);
    LDr2hist_9m(~isfinite(LDr2hist_9m)) = 0;
    
    %minr2 = 0.05; % If r^2<0.05, ignore --> noise. USUALLY HAVE USED THIS ONE.
    minr2 = 0.1;  % If r^2<0.05, ignore --> noise.
    minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    numSNPsInLDr2_gt_r2min_vec = sum(LDr2hist_9m(:,minr2bin:end),2);
    
    if 0
    minr2 = 0.05; minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    LDr2tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters, [size(LDr2hist_9m, 1), 1]);
    LDr4tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters.^2, [size(LDr2hist_9m, 1), 1]);
    end
    
    % calculate chi2 for two r2 bins ("low" and "high" r2)
    r2_bin_index = false(2, ldr2binsNum);
    r2_bin_index(1, minr2bin:(ldr2binsNum/2)) = true;
    r2_bin_index(2, ((ldr2binsNum/2) + 1):end) = true;
    
    alex_sum_r2=zeros(nsnps9m, 2);
    alex_sum_r4=zeros(nsnps9m, 2);
    for r2_bini=1:2
        %alex_sum_r2(:, r2_bini) = sum(LDr2tot_9m(:, r2_bin_index(r2_bini, :)), 2);
        %alex_sum_r4(:, r2_bini) = sum(LDr4tot_9m(:, r2_bin_index(r2_bini, :)), 2);
         alex_sum_r2(:, r2_bini) = sum(LDr2_9m(:, r2_bin_index(r2_bini, :)), 2);
        alex_sum_r4(:, r2_bini) = sum(LDr4_9m(:, r2_bin_index(r2_bini, :)), 2);
    end
    alex_shape_param = alex_sum_r4 ./ alex_sum_r2;

    TLD_MAX=600;
    LD_BLOCK_SIZE_MAX = 2000;
    MAF_THRESH = 0.01;
    mhcvec = chrnumvec==6 & posvec >= 25e6 & posvec <= 35e6;

    z = @(logpvec, zvec) -norminv(10.^-logpvec/2).*sign(zvec);

    %tasks = cell(10,4,3);
    tasks=cell(2, 1);
    for beta_idx=1:1 %10
    for pi_idx=1:1 %4
    for h2_idx=1:1 %3
    for use_1m=1:2
        task.zvec   = nan(nsnps9m,1);
        if ~USE_DATA_DEC_17, 
            zvec = z(-log10(p_HG80p3m.pvecs{beta_idx,pi_idx,h2_idx}), randn(size(p_HG80p3m.pvecs{beta_idx,pi_idx,h2_idx})));
            task.gwasParams = parameterResults.parameterResults{beta_idx, pi_idx,h2_idx};
        else
            zvec = dataDec17.zvec; 
            task.gwasParams = dataDec17.gwasParams;
        end

        task.zvec(i9m_9m_n_1kGHG1pc) = zvec(i1kGHG1pc_9m_n_1kGHG1pc);
        task.nvec   = 100000 * ones(size(task.zvec));

        % Save defvec
        defvec0 = isfinite(Hvec_9m + tldvec_9m) & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX) & (tldvec_9m <= TLD_MAX) & (mafvec_9m >= MAF_THRESH);
        task.defvec = defvec0 & ~mhcvec & isfinite(task.zvec); fprintf('sum(defvec_1M) = %i\n', sum(task.defvec));
        if use_1m==1
            task.defvec = task.defvec & mask_in_9m ;
        end

        nprune = 10;
        task.pruneidxmat = false(size(task.defvec,1), nprune);
        LDmat = LDr2_p8sparse;
        for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~task.defvec) = NaN; task.pruneidxmat(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end;
        sum(task.pruneidxmat)
                
        %tasks{beta_idx, pi_idx, h2_idx} = task;
        tasks{use_1m} = task;
        fprintf('|');
    end
    end
    end
        fprintf('\n');
    end
    
    save(workspace_data_file, '-v7.3', 'tasks', 'alex_shape_param', 'alex_sum_r2', 'alex_sum_r4', 'total_het',  'chrnumvec', 'posvec', 'mask_in_9m', 'index_to_9m', 'mafvec_9m', 'Hvec_9m', 'tldvec_9m', 'numSNPsInLDr2_gt_r2min_vec', 'nprune');
    keyboard
    clear
    return
else
    load(workspace_data_file);
    %ref_ld  = struct('sum_r2', alex_sum_r2, 'chi_r4', alex_shape_param);
    alex_sum_r4 = alex_sum_r2 .* alex_shape_param;
    ref_ld  = struct('sum_r2', alex_sum_r2, 'sum_r4_biased', alex_sum_r4, 'sum_r2_biased', alex_sum_r2);
end

for beta_idx=1 % 2:2
for pi_idx=1 % 3:34
for h2_idx=1 % 1:3

%task = tasks{beta_idx, pi_idx, h2_idx};
task = tasks{2}
disp_trait1_name=sprintf('syntheticDec17-%i-%.2e-%.2f', beta_idx, task.gwasParams.pi1True, task.gwasParams.h2True); zvec=task.zvec; zvec(~task.defvec) = nan; nvec=task.nvec; pruneidxmat=task.pruneidxmat;

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
print(figures.tot, sprintf('%s_%i.pdf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_HL.pdf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
end

params = struct('sig2_zero', 1, 'pi_vec', task.gwasParams.pi1True, 'sig2_beta', task.gwasParams.sig2betaEst); result_true_params_cdf = [];
options.title = sprintf('%s', disp_trait1_name);
options.use_poisson = 1;
[result_true_params_cdf, figures] = UGMG_qq_plot(params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_true_params_cdf);
print(figures.tot, sprintf('%s_%i_2018_01_22_poisson_without_r2bins_kmax49_true_params_with_Dominic_model.pdf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_2018_01_22_poisson_without_r2bins_kmax49_true_params_HL.pdf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')

if 0
true_sig2_zero = 1.088;
params = struct('sig2_zero', true_sig2_zero, 'pi_vec', task.gwasParams.pi1True, 'sig2_beta', task.gwasParams.sig2betaEst); result_true_model_sig2zero_params_cdf = [];
options.title = sprintf('%s', disp_trait1_name);
[result_true_model_sig2zero_params_cdf, figures] = UGMG_qq_plot(params, zvec, Hvec_9m, nvec, pruneidxmat, ref_ld, options, result_true_model_sig2zero_params_cdf);
print(figures.tot, sprintf('%s_%i_true_params_model_sig2zero.pdf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
set(figures.bin,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(figures.bin, sprintf('%s_%i_true_params_model_sig2zero_HL.pdf', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')
end

if 0
    %dominic_qq=load('sim_pi1_1em3_h2_p4_pruned_qq_values.mat')
    dominic_qq=load('2017_12_22_sim_pi1_1em3_h2_p4_pruned_qq_values (1).mat')
    figure(figures.tot)
    hold on
    plot(dominic_qq.qq_x_data, dominic_qq.qq_y_data, '-b', dominic_qq.qq_x_model, dominic_qq.qq_y_model, '-y')
    lgd=legend('Data - Alex', 'Model - Alex', 'Expected', 'Data - Dominic', 'Model - Dominic', 'Location', 'SouthEast'); lgd.FontSize = 19; 
end

if 0
    sig2betaEst=[]
    for i=1:10
        sig2betaEst=[sig2betaEst setUpParameterResults_HG80p3m.parameterResults{i, 4, 3}.sig2betaEst];
    end
    sig2betaEst
end

end
end
end