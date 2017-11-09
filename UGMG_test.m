addpath('DERIVESTsuite');
addpath('H:\NORSTORE\oleksanf\holland_genetics\GPSIM');  % for GenStats_meta
addpath('H:\NORSTORE\oleksanf\holland_matlab\qq');

nprune = 10;
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
    
    minr2 = 0.05; minr2bin = find(ldr2binsEdges < minr2,1,'last')+1;
    LDr2tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters, [size(LDr2hist_9m, 1), 1]);
    LDr4tot_9m = LDr2hist_9m .* repmat(ldr2binsCenters.^2, [size(LDr2hist_9m, 1), 1]);
    sum_r2 = sum(LDr2tot_9m(:, minr2bin:end), 2);
    sum_r4 = sum(LDr4tot_9m(:, minr2bin:end), 2);
    alex_shape_param = sum_r4 ./ sum_r2;

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
    scz.defvec_1M = mask_in_9m & isfinite(scz.zvec + Hvec_9m + tldvec_9m); fprintf('sum(scz.defvec_1M) = %i\n', sum(scz.defvec_1M));
    bip.defvec_1M = mask_in_9m & isfinite(bip.zvec + Hvec_9m + tldvec_9m); fprintf('sum(bip.defvec_1M) = %i\n', sum(bip.defvec_1M));
    scz.defvec_5M =               isfinite(scz.zvec + Hvec_9m + tldvec_9m); fprintf('sum(scz.defvec_5M) = %i\n', sum(scz.defvec_5M));
    bip.defvec_5M =               isfinite(bip.zvec + Hvec_9m + tldvec_9m); fprintf('sum(bip.defvec_5M) = %i\n', sum(bip.defvec_5M));
    
    if 0
        % Further filtering...
        scz.defvec_1M = scz.defvec_1M & (mafvec_9m >= MAF_THRESH) & (tldvec_9m <= TLD_MAX) & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX); % & mhcvec
        bip.defvec_1M = bip.defvec_1M & (mafvec_9m >= MAF_THRESH) & (tldvec_9m <= TLD_MAX) & (numSNPsInLDr2_gt_r2min_vec <= LD_BLOCK_SIZE_MAX); % & mhcvec
        fprintf('sum(scz.defvec_1M & further_filters) = %i\n', sum(scz.defvec_1M));
        fprintf('sum(bip.defvec_1M & further_filters) = %i\n', sum(bip.defvec_1M));
    end

    % Generate random pruning matrices
    load('H:\NORSTORE\oleksanf\holland_genetics\LD_1kGPhase3_hg19_compatibleWith_info9m_hg18\LD_r2_gt_p8_chrs_1_22_from_1kGPhase3_hg19_1kGHG1pc_compatibleWith_info9m_hg18.mat');
    LDmat = LDr2_p8sparse;   clear LDr2_p8sparse;   % 9279485x9279485
    scz.pruneidxmat_1M = false(size(scz.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~scz.defvec_1M) = NaN; scz.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    bip.pruneidxmat_1M = false(size(bip.defvec_1M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~bip.defvec_1M) = NaN; bip.pruneidxmat_1M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');

    scz.pruneidxmat_5M = false(size(scz.defvec_5M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~scz.defvec_5M) = NaN; scz.pruneidxmat_5M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');
    bip.pruneidxmat_5M = false(size(bip.defvec_5M,1), nprune);
    for repi=1:nprune, tmp=rand(size(chrnumvec,1),1); tmp(~bip.defvec_5M) = NaN; bip.pruneidxmat_5M(:,repi) = isfinite(FastPrune(tmp, LDmat)); fprintf('.'); end; fprintf('\n');

    %save('data_2017_11_08.mat', 'scz', 'bip', 'alex',                     'chrnumvec', 'posvec', 'mask_in_9m', 'index_to_9m', 'mafvec_9m', 'Hvec_9m', 'tldvec_9m', 'numSNPsInLDr2_gt_r2min_vec');
    save('data_2017_11_08.mat', 'scz', 'bip',          'alex_shape_param', 'chrnumvec', 'posvec', 'mask_in_9m', 'index_to_9m', 'mafvec_9m', 'Hvec_9m', 'tldvec_9m', 'numSNPsInLDr2_gt_r2min_vec');
    keyboard
    clear
    return
else
    load('data_2017_11_08.mat');
    ref_ld  = struct('sum_r2', tldvec_9m, 'chi_r4', alex_shape_param);
    %ref_ld  = struct('sum_r2', alex.ref_ld2, 'chi_r4', alex.biased_ref_ld4 ./ alex.biased_ref_ld2);
end


%task=scz; disp_trait1_name='SCZ'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_1M;clear('task');
%task=scz; disp_trait1_name='SCZ'; zvec=task.zvec; zvec(~task.defvec_5M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_5M;clear('task');
%task=bip; disp_trait1_name='BIP'; zvec=task.zvec; zvec(~task.defvec_1M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_1M;clear('task');
task=bip; disp_trait1_name='BIP'; zvec=task.zvec; zvec(~task.defvec_5M) = nan; nvec=task.Neff; pruneidxmat=task.pruneidxmat_5M;clear('task');

hits = sum(pruneidxmat, 2); w_ld = 1./hits; w_ld(hits==0) = nan;

if 1
    options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs
    options.verbose = true;
    options.ci_alpha = nan;
    
    result = BGMG_fit(zvec, Hvec_9m, nvec, w_ld, ref_ld, options);
    
    BGMG_util.result2str(1, result)
    result.univariate{1}.params
end

if 1
    options.calculate_z_cdf = true;
    options.calculate_z_cdf_limit = ceil(min(max(abs(zvec)), 15));
    options.calculate_z_cdf_step = 0.005;
    options.calculate_z_cdf_weights = hits ./ nansum(hits);
    [~, result_cdf] = BGMG_univariate_cost(result.univariate{1}.params, zvec, Hvec_9m, nvec, w_ld, ref_ld, options);
end

% data QQ plot
dz = 0.4; 
zextreme = 38.0;
zextreme = min(max(abs(zvec)),zextreme);
nzbins = 2*ceil( zextreme/dz ) + 1;   % Guaranteed odd.
zvals_disc = linspace(-zextreme,zextreme,nzbins); % zvals_disc(2)-zvals_disc(1) == dz;       Values are bin centers. Middle bin value is 0 (center of center bin).
nzvals = length(zvals_disc);                      % nzbins
hv_z = linspace(0,zvals_disc(end),10000);

data_logpvec = zeros(size(hv_z));
for repi = 1:nprune
    data_logpvecI = -log10(normcdf(-abs(zvec(pruneidxmat(:, repi))))*2);
    [data_logqvecI, hv_logp] = GenStats_QQ_plot_amd(data_logpvecI,hv_z);
    data_logpvec = data_logpvec + data_logqvecI;
    lamGC_empirical(repi) = nanmedian(zvec(pruneidxmat(:, repi)).^2)/chi2inv(   0.5,1);
end
data_logpvec = data_logpvec / nprune;

z_grid = result_cdf.cdf_z_grid;
model_logpvec = -log10(2*interp1(-z_grid(z_grid<=0), result_cdf.cdf(z_grid<=0), hv_z'));

% GenStats_QQ plot
clf;hold on
hData     = plot(data_logpvec, hv_logp,  '-','LineWidth',1); hold on;
hModel    = plot(model_logpvec,hv_logp, '-','LineWidth',1); hold on;

%model_logpvec0 = -log10(2*model_cdf(z_grid <= 0));
%model_hv_logp0 = -log10(2*normcdf(-abs(z_grid(z_grid <= 0))));
%hModel0    = plot(model_logpvec0,model_hv_logp0, '-','LineWidth',1); hold on;

lamGC_data = lamGCfromQQ(data_logpvec, hv_logp);
lamGC_model = lamGCfromQQ(model_logpvec, hv_logp);

if 0
    % Yet another data qq plot (directly from data points)
    zvec_qq = zvec; weights = hits;
    %weights=ones(size(hits));
    defvec=isfinite(zvec+weights);
    zvec_qq=zvec_qq(defvec); weights=weights(defvec); weights=weights/sum(weights);

    [data_logpvec2, si] = sort(-log10(2*normcdf(-abs(zvec_qq))));
    hv_logpvec2=-log10(cumsum(weights(si),1,'reverse'));
    plot(hv_logpvec2, data_logpvec2,'g'); 
end

qqlimy = 20; qqlimx = 7;

fontsize=19;
plot([0 qqlimy],[0 qqlimy], 'k--');
xlim([0 qqlimx]); ylim([0 qqlimy]);
lgd=legend('Data', 'Model', 'Expected', 'Location', 'SouthEast');
lgd.FontSize = fontsize;
xlabel('Empirical -log 10(q)','fontsize',fontsize)
ylabel('Nominal -log 10(p)','fontsize',fontsize)
title(disp_trait1_name,'fontsize',fontsize);
xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);
yt = get(gca, 'YTick');set(gca, 'FontSize', fontsize);
set(gca, 'box','off');
params=result.univariate{1}.params;
if isfield(result.univariate{1}, 'ci'), ci=result.univariate{1}.ci; end;
text(0.5,18,sprintf('$$ n_{snps} = %i $$', sum(isfinite(zvec))),'FontSize',fontsize,'Interpreter','latex');
text(0.5,16,sprintf('$$ \\hat\\sigma_0^2 = %.3f $$', params.sig2_zero),'FontSize',fontsize,'Interpreter','latex');
text(0.5,14,sprintf('$$ \\hat\\pi^u_1 = %.6f $$', params.pi_vec),'FontSize',fontsize,'Interpreter','latex');
text(0.5,12,sprintf('$$ \\hat\\sigma_{\\beta}^2 = %.6f $$', params.sig2_beta),'FontSize',fontsize,'Interpreter','latex');
text(0.5,10,sprintf('$$ \\hat\\lambda_{model} = %.3f $$', lamGC_model),'FontSize',fontsize,'Interpreter','latex');
text(0.5,8,sprintf('$$ \\hat\\lambda_{data} = %.3f $$', lamGC_data),'FontSize',fontsize,'Interpreter','latex');
text(0.5,6,sprintf('$$ \\hat h^2 = %.3f $$', params.pi_vec*params.sig2_beta*options.total_het),'FontSize',fontsize,'Interpreter','latex');
if exist('ci', 'var'), text(0.5,10,sprintf('$$ \\hat  h^2 = %.6f $$', ci.h2.point_estimate),'FontSize',fontsize,'Interpreter','latex'); end;

print(gcf, sprintf('%s_%i', disp_trait1_name, sum(isfinite(zvec))),'-dpdf')

% TBD: QQ plots with 5M SNPs
% TBD: run fully on Dominic's data, incl. TLD and shape param