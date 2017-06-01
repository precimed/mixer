if ~exist('LDr2_p1sparse', 'var')
    t = 'H:\NORSTORE\MMIL';                       if exist(t, 'dir'), mmil_gwas_path = t; end;
    t = '/work/users/oleksanf/NORSTORE/MMIL';     if exist(t, 'dir'), mmil_gwas_path = t; end;
    t = '/space/syn03/1/data/GWAS';               if exist(t, 'dir'), mmil_gwas_path = t; end;
    if ~exist('mmil_gwas_path', 'var'), error('Unable to locate data folder'); end;

    load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/infomat.mat')); 
    load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/annomat.mat'))
    load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p1.mat'));LDmat = LDr2_p1sparse;
    load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p8.mat'));
    %load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p8.mat'))
end

if 0
files = {};
files{1} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp2295c41c_8050_4408_8fdb_d9867cf525ab_pi=3e-02_gencorr=0.50_h2=0.90.mat';
files{2} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp895a3b23_5816_41fa_a170_0812b7ea602b_pi=1e-02_gencorr=0.50_h2=0.90.mat';
files{3} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp93fdb2c0_196c_4ef4_9410_4ee79fdc0677_pi=3e-03_gencorr=0.50_h2=0.90.mat';
files{4} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp20797da3_d059_47fa_bb0e_d045c0c9f32c_pi=1e-03_gencorr=0.50_h2=0.90.mat';
files{5} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp919ccf53_8917_498f_a7f8_9330add81adc_pi=3e-04_gencorr=0.50_h2=0.90.mat';
files{6} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp91181fa4_8365_4cd0_b1a5_2773f58470e4_pi=1e-04_gencorr=0.50_h2=0.90.mat';
files{7} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp7bbfb030_ac34_41df_b990_dbc61cc98129_pi=3e-05_gencorr=0.50_h2=0.90.mat';
files{8} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp8d706596_04ae_4407_898d_4432e3a74170_pi=1e-05_gencorr=0.50_h2=0.90.mat';
files{9} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp7fdc765b_8a1c_40ae_b259_b78eed64354e_pi=3e-06_gencorr=0.50_h2=0.90.mat';
files{10} = 'H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp_tp77b11cdb_3ea1_4c80_b32a_5c4019c37516_pi=1e-06_gencorr=0.50_h2=0.90.mat';
for i=1:length(files), frames{i}=load(files{i}); end;
end

resultStructs = {};

%files = dir('H:\NORSTORE\oleksanf\SynGP_GenCorr\2m\gwaslogp*.mat');

results = {};
for frame_index=1:length(frames)
    frame = frames{frame_index}.frame;
    fprintf('Processing frame %i, %s\n', frame_index, files{frame_index});
    if length(frame.gwasbeta) ~= 2558411, error('wrong template'); end;
    zvec = frame.gwaszscore; zmat=zvec;
    Hvec = 2 * mafvec .* (1-mafvec);
    Nvec = frame.gwassize; Nmat=Nvec;
    w_ld = annomat(:, end);
    
    fnUVT1 = @(x)GMM_ofrei_univariate_observed_cost(x, zmat(:, 1), Hvec, Nmat(:, 1), w_ld, @GMM_ofrei_univariate_mapparams);
    [cost, pdf, fdr] = fnUVT1(GMM_ofrei_univariate_mapparams(struct('sigma_beta', 0.123, 'sigma0', 1.05, 'pivec', 0.01)));

    fnUVT2 = @(x)GMM_ofrei_univariate_observed_cost(x, zmat(:, 2), Hvec, Nmat(:, 2), w_ld, @GMM_ofrei_univariate_mapparams);
    [cost, pdf, fdr] = fnUVT2(GMM_ofrei_univariate_mapparams(struct('sigma_beta', 0.123, 'sigma0', 1.05, 'pivec', 0.01)));

    fnBVT = @(x)GMM_ofrei_bivariate_observed_cost(x, zmat, Hvec, Nmat, w_ld, @GMM_ofrei_bivariate_mapparams);
    [cost, pdf, condfdr1, condfdr2, conjfdr] = fn2(GMM_ofrei_bivariate_mapparams(struct('sigma_beta', [0.123 0.234], 'rho_beta', 0.2, 'sigma0', [1.05 1.10], 'rho0', 0.02, 'pivec', [0.01 0.01 0.05])));
    
    
    break
end