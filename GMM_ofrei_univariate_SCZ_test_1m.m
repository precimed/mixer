outdir = 'figures'; if ~exist(outdir,'dir'), mkdir(outdir); end
if ~exist('baseline_ld', 'var')
    baseline_ld =  load('H:\NORSTORE\oleksanf\1m\1000G_Phase3_baseline_ldscores.mat');
    eur_w_ld = load('H:\NORSTORE\oleksanf\1m\eur_w_ld.mat');
    load('H:\NORSTORE\oleksanf\1m\ldmat_1m_p8.mat', 'mafvec'); % mafvec
    load('H:\NORSTORE\oleksanf\1m\infomat.mat');
end

if ~exist('files', 'var')
files = {};
norstore = 'H:\NORSTORE\oleksanf\SYNGP\EUR_100K_1M_frames';
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp28ba7c95_b93a_47f4_9059_e894f0a5ffe4_pi=1e-01_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp051a0937_314a_4cff_8add_c95988de8e26_pi=3e-02_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp1ea484d9_6d47_4f10_a478_1062f2e33b2b_pi=1e-02_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp111e0f23_6200_43a1_a975_6bb30854c4b1_pi=3e-03_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp052c47a6_fc4a_4bed_bc01_5e496129321f_pi=1e-03_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp1005c16b_11b3_4c40_b757_0348e6c663c3_pi=3e-04_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp29aaefc1_aa63_4d9b_a099_1e27c66fb840_pi=1e-04_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp1e52676d_6823_4eb2_b6f1_5a0babc9d035_pi=3e-05_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp2fa7af92_67ac_4d52_84fa_8079aa112f06_pi=1e-05_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp1b533cab_8e56_4968_a587_cca04264f68b_pi=3e-06_h2=0.90_AllSNPs.mat');
files{end+1, 1} = fullfile(norstore, 'gwaslogp_tp0f6d6409_06c1_4d28_98e3_2f0d5a0ce7ac_pi=1e-06_h2=0.90_AllSNPs.mat');
end

if ~exist('results', 'var'), results = {}; end;

for file_index = (1+length(results)):length(files)
fprintf('processing file %s...\n',files{file_index});
load(files{file_index});
zvec = frame.gwaszscore;
Hvec = 2 * mafvec .* (1-mafvec);
Nvec = frame.gwassize;
ref_ld = baseline_ld.annomat(:, 1);
w_ld = eur_w_ld.annomat;
mapparams = @GMM_ofrei_univariate_mapparams;
fn = @(x)GMM_ofrei_univariate_cost(x, zvec, Hvec, Nvec, ref_ld, w_ld, mapparams, struct('x_limit', 20, 'x_step', 0.1, 'convolution', true));

if 0
lb = mapparams(struct('sigma_beta', 1e-7, 'sigma0', 1, 'pivec', 7e-6));
ub = mapparams(struct('sigma_beta', 1e-0, 'sigma0', 1.5, 'pivec', 0.3));
xPS = particleswarm(fn, 3, lb, ub, optimoptions(@particleswarm, 'Display', 'iter', 'UseParallel', true)); sPS = mapparams(xPS)
fn(xPS)
end

% Check all points in a wide range of values across log space
sbt_vec = logspace(-3, 0, 25);
pi_vec = logspace(-7, 0, 29);
sbt_mat = repmat(sbt_vec, [length(pi_vec), 1]); sbt_mat=sbt_mat(:);
pi_mat = repmat(pi_vec, [1, length(sbt_vec)]); pi_mat=pi_mat(:);

cost = zeros(size(sbt_mat)); sigma0 = 1;
for i=1:length(sbt_mat), cost(i) = fn(mapparams(struct('sigma_beta', sbt_mat(i), 'sigma0', sigma0, 'pivec', pi_mat(i)))); end;
cost = reshape(cost, [length(pi_vec), length(sbt_vec)]);
%cost-min(cost(:)); image(cost-min(cost(:)));
imagesc(log(1 + cost(1:end-2, :)-min(cost(:)))); colorbar;
[i, j] = ind2sub(size(cost), find(cost(:) == min(cost(:))));
x0 = mapparams(struct('sigma_beta', sbt_vec(j), 'sigma0', 1.05, 'pivec', pi_vec(i))); mapparams(x0)
xFMS = []; sFMS = []; % xFMS = fminsearch(fn, x0, struct('Display', 'off')); sFMS = mapparams(xFMS); fn(xFMS)

result = struct();
result.filename   = files{file_index};
result.frame_opts = frame.opts;
result.costmatrix = cost;
result.xFMS = xFMS;
result.sFMS = sFMS;
results{file_index, 1} = result;
end


for ri=1:length(results)
    r = results{ri};
    cost = r.costmatrix; imagesc(log(1 + cost(1:end, :)-min(cost(:))), [0, 13]); colorbar;
    [i, j] = ind2sub(size(r.costmatrix), find(r.costmatrix(:) == min(r.costmatrix(:))));
    x0 = mapparams(struct('sigma_beta', sbt_vec(j), 'sigma0', 1.05, 'pivec', pi_vec(i))); s0=mapparams(x0); r.s0 = s0;
    yticks = find(floor(-log10(pi_vec)) == -log10(pi_vec)); set(gca, 'YTick', yticks, 'YTickLabels', {'10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}', '10^{-0}'}); ylabel('\pi_1');
    xticks = find(floor(-log10(sbt_vec)) == -log10(sbt_vec)); set(gca, 'XTick', xticks, 'XTickLabels', {'10^{-6}','10^{-4}','10^{-2}','10^{-1}'}); xlabel('\sigma^2_{\beta}');
    title(sprintf('$$\\pi_1$$=%.1e, $$\\hat{\\pi_1}$$=%.2e, $$\\hat{\\sigma^2_{\\beta}}$$=%.2e', r.frame_opts.make_truebeta.pivec, r.s0.pivec,r.s0.sigma_beta^2),'Interpreter','Latex')
    saveas(gca,sprintf('%s/fit_map_1m_%i.png', outdir, ri), 'png');
    fprintf('%s, pi1=%.2e, sbt=%.2e, sig0=%.3f\n', r.filename, s0.pivec, s0.sigma_beta, s0.sigma0);
    %waitforbuttonpress
end


return;








    %if ~exist('frame2M, 'var'),
    %cd 'C:\Users\Oleksandr\Documents\GitHub\SynGP_norment'
    %find_config
    %load(fullfile(mmil_gwas_path, '../oleksanf/SynGP_GenCorr', 'gwaslogp_tpf9ebd0e1_802c_49ba_97e5_819d7b15e0db_pi=3e-03_gencorr=0.75_h2=0.90.mat'));
    %frame = reframe(frame, config, 'EUR_100K_2558K_ref');
    %clear 'config'
    %cd 'C:\Users\Oleksandr\Documents\GitHub\GMM'
    %end

    if ~exist('LDr2_p1sparse', 'var')
        t = 'H:\NORSTORE\MMIL';                       if exist(t, 'dir'), mmil_gwas_path = t; end;
        t = '/work/users/oleksanf/NORSTORE/MMIL';     if exist(t, 'dir'), mmil_gwas_path = t; end;
        t = '/space/syn03/1/data/GWAS';               if exist(t, 'dir'), mmil_gwas_path = t; end;
        if ~exist('mmil_gwas_path', 'var'), error('Unable to locate data folder'); end;

        load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/infomat.mat')); 
        load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/annomat.mat'))
        load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p1.mat'));LDmat = LDr2_p1sparse;
        load(fullfile(mmil_gwas_path, 'old_25_SNPs_processed_summary/GWAS_Annot/ldmat_2m_p8.mat'));
        load('H:\NORSTORE\oleksanf\w_ld_2m.mat')
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
    zvec = frame.gwaszscore(:, 1);
    Hvec = 2 * mafvec .* (1-mafvec);
    Nvec = frame.gwassize(:, 1);
    ref_ld = annomat(:, end);
    w_ld = ref_ld;
    mapparams = @GMM_ofrei_univariate_mapparams;
    fn = @(x)GMM_ofrei_univariate_cost(x, zvec, Hvec, Nvec, ref_ld, w_ld, mapparams, struct('x_limit', 20, 'x_step', 0.1, 'convolution', true));
    
    tmp = rand(length(zvec), 1); tmp(~isfinite(zvec)) = NaN;
    tmp_pruned = GenStats_FastPrune(tmp,LDr2_p8sparse);
    zvec(~isfinite(tmp_pruned)) = nan;

    % Check all points in a wide range of values across log space
    sbt_vec = logspace(-3, 0, 25);
    pi_vec = logspace(-7, 0, 29);
    sbt_mat = repmat(sbt_vec, [length(pi_vec), 1]); sbt_mat=sbt_mat(:);
    pi_mat = repmat(pi_vec, [1, length(sbt_vec)]); pi_mat=pi_mat(:);

    cost = zeros(size(sbt_mat)); sigma0 = 1;
    for i=1:length(sbt_mat), cost(i) = GMM_ofrei_univariate_cost(mapparams(struct('sigma_beta', sbt_mat(i), 'sigma0', sigma0, 'pivec', pi_mat(i))), zvec, Hvec, Nvec, ref_ld, w_ld, mapparams); end;
    cost = reshape(cost, [length(pi_vec), length(sbt_vec)]);
    %cost-min(cost(:)); image(cost-min(cost(:)));
    imagesc(log(1 + cost(1:end, :)-min(cost(:)))); colorbar;

    [i, j] = ind2sub(size(cost), find(cost(:) == min(cost(:))));
    x0 = mapparams(struct('sigma_beta', sbt_vec(j), 'sigma0', 1.05, 'pivec', pi_vec(i))); mapparams(x0)
    xFMS = []; sFMS = []; % xFMS = fminsearch(fn, x0, struct('Display', 'off')); sFMS = mapparams(xFMS); fn(xFMS)

    result = struct();
    result.filename   = files{frame_index};
    result.frame_opts = frame.opts;
    result.costmatrix = cost;
    result.x0 = x0;
    result.s0 = mapparams(x0)
    result.xFMS = xFMS;
    result.sFMS = sFMS;
    results{frame_index, 1} = result;
end


for frame_index=1:length(frames)
    %frame = load(fullfile('H:\NORSTORE\oleksanf\SynGP_GenCorr\2m', files(i).name)); frame = frame.frame;
    frame = frames{frame_index}.frame;
    fprintf('Processing frame %i\n', frame_index);
    if length(frame.gwasbeta) ~= 2558411, error('wrong template'); end;


    zvec = frame.gwaszscore(:, 1);

    nsnp = length(zvec);
    ivec_mhc = ((chrnumvec==6)&(posvec>=28477797)&(posvec<=33448354)); % From http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/region.cgi?name=MHC&asm=GRCh37

    TLDvec = annomat(:,end); Lvec = TLDvec;
    Hvec = 2*mafvec.*(1-mafvec);
    genvar = Hvec;

    defvec = isfinite(sum(zvec,2));
    nprune = 2;
    imat_prune = false(nsnp,nprune);
    for prunei = 1:nprune
      tmp = rand(nsnp,1);
      tmp(~defvec) = NaN;
      tmp_pruned = GenStats_FastPrune(tmp,LDmat);
      imat_prune(:,prunei) = isfinite(tmp_pruned);
    end

    options = [];
    resultStruct = GMM_univariate_fit(zvec,Hvec,Lvec,imat_prune);
    GMM_univariate_mapparams(resultStruct.x_fit, [])

    resultStructs{i} = resultStruct;
    resultStructs{i}.file = files(i).name;
    
    x=GMM_univariate_mapparams(resultStruct.x_fit, []);
    resultStructs{i}.pi1_hat = x.pi1;
    resultStructs{i}.pi1 = frame.opts.make_truebeta.pivec;
    resultStructs{i}
    %[cost fitstruct] = GMM_univariate_cost(resultStruct.x_fit,zvec,Hvec,Lvec,imat_prune,options); % Make plots

    % ToDo
    %  Compare fit and observed tail distributions (e.g., QQ-plots)
end

%resultStructs01 = resultStructs;
fileID = fopen('exp.txt','w');
for i=1:length(resultStructs)
    x=GMM_univariate_mapparams(resultStructs{i}.x_fit, []);
    %fprintf('true pi1 = %.5e - estimated pi1 = %.5e\n', resultStructs{i}.pi1(1), resultStructs{i}.pi1_hat(1));
    fprintf(fileID, '%.5e\t%.5e\n', resultStructs{i}.pi1(1), resultStructs{i}.pi1_hat(1));
end
fclose(fileID);