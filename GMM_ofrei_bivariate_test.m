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
zvec = frame.gwaszscore; zmat=zvec;
Hvec = 2 * mafvec .* (1-mafvec);
Nvec = frame.gwassize; Nmat=Nvec;
ref_ld = baseline_ld.annomat(:, 1);
w_ld = eur_w_ld.annomat;
mapparams = @GMM2_bivariate_mapparams;
fn = @(x)GMM2_univariate_cost(x, zvec, Hvec, Nvec, ref_ld, w_ld, mapparams, struct('x_limit', 20, 'x_step', 0.1, 'convolution', true));




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
        load('H:\NORSTORE\oleksanf\tld_hist_2m.mat')
        
        ref_ld = struct();
        ref_ld.sum_r2 = annomat(:, end);
        ref_ld.sum_r4 = tld_hist * mean(tld_bins, 2) .^ 2 .* (annomat(:, end) ./ sum(tld_hist, 2));
        hist(ref_ld.sum_r4 ./ ref_ld.sum_r2, 25)
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

results={};
for i=1:length(frames)
%    params = struct('pivec', [0.001 0.930 0.017 ], 'rho_beta', 0.315, 'sigma_beta', sqrt([2.74e-05 1.24e-05 ]), 'rho0', 0.059, 'sigma0', sqrt([1.098 1.101 ]));
%    [~, result] = GMM2_bivariate_cost(params, frames{i}.frame.gwaszscore, Hvec, frames{i}.frame.gwassize, w_ld, ref_ld, [], struct('verbose', true, 'calculate_beta_hat', false, 'calculate_global_fdr', false));
        

    results{i} = GMM2_fit(frames{i}.frame.gwaszscore, Hvec, frames{i}.frame.gwassize, w_ld, ref_ld, options)
  %results{i}.first = GMM2_fit(frames{i}.frame.gwaszscore(:, 1), Hvec, frames{i}.frame.gwassize(:, 1), w_ld, ref_ld, options);
  %results{i}.second = GMM2_fit(frames{i}.frame.gwaszscore(:, 2), Hvec, frames{i}.frame.gwassize(:, 2), w_ld, ref_ld, options);
end
for i=1:length(frames)
    %p = results{i}.first.univariate{1}.params; fprintf('file=%s sigma0=%.3f sigma_beta=%.3e pi=%.3e\n', files{i}, p.sigma0, p.sigma_beta, p.pivec);
    %p = results{i}.second.univariate{1}.params; fprintf('file=%s sigma0=%.3f sigma_beta=%.3e pi=%.3e\n', files{i}, p.sigma0, p.sigma_beta, p.pivec);
    p = results{i}.univariate{1}.params; fprintf('file=%s sigma0=%.3f sigma_beta=%.3e pi=%.3e\n', files{i}, p.sigma0, p.sigma_beta, p.pivec);
    p = results{i}.univariate{2}.params; fprintf('file=%s sigma0=%.3f sigma_beta=%.3e pi=%.3e\n', files{i}, p.sigma0, p.sigma_beta, p.pivec);
    p = results{i}.bivariate.params; disp(p); %fprintf('file=%s sigma0=%.3f sigma_beta=%.3e pi=%.3e\n', files{i}, p.sigma0, p.sigma_beta, p.pivec);
end

if 0
    options = struct('verbose', true);
    z = @(logpvec, zvec) -norminv(10.^-logpvec/2).*sign(zvec);
    Hvec = 2*mafvec .* (1-mafvec); w_ld = 1./w_ld_p1; %ref_ld = annomat(:, end);
    scz = load('H:\NORSTORE\MMIL\old_25_SNPs_processed_summary\GWAS_Data\PGC5_SCZ.mat');  scz.zvec=z(scz.logpvec_pgc5_scz, scz.zvec_pgc5_scz);
    cog = load('H:\NORSTORE\MMIL\old_25_SNPs_processed_summary\GWAS_Data\COG_charge.mat');cog.zvec=z(cog.logpvec_cog, cog.zvec_cog);
    tg  = load('H:\NORSTORE\MMIL\old_25_SNPs_processed_summary\GWAS_Data\TG_2013.mat');   tg.zvec=z(tg.logpvec_tg_2013,tg.zvec_tg_2013);
    iq  = load('H:\NORSTORE\GWAS_SUMSTAT\MAT_2M\CTG_IQ_2017_2m.mat');iq.zvec=z(iq.logpvecIQ, iq.zvecIQ);
    Nmat_dummy = 50000 * ones(length(Hvec), 2);
    figure(1); results_scz_cog = GMM2_fit([scz.zvec, cog.zvec], Hvec, Nmat_dummy, w_ld, ref_ld, options)
    figure(2); results_scz_tg = GMM2_fit([scz.zvec, tg.zvec], Hvec, Nmat_dummy, w_ld, ref_ld, options)
    figure(3); results_cog_tg = GMM2_fit([cog.zvec, tg.zvec], Hvec, Nmat_dummy, w_ld, ref_ld, options)
    
    figure(4); results_scz_iq = GMM2_fit([scz.zvec, iq.zvec], Hvec, Nmat_dummy, w_ld, ref_ld, options)
    
    results_scz = GMM2_fit(scz.zvec, Hvec, Nmat_dummy(:, 1), w_ld, ref_ld, options); results_scz = results_scz.univariate{1};
    results_iq  = GMM2_fit(iq.zvec, Hvec, Nmat_dummy(:, 2), w_ld, ref_ld, options); results_iq = results_iq.univariate{1};

    results_scz.params
    results_iq.params
    figure(4); results_scz_iq = GMM2_fit([scz.zvec, iq.zvec], Hvec, Nmat_dummy, w_ld, ref_ld, struct('verbose', true, 'params1', results_scz.params, 'params2', results_iq.params))

    costmap = results_scz.univariate{1}.costmap; imagesc(log(1 + costmap(1:end, :)-min(costmap(:)))); colorbar; 
    
    zmat= [scz.zvec, cog.zvec];
    stat = results.bivariate.conjfdr; [~, idx]=sort(stat); idx = idx(1:10:300); [zmat(idx,:), stat(idx), -log10(stat(idx))]
    stat = results.bivariate.condfdr1; [~, idx]=sort(stat); idx = idx(1:100); [zmat(idx,:), stat(idx), -log10(stat(idx))]

    GMM2_bivariate_cost(struct('sigma0', [1.0645 1.0082], 'rho0', -0.0169, 'sigma_beta', [0.0085 0.0040], 'rho_beta', -0.3046, 'pivec', [0.0762 0.0721 0.1837]), [scz.zvec, cog.zvec], Hvec, Nmat_dummy, w_ld, ref_ld, @GMM2_bivariate_mapparams, options)
    
    GMM2_univariate_cost(struct('sigma0', [1.0645], 'sigma_beta', 0.0085 , 'pivec', 0.0762), scz.zvec, Hvec, Nmat_dummy(:, 1), w_ld, ref_ld, [], struct())
    
    results_scz = GMM2_fit(scz.zvec, Hvec, Nmat_dummy(:, 1), w_ld, ref_ld, struct('logging', true, 'plot_costlines', true, 'plot_costmaps', true));
    [~, result] = GMM2_univariate_cost(results_scz.univariate{1}.params, scz.zvec, Hvec, Nmat_dummy(:, 1), w_ld, ref_ld, [], struct('zmax', +Inf))
    idx =1:100:length(Hvec); a=[result.beta_hat_mean(idx) result.beta_hat_se(idx) result.beta_hat_median(idx) result.beta_hat_mean(idx) ./ result.beta_hat_se(idx) scz.zvec(idx)]; scatter(a(:, 5), a(:, 4), 10, Hvec(idx) * [2 0 0], 'filled')
    
    idx=1:1000;
    params= results_scz_iq.bivariate.params;
    params = struct('pivec', [5e-3, 2e-3, 1e-3], 'sigma_beta', [1e-3, 1e-3], 'sigma0', [1.17 1.05], 'rho0', 0.002, 'rho_beta', -0.3);
    [~, result] = GMM2_bivariate_cost(params, [scz.zvec(idx), iq.zvec(idx)], Hvec(idx), Nmat_dummy(idx, :), w_ld(idx), struct('sum_r2', ref_ld.sum_r2(idx), 'sum_r4', ref_ld.sum_r4(idx)), options);
    %idx=1:10; [result.beta_hat_mean(idx, :) result.beta_hat_var(idx, :) result.beta_hat_cov(idx, :)]
    plot(result.beta_hat_mean(:, :))
    plot(result.beta_hat_median(:, :))
    plot(-log10(result.condfdr_global(:, :)))hmnhgbvc =-0o9,mxiiiiiq
    plot(-log10(result.condfdr_local(:, :)))
    plot(-log10(result.conjfdr_global(:, :)))
    plot(-log10(result.conjfdr_local(:, :)))
    plot([scz.zvec(idx, :) iq.zvec(idx, :)])
    
    % SCZ vs IQ
    %Bivariate : pi=[0.003 0.001 0.001 ], rho_beta=-0.294, sigma_beta^2=[1.29e-04 6.60e-05 ], (eff. [1.21e+00 6.21e-01 ]), rho0=-0.018, sigma0^2=[1.188 1.076 ], cost=3.455e+05
    
    
    %idx=1:100000; scatter(scz.zvec(idx), result.beta_hat_mean(idx, 1) ./ sqrt(result.beta_hat_var(idx, 1)), 10, Hvec(idx) * [2 0 0], 'filled')
    %idx=1:100000; scatter(iq.zvec(idx), result.beta_hat_mean(idx, 2) ./ sqrt(result.beta_hat_var(idx, 2)), 10, abs(scz.zvec(idx)) * [0.3 0 0], 'filled')
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
    ref_ld = annomat(:, end);
    w_ld = ref_ld;
    
    fn = @(x)GMM2_bivariate_cost(x, zmat, Hvec, Nmat, w_ld, @GMM2_bivariate_mapparams);
    [cost, pdf, condfdr1, condfdr2, conjfdr] = fn(GMM2_bivariate_mapparams(struct('sigma_beta', [0.123 0.234], 'rho_beta', 0.2, 'sigma0', [1.05 1.10], 'rho0', 0.02, 'pivec', [0.01 0.01 0.05])));
    
    mapparams = @GMM2_bivariate_mapparams;
    fn = @(x)GMM_ofrei_bivariate_cost(x, zvec, Hvec, Nvec, ref_ld, w_ld, mapparams, struct('x_limit', 20, 'x_step', 0.1));
    
    tmp = rand(length(zvec), 1); tmp(~isfinite(sum(zvec, 2))) = NaN;
    tmp_pruned = GenStats_FastPrune(tmp,LDr2_p8sparse);
    zvec(~isfinite(tmp_pruned), :) = nan;

    if 0
        mapparams = @GMM2_bivariate_mapparams;
        params=(struct('sigma0', [1.05 1.10], 'rho0', 0.01, 'sigma_beta', [5e-5 5e-5], 'rho_beta', 0.7, 'pivec', [0.25 0.25 0.25]));
        options=struct('x_step', 0.5, 'x_limit', 20);
        
      %  GMM_ofrei_bivariate_cost(mapparams(params), zvec, Hvec, Nvec, ref_ld, w_ld, mapparams, options);
        tld=20; ld_bins=10; GMM_ofrei_bivariate_plot_callback([], struct('sigma0', [1.05 1.10], 'rho0', 0.01,...
            'sigma_beta', [1e-4, 1e-4], 'rho_beta', -0.95, 'pivec', [0.001 0.001 0.001]), tld, true, ld_bins);
    end
    
    % Check all points in a wide range of values across log space
    sbt_vec = logspace(-3, 0, 25);
    pi_vec = logspace(-7, 0, 29);
    sbt_mat = repmat(sbt_vec, [length(pi_vec), 1]); sbt_mat=sbt_mat(:);
    pi_mat = repmat(pi_vec, [1, length(sbt_vec)]); pi_mat=pi_mat(:);

    cost = zeros(size(sbt_mat)); sigma0 = 1;
    for i=1:length(sbt_mat), cost(i) = GMM2_univariate_cost(mapparams(struct('sigma_beta', sbt_mat(i), 'sigma0', sigma0, 'pivec', pi_mat(i))), zvec, Hvec, Nvec, ref_ld, w_ld, mapparams); end;
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