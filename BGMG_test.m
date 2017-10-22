addpath('DERIVESTsuite');

if ~exist('mafvec', 'var'),
    load(fullfile('test_data/ldmat_1m_p8.mat'),'mafvec');
    load('test_data/infomat.mat');

    if 0
        % old stuff
        ref_ld2 = load('test_data/1000G_EUR_Phase3_ref_ld2.mat');
        biased_ref_ld2 = load('test_data/1000G_EUR_Phase3_ref_ld2.mat');
        biased_ref_ld4 = load('test_data/1000G_EUR_Phase3_ref_ld4.mat');
        w_ld2 = load('test_data/1000G_EUR_Phase3_w_ld2.mat');
    else
        % re-generated on 2017.10.20
        ref_ld2 = load('H:/work/bgmg-annot/ref_l2.mat');
        biased_ref_ld4 = load('H:/work/bgmg-annot/biased_ref_l4.mat');
        biased_ref_ld2 = load('H:/work/bgmg-annot/biased_ref_l2.mat');
        w_ld2 = load('H:/work/bgmg-annot/w_ld.mat');
    end

    scz = load('test_data/PGC_SCZ_2014_noMHC.mat');
    bip = load('test_data/PGC_BIP_2012_lift_noMHC.mat');
    bip2 = load('test_data/PGC_BIP_2016_qc_noMHC.mat');
    cd  = load('test_data/IIBDGC_CD_2015_EUR_noMHC.mat');
    tg  = load('test_data/LIPIDS_TG_2013_noMHC.mat');
    
    mhc = (chrnumvec==6) & (posvec > 25e6) & (posvec < 35e6);
    scz.zvec(mhc)=nan;
    bip.zvec(mhc)=nan;
    bip2.zvec(mhc)=nan;
    tg.zvec(mhc)=nan;
    cd.zvec(mhc)=nan;
end

ref_ld = struct('sum_r2', biased_ref_ld2.annomat, 'chi_r4', biased_ref_ld4.annomat ./ biased_ref_ld2.annomat);
Hvec = 2*mafvec .* (1-mafvec);  w_ld  = w_ld2.annomat;

% mafvec_all = load('H:\work\bgmg-annot\mafvec_all_snps_1000G_Phase3_frq.mat')
% TotalHet =  2*sum(mafvec_all.mafvec .* (1-mafvec_all.mafvec))
TotalHET = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs

options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs
options.verbose = true;

if exist('do_SCZ', 'var') && do_SCZ
scz_result = BGMG_fit(scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, options);
BGMG_util.result2str(scz_result)
end

if exist('do_SCZ_BIP', 'var') && do_SCZ_BIP
options.use_legacy_impl = false; scz_bip_new = BGMG_fit([scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
%options.use_legacy_impl = true;  scz_bip_old = BGMG_fit([scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_bip_new.bivariate.params, [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
BGMG_util.result2str(scz_bip_new)
end

if exist('do_SCZ_BIP2', 'var') && do_SCZ_BIP2
options.use_legacy_impl = false; scz_bip2_new = BGMG_fit([scz.zvec bip2.zvec], Hvec, [scz.nvec bip2.nvec], w_ld, ref_ld, options);
%options.use_legacy_impl = true;  scz_bip2_old = BGMG_fit([scz.zvec bip2.zvec], Hvec, [scz.nvec bip2.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_bip2_new.bivariate.params, [scz.zvec bip2.zvec], Hvec, [scz.nvec bip2.nvec], w_ld, ref_ld, options);
BGMG_util.result2str(scz_bip2_new)
end

if exist('do_SCZ_TG', 'var') && do_SCZ_TG
options.use_legacy_impl = false; scz_tg_new = BGMG_fit([scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);
%options.use_legacy_impl = true;  scz_tg_old = BGMG_fit([scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_tg_new.bivariate.params, [scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);
BGMG_util.result2str(scz_tg_new)
end

if exist('do_BIP_TG', 'var') && do_BIP_TG
options.use_legacy_impl = false; bip_tg_new = BGMG_fit([bip.zvec tg.zvec], Hvec, [bip.nvec tg.nvec], w_ld, ref_ld, options);
%options.use_legacy_impl = true;  bip_tg_old = BGMG_fit([bip.zvec tg.zvec], Hvec, [bip.nvec tg.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(bip_tg_new.bivariate.params, [bip.zvec tg.zvec], Hvec, [bip.nvec tg.nvec], w_ld, ref_ld, options);
BGMG_util.result2str(bip_tg_new)
end

if exist('do_SCZ_CD', 'var') && do_SCZ_CD
options.use_legacy_impl = false; scz_cd_new = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
%options.use_legacy_impl = true;  scz_cd_old = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_cd_new.bivariate.params, [scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
BGMG_util.result2str(scz_cd_new)
end

if exist('do_TG_CD', 'var') && do_TG_CD
options.use_legacy_impl = false; tg_cd_new = BGMG_fit([tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);
%options.use_legacy_impl = true;  tg_cd_old = BGMG_fit([tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(tg_cd_new.bivariate.params, [tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);
BGMG_util.result2str(tg_cd_new)
end

return

o2 = options; o2.calculate_z_cdf = true; [cost, result_scz] = BGMG_univariate_cost(scz_bip_new.univariate{1}.params, scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, o2);
o2 = options; o2.calculate_z_cdf = true; [cost, result_bip] = BGMG_univariate_cost(scz_bip_new.univariate{2}.params, bip.zvec, Hvec, bip.nvec, w_ld, ref_ld, o2);
o2 = options; o2.calculate_z_cdf = true; [cost, result_cd] = BGMG_univariate_cost(scz_cd_new.univariate{2}.params, cd.zvec, Hvec, cd.nvec, w_ld, ref_ld, o2);
o2 = options; o2.calculate_z_cdf = true; [cost, result_tg] = BGMG_univariate_cost(bip_tg_new.univariate{2}.params, tg.zvec, Hvec, tg.nvec, w_ld, ref_ld, o2);


%scz_params = BGMG_fit(scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, options);
%bip_params = BGMG_fit(bip.zvec, Hvec, bip.nvec, w_ld, ref_ld, options);
%tg_params = BGMG_fit(tg.zvec, Hvec, tg.nvec, w_ld, ref_ld, options);
%cd_params = BGMG_fit(cd.zvec, Hvec, cd.nvec, w_ld, ref_ld, options);

%sbo = options; sbo.params1 = scz_params.univariate{1}.params; sbo.params2 = bip_params.univariate{1}.params;
scz_bip_params2 = BGMG_fit([scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, sbo);

BGMG_bivariate_cost(scz_bip_params.bivariate.params, [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_tg_newU.bivariate.params, [scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);

sbo = options; sbo.params1 = scz_params.univariate{1}.params; sbo.params2 = tg_params.univariate{1}.params;
scz_tg_params2 = BGMG_fit([scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, sbo);

sbo = options; sbo.params1 = scz_params.univariate{1}.params; sbo.params2 = cd_params.univariate{1}.params;
scz_cd_params2 = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, sbo);

sbo = options; sbo.params1 = tg_params.univariate{1}.params; sbo.params2 = cd_params.univariate{1}.params;
tg_cd_params = BGMG_fit([tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, sbo);

results.bivariate.params = struct('pi_vec', [1e-7 1e-7 0.003 ], ...
       'rho_beta', [0 0 0.835 ], ...
       'sig2_beta', [5.06e-05 0  5.06e-05; 0 1.02e-05 1.02e-05], ...
       'rho_zero', 0.181, ...
       'sig2_zero', [1.209; 1.051 ]);
%m=eye(9);m(3:end, 3:end)=inv(hess(3:end, 3:end));
%ci_sample = mvnrnd(BGMG_mapparams3(results.bivariate.params), m, options.ci_sample);


if 0
% univariate cost function
%BGMG_univariate_cost(struct('pivec', 0.1, 'sigma0', 1.05, 'sigma_beta', 1e-3), scz.zvec, Hvec, scz.nvec, w_ld, ref_ld)
%BGMG_bivariate_cost(struct('pivec', [0.1 0.2 0.3], 'sigma0', [1.05 1.1], 'sigma_beta', [1e-3, 1e-4], 'rho0', -0.1, 'rho_beta', 0.8), [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld)

% Univariate two-component mixture model (null + causal)
BGMG_univariate_cost(struct('pi_vec', 0.1, 'sig2_zero', 1.05, 'sig2_beta', 1e-3), scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, options)

% Univariate three-component mixture model (null + causal1 + causal2), example: AD APOE
BGMG_univariate_cost(struct('pi_vec', [0.1 0.2], 'sig2_zero', 1.05, 'sig2_beta', [1e-3 1e-4]), scz.zvec, Hvec, scz.nvec, w_ld, ref_ld)

% Bivariate two-component mixture model (null + pleio)
BGMG_bivariate_cost(struct('pi_vec', [0.1], 'sig2_zero', [1.05; 1.10], 'rho_zero', -0.1, 'sig2_beta', [1e-3; 1e-4], 'rho_beta', -0.5), ...
                    [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld)

% Bivariate three-component mixture model (null + indep1 + indep2)
BGMG_bivariate_cost(struct('pi_vec', [0.1 0.2], 'sig2_zero', [1.05; 1.10], 'rho_zero', -0.1, 'sig2_beta', [1e-3 0; 0 1e-4], 'rho_beta', [0 0]), ...
                    [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld)

% Bivariate four-component mixture model (null + indep1 + indep2)
BGMG_bivariate_cost(struct('pi_vec', [0.1 0.2 0.3], 'sig2_zero', [1.05; 1.10], 'rho_zero', -0.1, 'sig2_beta', [1e-3 0 1e-3; 0 1e-4 1e-4], 'rho_beta', [0 0 0.5]), [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options)


BGMG_bivariate_cost(struct('pi_vec', [2.649e-03 1.128e-06 4.255e-05 ], 'rho_beta', [0 0 0.050 ], 'sig2_beta', [5.92e-05 0 5.92e-05; 0 2.76e-03 2.76e-03], 'rho_zero', 0.067, 'sig2_zero', [1.181; 1.092 ]), [scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options)
end

if 0
    if ~isfield(options, 'calculate_z_cdf_limit'), options.calculate_z_cdf_limit = 15; end;
    if ~isfield(options, 'calculate_z_cdf_step'), options.calculate_z_cdf_step = 0.25; end;
    z_grid =  (-options.calculate_z_cdf_limit:options.calculate_z_cdf_step:options.calculate_z_cdf_limit);
    
    clf
    zvec = tg.zvec; result_cdf = result_tg.cdf;
    numstrata = 4;
    
    x = ref_ld.sum_r2;
    y = Hvec;
    % y = ref_ld.chi_r4;
    defvec=isfinite(zvec+x+y);
    zvec=zvec(defvec); x=x(defvec); y=y(defvec); result_cdf = result_cdf(defvec, :);
    
    strat = false(numstrata,numstrata,sum(defvec));
    xq = [-Inf, quantile(x,numstrata-2), +Inf];
    yq = [-Inf, quantile(y,numstrata-2), +Inf];
    for i=1:numstrata
        for j=1:numstrata
            idx = true(size(x));
            titl = '';
            if i ~= numstrata, idx = idx & ((x >= xq(i)) & (x <= xq(i+1))); titl = sprintf('%s%.1f<=TLD <=%.1f', titl, xq(i), xq(i+1)); end;
            if i ~= numstrata && j ~= numstrata, titl = sprintf('%s\n', titl); end;
            if j ~= numstrata, idx = idx & ((y >= yq(j)) & (y <= yq(j+1))); titl = sprintf('%s%.3f<=HVEC<=%.3f', titl, yq(j), yq(j+1)); end;
            subplot(numstrata,numstrata, (i-1)*numstrata+j);
            title(titl);
            
            hold on
            qqlim=6;
            plot(-log10(nanmean(result_cdf(idx, :))),-log10(normcdf(z_grid,0,1)),'b'); 

            zvecI = sort(zvec(idx));
            plot(-log10((1:length(zvecI)) / length(zvecI)),-log10(normcdf(zvecI, 0, 1)),'g'); 

            plot([0 qqlim],[0 qqlim], 'k');
            xlim([0 qqlim]); ylim([0 qqlim]);
        end
    end
    
    
    %legend('model-conv', 'model-amd', 'model-kl', 'null', 'Location', 'SouthEast');
end