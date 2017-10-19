addpath('DERIVESTsuite');

if ~exist('mafvec', 'var'),
    load(fullfile('test_data/ldmat_1m_p8.mat'),'mafvec');
    load('test_data/infomat.mat');
    ref_ld2 = load('test_data/1000G_EUR_Phase3_ref_ld2.mat');
    ref_ld4 = load('test_data/1000G_EUR_Phase3_ref_ld4.mat');
    w_ld2 = load('test_data/1000G_EUR_Phase3_w_ld2.mat');
    scz = load('test_data/PGC_SCZ_2014_noMHC.mat');
    bip = load('test_data/PGC_BIP_2012_lift_noMHC.mat');
    cd  = load('test_data/IIBDGC_CD_2015_EUR_noMHC.mat');
    tg  = load('test_data/LIPIDS_TG_2013_noMHC.mat');
    
    mhc = (chrnumvec==6) & (posvec > 25e6) & (posvec < 35e6);
    scz.zvec(mhc)=nan;
    bip.zvec(mhc)=nan;
    tg.zvec(mhc)=nan;
    cd.zvec(mhc)=nan;
end

ref_ld = struct('sum_r2', ref_ld2.annomat, 'sum_r4', ref_ld4.annomat);
Hvec = 2*mafvec .* (1-mafvec);  w_ld  = w_ld2.annomat;
TotalHET = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs

options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs
options.verbose = true;

if exist('do_SCZ_BIP', 'var') && do_SCZ_BIP
options.use_legacy_impl = false; scz_bip_new = BGMG_fit([scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
options.use_legacy_impl = true;  scz_bip_old = BGMG_fit([scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_bip_new.bivariate.params, [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld, options);
end

if exist('do_SCZ_TG', 'var') && do_SCZ_TG
options.use_legacy_impl = false; scz_tg_new = BGMG_fit([scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);
options.use_legacy_impl = true;  scz_tg_old = BGMG_fit([scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_tg_new.bivariate.params, [scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);
end

if exist('do_SCZ_CD', 'var') && do_SCZ_CD
options.use_legacy_impl = false; scz_cd_new = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
options.use_legacy_impl = true;  scz_cd_old = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(scz_cd_new.bivariate.params, [scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
end

if exist('do_TG_CD', 'var') && do_TG_CD
options.use_legacy_impl = false; tg_cd_new = BGMG_fit([tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);
options.use_legacy_impl = true;  tg_cd_old = BGMG_fit([tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);
BGMG_bivariate_cost(tg_cd_new.bivariate.params, [tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);
end

return 

%scz_params = BGMG_fit(scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, options);
%bip_params = BGMG_fit(bip.zvec, Hvec, bip.nvec, w_ld, ref_ld, options);
%tg_params = BGMG_fit(tg.zvec, Hvec, tg.nvec, w_ld, ref_ld, options);
%cd_params = BGMG_fit(cd.zvec, Hvec, cd.nvec, w_ld, ref_ld, options);

%BGMG_univariate_cost(scz_params.univariate{1}.params, scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, options);

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