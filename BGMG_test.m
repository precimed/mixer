addpath('DERIVESTsuite');

if ~exist('mafvec', 'var'),
    load(fullfile('test_data/ldmat_1m_p8.mat'),'mafvec');
    ref_ld2 = load('test_data/1000G_EUR_Phase3_ref_ld2.mat');
    ref_ld4 = load('test_data/1000G_EUR_Phase3_ref_ld4.mat');
    w_ld2 = load('test_data/1000G_EUR_Phase3_w_ld2.mat');
    scz = load('test_data/PGC_SCZ_2014_noMHC.mat');
    bip = load('test_data/PGC_BIP_2012_lift_noMHC.mat');
    cd  = load('test_data/IIBDGC_CD_2015_EUR_noMHC.mat');
    tg  = load('test_data/LIPIDS_TG_2013_noMHC.mat');
end

ref_ld = struct('sum_r2', ref_ld2.annomat, 'sum_r4', ref_ld4.annomat);
Hvec = 2*mafvec .* (1-mafvec);  w_ld  = w_ld2.annomat;
TotalHET = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs

% univariate cost function
%BGMG_univariate_cost(struct('pivec', 0.1, 'sigma0', 1.05, 'sigma_beta', 1e-3), scz.zvec, Hvec, scz.nvec, w_ld, ref_ld)
%BGMG_bivariate_cost(struct('pivec', [0.1 0.2 0.3], 'sigma0', [1.05 1.1], 'sigma_beta', [1e-3, 1e-4], 'rho0', -0.1, 'rho_beta', 0.8), [scz.zvec bip.zvec], Hvec, [scz.nvec bip.nvec], w_ld, ref_ld)

% Univariate two-component mixture model (null + causal)
BGMG_univariate_cost(struct('pi_vec', 0.1, 'sig2_zero', 1.05, 'sig2_beta', 1e-3), scz.zvec, Hvec, scz.nvec, w_ld, ref_ld)

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

options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs
options.verbose = true;

scz_params = BGMG_fit(scz.zvec, Hvec, scz.nvec, w_ld, ref_ld, options);
bip_params = BGMG_fit(bip.zvec, Hvec, bip.nvec, w_ld, ref_ld, options);
tg_params = BGMG_fit(tg.zvec, Hvec, tg.nvec, w_ld, ref_ld, options);
cd_params = BGMG_fit(cd.zvec, Hvec, cd.nvec, w_ld, ref_ld, options);

sbo = options;
sbo.params1 = scz_params.univariate{1}.params;
sbo.params2 = bip_params.univariate{1}.params;
scz_bip_params = BGMG_fit([scz.zvec bip.zvec], Hvec, [scz.nvec scz.nvec], w_ld, ref_ld, options);

scz_tg_params2 = BGMG_fit([scz.zvec tg.zvec], Hvec, [scz.nvec tg.nvec], w_ld, ref_ld, options);

scz_cd_params = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, options);
scz_cd_params = BGMG_fit([scz.zvec cd.zvec], Hvec, [scz.nvec cd.nvec], w_ld, ref_ld, sbo);

tg_cd_params = BGMG_fit([tg.zvec cd.zvec], Hvec, [tg.nvec cd.nvec], w_ld, ref_ld, options);

results.bivariate.params = struct('pi_vec', [1e-7 1e-7 0.003 ], ...
       'rho_beta', [0 0 0.835 ], ...
       'sig2_beta', [5.06e-05 0  5.06e-05; 0 1.02e-05 1.02e-05], ...
       'rho_zero', 0.181, ...
       'sig2_zero', [1.209; 1.051 ]);
%m=eye(9);m(3:end, 3:end)=inv(hess(3:end, 3:end));
%ci_sample = mvnrnd(BGMG_mapparams3(results.bivariate.params), m, options.ci_sample);
