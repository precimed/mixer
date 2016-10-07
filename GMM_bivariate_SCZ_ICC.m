if 1
%  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/LDmat_2558411_01_2M.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/LDmat_2558411_08_2M.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/infomat.mat','mafvec','chrnumvec','posvec','snpidlist');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/annomat.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/ICC_2015.mat'); % zvec_icc seems incorrectly scaled -- check with Yunpeng
  zvec_icc_corr = -norminv(10.^-logpvec_icc/2).*sign(zvec_icc);
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/PGC5_SCZ.mat');
  zmat = cat(2,zvec_icc_corr,zvec_pgc5_scz); nsnp = size(zmat,1);
else
%  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/LDmat_9m_01.mat');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/LDmat_9m_08.mat');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/infomat.mat','mafvec','chrnumvec','posvec','snpidlist');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/annomat.mat');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC_BIP32_9m.mat'); zvec_pgc_bip32 = -norminv(10.^-logpvec_pgc_bip32).*sign(zvec_pgc_bip32);
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC2_SCZ52.mat'); zvec_pgc2_scz = -norminv(10.^-logpvec_pgc2_scz).*sign(zvec_pgc2_scz);
  load('~/Downloads/BPvsCONT.txt.mat'); zvec_bip_noscz = -norminv(10.^-logpvec/2).*sign(zvec);
  load('~/Downloads/SCZvsCONT.txt.mat'); zvec_scz_nobip = -norminv(10.^-logpvec/2).*sign(zvec);
%  zmat = cat(2,zvec_pgc_bip32,zvec_pgc2_scz); nsnp = size(zmat,1);
  zmat = cat(2,zvec_bip_noscz,zvec_scz_nobip); nsnp = size(zmat,1);
end

var(zmat(isfinite(sum(zmat,2)),:),1)

Hvec = 2*mafvec.*(1-mafvec);
Lvec = annomat(:,end);

defvec = isfinite(sum(zmat,2));
nprune = 2;
imat_prune = false(nsnp,nprune);
for prunei = 1:nprune
  tmp = rand(nsnp,1);
  tmp(~defvec) = NaN;
  tmp_pruned = GenStats_FastPrune(tmp,LDmat);
  imat_prune(:,prunei) = isfinite(tmp_pruned);
end

% Plot observed vs. predicted univariate z-score distributions
% Should fit with constraits
%   sig0 = 1
%   pi1 = 1 => variance proportional to Lvec
% Should modify GMM_univariate_cost to allow for two non-zero causal SNPs (get back quadratic fitting code from Dominic), or many (variance proportional to Lvec)
% Should incorporate Neff for each trait
% Should map cost at different combos of pi1 and sigb (and heritabilities: pi1*sigb^2)
% Sample from posterior distribution (MCMC?)

% Fit bivariate model

options_fit = struct();
%options_fit = struct('rho0',0);
options_fit = struct('pi1',1e-6,'pi2',1e-6,'pi3',1-2*1e-6);

resultStruct = GMM_bivariate_fit(zmat,Hvec,Lvec,imat_prune,options_fit);

% ToDo
%   Check with height or other well powered, independent phenotype 

% ICC & SCZ 2.5m, w. rho0 = 0 constraint only:
      pi1: 9.1356e-04
      pi2: 1.3353e-10
      pi3: 0.0069
    sig01: 1.0051
    sig02: 1.1038
    sigb1: 0.5111
    sigb2: 2.1008
     rho0: 0
     rhob: 0.5019
      pi0: 0.9922

% ICC & SCZ 2.5m
      pi1: 9.6083e-11
      pi2: 1.2268e-08
      pi3: 0.0054
    sig01: 1.0098
    sig02: 1.1214
    sigb1: 0.5058
    sigb2: 2.1923
     rho0: 0.0202
     rhob: 0.3810
      pi0: 0.9946

% ICC & SCZ 2.5m, w. Pi3 = 1 constraint
      pi1: 1.0000e-06
      pi2: 1.0000e-06
      pi3: 1.0000
    sig01: 1.0063
    sig02: 1.0995
    sigb1: 0.2735
    sigb2: 1.0640
     rho0: 0.0169
     rhob: 0.3404
      pi0: 0


% BIP & SCZ 9m
      pi1: 1.2636e-04
      pi2: 1.8177e-05
      pi3: 0.0024
    sig01: 1.0033
    sig02: 1.0314
    sigb1: 1.6068
    sigb2: 2.9150
     rho0: 0.0975
     rhob: 0.6640
      pi0: 0.9974

% Non-overlapping BIP & SCZ
      pi1: 6.3764e-04
      pi2: -2.2204e-16
      pi3: 0.0052
    sig01: 1.0074
    sig02: 1.0383
    sigb1: 1.4208
    sigb2: 2.3265
     rho0: 0.0451
     rhob: 0.6502
      pi0: 0.9942

w. rho0==0:
      pi1: 0.0010
      pi2: -2.2204e-16
      pi3: 0.0064
    sig01: 1.0041
    sig02: 1.0339
    sigb1: 1.3635
    sigb2: 2.2267
     rho0: 0
     rhob: 0.7587
      pi0: 0.9925

