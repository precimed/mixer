if 1
%  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/LDmat_2558411_01_2M.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/LDmat_2558411_08_2M.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/infomat.mat','mafvec','chrnumvec','posvec','snpidlist');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/annomat.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/PGC5_SCZ.mat'); zvec_pgc5_scz = -norminv(10.^-logpvec_pgc5_scz).*sign(zvec_pgc5_scz);
%  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/PGC2_BIP.mat'); zvec_pgc2_bip = -norminv(10.^-logpvec_pgc2_bip).*sign(zvec_pgc2_bip);
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/PGC_BIP32_2.5m.mat'); zvec_pgc_bip32 = -norminv(10.^-logpvec_pgc_bip32).*sign(zvec_pgc_bip32);
  load('~/GWAS_Data/PGC_ADHD.mat'); zvec_pgc_adhd = -norminv(10.^-logpvec_pgc_adhd).*sign(zvec_pgc_adhd);
%  zmat = cat(2,zvec_pgc_adhd,zvec_pgc5_scz); nsnp = size(zmat,1);
  zmat = cat(2,zvec_pgc_adhd,zvec_pgc_bip32); nsnp = size(zmat,1);
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
  logpmat = -log10(normcdf(-abs(zmat))*2);
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

resultStruct = GMM_bivariate_fit(zmat,Hvec,Lvec,imat_prune,options_fit);

% 

% ADHD & BIP 2.5m
R      pi1: 0.0021
      pi2: 5.4589e-13
      pi3: 0.0040
    sig01: 0.9898
    sig02: 1.0230
    sigb1: 0.3388
    sigb2: 1.4973
     rho0: 0.0196
     rhob: 0.3426
      pi0: 0.9939


% ADHD & SCZ 2.5m
      pi1: -2.2204e-16 - 1.8216e-35i
      pi2: 1.6995e-13
      pi3: 0.0030
    sig01: 0.9914
    sig02: 1.0673
    sigb1: 0.3558
    sigb2: 2.4080
     rho0: 0.0217
     rhob: 0.1136
      pi0: 0.9970 + 0.0000i

