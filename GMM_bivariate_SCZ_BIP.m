if 0
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

% Fit bivariate model

%options_fit = struct();
options_fit = struct('rho0',0,'ngm',2);

x0 = GMM_bivariate_mapparams(struct('sig01',1.0,'sig02',1.0,'rho0',0,'pivec',[1e-1 2e-2],'sig1vec',[1.0 0.5],'sig2vec',[1.2 2.3],'rhovec',[0.4 0.95]),options_fit)

resultStruct = GMM_bivariate_histfit(sumstats,options_fit,x0);

% ToDo
%   Should use more aggressive pruning LD R^2 threshold of 0.1
%   Should incorporate Neff for each trait

%   Weight fit (cost function) towards tails?
 
%   Try 2-D smoothing
%     Of counts
%     Of log-transformed counts
%     Of counts divided by (smoothed?) standard deviation, to make noise similar across space
%     Look for code for 2-D smoothing of histograms / count data
%     Look for polar coordinate smoothing approach that enforces smoothness as function of angle (Fourier basis?) and monotonic and smooth log transformed radial function

%   Re-scale estimates, to be consistent with causal betas (show predicted and observed z-score histograms for different ranges of L and H)

%   Also look at kernelcmi (kernel conditional mutual information)
%   Revisit copula fitting -- see if it provides good fit to bivariate data
% 

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

% Should also try with overlapping PGC2 sub-study data, check discovery / replication z-scores

