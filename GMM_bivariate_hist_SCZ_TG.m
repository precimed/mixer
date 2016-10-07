
% Figure out actual Neff for PGC2 SCZ and BIP
if 0
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC2_BIP32_multisites_9m.mat');
  Neff_bip = sum(2./(1./Nvec_A+1./Nvec_U));
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC2_SCZ52_multisites.mat');
  Neff_scz = sum(2./(1./N_A_vec+1./N_U_vec));
end

if 1
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/LDmat_2558411_01_2M.mat'); Dmat_01 = LDmat; % Need to get new LDmats from Dominic; coordinate directory location w. Oslo
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/LDmat_2558411_08_2M.mat'); % Need to get new LDmats from Dominic; coordinate directory location w. Oslo
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/infomat.mat','mafvec','chrnumvec','posvec','snpidlist');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Annot/annomat.mat');
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/ICC_2015.mat'); Neff_ICC = 12684; % zvec_icc seems incorrectly scaled -- check with Yunpeng
  zvec_icc_corr = -norminv(10.^-logpvec_icc/2).*sign(zvec_icc);
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/PGC5_SCZ.mat');
  load('~/GWAS_Data/COG_charge.mat'); Neff_cog = 53949/2; % Adjust for "per arm"
  zvec_cog_corr = -norminv(10.^-logpvec_cog/2).*sign(zvec_cog);
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/GIANT_Height_2014.mat'); Neff_height = NaN; zvec_giant_height_corr = -norminv(10.^-logpvec_gaint_height_2014/2).*sign(zvec_gaint_height_2014);
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/BLIP_TG.mat'); Neff_tg = NaN; 
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/23andme_neuroticism.mat'); Neff_neuroticism = 59225/2; 
  zvec_23andme_neuroticism_corr = -norminv(10.^-logpvec_23andme_neuroticism/2).*sign(zvec_23andme_neuroticism);
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/ENIGMA2_Putamen.mat'); Neff_putamen = 12500/2;
  load('/space/syn03/1/data/GWAS/old_25_SNPs_processed_summary/GWAS_Data/MANU_PD.mat'); Neff_pd = NaN;
%  zmat = cat(2,zvec_icc_corr,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_ICC Neff_pgc2_scz]; traitlist = {'ICC' 'SCZ'};
%  zmat = cat(2,zvec_cog_corr,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_cog Neff_pgc2_scz]; traitlist = {'COG' 'SCZ'}; 
  zmat = cat(2,zvec_blip_tg,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_tg Neff_pgc2_scz]; traitlist = {'TG' 'SCZ'}; % Need to locate new version of TG summary stats
%  zmat = cat(2,zvec_giant_height_corr,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_height Neff_pgc2_scz]; traitlist = {'Height' 'SCZ'}; 
%  zmat = cat(2,zvec_23andme_neuroticism_corr,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_neuroticism Neff_pgc2_scz]; traitlist = {'Neuroticism' 'SCZ'};
%  zmat = cat(2,zvec_putamen,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_putamen Neff_pgc2_scz]; traitlist = {'Putamen' 'SCZ'}; % Overlapping samples?
%  zmat = cat(2,zvec_pd,zvec_pgc5_scz); nsnp = size(zmat,1); Neffvec = [Neff_putamen Neff_pgc2_scz]; traitlist = {'Putamen' 'SCZ'}; % Overlapping samples?
else
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/LDmat_9m_01.mat'); LDmat_01 = LDmat; % Should update these with Dominic's new ones
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/LDmat_9m_08.mat');
%  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/infomat.mat','mafvec','chrnumvec','posvec','snpidlist');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/infomat.mat');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Annot/annomat.mat');
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC_BIP32_9m.mat'); zvec_pgc_bip32 = -norminv(10.^-logpvec_pgc_bip32/2).*sign(zvec_pgc_bip32); Neff_pgc_bip32 = 23291; % Neff should be saved in file
  load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC2_SCZ52.mat'); zvec_pgc2_scz = -norminv(10.^-logpvec_pgc2_scz/2).*sign(zvec_pgc2_scz); Neff_pgc2_scz = 38163; % Neff should be saved in file
%  load('/space/syn03/1/data/cfan/PGC/BPvsCONT.txt.mat'); zvec_bip_noscz = -norminv(10.^-logpvec/2).*sign(zvec); % Need Neff for this! Should also be moved to standard location
%  load('/space/syn03/1/data/cfan/PGC/SCZvsCONT.txt.mat'); zvec_scz_nobip = -norminv(10.^-logpvec/2).*sign(zvec); % Need Neff for this! Should also be moved to standard location
  zmat = cat(2,zvec_pgc_bip32,zvec_pgc2_scz); nsnp = size(zmat,1); Neffvec = [Neff_pgc_bip32 Neff_pgc2_scz]; traitlist = {'BIP' 'SCZ'}; % Overlapping samples
%  zmat = cat(2,zvec_bip_noscz,zvec_scz_nobip); nsnp = size(LDmat,1); traitlist = {'BIP' 'SCZ'}; % Non-overlapping samples
end

ivec_mhc = ((chrnumvec==6)&(posvec>=28477797)&(posvec<=33448354)); % From http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/region.cgi?name=MHC&asm=GRCh37

TLDvec = annomat(:,end);
Hvec = 2*mafvec.*(1-mafvec);
genvar = Hvec;

zmax = -norminv(5e-8/2); % Consider only non-GWAS significant SNPs
%zmax = -norminv(1e-10/2);
zedgevec = linspace(-zmax,zmax,30);

TLDmax = 200;
figure(1); clf; plot(cat(1,-zmat(TLDvec<=TLDmax,2),zmat(TLDvec<=TLDmax,2)),cat(1,-zmat(TLDvec<=TLDmax,1),zmat(TLDvec<=TLDmax,1)),'o'); xlim(2*zmax*[-1 1]); ylim(2*zmax*[-1 1]); axis equal;
h=xlabel(sprintf('z-score in %s',traitlist{2})); set(h,'FontSize',18,'FontWeight','bold');
h=ylabel(sprintf('z-score in %s',traitlist{1})); set(h,'FontSize',18,'FontWeight','bold');


nrep = 1;
gvlist = linspace(0,0.5,6); ngvbins = length(gvlist)-1; 
ldlist = [0 2 5 20 50 100 200]; nldbins = length(ldlist)-1;
pruneflag = true;
excludeMHCflag = false;
params_fit = struct('edgevec_z',zedgevec,'gvlist',gvlist,'ldlist',ldlist,'pruneflag',pruneflag,'excludeMHCflag',excludeMHCflag);
sumstats = cell(1,nrep);
for repi = 1:nrep
  zvec1 = zmat(:,1); % Neff1 = Neffvec_fit(1);
  zvec2 = zmat(:,2); % Neff2 = Neffvec_fit(2);
  ld_gv_freq = NaN(length(ldlist)-1,length(gvlist)-1);
  hcmats = cell(nldbins,ngvbins);
  if excludeMHCflag, zvec1(ivec_mhc) = NaN; zvec2(ivec_mhc) = NaN; end
  if pruneflag, tmp = NaN(size(zvec1)); tmp(isfinite(zvec1+zvec2)) = rand(1,sum(isfinite(zvec1+zvec2))); tmp = GenStats_FastPrune(tmp,LDmat); zvec1(~isfinite(tmp)) = NaN; zvec2(~isfinite(tmp)) = NaN; end
  zmat_tmp = cat(2,zvec1,zvec2);
  defvec = isfinite(sum(zmat_tmp,2));
  for ldi = 1:nldbins
    for gvi = 1:ngvbins
      ivec = find(genvar>=gvlist(gvi)&genvar<=gvlist(gvi+1) & TLDvec>=ldlist(ldi)&TLDvec<ldlist(ldi+1) & defvec);
      gvfreq(gvi) = sum(ivec)/sum(defvec);
      hcmat = hist3(zmat_tmp(ivec,:),'Edges',{zedgevec,zedgevec}); hcmat = hcmat(1:end-1,1:end-1);
      hcmats{ldi,gvi} = hcmat;
    end
  end
  sumstats{repi} = struct('hcmats',{hcmats});
end

% ToDo
%   Make general form:
%     Fit each trait separately
%     Then, fit pleiotropic component 
 
% Fit data to model (Each trait separately) -- appears to have convergence / local minimum issues; try to initialize w. indep. fits for each trait;

options_fit = struct('mixtypevec',[1 2]);
%x0_struct = struct('sig01',1.0,'sig02',1.0,'rho0',0.0,'pivec',[0.01 0.01],'sig1vec',[1 0],'sig2vec',[0 1],'rhovec',[0.0 0.0]);
x0_struct = struct('sig01',0.9727,'sig02',1.1040,'rho0',-0.0088,'pivec',[0.0010 0.0086],'sig1vec',[2.4024 0],'sig2vec',[0 1.8622],'rhovec',[0 0]); % TG / SCZ solution
options_fit.ngm = length(x0_struct.pivec);
x0 = GMM_bivariate_mapparams(x0_struct,options_fit);
resultStruct1 = GMM_bivariate_histfit(sumstats,params_fit,options_fit,x0); x_fit1_struct = GMM_bivariate_mapparams(resultStruct1.x_fit,options_fit);; % Should perhaps try improved optimization methods; multistart?
[cost fitstruct] = GMM_bivariate_histcost(resultStruct1.x_fit,sumstats,params_fit,options_fit,true,zmat,TLDvec,Hvec);

resultStruct1_indep = resultStruct1; x_fit_struct_indep = GMM_bivariate_mapparams(resultStruct1_indep.x_fit,options_fit);


% Fit data to model (including pleiotropic component)

options_fit = struct('mixtypevec',[1 2 3]);
%x0_struct = struct('sig01',x_fit_struct_indep.sig01,'sig02',x_fit_struct_indep.sig02,'rho0',x_fit_struct_indep.rho0,'pivec',[x_fit_struct_indep.pivec 0.0001],'sig1vec',[x_fit_struct_indep.sig1vec 1.0],'sig2vec',[x_fit_struct_indep.sig2vec 1.0],'rhovec',[0 0 0]);
x0_struct = struct('sig01',0.9635,'sig02',1.1004,'rho0',-0.0054,'pivec',[8.2050e-05 0.0072 0.0023],'sig1vec',[9.2956 0 1.6247],'sig2vec',[0 1.6556 2.1565],'rhovec',[0 0 -0.0678]); % Best fit for TG / SCZ
options_fit.ngm = length(x0_struct.pivec)
x0 = GMM_bivariate_mapparams(x0_struct,options_fit);
resultStruct1 = GMM_bivariate_histfit(sumstats,params_fit,options_fit,x0); x_fit1_struct = GMM_bivariate_mapparams(resultStruct1.x_fit,options_fit); % Why is this not updating sig1vec and sigg2vec values exactly equal to 1.0000??
[cost fitstruct] = GMM_bivariate_histcost(resultStruct1.x_fit,sumstats,params_fit,options_fit,true,zmat,TLDvec,Hvec); 

% Should perhaps ultimately allow for each component to be pleiotropic? (all mixtypevec == 3)

% Could also empirically highlight SNPs / hist counts inconsistent with independence (a la copula, but based on pdfs) 


% Fit data to model (pleiotropic component only)

options_fit = struct();
%options_fit = struct('rho0',0);

%x0_struct = struct('sig01',1.0,'sig02',1.0,'rho0',0,'pivec',[1e-3],'sig1vec',[1.0],'sig2vec',[1.0],'rhovec',[0.0]);
%x0_struct = struct('sig01',1.0076,'sig02',1.0249,'rho0',0.0871,'pivec',0.0027,'sig1vec',1.3537,'sig2vec',2.0187,'rhovec',0.7203); % Including MHC
x0_struct = struct('sig01',1.0115,'sig02',1.0452,'rho0',0.2267,'pivec',0.0105,'sig1vec',1.3143,'sig2vec',1.8142,'rhovec',0.6931); % Excluding MHC

options_fit.ngm = length(x0_struct.pivec);
x0 = GMM_bivariate_mapparams(x0_struct,options_fit);
resultStruct1 = GMM_bivariate_histfit(sumstats,params_fit,options_fit,x0); x_fit1_struct = GMM_bivariate_mapparams(resultStruct1.x_fit,options_fit);;

[cost fitstruct pmat_post dmat_post] = GMM_bivariate_histcost(resultStruct1.x_fit,sumstats,params_fit,options_fit,true,zmat,TLDvec,Hvec);

lookuptables = GMM_bivariate_maketables(resultStruct1.x_fit,params_fit,options_fit,1,1);
[pmat_post2 dmat_post2] = GMM_bivariate_lookup(lookuptables,params_fit,options_fit,zmat,TLDvec,Hvec);

%figure(100); clf; plot(-log10(pmat_post(:,1)),dmat_post.^2,'.');

logpvec_post = -log10(pmat_post(:,1));
logpvec_post_pruned = GenStats_FastPrune(logpvec_post,LDmat_01);

sum(logpvec_post_pruned>-log10(0.01))

genenamelist = unique(geneNames(logpvec_post>-log10(0.05)));

d2mat_post = [dmat_post(:,1).^2/Neffvec(1) dmat_post(:,2).^2/Neffvec(2)]; % Adjust effect sizes for sample size
d2vec_post = max(d2mat_post,[],2);
colmat = [d2mat_post(:,1)./d2vec_post d2mat_post(:,2)./d2vec_post zeros(size(d2vec_post))];
ivec = find(logpvec_post>-log10(0.05));
figure(101); clf; scatter(ivec,logpvec_post(ivec),[],colmat(ivec,:)); axis tight;

% ToDo
%   Compute replication probability in finite sample
%   Test in SCZ & BIP substudies
%   Add covariates (replace TLD w. RES)
%   Always include non-pleiotropic mixture components
%   Test polygenic prediction, using a posteriori effect size estimates


% Test in BIP substudies

load('/space/syn03/1/data/GWAS/new_9m_SNPs_processed_summary/GWAS_Data/PGC2_BIP32_multisites_9m.mat'); logpmat_bip = logpmat; zmat_bip = -norminv(10.^-logpmat/2).*sign(betamat); Neffvec_bip = 2./(1./Nvec_A+1./Nvec_U); Neff_bip = sum(Neffvec_bip);

Neffvec = Neffvec_bip; 
zmat = zmat_bip; logpmat = logpmat_bip;
nsnp = size(zmat,1); nsamp = size(zmat,2); studyindlist_disc = [1:nsamp];
Neff = sum(Neffvec);

% Make stratified z2z plots

pthresh = 0.05;
zvals_disc = linspace(-10.0,10.0,201);

dflist = [2 3 4 5]/10; ndf = length(dflist); nr = round(sqrt(ndf)); nc = ceil(ndf/nr);
pruneflag = true;
nrep = 100;
datastructmat = cell(ndf,nrep);
hist_logpvals = linspace(0,50,501);
figure(101); clf; figure(102); clf; figure(103); clf;
z2meanmat = NaN(ndf,nrep); lamgcmat = NaN(ndf,nrep); Neffmat = NaN(ndf,nrep);
time0 = now;
for dfi = 1:ndf
  dfrac = dflist(dfi);
  dmat = false(nrep,nsamp); rmat = false(nrep,nsamp);
  for repi = 1:nrep
    done = false;
    while ~done
      ivec = randperm(nsamp);
      mi = floor(dfrac*nsamp+rand);
      dmat(repi,:) = false; dmat(repi,ivec(1:mi)) = true; rmat(repi,:) = ~dmat(repi,:);
      done = abs((sum(Neffvec(dmat(repi,:)))/sum(Neffvec))-dfrac) < 0.10;
    end
  end
  zmat_repl_sum = 0; z2mat_repl_sum = 0; repratemat_sum = 0; cntmat_sum = 0;
  time1 = now;
  for repi = 1:nrep
    Neffvec_disc = rowvec(double(Neffvec(dmat(repi,:)))); Neffvec_repl = rowvec(double(Neffvec(rmat(repi,:))));
    Neff_disc = sum(Neffvec_disc); Neff_repl = sum(Neffvec_repl);
    zvec_disc = sum(zmat(:,dmat(repi,:)).*repmat(sqrt(Neffvec_disc),[size(zmat,1) 1]),2)/sqrt(sum(Neffvec_disc));
    zvec_repl = sum(zmat(:,rmat(repi,:)).*repmat(sqrt(Neffvec_repl),[size(zmat,1) 1]),2)/sqrt(sum(Neffvec_repl));
    z2vec_repl = zvec_repl.^2;
    pvec_disc = 2*normcdf(-abs(zvec_disc),0,1); logphistcount = hist(-log10(pvec_disc),hist_logpvals);
    pvec_repl = normcdf(-sign(zvec_disc).*zvec_repl,0,1); repvec = double(pvec_repl<pthresh);
    if pruneflag, tmp = NaN(size(zvec_disc)); tmp(isfinite(zvec_disc+zvec_repl)) = rand(1,sum(isfinite(zvec_disc+zvec_repl))); tmp = GenStats_FastPrune(tmp,LDmat); zvec_disc(~isfinite(tmp)) = NaN; prunevec = isfinite(zvec_disc); else prunevec = []; end
    pvec_disc = 2*normcdf(-abs(zvec_disc),0,1); logphistcount_pruned = hist(-log10(pvec_disc),hist_logpvals);
keyboard
% Compute predicted zvec_repl, given z-score for BIP
% Compute lookup table, given Neff_disc/Neff
    [dummy dummy2 pmat_post_tmp dmat_post_tmp] = GMM_bivariate_histcost(resultStruct1.x_fit,sumstats,options_fit,[zvec_disc zvec_trait2],TLDvec,Hvec);
    for bini = 1:size(iimat,2)
      defvec = isfinite(zvec_disc+zvec_repl) & iimat(:,bini);
      [zvals_repl zvals_repl_cnt] = histmean([-zvec_disc(defvec); zvec_disc(defvec)],[-zvec_repl(defvec); zvec_repl(defvec)],zvals_disc);
      [z2vals_repl z2vals_repl_cnt] = histmean([-zvec_disc(defvec); zvec_disc(defvec)],[z2vec_repl(defvec); z2vec_repl(defvec)],zvals_disc);
      [repvals repvals_cnt] = histmean([-zvec_disc(defvec); zvec_disc(defvec)],[repvec(defvec); repvec(defvec)],zvals_disc);
      zmat_repl(:,bini) = zvals_repl;
      z2mat_repl(:,bini) = z2vals_repl;
      repratemat(:,bini) = repvals;
      cntmat(:,bini) = zvals_repl_cnt;
    end
    zmat_repl_sum = zmat_repl_sum + zmat_repl.*cntmat;
    z2mat_repl_sum = z2mat_repl_sum + z2mat_repl.*cntmat;
    repratemat_sum = repratemat_sum + repratemat.*cntmat;
    cntmat_sum = cntmat_sum + cntmat;
    datastructmat{dfi,repi} = struct('Neff_disc',Neff_disc,'Neff_repl',Neff_repl,'zmat_repl',zmat_repl,'cntmat',cntmat,'z2mat_repl',z2mat_repl,'repratemat',repratemat);
    figure(101); subplot(nr,nc,dfi); plot(zvals_disc,zmat_repl_sum(:,1:nbins)./cntmat_sum(:,1:nbins),zvals_disc,zmat_repl_sum(:,end)./cntmat_sum(:,end),'k-','LineWidt',2);
    h=xlabel('z_d_i_s_c'); set(h,'FontSize',22); h=ylabel('E(z_r_e_p_l | z_d_i_s_c)'); set(h,'FontSize',22); xlim([-6 6]); %ylim([-2 2]);
    figure(102); subplot(nr,nc,dfi); plot(zvals_disc,z2mat_repl_sum(:,1:nbins)./cntmat_sum(:,1:nbins),zvals_disc,z2mat_repl_sum(:,end)./cntmat_sum(:,end),'k-','LineWidt',2);
    h=xlabel('z_d_i_s_c'); set(h,'FontSize',22); h=ylabel('E(z_r_e_p_l^2 | z_d_i_s_c)'); set(h,'FontSize',22); xlim([-6 6]); %ylim([0 6]);
    figure(103); subplot(nr,nc,dfi); plot(zvals_disc,repratemat_sum(:,1:nbins)./cntmat_sum(:,1:nbins),zvals_disc,repratemat_sum(:,end)./cntmat_sum(:,end),'k-','LineWidt',2);
    h=xlabel('z_d_i_s_c'); set(h,'FontSize',22); h=ylabel('Expected Replication Rate'); set(h,'FontSize',22); ylim([0 1.1]); xlim([-6 6]);
    drawnow;
    fprintf(1,'rep %d of %d (now:%s done:%s)\n',repi,nrep,datestr(now,'HH:MM:SS'),datestr(time1+(now-time1)/repi*nrep,'HH:MM:SS'));
  end
  fprintf(1,'dfi %d of %d (now:%s done:%s)\n',dfi,ndf,datestr(now,'HH:MM:SS'),datestr(time0+(now-time0)/dfi*ndf,'HH:MM:SS'));
end

