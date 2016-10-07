
% Note that this needs updated TLDvec, computed by interpolating the 1kG TLDvec by chromosomal location -- use SCZ & BIP for now 

load('~/GWAS_Data/immune_sumstats.mat'); % Where did this come from? Should make sure each rep of sumstats_dfrac for different conditions have non-overlapping controls 

load('~/GWAS_Data/AutoimmuneMatFiles/maf.mat'); % Why the hell are MHC_start and MHC_end in here??
%load('~/GWAS_Annot/Autoimmune_annotmat.mat'); colnames_annot = annonames;

load('~/GWAS_Annot/Autoimmune_LDmat.mat');
LDmat = CS; clear CS; % Fix Yunpeng's naming scheme

% Should compute TLD for immunochip
% Double-check Yunpeng's LDmat

% Read in temporary data saved on ip65
tmp = load('~/tmp/AutoImmune/Control_covmat.mat');
load('~/tmp/AutoImmune/Control_covmat_snap.mat');
[C,IA] = unique([Icat Jcat],'rows');
Icat = Icat(IA);
Jcat = Jcat(IA);
Rcat = Rcat(IA);
nsnp = length(muvec);
fvec = colvec(muvec)/2; mafvec = min(fvec,1-fvec); Hvec = 2*mafvec.*(1-mafvec);
LDmat_r2 = (sparse(Icat,Jcat,double(Rcat),nsnp,nsnp)).^2;
TLDvec = sum(LDmat_r2,2); % Should perhaps generate this by interpolating 1kG TLD by chromosomal position (from same HG build)
LDmat = (LDmat_r2>=0.1); % Create binary LDmat for pruning
genvar = Hvec;

MHC_start = 52451; MHC_end = 59450;
ivec_MHC = false(1,nsnp); ivec_MHC(MHC_start:MHC_end) = true;

condlist = {'CD' 'UC' 'PSC' 'PS' 'AS'}; ncond = length(condlist);

seed = 1; iter = 1;
load(sprintf('~/immune_sumstats_%d_%d.mat',seed,iter));

condi1 = 1; condi2 = 2; traitlist = condlist([condi1,condi2]);
%condi1 = 1; condi2 = 5; traitlist = condlist([condi1,condi2]);
%condi1 = 1; condi2 = 3; traitlist = condlist([condi1,condi2]);
%condi1 = 2; condi2 = 3; traitlist = condlist([condi1,condi2]);

zmax = -norminv(5e-8/2);
zedgevec = linspace(-zmax,zmax,15);

nrep = size(sumstats_dfrac,2);
ndf = size(sumstats_dfrac,3);
gvlist = linspace(0,0.5,6); ngvbins = length(gvlist)-1; 
ldlist = [0 2 5 100]; nldbins = length(ldlist)-1;
pruneflag = 1;
excludeMHCflag = 1;
sumstats = cell(ndf,nrep);
for dfi = 1:ndf
  for repi = 1:nrep
    zvec1 = sumstats_dfrac{condi1,repi,dfi}.zvec_disc; Neff1 = sumstats_dfrac{condi1,repi,dfi}.Neff_disc;
    zvec2 = sumstats_dfrac{condi2,repi,dfi}.zvec_disc; Neff2 = sumstats_dfrac{condi2,repi,dfi}.Neff_disc;
    ld_gv_freq = NaN(length(ldlist)-1,length(gvlist)-1);
    hcmats = cell(nldbins,ngvbins);
    if excludeMHCflag, zvec1(ivec_MHC) = NaN; zvec2(ivec_MHC) = NaN; end
    if pruneflag, tmp = NaN(size(zvec1)); tmp(isfinite(zvec1+zvec2)) = rand(1,sum(isfinite(zvec1+zvec2))); tmp = GenStats_FastPrune(tmp,LDmat); zvec1(~isfinite(tmp)) = NaN; zvec2(~isfinite(tmp)) = NaN; end
    zmat = cat(2,zvec1,zvec2);
    defvec = isfinite(sum(zmat,2));
    for ldi = 1:nldbins
      for gvi = 1:ngvbins
        ivec = find(genvar>=gvlist(gvi)&genvar<=gvlist(gvi+1) & TLDvec>=ldlist(ldi)&TLDvec<ldlist(ldi+1) & defvec);
        gvfreq(gvi) = sum(ivec)/sum(defvec);
        hcmat = hist3(zmat(ivec,:),'Edges',{zedgevec,zedgevec}); hcmat = hcmat(1:end-1,1:end-1);
        hcmats{ldi,gvi} = hcmat;
      end
    end
    sumstats{dfi,repi} = struct('zmat',zmat,'genvar',genvar,'Neff1',Neff1,'Neff2',Neff2,'hcmats',{hcmats},'traitlist',{traitlist},'edgevec',zedgevec,'gvlist',gvlist,'ldlist',ldlist,'ld_gv_freq',ld_gv_freq);
  end
end

figure(10); clf; cnt = 0;
for ldi = 1:nldbins
  for gvi = 1:ngvbins
     hcmat = 0;
     for repi = 1:nrep
       hcmat = hcmat + sumstats{end,1}.hcmats{ldi,gvi};
     end
     hcmat = (hcmat+hcmat')/2; % Make symmetric
     cnt = cnt+1; subplot(nldbins,ngvbins,cnt); 
     imagesc(log10(hcmat./max(hcmat(:))),[-5 0]); colormap(hot); axis equal; axis xy; axis tight;
  end
end

hcmat = sumstats{1,1}.hcmats{1,end};
figure(1); clf; imagesc(log10(hcmat./max(hcmat(:))),[-5 0]); colormap(hot); axis equal; axis xy; axis tight; % set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels);

hcmat = sumstats{end,1}.hcmats{end-2,end};
figure(2); clf; imagesc(log10(hcmat./max(hcmat(:))),[-5 0]); colormap(hot); axis equal; axis xy; axis tight; % set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks,'YTickLabel',ticklabels);

% Implement fitting with mixture of Gaussians
%   Start with single non-null, then allow for two or three -- compare fits, and generalization performance

% ToDos
%   Compute TLDvec by interpolating TLDvec from 1kG by chr pos (ignore TLD for now)
%   Re-generate immune_sumstats.mat / sumstats_dfrac (originally from Yunpeng?)
%   Should generate sumstats_dfrac with overlapping controls across conditions, to test if method works in this case

% Run on SCZ / BIP data until valid TLDvec has been computed for ImmunoChip
% Generate summary stats for different random subsets of SCZ and BIP substudies, put into same workflow as AutoImmune, above
% Use full set of overlapping and non-overlapping SCZ & BIP summary stats (repi=1 and dfi=1) 


