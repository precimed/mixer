folder = 'C:\Users\oleksanf\Dropbox\analysis\2016_09_September_19_LDScoreRegression\1m';

if 0
    addpath('BVNcdf');
    SCZ_vs_CONT=load(fullfile(folder,  'sumstats\PGC_BIP_2016_SCZvsCONT.mat')); 
    BIP_vs_CONT=load(fullfile(folder,  'sumstats\PGC_BIP_2016_BIPvsCONT.mat'));
    SCZ        =load(fullfile(folder,  'sumstats\PGC_SCZ_2014.mat'));
    BIP        =load(fullfile(folder,  'sumstats\PGC_BIP_2016.mat'));
    TG         =load(fullfile(folder,  'sumstats\LIPIDS_TG_2013.mat'));

    eur_w_ld   =load(fullfile(folder,  'eur_w_ld.mat'));annomat=eur_w_ld.annomat;
                load(fullfile(folder,  'infomat.mat'));
                load(fullfile(folder,  'ldmat_1m_p1.mat'));
                load(fullfile(folder,  'ldmat_1m_p8.mat'));
                
    SCZ_vs_CONT.logpvec = -log10(2*normcdf(-abs(SCZ_vs_CONT.zvec)));
    BIP_vs_CONT.logpvec = -log10(2*normcdf(-abs(BIP_vs_CONT.zvec)));
    SCZ.logpvec = -log10(2*normcdf(-abs(SCZ.zvec)));
    BIP.logpvec = -log10(2*normcdf(-abs(BIP.zvec)));
    TG.logpvec = -log10(2*normcdf(-abs(TG.zvec)));


    LDmat = LDr2_p1sparse;
end


p1 = SCZ; p2 = TG;
%p1 = SCZ_vs_CONT; p2 = BIP_vs_CONT;
zmat = cat(2,p1.zvec,p2.zvec);
diary 'SCZ_vs_TG'
gmm_SCZ_TG = GMM_bivariate_hist_fit_many(chrnumvec, posvec, annomat, mafvec, zmat, LDmat);
diary off