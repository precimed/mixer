cd 'C:\Users\Oleksandr\Documents\GitHub\GMM'
addpath('H:\Dropbox\analysis\2016_09_September_22_Bivariate');

%folder = 'H:\Dropbox\analysis\2017_07_July_31_TestingGMM_vs_LDSR\bgmg_real_fromABEL_mat';
folder = 'H:\Dropbox\analysis\2017_07_July_31_TestingGMM_vs_LDSR\loci_discovery';
files = dir(fullfile(folder, '*.mat'));

if 0
data = {}
for i=1:length(files)
data{i} = load(fullfile(folder, files(i).name))
end
end


for i=1:length(files)
result = data{i}.result;
data1  = load(fullfile(fullfile(ldsr_path, 'MATLAB_Data', result.trait1)));
data2  = load(fullfile(fullfile(ldsr_path, 'MATLAB_Data', result.trait2)));

o = options; o.calculate_global_fdr=false; o.calculate_beta_hat=true; 
idx=1:10000;
[~, r] = GMM2_bivariate_cost(result.bivariate.params, [data1.zvec(idx) data2.zvec(idx)], Hvec(idx), [data1.nvec(idx) data2.nvec(idx)], w_ld(idx), struct('sum_r2', ref_ld.sum_r2(idx), 'sum_r4', ref_ld.sum_r4(idx)), [], o);
p=result.bivariate.params;

zscore = r.beta_hat_mean ./ sqrt(r.beta_hat_var);

a11 = mean(r.params_per_snp.a11);
a12 = mean(r.params_per_snp.a12);
a22 = mean(r.params_per_snp.a22);

mixture.pivec = [1-sum(p.pivec) p.pivec];
mixture.pivec(mixture.pivec < 1e-7) = 0;
if result.bivariate.params.uncertainty.pi3.pvalue > (0.05/15), mixture.pivec(4) = 0; end
mixture.sigma{1} = [p.sigma0(1).^2, p.sigma0(1)*p.sigma0(2)*p.rho0; p.sigma0(1)*p.sigma0(2)*p.rho0, p.sigma0(2).^2];
mixture.sigma{2} = [a11, 0; 0, 0];
mixture.sigma{3} = [0, 0; 0, a22];
mixture.sigma{4} = [a11, a12; a12, a22];

%subplot(3,5,i);
plot_snps([data1.zvec, data2.zvec], 'how', 'density', 'ticks', true)
plot_gmm(mixture, 'how', 'density')
titlestr = sprintf('%s vs %s', result.trait1, result.trait2);
titlestr = strrep(titlestr, '_noMHC', '');
titlestr = strrep(titlestr, '_qc', '');
titlestr = strrep(titlestr, '_lift', '');
titlestr = strrep(titlestr, '_', ' ');
title(titlestr)

end


if 0
    p=result.univariate{1}.params;
    o=options; o.calculate_beta_hat = true;
    [~, r]=GMM2_univariate_cost(p, data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, [], o);
    r.zvec = r.beta_hat_mean ./ r.beta_hat_se;

    binlimits = [-2 2];
    histogram(normrnd(0,1, size(zscore, 1), 1), 'BinLimits', binlimits)
    histogram(zscore(:, 1), 'BinLimits', binlimits)

    for i=1:2
        %if i==1, z = data1.zvec; col='g';  else z = r.zvec; col='b'; end;
        z = data1.zvec;
        qqlim=16;
        z=z(~isnan(z));
        
        sig0 = sqrt(nanmedian(zscore.^2) ./ chi2inv(0.5,1));
        sig0 = sqrt(median(z.^2) ./ chi2inv(0.5,1));

        if i==1, z = z/sig0; col='g';  else col='b'; end;

        logpvec = sort(-log10(normcdf(z, 0, 1)));
        a=1:(-1/length(logpvec)):0; a=a(1:length(logpvec)); a=-log10(a);
        hold on
        plot(a, logpvec,col,'LineWidth',2);
        plot([0 qqlim],[0 qqlim], 'k');
        xlim([0 6]); ylim([0 7]);
    end

end