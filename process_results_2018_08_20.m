% =========================== SIMULATIONS ===========================
folder = 'H:\work\simu_ugmg120_run12a';
files = dir([folder, '\*run12a.fit.test.mat']);
h2_vec={'0.1', '0.4', '0.7'};
pi_vec = {'1e-05', '0.0001', '0.001', '0.01'};

failed_count = 0; pi_thresh = 2.6e-02;


results.pi_mat = cell(length(h2_vec), length(pi_vec));
results.sig2_beta_mat = cell(length(h2_vec), length(pi_vec));
results.sig2_zero_mat = cell(length(h2_vec), length(pi_vec));
results.h2_mat = cell(length(h2_vec), length(pi_vec));
figure(123);
for h2_index = 1:length(h2_vec)
    for pi_index = 1:length(pi_vec)
        results.last_pi_vec = [];
        results.last_sig2_beta = [];
        results.last_sig2_zero = [];
        results.last_h2_vec = [];
        subplot(3,4,4*(h2_index-1)+pi_index); hold on;
        for repi=1:1  %10
            %continue
            fname = sprintf('simu_h2=%s_pi1u=%s_rep=%i.run12a.fit.test.mat', h2_vec{h2_index}, pi_vec{pi_index}, repi);
            fprintf('%s...\n', fname);
            r = load(fullfile(folder, fname));
            params = r.result.params.univariate{1};
            results.last_pi_vec(end+1, 1) = params.pi_vec;
            results.last_sig2_beta(end+1, 1) = params.sig2_beta;
            results.last_sig2_zero(end+1, 1) = params.sig2_zero;
            results.last_h2_vec(end+1, 1) = params.sig2_beta * params.pi_vec * r.result.options.total_het;
            
            pd = r.result.univariate{1}.qq_plot_data;
            
            if params.pi_vec > pi_thresh || (h2_index == 2 && pi_index == 1 && repi == 9), continue; end
            plot(pd.data_logpvec, pd.hv_logp, '-b', pd.model_logpvec, pd.hv_logp, '-r', 0:7, 0:7, '--k')
            xlim([0 7]); ylim('auto');
        end
        ylim auto
        legend('data', 'model', 'null', 'Location', 'NorthWest');
        results.pi_mat{h2_index, pi_index} = results.last_pi_vec;
        results.sig2_beta_mat{h2_index, pi_index} = results.last_sig2_beta;
        results.sig2_zero_mat{h2_index, pi_index} = results.last_sig2_zero;
        results.h2_mat{h2_index, pi_index} = results.last_h2_vec;
    end
end

figure(1); plot_index = 1;
for h2_index = 1:length(h2_vec)
    for pi_index = 1:length(pi_vec)
        last_pi_vec = results.pi_mat{h2_index, pi_index};
        last_pi_vec(last_pi_vec>=pi_thresh) = nan;
        if h2_index == 2 && pi_index == 1, last_pi_vec(9) = nan; end;
        failed_count = failed_count + sum(~isfinite(last_pi_vec));
        
        values_of_interest = results.h2_mat{h2_index, pi_index}; param_name = 'h2';
        values_of_interest(~isfinite(last_pi_vec)) = nan; 
        
        fprintf('h2=%s, pi=%s, median(%s)=%.3f 90%%ci=[%.3f, %.3f]\n', h2_vec{h2_index}, pi_vec{pi_index}, param_name, nanmedian(values_of_interest), min(values_of_interest), max(values_of_interest));
        subplot(3,4,plot_index); histogram(values_of_interest, 10);plot_index=plot_index+1;
        
    end
end
set(123,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(123, sprintf('synthetic.test.hapmap3_hardprune@p1_run12a.pdf'), '-dpdf')

fprintf('done, %i out of 120 runs failed\n', failed_count);
% ============================ REAL DATA ============================
if 0
folder = 'H:\work\SUMSTAT_11M_result';
files1 = dir([folder, '\*run7.*test.mat']);
files2 = dir([folder, '\*run8.*fit.mat']);
%files = [files1; files2]; indices=reshape(1:length(files), [12 2])';indices=indices(:);
files = files1; indices=1:length(files); indices=indices';
fprintf('trait\tsnps\tpi\t(se)\t\tsig2beta\t(se)\t\tsig2zero\t(se)\t\th2\t(se)\n');

sample_sizes.CTG_COG_2018 = 269867	;
sample_sizes.ENIGMA_PUT_2016_noPGC_noGC=11598	;
sample_sizes.GIANT_HEIGHT_2018_UKB = 709706	;
sample_sizes.IIBDGC_CD_2017_lift = 	4/(1/12194	+1/34915);
sample_sizes.IIBDGC_UC_2017_lift = 	4/(1/12366	+1/34915);
sample_sizes.LIPIDS_HDL_2013=187167	;
sample_sizes.LIPIDS_LDL_2013 = 173082;
sample_sizes.PGC_AD_2018_biorxiv_feb=4/(1/71880+1/383378);
sample_sizes.PGC_BIP_2016_qc = 4/(1/20352+1/31358);
sample_sizes.PGC_MDD_2018 = 4/(1/83836+1/177093);
sample_sizes.PGC_SCZ_2014 = 4/(1/35476+1/46839);
sample_sizes.SSGAC_EDU_2016 = 293723;

nCurrent=[]; sCurrent=[];
indices=[4, 6, 5, 2, 7, 3, 8, 11, 9, 1, 12, 10]';
for i=1:length(indices);
    if files(indices(i)).bytes <= 0, continue; end;
    file_name=files(indices(i)).name;
    r = load(fullfile(folder, file_name));
    params = r.result.params.univariate{1};
    if isfield('ci', r.result.univariate{1})
        params_ci = r.result.univariate{1}.ci;
    else
        params_ci.pi_vec.se=nan;
        params_ci.sig2_beta.se=nan;
        params_ci.sig2_zero.se=nan;
        params_ci.h2.se=nan;
    end
    h2 = r.result.options.total_het * params.pi_vec * params.sig2_beta;
    fit_hapmap = findstr(file_name, 'run5');
    fit_all = findstr(file_name, 'run6');
    if fit_hapmap, snps='hapmap3'; fi=2; else snps='allSNPs'; fi=1; end;
    fprintf('%s\t%s\t%.3e\t%.3e\t\t%.3e\t%.3e\t\t%.3f\t%.3f\t\t%.3f\t%.3e\n', file_name(1:end-4), snps, params.pi_vec, params_ci.pi_vec.se, params.sig2_beta, params_ci.sig2_beta.se, params.sig2_zero, params_ci.sig2_zero.se,  h2, params_ci.h2.se)
    if length(indices)~=12, continue; end;

    subplot_index = floor((i+1)/2); 
    subplot_index = i;
    title_str = file_name(1:end-4); title_str(title_str=='_') = '-';
        
    if isfield(r.result.univariate{1}, 'qq_plot_data')
        pd = r.result.univariate{1}.qq_plot_data;

        figure(fi); subplot(4,3, subplot_index);
        plot(pd.data_logpvec, pd.hv_logp, '-', pd.model_logpvec, pd.hv_logp, '--', 0:7, 0:7, '--k')
        title(title_str);
        xlim([0, 7]);% ylim([0, 20]);
    end
    
    if 0 && isfield(r.result.univariate{1}, 'loglike_adj_trajectory')
        lat=r.result.univariate{1}.loglike_adj_trajectory;
    
        figure(3); if i==1, clf; end; hold on;  subplot(4,3, subplot_index);
        plot(lat.pivec, lat.cost, '.-', r.result.univariate{1}.params.pi_vec, min(lat.cost), '*');
        title(title_str);

        x = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(r.result.univariate{1}.params);
        figure(4); if i==1, clf; end; hold on; subplot(4,3, subplot_index);
        plot(lat.xgrid, lat.cost, '.-', x(end), min(lat.cost), '*');
        title(title_str);
    end
    
    if isfield(r.result.univariate{1}, 'power_plot_data')
        ppd = r.result.univariate{1}.power_plot_data;
        figure(5); if i==1, clf; end; hold on; 
        x=findstr(title_str, '.');x=x(1);
        if i<=7, plot_marker='-', else plot_marker='--'; end;

        nCurrent(end+1, 1)= sample_sizes.(file_name(1:x-1));
        sCurrent(end+1, 1) = interp1(log10(ppd.power_nvec), ppd.power_svec, log10(nCurrent(end)));
        
        plot(log10(ppd.power_nvec), ppd.power_svec, plot_marker, 'DisplayName',sprintf('%s (%.1f%%)', title_str(1:x-1), 100*sCurrent(end)), 'LineWidth',2);
        show_legend_at_the_end = 1;
    end
end

if exist('show_legend_at_the_end', 'var') && show_legend_at_the_end, 
    figure(5); legend('off'); ax=legend('show', 'location', 'NorthWest'); ax.FontSize=15;
    ax=gca;
    ax.ColorOrderIndex = 1;
    for i=1:length(nCurrent)
        plot(log10(nCurrent(i)), sCurrent(i), '*', 'LineWidth',2)
    end
end;
set(5,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(5, sprintf('run7.test.power.pdf'), '-dpdf')




set(1,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(1, sprintf('run6.test.allSNPs.pdf'), '-dpdf')

set(3,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(3, sprintf('run6.test.loglike.pdf'), '-dpdf')

set(4,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(4, sprintf('run6.test.loglike_xgrid.pdf'), '-dpdf')

end

if 0
    point_estimates.pi_vec = []
    point_estimates.sig2_beta = [];
    point_estimates.sig2_zero = [];
    point_estimates.h2 = [];
    for repi=1:10
        r=load(['H:\work\simu_ugmg120_run14\', sprintf('simu_h2=0.4_pi1u=0.001_rep=%i.run15.fit.mat', repi)]);
        if all(isfinite(r.result.univariate{1}.ci_hess(:))) && isfinite(r.result.univariate{1}.ci.sig2_zero.mean)
            ci = r.result.univariate{1}.ci;
            fname='h2'; fprintf('%i %s: %.3f, %.3f, %.3f\n', repi, fname, ci.(fname).point_estimate, ci.(fname).mean, ci.(fname).se); point_estimates.(fname)(end+1, 1) = ci.(fname).point_estimate;
            fname='pi_vec'; fprintf('%i %s: %.3e, %.3e, %.3e\n', repi, fname, ci.(fname).point_estimate, ci.(fname).mean, ci.(fname).se); point_estimates.(fname)(end+1, 1) = ci.(fname).point_estimate;
            fname='sig2_beta'; fprintf('%i %s: %.3e, %.3e, %.3e\n', repi, fname, ci.(fname).point_estimate, ci.(fname).mean, ci.(fname).se); point_estimates.(fname)(end+1, 1) = ci.(fname).point_estimate;
            fname='sig2_zero'; fprintf('%i %s: %.3f, %.3f, %.3f\n', repi, fname, ci.(fname).point_estimate, ci.(fname).mean, ci.(fname).se); point_estimates.(fname)(end+1, 1) = ci.(fname).point_estimate;
            fprintf('\n');
        end
    end
    fnames = {'pi_vec', 'sig2_beta', 'sig2_zero', 'h2'};
    for i=1:4
        fprintf('%s - %.3e\n', fnames{i}, std(point_estimates.(fnames{i})));
    end
    
    result = r.result;
sig2_zero_ci = zeros(10000, 1);
pi_vec_ci = zeros(10000, 1);
sig2_beta_ci = zeros(10000, 1);
for i=1:10000,
     sig2_zero_ci(i) = ci_params{i}.sig2_zero;
     pi_vec_ci(i) = ci_params{i}.pi_vec;
     sig2_beta_ci(i) = ci_params{i}.sig2_beta;
end
figure(43); 
subplot(2,2,1); f=pi_vec_ci; histogram(f); title(sprintf('pi, %.3e ( %.3e )', mean(f), std(f)));
subplot(2,2,2); f=sig2_beta_ci; histogram(f); title(sprintf('sig2-beta, %.3e ( %.3e )', mean(f), std(f)));
subplot(2,2,3); f=sig2_zero_ci; histogram(f); title(sprintf('sig2-zero, %.3f ( %.3f )', mean(f), std(f)));
subplot(2,2,4); f=result.options.total_het * sig2_beta_ci .* pi_vec_ci; histogram(f); title(sprintf('h2, %.3f ( %.3f )', mean(f), std(f)));
figure(44);
subplot(2,2,1); pcolor(hist3([pi_vec_ci, sig2_beta_ci], [30 30])); title('pi vs sig2-beta');
subplot(2,2,2); pcolor(hist3([pi_vec_ci, sig2_zero_ci], [30 30])); title('pi vs sig2-zero');
subplot(2,2,3); pcolor(hist3([sig2_beta_ci, sig2_zero_ci], [30 30])); title('sig2-beta vs sig2-zero');
%set(total_fig,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
 print(43, sprintf('Uncertainty histograms'), '-dpdf')
 print(44, sprintf('Uncertainty heatmaps'), '-dpdf')
end