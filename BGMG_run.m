% This script can be run univariate or bivariate mixture analysis from command line like this:
%
%   matlab -nodisplay -nosplash -nodesktop -r "trait1='PGC_SCZ_2014'; trait2='LIPIDS_TG_2013'; BGMG_run; exit;"
%
% You may use this script with one trait (univariate analysis) or with two
% traits (bivariate analysis).
% The results are saved into a text file as well as into .mat file.
% Currently you may need to modify this file to pass additional
% parameters, later all parameters will be exposed as a text config file.
%
% You may also set trait1 to a file containing Nx2 zvec and nvec,
% this will trigger bivariate analysis even though trait2 is not set.
% This is handfull, for example, in simulations when two GWAS results
% are saved in one file.
%
% For information about various parameters look into BGMG_fit.m.
%try
fprintf('trait1: %s\n', trait1);
if exist('trait2', 'var'), fprintf('trait2: %s\n', trait2); end;

%DO_RANDPRUNE_p8 = 1;
%DO_RANDPRUNE_p2 = 1;

[~, trait1_name, ~] = fileparts(trait1);
if exist('trait2', 'var'), [~, trait2_name, ~] = fileparts(trait2); end;

if ~exist('data_path', 'var')
t = 'C:\Users\Oleksandr\Documents\GitHub\BGMG\reference_data';                       if exist(t, 'dir'), data_path = t; end;
%t = 'C:\Users\Oleksandr\Documents\GitHub\BGMG\reference_data_10m';                       if exist(t, 'dir'), data_path = t; end;
t = '/work/users/oleksanf/NORSTORE/GWAS_SUMSTAT/LDSR/MATLAB_Data';     if exist(t, 'dir'), data_path = t; end;
t = '/space/syn03/1/data/GWAS_SUMSTAT/LDSR/MATLAB_Data';               if exist(t, 'dir'), data_path = t; end;
end
if ~exist('data_path', 'var'), error('Unable to locate data folder'); end;

addpath('DERIVESTsuite');

%if ~exist('reference_data', 'var'), reference_data = 'reference_data_10m'; end;
if ~exist('reference_data', 'var'), reference_data = 'reference_data'; end;
fprintf('Loading reference data from %s...\n', reference_data);
load(fullfile(reference_data, 'mafvec.mat'),'mafvec');
load(fullfile(reference_data, 'chrnumvec.mat'),'chrnumvec');
load(fullfile(reference_data, 'posvec.mat'),'posvec');
ref_ld2 = load(fullfile(reference_data, 'ref_l2.mat'));
biased_ref_ld4 = load(fullfile(reference_data, 'biased_ref_l4.mat'));
biased_ref_ld2 = load(fullfile(reference_data, 'biased_ref_l2.mat'));
w_ld2 = load(fullfile(reference_data, 'w_ld.mat'));

ref_ld = struct('sum_r2', biased_ref_ld2.annomat, 'chi_r4', biased_ref_ld4.annomat ./ biased_ref_ld2.annomat);
Hvec = 2*mafvec .* (1-mafvec);  w_ld  = w_ld2.annomat;

% Remember to exclude MHC, also for summary stats called ending with
% no_MHC.mat. no_MHC suffix imply exclusion of 26-34, as done by LD score
% regresison. In some cases this seems to be not enough. 25-35 is the 
% standard from Ricopilli pipeline.
mhc = (chrnumvec==6) & (posvec > 25e6) & (posvec < 35e6);

if exist('DO_RANDPRUNE_p8', 'var') && DO_RANDPRUNE_p8 && ~exist('LDmat', 'var')
    load(fullfile(reference_data, 'binary_ld_p8.mat'))  % => LDmat
end
if exist('DO_RANDPRUNE_p2', 'var') && DO_RANDPRUNE_p2 && ~exist('LDmat', 'var')
    load(fullfile(reference_data, 'binary_ld_p2.mat'))  % => LDmat
end

% mafvec_all = load('mafvec_all_snps_1000G_Phase3_frq.mat')
% options.total_het = 2*sum(mafvec_all.mafvec .* (1-mafvec_all.mafvec))
options.total_het = 2 * 1037117.5140529468;  % Total heterozigosity across all SNPs
options.verbose = true;

data1  = load(fullfile(data_path, trait1)); data1.zvec(mhc, :) = nan; if isfield(data1, 'infovec'), data1.zvec(data1.infovec < 0.9, :) = nan; defvec = isfinite(data1.zvec); end;
if exist('data1_neff', 'var'), data1.nvec = ones(size(data1.zvec)) * data1_neff; end;
if exist('trait2', 'var'),
    data2 = load(fullfile(data_path, trait2)); data2.zvec(mhc, :) = nan;  if isfield(data2, 'infovec'), data2.zvec(data2.infovec < 0.9, :) = nan; defvec = isfinite(data1.zvec + data2.zvec); end;
    if exist('data2_neff', 'var'), data2.nvec = ones(size(data2.zvec)) * data2_neff; end;
end

if exist('LDmat', 'var')
    % do random prunin
    nprune = 100;
    fprintf('Generate random-pruning indiced weights for %i SNPs, %i iterations\n', sum(defvec), nprune);
    hits = zeros(size(defvec));
    defidx = find(defvec);
    for prune_index = 1:nprune
        si = defidx(randperm(length(defidx)));
        prunevec = false(size(defvec));
        for i = 1:length(si)
          if ~prunevec(si(i))
             prunevec(find(LDmat(:,si(i)))) = true; 
             prunevec(si(i)) = false;
          end
        end
        hits(~prunevec) = hits(~prunevec) + 1;
        fprintf('.');
    end
    w_ld = hits;
    w_ld(hits==0) = nan;
end

if exist('trait2', 'var'),
    result = BGMG_fit([data1.zvec, data2.zvec], Hvec, [data1.nvec, data2.nvec], w_ld, ref_ld, options);
    result.trait1 = trait1;
    result.trait2 = trait2;

    % Save the result in .mat file
    fname = sprintf('BGMG_run_%s-%s', trait1_name, trait2_name);
    save([fname '.mat'], 'result');
else
    result = BGMG_fit(data1.zvec, Hvec, data1.nvec, w_ld, ref_ld, options);
    result.trait1 = trait1;
    
    % Save the result in .mat file
    if size(data1.zvec, 2) == 1
        fname = sprintf('UGMG_run_%s', trait1_name);
    else
        fname = sprintf('BGMG_run_%s', trait1_name);
    end
    save([fname '.mat'], 'result');
end

fileID = fopen([fname '.txt'], 'w');
BGMG_util.result2str(1,      result)
BGMG_util.result2str(fileID, result)
fclose(fileID);
%catch err
%end

if ~exist('trait2', 'var')
    o2=options; o2.calculate_z_cdf = true;
    if ~isfield(o2, 'calculate_z_cdf_limit'), o2.calculate_z_cdf_limit = 15; end;
    if ~isfield(o2, 'calculate_z_cdf_step'), o2.calculate_z_cdf_step = 0.25; end;
    z_grid =  (-o2.calculate_z_cdf_limit:o2.calculate_z_cdf_step:o2.calculate_z_cdf_limit);
    [~, result_with_cdf] = BGMG_univariate_cost(result.univariate{1}.params, data1.zvec(:, 1), Hvec, data1.nvec(:, 1), w_ld, ref_ld, o2);
  
    fig=figure;
    zvec = data1.zvec; result_cdf = result_with_cdf.cdf;
    % weights = 1 ./ w_ld;
    %weights(w_ld < 1) = 1;
    weights=ones(size(data1.zvec));
    numstrata = 1;
    
    x = ref_ld.sum_r2;
    y = Hvec;
    % y = ref_ld.chi_r4;
    defvec=isfinite(zvec+x+y+weights);
           
    zvec=zvec(defvec); x=x(defvec); y=y(defvec); result_cdf = result_cdf(defvec, :); weights=weights(defvec);
    
    strat = false(numstrata,numstrata,sum(defvec));
    if (numstrata >= 3), xq = [-Inf, quantile(x,numstrata-2), +Inf]; end;
    if (numstrata >= 3), yq = [-Inf, quantile(y,numstrata-2), +Inf]; end;
    for i=1:numstrata
        for j=1:numstrata
            idx = true(size(x));
            titl = '';
            if i ~= numstrata, idx = idx & ((x >= xq(i)) & (x <= xq(i+1))); titl = sprintf('%s%.1f<=TLD <=%.1f', titl, xq(i), xq(i+1)); end;
            if i ~= numstrata && j ~= numstrata, titl = sprintf('%s\n', titl); end;
            if j ~= numstrata, idx = idx & ((y >= yq(j)) & (y <= yq(j+1))); titl = sprintf('%s%.3f<=HVEC<=%.3f', titl, yq(j), yq(j+1)); end;
            subplot(numstrata,numstrata, (i-1)*numstrata+j);
            title(titl);

            hold on
            qqlimy=20;
            qqlimx=7;
            
            data_cdf = 2*normcdf(-abs(zvec(idx)));
            weights_idx=weights(idx)/sum(weights(idx));
            
            [logpvec, si] = sort(-log10(data_cdf));
            a=1:(-1/length(logpvec)):0; a=a(1:length(logpvec)); a=-log10(a);
            %a=-log10(cumsum(weights_idx(si),1,'reverse'));
            plot(a,logpvec,'g'); 
            t=nanmean(result_cdf(idx, :));
            plot(-log10(2*min(t,1-t)),-log10(2*normcdf(-abs(z_grid),0,1)),'b'); 

            fontsize=19;
            plot([0 qqlimy],[0 qqlimy], 'k--');
            xlim([0 qqlimx]); ylim([0 qqlimy]);
            lgd=legend('Data', 'Model', 'Expected', 'Location', 'SouthEast');
            lgd.FontSize = fontsize;
            xlabel('Empirical -log 10(q)','fontsize',fontsize)
            ylabel('Nominal -log 10(p)','fontsize',fontsize)
            disp_trait1_name=trait1_name;disp_trait1_name(disp_trait1_name=='_') = '-';
            title(disp_trait1_name,'fontsize',fontsize);
            xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);
            yt = get(gca, 'YTick');set(gca, 'FontSize', fontsize);
            params=result.univariate{1}.params;
            ci=result.univariate{1}.ci;
            text(0.5,18,sprintf('$$ n_{total} = %i $$', sum(isfinite(data1.zvec))),'FontSize',fontsize,'Interpreter','latex');
            text(0.5,16,sprintf('$$ \\hat\\sigma_0^2 = %.3f $$', params.sig2_zero),'FontSize',fontsize,'Interpreter','latex');
            text(0.5,14,sprintf('$$ \\hat\\pi^u_1 = %.6f $$', params.pi_vec),'FontSize',fontsize,'Interpreter','latex');
            text(0.5,12,sprintf('$$ \\hat\\sigma_{\\beta}^2 = %.6f $$', params.sig2_beta),'FontSize',fontsize,'Interpreter','latex');
            text(0.5,10,sprintf('$$ \\hat  h^2 = %.6f $$', ci.h2.point_estimate),'FontSize',fontsize,'Interpreter','latex');
    
            

        end
    end
    print(fig, fname,'-dpdf')
end
