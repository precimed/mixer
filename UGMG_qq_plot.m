function [ugmg_cdf, fig_tot, fig_bin] = UGMG_qq_plot(params, zvec, Hvec, Nvec, pruneidxmat, ref_ld, options, ugmg_cdf)
    if ~isfield(options, 'title'), options.title = 'UNKNOWN TRAIT'; end;

    fig_tot = figure; hold on;
    fig_bin = figure('units','normalized','outerposition',[0 0 1 1]); hold on;

    % QQ plot for data and model
    nprune = size(pruneidxmat, 2);

    % Create a grid of HxL bins. Add last column for "all SNPs" bin.
    numstrata = 4;

    cdf_idx     = false(length(zvec), numstrata*numstrata + 1);
    x_bin = Hvec; y_bin = ref_ld.sum_r2;
    defvec = isfinite(zvec + x_bin + y_bin);
    xq = [-Inf, quantile(x_bin(defvec),numstrata-1), +Inf];
    yq = [-Inf, quantile(y_bin(defvec),numstrata-1), +Inf];

    % Find SNP indices for each HxL bin.
    fprintf('Distribution of SNPs into HxL bins:\n');
    for i=1:numstrata,
        for j=1:numstrata
            idx = isfinite(zvec) & ((x_bin >= xq(i)) & (x_bin <= xq(i+1))) & ((y_bin >= yq(j)) & (y_bin <= yq(j+1)));
            cdf_idx(:, (i-1)*numstrata + j) = idx;
            fprintf('%i ', sum(idx));
        end
        fprintf('\n');
    end
	cdf_idx(:, end) = isfinite(zvec);

    if ~exist('ugmg_cdf', 'var') || isempty(ugmg_cdf)
        % For each bin, calculate weights based on random-pruning
        cdf_weights = zeros(length(zvec), numstrata*numstrata + 1);
        for i=1:size(cdf_idx, 2)
            hits = sum(pruneidxmat & repmat(cdf_idx(:, i), [1 nprune]), 2);
            cdf_weights(:, i) = hits ./ nansum(hits);
        end

        options.verbose = true;
        options.calculate_z_cdf = true;
        options.calculate_z_cdf_limit = ceil(min(max(abs(zvec)), 15));
        options.calculate_z_cdf_step = 0.05;
        options.calculate_z_cdf_weights = cdf_weights;
        w_ld = ones(size(zvec)); % dummy parameter --- it is used in cost calculation, which we are not interested in.
        [~, ugmg_cdf] = BGMG_univariate_cost(params, zvec, Hvec, Nvec, w_ld, ref_ld, options);
        return
    end

    hv_z = linspace(0, min(max(abs(zvec)), 38.0), 10000);

    for ploti=1:size(cdf_idx, 2)
        stratj = mod(ploti, numstrata);
        if stratj==0, stratj=stratj+numstrata; end;
        strati = (ploti-stratj)/numstrata+1;
        is_bin_plot = (ploti ~= size(cdf_idx, 2));
        if is_bin_plot,
            figure(fig_bin);
            subplot(numstrata,numstrata, ploti);
        else
            figure(fig_tot);
        end

        data_logpvec = zeros(size(hv_z));
        for repi = 1:nprune
            data_logpvecI = -log10(normcdf(-abs(zvec(cdf_idx(:, ploti) & pruneidxmat(:, repi))))*2);
            [data_logqvecI, hv_logp] = GenStats_QQ_plot_amd(data_logpvecI,hv_z);
            data_logpvec = data_logpvec + data_logqvecI;
            %lamGC_empirical(repi) = nanmedian(zvec(pruneidxmat(:, repi)).^2)/chi2inv(   0.5,1);
        end
        data_logpvec = data_logpvec / nprune;

        z_grid = ugmg_cdf.cdf_z_grid;
        model_logpvec = -log10(2*interp1(-z_grid(z_grid<=0), ugmg_cdf.cdf(ploti, z_grid<=0), hv_z')); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)

        hData     = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
        hModel    = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;

        qq_options = [];
        if is_bin_plot
            qq_options.title = sprintf('$$ H \\in [%.3f,%.3f) $$ \n $$ L \\in [%.3f,%.3f) $$', xq(strati), xq(strati+1), yq(stratj), yq(stratj+1));
            qq_options.fontsize=15;
            qq_options.legend=false;
            qq_options.xlabel=false;
            qq_options.ylabel=false;
            qq_options.qqlimy = 16;
        else
            qq_options = params; % must be on the top of other lines, otherwise this assigment overwrites all qq_options
            qq_options.title = options.title;
            qq_options.h2 = sum(qq_options.pi_vec.*qq_options.sig2_beta)*options.total_het;
            if length(qq_options.pi_vec) > 1
                qq_options.h2vec = (qq_options.pi_vec.*qq_options.sig2_beta) *options.total_het;
            end
        end
        qq_options.lamGC_data = lamGCfromQQ(data_logpvec, hv_logp);
        qq_options.lamGC_model = lamGCfromQQ(model_logpvec, hv_logp);
        qq_options.n_snps = sum(isfinite(zvec) & cdf_idx(:, ploti));

        annotate_qq_plot(qq_options);
    end
end

function annotate_qq_plot(qq_options)
    % Tune appearance of the QQ plot, put annotations, title, etc...
    if ~isfield(qq_options, 'qqlimy'), qq_options.qqlimy = 20; end;
    if ~isfield(qq_options, 'qqlimx'), qq_options.qqlimx = 7; end;
    if ~isfield(qq_options, 'fontsize'), qq_options.fontsize = 19; end;
    if ~isfield(qq_options, 'legend'), qq_options.legend = true; end;
    if ~isfield(qq_options, 'xlabel'), qq_options.xlabel = true; end;
    if ~isfield(qq_options, 'ylabel'), qq_options.ylabel = true; end;

    has_opt = @(opt)(isfield(qq_options, opt));

    plot([0 qq_options.qqlimy],[0 qq_options.qqlimy], 'k--');
    xlim([0 qq_options.qqlimx]); ylim([0 qq_options.qqlimy]);
    if has_opt('legend') && qq_options.legend, lgd=legend('Data', 'Model', 'Expected', 'Location', 'SouthEast'); lgd.FontSize = qq_options.fontsize; end;
    if has_opt('xlabel') && qq_options.xlabel, xlabel('Empirical -log 10(q)','fontsize',qq_options.fontsize); end;
    if has_opt('ylabel') && qq_options.ylabel, ylabel('Nominal -log 10(p)','fontsize',qq_options.fontsize); end;
    if has_opt('title'), title(qq_options.title,'fontsize',qq_options.fontsize,'Interpreter','latex'); end;
    xt = get(gca, 'XTick');set(gca, 'FontSize', qq_options.fontsize);
    yt = get(gca, 'YTick');set(gca, 'FontSize', qq_options.fontsize);
    set(gca, 'box','off');

    loc = qq_options.qqlimy-2;
    if has_opt('n_snps'), text(0.5,loc,sprintf('$$ n_{snps} = %i $$', qq_options.n_snps),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('sig2_zero'), text(0.5,loc,sprintf('$$ \\hat\\sigma_0^2 = %.3f $$', qq_options.sig2_zero),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('pi_vec'), text(0.5,loc,sprintf('$$ \\hat\\pi^u_1 = %s $$', vec2str(qq_options.pi_vec)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('sig2_beta'), text(0.5,loc,sprintf('$$ \\hat\\sigma_{\\beta}^2 = %s $$', vec2str(qq_options.sig2_beta)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    h2vec = ''; if has_opt('h2vec'), h2vec = ['$$ \\; $$' vec2str(qq_options.h2vec, 3)]; end;
    if has_opt('h2'), text(0.5,loc,sprintf('$$ \\hat h^2 = %.3f%s$$', qq_options.h2, h2vec),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('lamGC_model'), text(0.5,loc,sprintf('$$ \\hat\\lambda_{model} = %.3f $$', qq_options.lamGC_model),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('lamGC_data'), text(0.5,loc,sprintf('$$ \\hat\\lambda_{data} = %.3f $$', qq_options.lamGC_data),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
end

function s = vec2str(vec, digits)
    if ~exist('digits', 'var'), digits = 6; end;
    format = ['%.', int2str(digits), 'f'];
    if length(vec) == 1
        s=sprintf(format, vec);
    else
        s=['[', sprintf(format, vec(1))];
        for i=2:length(vec)
            s=[s, ', ', sprintf(format, vec(i))];
        end
        s = [s ']'];
    end
end

function other_qq_plots
    % original code for hv_z
    dz = 0.4; zextreme = 38.0;
    zextreme = min(max(abs(zvec)),zextreme);
    nzbins = 2*ceil( zextreme/dz ) + 1;   % Guaranteed odd.
    zvals_disc = linspace(-zextreme,zextreme,nzbins); % zvals_disc(2)-zvals_disc(1) == dz;       Values are bin centers. Middle bin value is 0 (center of center bin).
    hv_z = linspace(0,zvals_disc(end),10000);

    % Yet another data qq plot (directly from data points)
    zvec_qq = zvec; weights = hits; %weights=ones(size(hits));
    defvec=isfinite(zvec+weights);
    zvec_qq=zvec_qq(defvec); weights=weights(defvec); weights=weights/sum(weights);

    [data_logpvec2, si] = sort(-log10(2*normcdf(-abs(zvec_qq))));
    hv_logpvec2=-log10(cumsum(weights(si),1,'reverse'));
    plot(hv_logpvec2, data_logpvec2,'g');

    % Yet another model qq plot (directly on z_grid - no interpolation to hv_z)
    model_logpvec0 = -log10(2*model_cdf(z_grid <= 0));
    model_hv_logp0 = -log10(2*normcdf(-abs(z_grid(z_grid <= 0))));
    hModel0    = plot(model_logpvec0,model_hv_logp0, '-','LineWidth',1); hold on;
end
