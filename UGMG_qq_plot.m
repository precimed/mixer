function [ugmg_cdf, figures] = UGMG_qq_plot(params, zvec, Hvec, Nvec, pruneidxmat_or_w_ld, ref_ld, options, ugmg_cdf)
    % QQ plot for data and model

    if ~isfield(options, 'title'), options.title = 'UNKNOWN TRAIT'; end;
    if ~isfield(options, 'plot_HL_bins'), options.plot_HL_bins = true; end;

    figures.tot = figure; hold on;
    if options.plot_HL_bins, figures.bin = figure('units','normalized','outerposition',[0 0 1 1]); hold on; end;

    if size(pruneidxmat_or_w_ld, 2) == 1
        weights = 1 ./ pruneidxmat_or_w_ld;
        weights(pruneidxmat_or_w_ld < 1) = 1;
    else
        pruneidxmat = pruneidxmat_or_w_ld;
        nprune = size(pruneidxmat, 2);
    end

    if options.plot_HL_bins
        % Create a grid of HxL bins. Add last column for "all SNPs" bin.
        numstrata = 4;

        cdf_idx     = false(length(zvec), numstrata*numstrata + 1);

        x_bin = Hvec; y_bin = sum(ref_ld.sum_r2, 2);
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
    else
        cdf_idx         = isfinite(zvec);
    end

    if ~exist('ugmg_cdf', 'var') || isempty(ugmg_cdf)
        % For each bin, calculate weights based on random-pruning
        cdf_weights = zeros(length(zvec), size(cdf_idx, 2));

        for i=1:size(cdf_idx, 2)
            if exist('nprune', 'var')
                hits = sum(pruneidxmat & repmat(cdf_idx(:, i), [1 nprune]), 2);
                cdf_weights(:, i) = hits ./ nansum(hits);
            else
                weights_idx = weights; weights_idx(~cdf_idx(:, i)) = 0;
                cdf_weights(:, i) =  weights_idx ./ nansum(weights_idx);
            end
        end

        options.verbose = true;
        options.calculate_z_cdf = true;
        options.calculate_z_cdf_limit = ceil(min(max(abs(zvec)), 15));
        options.calculate_z_cdf_step = 0.05;
        options.calculate_z_cdf_weights = cdf_weights;
        w_ld_dummy = ones(size(zvec)); % dummy parameter --- it is used in cost calculation, which we are not interested in.
        [~, ugmg_cdf] = BGMG_univariate_cost(params, zvec, Hvec, Nvec, w_ld_dummy, ref_ld, options);
    end

    hv_z = linspace(0, min(max(abs(zvec)), 38.0), 10000);

    for ploti=1:size(cdf_idx, 2)
        if options.plot_HL_bins && (ploti ~= size(cdf_idx, 2))
            stratj = mod(ploti, numstrata);
            if stratj==0, stratj=stratj+numstrata; end;
            strati = (ploti-stratj)/numstrata+1;
            is_bin_plot = true;
            figure(figures.bin);
            subplot(numstrata,numstrata, ploti);
        else
            is_bin_plot = false;
            figure(figures.tot);
        end

        data_logpvec = zeros(size(hv_z));
        if exist('nprune', 'var')
            for repi = 1:nprune
                data_logpvecI = -log10(normcdf(-abs(zvec(cdf_idx(:, ploti) & pruneidxmat(:, repi))))*2);
                [data_logqvecI, hv_logp] = GenStats_QQ_plot_amd(data_logpvecI,hv_z);
                data_logpvec = data_logpvec + data_logqvecI;
                %lamGC_empirical(repi) = nanmedian(zvec(pruneidxmat(:, repi)).^2)/chi2inv(   0.5,1);
            end
            data_logpvec = data_logpvec / nprune;
        else
            zvec_idx = zvec(cdf_idx(:, ploti)); weights_idx = weights(cdf_idx(:, ploti));
            defvec=isfinite(zvec_idx+weights_idx);
            zvec_idx=zvec_idx(defvec); weights_idx=weights_idx(defvec); weights_idx=weights_idx/sum(weights_idx);

            [data_y, si] = sort(-log10(2*normcdf(-abs(zvec_idx))));
            data_x=-log10(cumsum(weights_idx(si),1,'reverse'));
            data_idx = ([data_y(2:end); +Inf] ~= data_y);
            hv_logp = -log10(2*normcdf(-hv_z));
            data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);
        end

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
        qq_options.lamGC_data = BGMG_util.lamGCfromQQ(data_logpvec, hv_logp);
        qq_options.lamGC_model = BGMG_util.lamGCfromQQ(model_logpvec, hv_logp);
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
