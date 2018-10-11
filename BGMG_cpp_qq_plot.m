function plot_data = BGMG_cpp_qq_plot(params, trait_index, options)
    % QQ plot for data and model

    if ~isfield(options, 'title'), options.title = 'UNKNOWN TRAIT'; end;
    if ~isfield(options, 'downscale'), options.downscale = 10; end;
    if ~isfield(options, 'full_annotation'), options.full_annotation = true; end;
    plot_data = []; hold on;

    bgmglib = BGMG_cpp();
    weights_bgmg = bgmglib.weights;
    zvec = bgmglib.get_zvec(trait_index);
    
    % Check that weights and zvec are all defined
    if any(~isfinite(zvec)), error('all values must be defined'); end;
    if any(~isfinite(weights_bgmg)), error('all values must be defined'); end;

    % We may limit QQ plots to certain group of SNPs (let's say a bin of mafvec and TLD)
    % Such subset should be defined via mask (same convention as defvec, 1 in mask means "use for QQ plot)
    if ~isfield(options, 'mask'), options.mask = true(size(zvec)); end;

    % To speedup calculations we set temporary weights where many elements
    % are zeroed. For data QQ plots all elements are set to 1.
    % At the end of this function we restore weights_bgmg.
    downscale_mask_indices = false(sum(options.mask), 1); 
    downscale_mask_indices(1:options.downscale:end)=true;
    mask_indices = find(options.mask); downscale_indices = mask_indices(downscale_mask_indices);
    
    model_weights = zeros(size(weights_bgmg)); 
    model_weights(downscale_indices) = weights_bgmg(downscale_indices); 
    model_weights = model_weights ./ sum(model_weights);
    bgmglib.weights = model_weights;
    
    % Calculate data_logpvec
    data_weights = weights_bgmg(options.mask) / sum(weights_bgmg(options.mask));
    hv_z = linspace(0, min(max(abs(zvec)), 38.0), 10000);
    [data_y, si] = sort(-log10(2*normcdf(-abs(zvec(options.mask)))));
    data_x=-log10(cumsum(data_weights(si),1,'reverse'));
    data_idx = ([data_y(2:end); +Inf] ~= data_y);
    hv_logp = -log10(2*normcdf(-hv_z));
    data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

    % Calculate model_logpvec
    zgrid = single(0:0.05:38.0); 
    pdf = bgmglib.calc_univariate_pdf(trait_index, params.pi_vec, params.sig2_zero, params.sig2_beta, zgrid);
    pdf = pdf / sum(model_weights);
    if (zgrid(1) == 0), zgrid = [-fliplr(zgrid(2:end)) zgrid];pdf = [fliplr(pdf(2:end)) pdf]; end
    model_cdf = cumsum(pdf)  * (zgrid(2) - zgrid(1)) ;
    X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
    model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
    model_logpvec = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), hv_z')); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)

    % Restore original weights
    bgmglib.weights = weights_bgmg;

    hData  = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
    hModel = plot(model_logpvec,hv_logp, '-', 'LineWidth',1); hold on;

    plot_data.hv_logp = hv_logp;
    plot_data.data_logpvec = data_logpvec;
    plot_data.model_logpvec = model_logpvec;
    plot_data.data_x = data_x;
    plot_data.data_y = data_y;
    plot_data.params = params;
    plot_data.options = options;
    plot_data.options.calculate_z_cdf_weights='removed';
    plot_data.options.mafvec='removed';
    plot_data.options.mask='removed';
    plot_data.pdf = pdf;
    plot_data.pdf_zgrid = zgrid;
    plot_data.data_pval = -log10(2*normcdf(-abs(zvec(mask_modified))));
    plot_data.data_weights = data_weights;

    qq_options = params; % must be on the top of other lines, otherwise this assigment overwrites all qq_options
    qq_options.title = options.title;
    if isfield(params, 'pi_vec')
        qq_options.h2 = sum(qq_options.pi_vec.*qq_options.sig2_beta)*options.total_het;
        if length(qq_options.pi_vec) > 1
            qq_options.h2vec = (qq_options.pi_vec.*qq_options.sig2_beta) *options.total_het;
        end
    end
    
    qq_options.full_annotation = options.full_annotation;
    qq_options.lamGC_data = BGMG_util.lamGCfromQQ(data_logpvec, hv_logp);
    qq_options.lamGC_model = BGMG_util.lamGCfromQQ(model_logpvec, hv_logp);
    qq_options.n_snps = sum(options.mask);
    qq_options.qqlimy = hv_logp(find(isfinite(data_logpvec), 1, 'last' ))*1.05;

    plot_data.qq_options = qq_options;

    annotate_qq_plot(qq_options);
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
    if qq_options.full_annotation
        if has_opt('legend') && qq_options.legend, lgd=legend('Data', 'Model', 'Expected', 'Location', 'SouthEast'); lgd.FontSize = qq_options.fontsize; end;
        if has_opt('xlabel') && qq_options.xlabel, xlabel('Empirical -log 10(q)','fontsize',qq_options.fontsize); end;
        if has_opt('ylabel') && qq_options.ylabel, ylabel('Nominal -log 10(p)','fontsize',qq_options.fontsize); end;
    end
    if has_opt('title'), title(qq_options.title,'fontsize',qq_options.fontsize,'Interpreter','latex'); end;
    xt = get(gca, 'XTick');set(gca, 'FontSize', qq_options.fontsize);
    yt = get(gca, 'YTick');set(gca, 'FontSize', qq_options.fontsize);
    set(gca, 'box','off');

    loc_delta = 0.1 * qq_options.qqlimy;
    loc = qq_options.qqlimy-loc_delta;
    if has_opt('n_snps'), text(0.5,loc,sprintf('$$ n_{snps} = %i $$', qq_options.n_snps),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
    if qq_options.full_annotation
    if has_opt('sig2_zero'), text(0.5,loc,sprintf('$$ \\hat\\sigma_0^2 = %.3f $$', qq_options.sig2_zero),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
    if has_opt('pi_vec'), text(0.5,loc,sprintf('$$ \\hat\\pi^u_1 = %s $$', vec2str(qq_options.pi_vec)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
    if has_opt('sig2_beta'), text(0.5,loc,sprintf('$$ \\hat\\sigma_{\\beta}^2 = %s $$', vec2str(qq_options.sig2_beta)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
    h2vec = ''; if has_opt('h2vec'), h2vec = ['$$ \\; $$' vec2str(qq_options.h2vec, 3)]; end;
    if has_opt('h2'), text(0.5,loc,sprintf('$$ \\hat h^2 = %.3f%s$$', qq_options.h2, h2vec),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
    end
    if has_opt('lamGC_model'), text(0.5,loc,sprintf('$$ \\hat\\lambda_{model} = %.3f $$', qq_options.lamGC_model),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
    if has_opt('lamGC_data'), text(0.5,loc,sprintf('$$ \\hat\\lambda_{data} = %.3f $$', qq_options.lamGC_data),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - loc_delta; end;
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
