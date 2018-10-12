function [figures, plot_data] = BGMG_cpp_stratified_qq_plot(params, options)
    % QQ plot for data and model

    if ~isfield(options, 'title'), options.title = 'UNKNOWN TRAIT'; end;
    if ~isfield(options, 'downscale'), options.downscale = 10; end;
    plot_data = {}; figures.tot = figure; hold on;

    bgmglib = BGMG_cpp();
    weights_bgmg = bgmglib.weights;
    zmat = [bgmglib.zvec1, bgmglib.zvec2];

    % Check that weights and zvec are all defined
    if any(~isfinite(zmat(:))), error('all values must be defined'); end;
    if any(~isfinite(weights_bgmg)), error('all values must be defined'); end;

    % To speedup calculations we set temporary weights where many elements
    % are zeroed. The remaining elements are set to 1 (e.i. there is no
    % weighting). For data QQ plots all elements are set to 1.
    % At the end of this function we restore weights_bgmg.
    downscale_indices = false(size(weights_bgmg)); downscale_indices(1:options.downscale:end)=true;
    model_weights = zeros(size(weights_bgmg));  model_weights(downscale_indices) = weights_bgmg(downscale_indices); model_weights = model_weights ./ sum(model_weights);
    bgmglib.weights = model_weights;
    data_weights = weights_bgmg / sum(weights_bgmg);

    % Calculate bivariate pdf on a regular 2D grid (z1, z2)
    zgrid = single(-38:0.25:38);
    [zgrid1, zgrid2] = meshgrid(zgrid, zgrid);
    zgrid1 = zgrid1(zgrid >= 0, :); zgrid2 = zgrid2(zgrid >= 0, :);  % taking symmetry into account to speedup by a factor of 2
    pdf = bgmglib.calc_bivariate_pdf(params.pi_vec, params.sig2_beta(:, end),  params.rho_beta(end), params.sig2_zero, params.rho_zero, zgrid1(:), zgrid2(:));
    pdf = reshape(pdf, size(zgrid1));
    clear('zgrid1', 'zgrid2');
    pdf = [fliplr(flipud(pdf(2:end, :))); pdf];
    pdf = pdf / sum(model_weights);
    
    % Restore original weights
    bgmglib.weights = weights_bgmg;

    hv_z = linspace(0, min(max(abs(zgrid)), 38.0), 10000);
    hv_logp = -log10(2*normcdf(-hv_z));

    zthresh_vec = [0 1 2 3];
    for zthresh_index = 1:length(zthresh_vec)
        zthresh = zthresh_vec(zthresh_index);
        % Calculate data_logpvec
        zvec1 = zmat(:, 2);
        zvec2 = zmat(:, 1);
        zvec = zvec1(abs(zvec2)>=zthresh);
        weights = data_weights(abs(zvec2)>=zthresh); weights = weights ./ sum(weights);
        [data_y, si] = sort(-log10(2*normcdf(-abs(zvec))));
        data_x=-log10(cumsum(weights(si),1,'reverse'));
        data_idx = ([data_y(2:end); +Inf] ~= data_y);
        data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

        pdf_cond = sum(pdf(:, abs(zgrid)>=zthresh), 2);
        pdf_cond = pdf_cond ./ sum(pdf_cond);
        model_cdf = cumsum(pdf_cond');
        X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
        model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
        model_logpvec = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), hv_z')); 

        ax = gca;ax.ColorOrderIndex = zthresh_index;
        hData  = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
        ax = gca;ax.ColorOrderIndex = zthresh_index;
        hModel = plot(model_logpvec,hv_logp, '-.', 'LineWidth',1); hold on;
        
        plot_data{zthresh_index}.hv_logp = hv_logp;
        plot_data{zthresh_index}.data_logpvec = data_logpvec;
        plot_data{zthresh_index}.model_logpvec = model_logpvec;
        plot_data{zthresh_index}.params = params;
        plot_data{zthresh_index}.zthresh = zthresh;
        plot_data{zthresh_index}.pdf = pdf;
        plot_data{zthresh_index}.pdf_zgrid = zgrid;
    end  
    
    qq_options=[];
    qq_options.pi_vec = params.pi_vec;
    qq_options.sig2_zero = params.sig2_zero;
    qq_options.sig2_beta = params.sig2_beta(:, end);
    qq_options.rho_beta = params.rho_beta(:, end);
    qq_options.rho_zero = params.rho_zero;
    qq_options.title = options.title;

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
    if has_opt('legend') && qq_options.legend, lgd=legend('Data(z_1)', 'Model(z_1)', 'Data(z_1 : |z_2| \geq 1)', 'Model(z_1 : |z_2| \geq 1)', 'Data(z_1 : |z_2| \geq 2)', 'Model(z_1 : |z_2| \geq 2)', 'Data(z_1 : |z_2| \geq 3)', 'Model(z_1 : |z_2| \geq 3)', 'Expected', 'Location', 'SouthEast'); lgd.FontSize = qq_options.fontsize/2; end;
    if has_opt('xlabel') && qq_options.xlabel, xlabel('Empirical -log 10(q)','fontsize',qq_options.fontsize); end;
    if has_opt('ylabel') && qq_options.ylabel, ylabel('Nominal -log 10(p)','fontsize',qq_options.fontsize); end;
    if has_opt('title'), title(qq_options.title,'fontsize',qq_options.fontsize,'Interpreter','latex'); end;
    xt = get(gca, 'XTick');set(gca, 'FontSize', qq_options.fontsize);
    yt = get(gca, 'YTick');set(gca, 'FontSize', qq_options.fontsize);
    set(gca, 'box','off');

    loc = qq_options.qqlimy-2;
    if has_opt('n_snps'), text(0.5,loc,sprintf('$$ n_{snps} = %i $$', qq_options.n_snps),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('sig2_zero'), text(0.5,loc,sprintf('$$ \\hat\\sigma_0^2 = %s $$', vec2str(qq_options.sig2_zero)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('pi_vec'), text(0.5,loc,sprintf('$$ \\hat\\pi^u_1 = %s $$', vec2str(qq_options.pi_vec)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('sig2_beta'), text(0.5,loc,sprintf('$$ \\hat\\sigma_{\\beta}^2 = %s $$', vec2str(qq_options.sig2_beta)),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('rho_beta'), text(0.5,loc,sprintf('$$ \\rho_\\beta = %.3f $$', qq_options.rho_beta),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
    if has_opt('rho_zero'), text(0.5,loc,sprintf('$$ \\rho_0 = %.3f $$', qq_options.rho_zero),'FontSize',qq_options.fontsize,'Interpreter','latex'); loc = loc - 2; end;
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
