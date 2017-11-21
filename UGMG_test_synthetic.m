addpath('DERIVESTsuite');
addpath('from_Anders');

nobs=1e6;
%nobs=1e4;

figure(11); clf;
for iterHIST = 1:3
for iterPARAMS = 1:3
fprintf('.');   
% Simulate bivariate associateion stats, given counts of potential causal SNPs across different LD r^2 bins

min_r2 = 0.05^2; nbins_r2 = 20;
%min_r2 = 0.1; nbins_r2 = 20; % Try higher minimum -- seems to reduce ringing / "foothills"
edges_r2 = [1 linspace(1,min_r2,nbins_r2)]; vals_r2 = (edges_r2(1:end-1)+edges_r2(2:end))/2;

if iterHIST==1, hist_r2_type='hist: perfect LD'; hist_r2 = zeros(size(vals_r2)); hist_r2(1) = 1; end    % Assume perfect LD with all causal SNPs
if iterHIST==3, hist_r2_type='hist: weak LD'; hist_r2 = zeros(size(vals_r2)); hist_r2(end) = 1;  end  % Assume same weak LD with all causal SNPs
if iterHIST==2, hist_r2_type='hist: exp.decay LD r'; hist_r2 = -diff(expcdf(sqrt(edges_r2)));  end;           % Assume exponential decay of LD r, with uniform sampling of SNPs by position
hist_r2 = hist_r2 / sum(hist_r2);

nsnp = 200;
count_r2 = round(nsnp*hist_r2);

clear params
if iterPARAMS == 1, params.pi_1 = 1e-4; params.pi_2 = 1e-4; params.pi_3 = 0; params.sigb_1 = 3.0; params.sigb_2 = 3.0; params.sig0_1 = 1; params.sig0_2 = 1; params.rhob = 0; params.rho0 = 0; end; % Low polygenisity for both
if iterPARAMS == 2, params.pi_1 = 1e-3; params.pi_2 = 1e-3; params.pi_3 = 0; params.sigb_1 = 3.0; params.sigb_2 = 3.0; params.sig0_1 = 1; params.sig0_2 = 1; params.rhob = 0; params.rho0 = 0; end; % Uneven polygenisities
if iterPARAMS == 3, params.pi_1 = 1e-2; params.pi_2 = 1e-2; params.pi_3 = 0; params.sigb_1 = 3.0; params.sigb_2 = 3.0; params.sig0_1 = 1; params.sig0_2 = 1; params.rhob = 0.7; params.rho0 = 0; end; % Including pleiotropic comp

delvals = linspace(-21,21,2^14+1); delstep = delvals(2)-delvals(1); 

[zmat delmat_total delmat_indep delmat_pleio] = BGMG_generate_random_stats(params,vals_r2,count_r2,nobs);
hc_z = hist3(zmat,{delvals delvals}); 
pdfmat_z_obs = hc_z/(sum(hc_z(:))*delstep^2);
[qqvec1_obs qqvec2_obs logpvals_qq zvals_qq] = BGMG_compute_qq(pdfmat_z_obs,delvals);

% univariate kurtosis approximation
sum_r2 = sum(vals_r2    .* count_r2);
sum_r4 = sum(vals_r2.^2 .* count_r2);
pi1u = params.pi_1 + params.pi_3;
eta_factor = pi1u .* sum_r2 + (1-pi1u) * (sum_r4 ./ sum_r2);
pi1u_delta = pi1u .* sum_r2 ./ eta_factor;
sig1_delta = params.sigb_1.^2 * eta_factor;

Nvec = 1; Hvec = 1;
model_cdf = (1-pi1u_delta) .* normcdf(delvals, 0, params.sig0_1) + pi1u_delta .* normcdf(delvals, 0, sqrt(sig1_delta .* Hvec .* Hvec + params.sig0_1.^2));
model_pdf = (1-pi1u_delta) .* normpdf(delvals, 0, params.sig0_1) + pi1u_delta .* normpdf(delvals, 0, sqrt(sig1_delta .* Hvec .* Hvec + params.sig0_1.^2));

%hv_logp = fliplr(-log10(2*normcdf(delvals(delvals <= 0))));  
model_logpvec = fliplr(-log10(2 * model_cdf(delvals <= 0)));

subplot(3, 3, 3*(iterHIST-1) + iterPARAMS);hold on
plot(logpvals_qq,logpvals_qq,'k--',-log10(qqvec1_obs),logpvals_qq,model_logpvec,logpvals_qq,'LineWidth',1); xlim([0 6]);
ylim([0 20])
%h=title('Q-Q Plot Trait 1'); set(h,'FontSize',18);
h=ylabel(sprintf('Nominal -log_1_0(p)')); set(h,'FontSize',10);
h=xlabel(sprintf('Empirical -log_1_0(q)')); set(h,'FontSize',10);
simulated_label = sprintf('DATA : std=%.3f, kurt=%.3f', std(zmat(:, 1)), kurtosis(zmat(:, 1)));
model_label = sprintf('MODEL: std=%.3f kurt=%.3f',sqrt(sum(model_pdf*delstep .* delvals.^2)), sum(model_pdf*delstep .* (delvals.^4)) ./ sum(model_pdf*delstep .* (delvals.^2)).^2);
h=legend({'Expected' 'Simulated' 'Model'},'Location','SE'); set(h,'FontSize',12);

for extra_lines=2:3
    [zmat delmat_total delmat_indep delmat_pleio] = BGMG_generate_random_stats(params,vals_r2,count_r2,nobs);
    hc_z = hist3(zmat,{delvals delvals}); 
    pdfmat_z_obs = hc_z/(sum(hc_z(:))*delstep^2);
    [qqvec1_obs qqvec2_obs logpvals_qq zvals_qq] = BGMG_compute_qq(pdfmat_z_obs,delvals);
    
    ax = gca; ax.ColorOrderIndex = 1;
    plot(logpvals_qq,logpvals_qq,'k--', -log10(qqvec1_obs),logpvals_qq,model_logpvec,logpvals_qq,'LineWidth',1); xlim([0 6]);
    
    ax = gca; ax.ColorOrderIndex = 1;
    plot(logpvals_qq,logpvals_qq,'k--', -log10(qqvec2_obs),logpvals_qq,model_logpvec,logpvals_qq,'LineWidth',1); xlim([0 6]);
end

qq_options.fontsize = 12;
%text_yloc = max(logpvals_qq(isfinite(-log10(qqvec1_obs)))); text_step = text_yloc / 12; text_yloc = text_yloc-1.5*text_step;
text_yloc = ylim; text_yloc=text_yloc(2); text_step = text_yloc / 12; text_yloc = text_yloc-1.5*text_step;
text(0.5,text_yloc,model_label,'FontSize',qq_options.fontsize,'Interpreter','latex'); text_yloc = text_yloc - text_step; 
text(0.5,text_yloc,simulated_label,'FontSize',qq_options.fontsize,'Interpreter','latex'); text_yloc = text_yloc - text_step; 
text(0.5,text_yloc,sprintf('$$ \\pi^u_1 = %.5f $$', params.pi_1+params.pi_3),'FontSize',qq_options.fontsize,'Interpreter','latex'); text_yloc = text_yloc - text_step; 
text(0.5,text_yloc,sprintf('$$ N H \\sigma^2_{\\beta_1} = %.5f $$', params.sigb_1.^2),'FontSize',qq_options.fontsize,'Interpreter','latex'); text_yloc = text_yloc - text_step; 
text(0.5,text_yloc,sprintf('$$ TLD = %.5f $$', sum(count_r2 .* vals_r2)),'FontSize',qq_options.fontsize,'Interpreter','latex'); text_yloc = text_yloc - text_step; 
text(0.5,text_yloc,sprintf(hist_r2_type),'FontSize',qq_options.fontsize,'Interpreter','latex'); text_yloc = text_yloc - text_step; 

end
end

set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','a3'); % https://se.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html
print(gcf, 'simulated.pdf','-dpdf')
fprintf(' Done.\n');
