if libisloaded('bgmg'), unloadlibrary('bgmg'); end;
if ~libisloaded('bgmg'), loadlibrary('H:\GitHub\BGMG\src\build_win\bin\RelWithDebInfo\bgmg.dll', 'H:\GitHub\BGMG\src\bgmg_matlab.h'); end;

%loadlibrary('H:\GitHub\BGMG\src\build_win\bin\Debug\bgmg.dll', 'H:\GitHub\BGMG\src\bgmg_matlab.h');
%libfunctions('bgmg')
%unloadlibrary('bgmg');

if ~exist('ref', 'var')
ref = load('H:\Dropbox\shared\BGMG\HAPGEN_EUR_100K_11015883_reference_bfile_merged_ldmat_p01_SNPwind50k_per_allele_4bins.mat');
end

%hits = sum(ref.pruneidxmat, 2); weights = hits / size(ref.pruneidxmat, 2);

def2= load('H:\Dropbox\shared\BGMG\defvec_hapmap3.mat');

defvec = true(size(ref.defvec)) ;
defvec = def2.defvec;
defvec = defvec(ref.chrnumvec == 1);

%defvec = defvec & (ref.chrnumvec == 1);  
%defvec = defvec & (weights > 0) &
%defvec = defvec & def2.defvec;

weights = ones(size(defvec));
tag_indices = find(defvec);

m2c = @(x)(x-1); % convert matlab to c indices
check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
context=0;kmax = 5000; num_components=3;  kmax=10;
calllib('bgmg', 'bgmg_set_tag_indices', 0, length(defvec), length(tag_indices), m2c(tag_indices));  check();
calllib('bgmg', 'bgmg_set_option', 0,  'r2min', 0.01); check();
calllib('bgmg', 'bgmg_set_option', 0,  'kmax', kmax); check();
calllib('bgmg', 'bgmg_set_option', 0,  'max_causals', floor(0.02 * length(defvec))); check();  
calllib('bgmg', 'bgmg_set_option', 0,  'num_components', num_components); check();

fprintf('Loading LD structure...\n');
for c= 1 % 22:-1:1
    fprintf('chr %i...\n', c);
    if ~exist('last_c', 'var') || (last_c ~= c)
        c1 = load(['H:\work\hapgen_ldmat2_plink\', sprintf('bfile_merged_10K_ldmat_p01_SNPwind50k_chr%i.ld.mat', c)]); last_c = c;
        c1.index_A = c1.index_A + 1; 
        c1.index_B = c1.index_B + 1; % 0-based, comming from python
    end
    calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, length(c1.r2), m2c(c1.index_A), m2c(c1.index_B), c1.r2);  check();
%    clear('c1');
end

tic
calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();
calllib('bgmg', 'bgmg_set_weights_randprune', 0, 1000, 0.8);  check();  % standard randprune: 10 iterations at 0.8
hvec = ref.mafvec .* (1-ref.mafvec) * 2; hvec = hvec(ref.chrnumvec == 1);
calllib('bgmg', 'bgmg_set_hvec', 0, length(hvec), hvec);  check();
toc

if 0
    pBuffer = libpointer('singlePtr', zeros(sum(defvec), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', 0, sum(defvec), pBuffer);  check(); 
    weights = pBuffer.Value;
    clear pBuffer
end





trait=1; 

%calllib('bgmg', 'bgmg_set_zvec', 0, trait, sum(defvec), zvec(defvec));  check();
%calllib('bgmg', 'bgmg_set_nvec', 0, trait, sum(defvec), nvec(defvec));  check();
calllib('bgmg', 'bgmg_set_weights', 0, sum(defvec), weights(defvec));  check();

% weights(defvec)=1;calllib('bgmg', 'bgmg_set_weights', 0, sum(defvec), weights(defvec));  check();



%dat = load('H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_traits\simu_h2=0.1_pi1u=0.0001_rep=3.trait1.mat');zvec = dat.zvec;
%nvec = ones(size(zvec)) * 100000;

%def1= load('H:\Dropbox\shared\BGMG\defvec_HAPGEN_EUR_100K.mat'); def1=def1.defvec; def1 = def1(1:num_snp);
%def2= load('H:\Dropbox\shared\BGMG\defvec_hapmap3.mat'); def2=def2.defvec; def2 = def2(1:num_snp);
%defvec = def1 & def2 & isfinite(zvec) & (weights>=1); clear('def1', 'def2');
%hvec = ref.hvec; %2 * ref.mafvec(1:num_snp) .* ref.mafvec(1:num_snp);

% idea1 - [DONE] smooth num_causals (float-point) to avoid jumps in cost functions - implemented
% idea2 - [DONE] ignore low z scores , e.i. fit tail only - didn't help
% idea3 - TBD: boost kmax (but keep little info about LD structure by permuting snp_order
% idea4 - make snp_order such that each SNP in the template became causal at least once.

if 0
figure(2)
num_snp = length(zvec); 
pBuffer = libpointer('singlePtr', zeros(kmax*sum(defvec), 1, 'single'));
calllib('bgmg', 'bgmg_retrieve_tag_r2_sum', 0, 0, round(dat.causal_pi * num_snp), kmax*sum(defvec), pBuffer);  check(); 
tag_r2_sum = pBuffer.Value;
tag_r2_sum=reshape(tag_r2_sum, [kmax, sum(defvec)])';
clear pBuffer
end

clf;figure(1); hold on;

h2vec_str = {'0.1', '0.4', '0.7'};
pivec_str = {'1e-05', '0.0001', '0.001', '0.01'};

if 1
    % test bivariate stuff 
    trait_name_vec={};
    %trait_name_vec{end+1, 1} = 'simu_h2=0.1_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=1e-03_rep=1_tag1=completePolygenicOverlap_tag2=evenPolygenicity';
    %trait_name_vec{end+1, 1} = 'simu_h2=0.1_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=3e-04_rep=1_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity';
    %trait_name_vec{end+1, 1} = 'simu_h2=0.1_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=1e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity'; 
    trait_name_vec{end+1, 1} = 'simu_h2=0.4_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=1e-03_rep=1_tag1=completePolygenicOverlap_tag2=evenPolygenicity';
    trait_name_vec{end+1, 1} = 'simu_h2=0.4_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=3e-04_rep=1_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity';
    trait_name_vec{end+1, 1} = 'simu_h2=0.4_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=1e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity'; 
    %trait_name_vec{end+1, 1} = 'simu_h2=0.7_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=1e-03_rep=1_tag1=completePolygenicOverlap_tag2=evenPolygenicity';
    %trait_name_vec{end+1, 1} = 'simu_h2=0.7_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=3e-04_rep=1_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity';
    %trait_name_vec{end+1, 1} = 'simu_h2=0.7_rg=0.0_pi1u=1e-03_pi2u=1e-03_pi12=1e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity'; 
    
    %trait_name_vec{1} = 'simu_h2=0.4_rg=0.0_pi1u=1e-02_pi2u=1e-02_pi12=1e-02_rep=1_tag1=completePolygenicOverlap_tag2=evenPolygenicity';fig=1;
    %trait_name_vec{2} = 'simu_h2=0.4_rg=0.0_pi1u=1e-02_pi2u=1e-02_pi12=3e-03_rep=1_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity';fig=2;
    %trait_name_vec{3} ='simu_h2=0.4_rg=0.0_pi1u=1e-02_pi2u=1e-02_pi12=1e-04_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity'; fig=3;
    
    %trait_name_vec{1} = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=3e-03_rep=1_tag1=completePolygenicOverlap_tag2=evenPolygenicity';fig=1;
    %trait_name_vec{2} = 'simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=8e-04_rep=1_tag1=partial25PolygenicOverlap_tag2=evenPolygenicity';fig=2;
    %trait_name_vec{3} ='simu_h2=0.4_rg=0.0_pi1u=3e-03_pi2u=3e-03_pi12=9e-06_rep=1_tag1=randomPolygenicOverlap_tag2=evenPolygenicity'; fig=3;
    
    dat_trait1 = {}; dat_trait2 = {};
    for fig=1:length(trait_name_vec);
        trait_name= trait_name_vec{fig};
        folder = 'H:\NORSTORE\oleksanf\11015833\simu_BGMG_chr1_pqexp0\';
        %folder = 'H:\NORSTORE\oleksanf\11015833\SIMU_BGMG2\';
        dat_trait1{fig} = load([folder, trait_name, '.trait1.mat']);
        dat_trait2{fig} = load([folder, trait_name, '.trait2.mat']);

        dat_trait1{fig}.nvec = ones(size(dat_trait1{fig}.zvec)) * 10000;
        dat_trait2{fig}.nvec = ones(size(dat_trait2{fig}.zvec)) * 10000;
    end
    
    % type of plots
    plots = {};
    %plots.pi12frac -     at correct pi1u, pi2u, varying pi12 fraction from 0 to 1
    %plots.pivec -        varying pi1, pi2, pi12 from 0.1 to 10 of the correct value
    %plots.sig2beta -     varying sig2beta from 0.1 to 10 of the correct value
    %plots.const_h2 -     simultaneously varying pi and sig2beta at correct h2=pi*sig2beta
    %plots.sig2zero -     varying sig2zero from 0.1 to 10 of the correct value
    %plots.rhozero -      varying rhozero from -1 to 1
    %plots.rhobeta -      varying rhobeta from -1 to 1

    set_data = @(dt1, dt2) calllib('bgmg', 'bgmg_set_zvec', 0, 1, sum(defvec), dt1.zvec(defvec)) + ...
                           calllib('bgmg', 'bgmg_set_nvec', 0, 1, sum(defvec), dt1.nvec(defvec)) + ...
                           calllib('bgmg', 'bgmg_set_zvec', 0, 2, sum(defvec), dt2.zvec(defvec)) + ...
                           calllib('bgmg', 'bgmg_set_nvec', 0, 2, sum(defvec), dt2.nvec(defvec));
         
	logspace_43 = logspace(log10(3/4), log10(4/3), 11);
    for fig=1:length(dat_trait1)
        set_data( dat_trait1{fig},  dat_trait2{fig} );
        pi1u = dat_trait1{fig}.causal_pi;
        pi2u = dat_trait2{fig}.causal_pi;
        pi12 = sum(pi1u(pi1u>0 & pi2u>0)); pi1u=sum(pi1u); pi2u=sum(pi2u);
        
        plot_name = 'pi12_frac';
        plots{fig}.(plot_name).x = 0:0.05:1; plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [(1-x), (1-x), x] * pi1u;
            sig2_beta = [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0]; rho_beta = 0; rho_zero = 0;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;
        
        plot_name = 'pivec_scale';
        plots{fig}.(plot_name).x = logspace_43; plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [x, x, x] .* [pi1u-pi12, pi2u-pi12, pi12];
            sig2_beta = [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0]; rho_beta = 0; rho_zero = 0;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;

        plot_name = 'sig2beta_scale';
        plots{fig}.(plot_name).x = logspace_43; plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
            sig2_beta = [x, x] .* [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0]; rho_beta = 0; rho_zero = 0;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;
        
        plot_name = 'pivec_scale_const_h2';
        plots{fig}.(plot_name).x = logspace_43; plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [x, x, x] .* [pi1u-pi12, pi2u-pi12, pi12];
            sig2_beta = [1./x, 1./x] .* [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0]; rho_beta = 0; rho_zero = 0;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;
        
        plot_name = 'sig2zero';
        plots{fig}.(plot_name).x = logspace(log10(0.5), log10(2), 11); plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
            sig2_beta = [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0] * x; rho_beta = 0; rho_zero = 0;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;
        
        plot_name = 'rhozero';
        plots{fig}.(plot_name).x = -0.1:0.01:0.1; plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
            sig2_beta = [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0]; rho_beta = 0; rho_zero = x;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;
        
        plot_name = 'rhobeta';
        plots{fig}.(plot_name).x = -0.2:0.02:0.2; plots{fig}.(plot_name).y = [];
        for x = plots{fig}.(plot_name).x, 
            pi_vec = [pi1u-pi12, pi2u-pi12, pi12];
            sig2_beta = [dat_trait1{fig}.sigsq, dat_trait2{fig}.sigsq];
            sig2_zero = [1.0, 1.0]; rho_beta = x; rho_zero = 0;
            plots{fig}.(plot_name).y(end+1, 1) = calllib('bgmg', 'bgmg_calc_bivariate_cost', 0, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);  check(); 
        end;
    %end
    
    figure(1); hold on; 
    %for fig=1:3
        fnames = fieldnames(plots{fig});
        for fni=1:length(fnames)
            subplot(3,3,fni);hold on; 
            plot_name=fnames{fni};
            y = plots{fig}.(plot_name).y; y(y>1e99)=nan;
            plot(plots{fig}.(plot_name).x, y-min(y), '-*');
            plot_name(plot_name=='_')=' ';title(plot_name);
            %if fig==3
            %    legend({'complete', 'partial25', 'random'})
            %end
        end
        drawnow;
    end
    
    calllib('bgmg', 'bgmg_set_option', 0,  'diag', 0); check();
end
return

for repi=1:10
for h2_index = [2 3 1]
for pi_index = [3 2 4 1]
   subplot(3, 4, (h2_index-1)*4 + pi_index);
   title(sprintf('h2=%s, pi=%s', h2vec_str{h2_index}, pivec_str{pi_index}));
%end;end;end;

hold on
ax = gca;
ax.ColorOrderIndex = 1;
dat= load(['H:\NORSTORE\oleksanf\11015833\simu_ugmg_120_traits\', sprintf('simu_h2=%s_pi1u=%s_rep=%i.trait1.mat', h2vec_str{h2_index}, pivec_str{pi_index}, repi)]);

dat.sig2zero = 1;
zvec = dat.zvec; nvec = ones(size(zvec)) * 100000;
zvec(~isfinite(zvec)) = 100;
%sum(~isfinite(zvec(defvec)))

zvec_test = zvec; %zvec_test(abs(zvec_test) < 0.5) = nan;    % <-------------------- try idea with not fitting bad z scores
calllib('bgmg', 'bgmg_set_zvec', 0, trait, sum(defvec), zvec_test(defvec));  check();
calllib('bgmg', 'bgmg_set_nvec', 0, trait, sum(defvec), nvec(defvec));  check();

fprintf('Fitting params...\n');
model_params{1} = struct('pi_vec', dat.causal_pi, 'sig2_beta', dat.sigsq, 'sig2_zero', dat.sig2zero);
model_params{2} = BGMG_util.UGMG_mapparams1(fminsearch(@BGMG_util.UGMG_CPP_fminsearch_cost, BGMG_util.UGMG_mapparams1(struct('pi_vec', dat.causal_pi, 'sig2_beta', dat.sigsq, 'sig2_zero', 1.05)), struct('Display', 'off')));

% calculate model cdf
zgrid = single(0:0.05:15); 

model_logpvec={};
for mi = 1:2
    m=model_params{mi};
    pBuffer = libpointer('singlePtr', zeros(length(zgrid), 1, 'single'));
    %cost=calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, dat.sig2zero, dat.sigsq);  check(); 
    calllib('bgmg', 'bgmg_calc_univariate_pdf', 0, m.pi_vec, m.sig2_zero, m.sig2_beta, length(zgrid), zgrid, pBuffer);  check(); 
    pdf = pBuffer.Value'; clear pBuffer
    pdf = pdf / sum(weights(defvec));
    if (zgrid(1) == 0), zgrid = [-fliplr(zgrid(2:end)) zgrid];pdf = [fliplr(pdf(2:end)) pdf]; end
    model_cdf = cumsum(pdf)  * (zgrid(2) - zgrid(1)) ;
    X = model_cdf;X1 = ones(size(X, 1), 1); X0 = zeros(size(X, 1), 1);
    model_cdf = 0.5 * ([X0, X(:, 1:(end-1))] + [X(:, 1:(end-1)), X1]);
    model_logpvec{mi} = -log10(2*interp1(-zgrid(zgrid<=0), model_cdf(zgrid<=0), zgrid(zgrid>=0))); % hv_z is to fine, can't afford calculation on it - do interpolation instead; don't matter for QQ plot (no visual difference), but lamGCfromQQ doesn't work for z_grid (to coarse)
end

% calculate data cdf        
zvec_idx = zvec(defvec); weights_idx = weights(defvec); 
weights_idx=weights_idx/sum(weights_idx);
[data_y, si] = sort(-log10(2*normcdf(-abs(zvec_idx))));
data_x=-log10(cumsum(weights_idx(si),1,'reverse'));
data_idx = ([data_y(2:end); +Inf] ~= data_y);
hv_logp = -log10(2*normcdf(-zgrid(zgrid >= 0)));
data_logpvec = interp1(data_y(data_idx), data_x(data_idx), hv_logp);

% plot QQ plots
hData        = plot(data_logpvec, hv_logp, '-', 'LineWidth',1); hold on;
hModel_true  = plot(model_logpvec{1},hv_logp, '-', 'LineWidth',1); hold on;
hModel_fit   = plot(model_logpvec{2},hv_logp, '-', 'LineWidth',1); hold on;

qq_options=[];
if ~isfield(qq_options, 'qqlimy'), qq_options.qqlimy = 20; end;
if ~isfield(qq_options, 'qqlimx'), qq_options.qqlimx = 7; end;
plot([0 qq_options.qqlimy],[0 qq_options.qqlimy], 'k--');
xlim([0 qq_options.qqlimx]); ylim([0 qq_options.qqlimy]);
drawnow
end
end
end

return

koef_vec = logspace(-1, 1, 31);
cost_pi = [];cost_sig2beta = []; cost_sig2zero = [];
for koef = koef_vec
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi * koef, dat.sig2zero, dat.sigsq);  check(); fprintf('%.3f\t', cost); cost_pi(end+1, 1) = cost;
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, dat.sig2zero * koef, dat.sigsq);  check(); fprintf('%.3f\t', cost); cost_sig2zero(end+1, 1) = cost;
    cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, dat.causal_pi, dat.sig2zero, dat.sigsq * koef);  check(); fprintf('%.3f\n', cost); cost_sig2beta(end+1, 1) = cost;
end
cost_pi(cost_pi > 1e99) = nan;
figure(2);hold on; subplot(1,3,1); plot(log10(koef_vec), cost_pi)
figure(2);hold on; subplot(1,3,2); plot(log10(koef_vec), cost_sig2beta)
figure(2);hold on; subplot(1,3,3); plot(log10(koef_vec), cost_sig2zero)


dat.causal_pi = exp(f_opt(1));
dat.sig2zero = f_opt(2);
dat.sigsq = exp(f_opt(3));
