function [zvec, Fpdfmat_z, Fpdfmat_total, Fpdfmat_indep, Fpdfmat_pleio, pdfmat_z, pdfmat_total, pdfmat_indep, pdfmat_pleio] = BGMG_predict_pdfs_poisson(...
    params,vals_r2,count_r2,delvals,Neff,H,dispflag)

Fpdfmat_z = 0; Fpdfmat_total = 0; Fpdfmat_indep = 0; Fpdfmat_pleio = 0; % irrelevant
% zvec - grid in the space of z values; 
% pdfmat_z - total bivariate probability density function (z score)
% pdfmat_total - total bivariate probability density function (delta)
% pdfmat_indep - bivariate probability density function (delta, independent components i.e. 1 and 2)
% pdfmat_pleio - bivariate probability density function (delta, only pleiotropic component, i.e. 3)

%params = 
%      pi_1: 0.0100
%      pi_2: 1.0000e-03
%      pi_3: 0.0100
%    sigb_1: 3
%    sigb_2: 1.5000
%    sig0_1: 1
%    sig0_2: 1
%      rhob: 0.7000
%      rho0: 0

% vals_r2 =  [1.0000    0.9738    0.9213    0.8688    0.8163    0.7638    0.7113    0.6588    0.6062    0.5537    0.5012    0.4487    0.3962    0.3437    0.2912    0.2388    0.1862    0.1337    0.0812    0.0287      ]
% count_r2 =   [     0     2     2     2     2     2     2     2     3     3     3     3     4     4     5     6     7     9    12    27 ]
% delvals  - equally spaced grid with 512 nodes, both negative and positive  values

% If Neff or H not specified, assume sigb is in units of del (already scaled)
if ~exist('Neff','var') || isempty(Neff), Neff = 1; end 
if ~exist('H','var') || isempty(H), H = 1; end
if ~exist('dispflag','var') || isempty(dispflag), dispflag = false; end

nsnp = sum(count_r2); nbins_r2 = length(vals_r2);

%upsampfact = 16; % This is needed to handle tiny sig_b * r^2 -- should find more computationally efficient solution
zvec = linspace(delvals(1),-delvals(1),length(delvals)+1); zvec = zvec(1:end-1);
zstep = zvec(2)-zvec(1); 

sigd1 = sqrt(Neff*H)*params.sigb_1;
sigd2 = sqrt(Neff*H)*params.sigb_2;

% For each component let's build sigma grid and calculate weights on that grid
nbinsbins=10;
idx = false(nbins_r2, nbinsbins);
for t=1:nbinsbins, idx((1:nbins_r2/nbinsbins) + (t-1)*(nbins_r2/nbinsbins), t) = true; end

pivec = [params.pi_1, params.pi_2, params.pi_3];
grid_step = [1/3 1/3 1/3];  % w.r.t. Sigma
grid_size = 12;
kmax = 6;

grid_coef = zeros(length(pivec), grid_size); grid_coef(:, 1) = 1;
grid_coefK= zeros(length(pivec), kmax+1); grid_coefK(:, 1)= 1;
grid_stepK= zeros(length(pivec), 1);
for pi_index = 1:length(pivec)
    lambda=[];sig1_delta=[];
    for r2bini=1:nbinsbins
        i=idx(:, r2bini);
        sum_r2 = sum(vals_r2(:, i)    .* count_r2(:, i));
        sum_r4 = sum(vals_r2(:, i).^2 .* count_r2(:, i));
        chi2 = (sum_r4 ./ sum_r2);
        lambda(r2bini) = pivec(pi_index) * sum_r2 ./ ((1-pivec(pi_index)) * chi2);
        sig1_delta(r2bini) = (1-pivec(pi_index)) * chi2;  % w.r.t. Sigma(pi_index)
    end

    for i=1:length(lambda)
        k              = (0:kmax)';
        p              = lambda(i).^ k * exp(-lambda(i)) ./ factorial(k);
        grid_idx       = k * sig1_delta(i) / grid_step(pi_index);
        grid_idx       = min(grid_idx, grid_size - 2);
        grid_idx_floor = floor(grid_idx);
        grid_coef_tmp  = accumarray(1+grid_idx_floor, p .* (1 - grid_idx + grid_idx_floor), [grid_size 1]) + ...
                         accumarray(2+grid_idx_floor, p .* (0 + grid_idx - grid_idx_floor), [grid_size, 1]);
        grid_coef_tmp  = conv(grid_coef(pi_index, :), grid_coef_tmp);
        grid_coef(pi_index, :)=grid_coef_tmp(1:grid_size);
    end
    
    if length(lambda) == 1  % same as nbinsbins==1
        grid_coefK(pi_index, :) = p;
        grid_stepK(pi_index) = sig1_delta;
    end
end
assert(all(sum(grid_coef, 2) > 0.999));

if 0 % explicitly use poisson grid (NOT the projection to fixed sigma grid)
    assert(nbinsbins==1);
    grid_size = (kmax + 1);
    grid_coef = grid_coefK;
    grid_step = grid_stepK;
end

C0_delta = [1e-3 0; 0 1e-3];
C0 = [params.sig0_1^2 params.sig0_1*params.sig0_2*params.rho0; params.sig0_1*params.sig0_2*params.rho0 params.sig0_2^2];
C1 = grid_step(1) * [sigd1^2 0; 0 0];
C2 = grid_step(2) * [0 0; 0 sigd2^2];
C3 = grid_step(3) * [sigd1^2 sigd1*sigd2*params.rhob; sigd1*sigd2*params.rhob sigd2^2];

%model_cdf = zeros(length(delvals));
pdfmat_z = zeros(length(delvals));
pdfmat_total = zeros(length(delvals));
pdfmat_indep = zeros(length(delvals));
pdfmat_pleio = zeros(length(delvals));
[zmatX, zmatY] = meshgrid(delvals);

for k1=1:grid_size
    for k2=1:grid_size
        for k3=1:grid_size
            fprintf('[%i %i %i]\n', k1, k2, k3);
            weight = grid_coef(1, k1) .* grid_coef(2, k2) .* grid_coef(3, k3);
            C = (k1-1) * C1 + (k2-1) * C2 + (k3-1) * C3;
            pdfmat_z     = pdfmat_z     + weight * reshape(mvnpdf([zmatX(:) zmatY(:)], [0 0], C + C0), size(pdfmat_z));
            pdfmat_total = pdfmat_total + weight * reshape(mvnpdf([zmatX(:) zmatY(:)], [0 0], C + C0_delta), size(pdfmat_z));
        end
    end
end
for k1=1:grid_size
    for k2=1:grid_size
        fprintf('[%i %i -]\n', k1, k2);
        weight = grid_coef(1, k1) .* grid_coef(2, k2);
        C = (k1-1) * C1 + (k2-1) * C2;
        pdfmat_indep = pdfmat_indep + weight * reshape(mvnpdf([zmatX(:) zmatY(:)], [0 0], C + C0_delta), size(pdfmat_z));
    end
end
for k3=1:grid_size
    fprintf('[- - %i]\n', k1, k2);
    weight = grid_coef(3, k3);
    C = (k3-1) * C3;
    pdfmat_pleio = pdfmat_pleio + weight * reshape(mvnpdf([zmatX(:) zmatY(:)], [0 0], C + C0_delta), size(pdfmat_z));
end

pdfmat_z=pdfmat_z';
pdfmat_total=pdfmat_total';
pdfmat_indep=pdfmat_indep';
pdfmat_pleio=pdfmat_pleio';
