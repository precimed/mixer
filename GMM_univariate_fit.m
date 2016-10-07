function resultStruct = GMM_univariate_fit(zvec,Hvec,Lvec,imat_prune,options,x0,fignumoffset,bini)

if ~exist('options','var'), options = []; end
if ~exist('x0','var') | isempty(x0)
  x0 = GMM_univariate_mapparams(struct('pi1',0.001,'sig0',1.0,'sigb',5e-3),options); % Initial parameter guess
%  x0 = GMM_univariate_mapparams(struct('pi1',0.001,'sig0',1.0,'sigb',1.1),options); % Initial parameter guess
%  x0 = GMM_univariate_mapparams(struct('pi1',0.0062,'sig0',1.1071,'sigb',2.1266),options); % Initial parameter guess
end
if ~exist('fignumoffset','var'), fignumoffset = 0; end
if ~exist('bini','var'), bini = []; end


costfun = @(x) GMM_univariate_cost(x,zvec,Hvec,Lvec,imat_prune,options,fignumoffset,bini);

x_fit = x0;
[x_fit cost] = fminsearch_stochastic(costfun,x_fit); % First do stochastic optimization
[x_fit cost] = fminsearch(costfun,x_fit,statset('MaxIter',1000,'MaxFunEvals',1000,'Display','iter'));

GMM_univariate_mapparams(x_fit,options) % Print out x_fit as struct

%[cost fitstruct] = GMM_univariate_cost(x_fit,zvec,Hvec,Lvec,imat_prune,options,fignumoffset,bini); % Make plots

resultStruct = struct('x_fit',x_fit);

