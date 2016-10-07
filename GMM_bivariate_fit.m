function resultStruct = GMM_bivariate_fit(zmat,Hvec,Lvec,imat_prune,options,x0)

ngm = options.ngm;

if ~exist('x0','var') | isempty(x0)
  x0 = GMM_bivariate_mapparams(struct('sig01',1.0,'sig02',1.0,'rho0',0,'pivec',1e-2*ones(1,ngm),'sig1vec',1.0*rand(1,ngm),'sig2vec',1.0*rand(1,ngm),'rhovec',zeros(1,ngm)),options); % Initialize assuming fitting z-scores directly
%  x0 = GMM_bivariate_mapparams(struct('sig01',1.0,'sig02',1.0,'rho0',0,'pivec',[1e-2 1e-4],'sig1vec',[0.4 0.5],'sig2vec',[0.6 2.5],'rhovec',[0.4 0.95]),options);
%  x0 = GMM_bivariate_mapparams(struct('sig01',1.0,'sig02',1.0,'rho0',0,'pivec',[1e-1 2e-2],'sig1vec',[1.0 0.5],'sig2vec',[1.2 2.3],'rhovec',[0.4 0.95]),options);
end

[cost fitstruct] = GMM_bivariate_cost(x0,zmat,Hvec,Lvec,imat_prune,options); % Make plots before fitting

costfun = @(x)GMM_bivariate_cost(x,zmat,Hvec,Lvec,imat_prune,options);
[x_fit cost] = fminsearch(costfun,x0,statset('MaxIter',1000,'MaxFunEvals',1000,'Display','iter'));

[cost fitstruct] = GMM_bivariate_cost(x_fit,zmat,Hvec,Lvec,imat_prune,options); % Make plots after fitting

GMM_bivariate_mapparams(x_fit,options) % Print out x_fit as struct

if ~exist('fistruct','var'), fitstruct = struct(); end

resultStruct = struct('cost',cost,'x_fit',x_fit,'fitstruct',fitstruct);

keyboard

% ToDo
%   compute conditional & conjunctional tdr, etc.

