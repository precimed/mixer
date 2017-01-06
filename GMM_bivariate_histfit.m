function resultStruct = GMM_bivariate_histfit(sumstats,params,options,x0)

if isfield(options, 'MaxIter1'), MaxIter1 = options.MaxIter1; else MaxIter1 = 1000; end;
if isfield(options, 'MaxIter2'), MaxIter2 = options.MaxIter2; else MaxIter2 = 5000; end;

ngm = options.ngm;

%[cost fitstruct] = GMM_bivariate_histcost(x0,sumstats,options); % Make plots before fitting

costfun = @(x)GMM_bivariate_histcost(x,sumstats,params,options);

fprintf('fminsearch_stochastic\n');
[x_fit_stochastic cost] = fminsearch_stochastic(costfun,x0,statset('MaxIter',MaxIter1,'MaxFunEvals',MaxIter1,'Display','iter'));
GMM_bivariate_mapparams(x_fit_stochastic,options) % Print out x_fit as struct

fprintf('fminsearch\n');
[x_fit cost] = fminsearch(costfun,x_fit_stochastic,statset('MaxIter',MaxIter2,'MaxFunEvals',MaxIter2,'Display','off'));
%[cost fitstruct] = GMM_bivariate_histcost(x_fit,sumstats,params,options); % Make plots after fitting
GMM_bivariate_mapparams(x_fit,options) % Print out x_fit as struct

%resultStruct = struct('cost',cost,'x_fit',x_fit,'fitstruct',fitstruct);
resultStruct = struct('cost',cost,'x_fit',x_fit);

% ToDo
%   compute conditional & conjunctional tdr, etc.

