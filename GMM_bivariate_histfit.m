function resultStruct = GMM_bivariate_histfit(sumstats,params,options,x0)

ngm = options.ngm;

%[cost fitstruct] = GMM_bivariate_histcost(x0,sumstats,options); % Make plots before fitting

costfun = @(x)GMM_bivariate_histcost(x,sumstats,params,options);

[x_fit_stochastic cost] = fminsearch_stochastic(costfun,x0,statset('MaxIter',1000,'MaxFunEvals',1000,'Display','iter'));
GMM_bivariate_mapparams(x_fit_stochastic,options) % Print out x_fit as struct

[x_fit cost] = fminsearch(costfun,x_fit_stochastic,statset('MaxIter',5000,'MaxFunEvals',5000,'Display','iter'));
%[cost fitstruct] = GMM_bivariate_histcost(x_fit,sumstats,params,options); % Make plots after fitting
GMM_bivariate_mapparams(x_fit,options) % Print out x_fit as struct

%resultStruct = struct('cost',cost,'x_fit',x_fit,'fitstruct',fitstruct);
resultStruct = struct('cost',cost,'x_fit',x_fit);

% ToDo
%   compute conditional & conjunctional tdr, etc.

