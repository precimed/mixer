function [x_fit cost] = fminsearch_stochastic(costfun, x0, options)

if ~exist('options','var'), options = struct(); end
niter = getfielddefault(options,'MaxIter',1000);
randscale = getfielddefault(options,'randscale',0.1);

x_fit = x0;
mincost = inf;
for iter = 0:niter
%  x = x_fit + (iter>0&iter<niter)*randscale*randn(size(x_fit));
  x = x_fit + (iter>0&iter<niter)*exprnd(randscale,size(x_fit)).*sign(rand(size(x_fit))-0.5); % Use exponential distribution to escape local minima
  cost = costfun(x);
  if cost<mincost
    mincost = cost;
    x_fit = x;
    %fprintf(1,'*** iter=%d of %d: cost=%f\n',iter,niter,cost);
  end
end

% ToDo
%   use random distribution with heavier tails (e.g., dual exponential)
%   adapt step size based on history of decrease in cost as function of step size
