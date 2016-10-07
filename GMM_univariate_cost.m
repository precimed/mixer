function [cost statstruct] = GMM_univariate_cost(x,zvec,Hvec,Lvec,imat_prune,options,fignumoffset,bini)

if ~exist('options','var'), options = []; end
if ~exist('fignumoffset','var'), fignumoffset = 0; end
if ~exist('bini','var'), bini = []; end

x_struct = GMM_univariate_mapparams(x,options); 

v2struct(x_struct);

nsnp = length(zvec);
nprune = size(imat_prune,2);
pi0effvec = (1-pi1).^Lvec;
cost = 0;
for prunei = 1:nprune
  if pi1 < 1
    pi0effvec_tmp = pi0effvec(imat_prune(:,prunei));
    pi1effvec_tmp = 1-pi0effvec_tmp;
    pvec_tmp = pi0effvec_tmp.*normpdf(zvec(imat_prune(:,prunei)),0,sig0)+pi1effvec_tmp.*normpdf(zvec(imat_prune(:,prunei)),0,sqrt(sig0^2+Hvec(imat_prune(:,prunei))*sigb^2));
  else
    pvec_tmp = normpdf(zvec(imat_prune(:,prunei)),0,sqrt(sig0^2+Lvec(imat_prune(:,prunei)).*Hvec(imat_prune(:,prunei))*sigb^2));
  end
  cost = cost + sum(-log(pvec_tmp));
end

if nargout==1, return; end

fitstruct = struct();

% Add code to plot fit

