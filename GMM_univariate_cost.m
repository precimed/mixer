function [cost statstruct] = GMM_univariate_cost(x,zvec,Hvec,Lvec,imat_prune,options)

if ~exist('options','var'), options = []; end

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

keyboard

% Add code to plot fit

zvals_edges =  linspace(-15,15,31); zvals_mean = (zvals_edges(2:end)+zvals_edges(1:end-1))/2;

nsnp = length(zvec);
nprune = size(imat_prune,2);
pi0effvec = (1-pi1).^Lvec;
cdfvec_sum = 0; cntvec_sum = 0;
for prunei = 1:nprune
  pi0effvec_tmp = pi0effvec(imat_prune(:,prunei));
  pi1effvec_tmp = 1-pi0effvec_tmp;
  cdfmat_tmp = NaN(length(pi0effvec_tmp),length(zvals_edges));
  for  zi =  1:length(zvals_edges)
    cdfmat_tmp(:,zi) = pi0effvec_tmp.*normcdf(zvals_edges(zi),0,sig0)+pi1effvec_tmp.*normcdf(zvals_edges(zi),0,sqrt(sig0^2+Hvec(imat_prune(:,prunei))*sigb^2));
  end
  cdfvec_sum = cdfvec_sum  + mean(cdfmat_tmp,1);
  cntvec_tmp = histc([-zvec; zvec],zvals_edges); cntvec_tmp = [cntvec_tmp(1:end-1)]; 
  cntvec_sum = cntvec_sum + cntvec_tmp;
end
cntvec_pred = diff(cdfvec_sum)/cdfvec_sum(end);

figure; plot(zvals_mean,cntvec_sum/sum(cntvec_sum),zvals_mean,cntvec_pred/sum(cntvec_pred));
figure; semilogy(zvals_mean,cntvec_sum/sum(cntvec_sum),zvals_mean,cntvec_pred/sum(cntvec_pred)); % Doesn't seem to capture tails

