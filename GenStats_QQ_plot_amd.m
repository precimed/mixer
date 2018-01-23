function [qqvec hv_logp] = GenStats_QQ_plot(logpvec,hv_z)

if exist('hv_z','var') 
  hv_logp = -log10(2*normcdf(-abs(hv_z)));
else
  hv_logp = linspace(0,200,10000);
end
%qqvec = NaN(size(logpvec));
qqvec = NaN(size(hv_z));
hc = hist(logpvec,hv_logp);

chc = cumsum(hc)/sum(hc);  % Goes from 0 to 1. Flip it around: cdf_arr==(1-chc) goes from 1 to 0. cdf_arr(1) gives the proportion or fraction of p-values less than or equal to 1: cdf_arr(1)=1.
%qqvecOld = -log10(1-chc); % The idea is that qqvec(i) = -log10(of the proportion of SNPs with p-values as significant as or more significant than the threshold given by pthresh(i))
                           % where pthresh(i) is given by  hv_logp(i) = -log10(pthresh(i)),  i.e.,  pthresh(i) = 10^-hv_logp(i);  this is the basic definition of a qq plot.
                           % Note: hv_logp(1) = 0,  i.e,  pthresh(1) = 1;
                           % Define prop(i)==1-chc(i) = the proportion of SNPs with p-values as significant as or more significant than the threshold given by pthresh(i);
                           % Since pthresh(1) = 1, and no SNP can have a p-value greater than 1 (!),
                           % then prop(1)=1-chc(1) = 1,  i.e., qqvec(1) = 0,  which of course requires  chc(1) = 0;
                           % However, "chc(1) = 0" will not be true if there are SNPs with z-scores identically equal to 0. So,...

prop = zeros(size(hv_z));
shc = sum(hc);
prop(1) = (shc)/shc;
%prop(2) = (shc-hc(1))/shc;        % etc...
%prop(3) = (shc-hc(1)-hc(2))/shc;  % etc...
prop(2:end) = (shc-cumsum(hc(1:end-1)))/shc;
qqvec = -log10(prop);   % NOTE: qqvec(1)==0 and  qqvec(2:end) = qqvecOld(1:end-1);   Just a shift.  But correct!

if nargout==0
  plot([0 7],[0 7],'k:',qqvec,hv_logp,'LineWidth',2);
  xlim([0 7]);
  ylim([0 10]);
  h=xlabel('Empirical -log_1_0(q)');
  set(h,'FontSize',20);
  h=ylabel('Nominal -log_1_0(p)');
  set(h,'FontSize',20);
end

% chc   = cumsum(hc)/sum(hc) goes from 0 to 1. Flip it around:
% cdf_arr = (1-chc)            goes from 1 to 0.
% qqvec = -log10(cdf_arr)      goes from 0 to infinity (approximately 200).
% Note that cdf_arr is indeed a cumulative distribution function! Therefore, cdf_arr(i) can be interpreted as a p-value!
% cdf_arr(1) gives the fraction or proportion of data p-values LESS THAN or equal to 1: cdf_arr(1)=1.
% cdf_arr(i) gives the fraction or proportion of data p-values LESS THAN or equal to 10^(-hv_logp(i)).
% 
% Let p(i)=2*normcdf(-abs(hv_z(i))),  the normal two-tailed p-value corresponding to |z|=hv_z(i).
% 
% IF SNP z-scores were distributed normally,
% then the PROPORTION of SNPs with -abs(z) more negative than -abs(hv_z(i)) would be just pi, i.e.,
% then the PROPORTION of SNPs with p-value LESS THAN pi=2*normcdf(-abs(hv_z(i))) would be just pi,
% i.e., then would have  cdf_arr(i) = 2*normcdf(-abs(hv_z(i)))  i.e.,  qqvec(i) = hv_logp(i)  for all i -- a straight line!
% Plot  x=-log10(cdf_arr)  versus  y=hv_logp.
% x(1)=0; y(1)=0.
