function y = erf_amd(x,s)

if ~exist('s'), s = 1; end

minval = eps;
if s == 0
  x = x*(1-minval);
  y = erfinv(x);
else
  y = erf(x);
  y = y/(1-minval);
end

