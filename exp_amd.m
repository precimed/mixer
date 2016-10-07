function y = exp_amd(x,s)

if ~exist('s'), s = 1; end

minval = eps;
if s == 0
  x = minval+x;
  y = log(x);
else
  y = exp(x);
  y = (y-exp(log(minval)));
end

