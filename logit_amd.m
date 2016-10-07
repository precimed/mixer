function y = logit_amd(x,s)

if ~exist('s'), s = 0; end

minval = eps;
if s == 0
  x = minval+x*(1-2*minval);
  y = logit(x,s);
else
  y = logit(x,s);
  y = (y-logit(logit(minval,0),1))/(1-2*minval);
end

