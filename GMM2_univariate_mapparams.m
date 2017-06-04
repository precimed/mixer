function ov = GMM2_univariate_mapparams(iv, options) 

if ~exist('options', 'var'), options=struct(); end;
if isstruct(iv)
  ov = [];
  if ~isfield(options,'sigma0'), ov = cat(2,ov,exp_amd(iv.sigma0,0)); end
  if ~isfield(options,'sigma_beta'), ov = cat(2,ov,exp_amd(iv.sigma_beta,0)); end
  if ~isfield(options,'pivec'), ov = cat(2,ov,logit_amd(iv.pivec,0)); end
else
  ov = struct(); cnti = 1;
  if isfield(options,'sigma0'), ov.sigma0=options.sigma0; else ov.sigma0=exp_amd(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'sigma_beta'), ov.sigma_beta=options.sigma_beta; else ov.sigma_beta=exp_amd(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'pivec'), ov.pivec=options.pivec; else ov.pivec=logit_amd(iv(cnti),1); cnti=cnti+1; end
end

end

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

end

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

end

function y = logit(x,invflag)

    % logit transformation
    %
    % Input:
    % ------
    % x,            input value
    % invflag,      logical, inverse?
    %
    % Return:
    % ------
    % y,            transformed value

    if exist('invflag','var') & invflag
      y = exp(x)./(1+exp(x));
    else
      y = log(x./(1-x));
    end

end