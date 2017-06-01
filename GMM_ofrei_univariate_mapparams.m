function ov = GMM_ofrei_univariate_mapparams(iv, options) 

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
