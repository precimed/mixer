function ov = GMM_univariate_mapparams(iv,options) 

if isstruct(iv)
  ov = [];
  if ~isfield(options,'pi1'), ov = cat(2,ov,logit_amd(iv.pi1,0)); end 
  if ~isfield(options,'sig0'), ov = cat(2,ov,exp_amd(iv.sig0,0)); end 
  if ~isfield(options,'sigb'), ov = cat(2,ov,exp_amd(iv.sigb,0)); end 
else
  ov = struct(); cnti = 1;
  if isfield(options,'pi1'), ov.pi1=options.pi1; else ov.pi1=logit(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'sig0'), ov.sig0=options.sig0; else ov.sig0=exp_amd(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'sigb'), ov.sigb=options.sigb; else ov.sigb=exp_amd(iv(cnti),1); cnti=cnti+1; end
end
