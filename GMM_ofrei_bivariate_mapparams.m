function ov = GMM_ofrei_bivariate_mapparams(iv, options) 

if ~exist('options', 'var'), options=struct(); end;
if isstruct(iv)
  ov = [];
  if ~isfield(options,'sigma0'), ov = cat(2,ov,exp_amd(iv.sigma0,0)); end           % vector, 2 components - one per phenotype
  if ~isfield(options,'rho0'), ov = cat(2,ov,erf_amd(iv.rho0,0)); end;              % scalar
  if ~isfield(options,'sigma_beta'), ov = cat(2,ov,exp_amd(iv.sigma_beta,0)); end   % vector, 2 components - one per phenotype
  if ~isfield(options,'rho_beta'), ov = cat(2,ov,erf_amd(iv.rho_beta,0)); end;          % scalar
  if ~isfield(options,'pivec'), ov = cat(2,ov,logit_amd(iv.pivec,0)); end           % vector, 3 components - 1, 2, both.
else
  ov = struct(); cnti = 1;
  if isfield(options,'sigma0'), ov.sigma0=options.sigma0; else ov.sigma0=exp_amd(iv(cnti:(cnti+1)),1); cnti=cnti+2; end
  if isfield(options,'rho0'), ov.rho0=options.rho0; else ov.rho0 = erf_amd(iv(cnti), 1); cnti=cnti+1; end;              
  if isfield(options,'sigma_beta'), ov.sigma_beta=options.sigma_beta; else ov.sigma_beta=exp_amd(iv(cnti:(cnti+1)),1); cnti=cnti+2; end
  if isfield(options,'rho_beta'), ov.rho_beta=options.rho_beta; else ov.rho_beta = erf_amd(iv(cnti), 1); cnti=cnti+1; end;              
  if isfield(options,'pivec'), ov.pivec=options.pivec; else ov.pivec=logit_amd(iv(cnti:(cnti+2)),1); cnti=cnti+3; end
end

end