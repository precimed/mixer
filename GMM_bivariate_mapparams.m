function ov = cm3_bivariate_mapparams(iv,options) 

ngm = options.ngm;
if isfield(options,'mixtypevec')
  mixtypevec = options.mixtypevec;
  ivec_trait1 = (mixtypevec==1)|(mixtypevec==3);
  ivec_trait2 = (mixtypevec==2)|(mixtypevec==3);
else
  ivec_trait1 = true(1,ngm);
  ivec_trait2 = true(1,ngm);
end
ivec_trait12 = ivec_trait1&ivec_trait2;;

if isstruct(iv)
  ov = [];
  if ~isfield(options,'sig01'), ov = cat(2,ov,exp_amd(iv.sig01,0)); end
  if ~isfield(options,'sig02'), ov = cat(2,ov,exp_amd(iv.sig02,0)); end
  if ~isfield(options,'rho0'), ov = cat(2,ov,erf_amd(iv.rho0,0)); end
  if ~isfield(options,'pivec'), ov = cat(2,ov,logit_amd(iv.pivec,0)); end
  if ~isfield(options,'sig1vec'), ov = cat(2,ov,exp_amd(iv.sig1vec(ivec_trait1),0)); end
  if ~isfield(options,'sig2vec'), ov = cat(2,ov,exp_amd(iv.sig2vec(ivec_trait2),0)); end
  if ~isfield(options,'rhovec'), ov = cat(2,ov,erf_amd(iv.rhovec(ivec_trait12),0)); end
else
  ov = struct(); cnti = 1;
  if isfield(options,'sig01'), ov.sig01=options.sig01; else ov.sig01=exp_amd(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'sig02'), ov.sig02=options.sig02; else ov.sig02=exp_amd(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'rho0'), ov.rho0=options.rho0; else ov.rho0=erf_amd(iv(cnti),1); cnti=cnti+1; end
  if isfield(options,'pivec'), ov.pivec=options.pivec; else ov.pivec=logit_amd(iv(cnti-1+[1:ngm]),1); cnti=cnti+ngm; end
  if isfield(options,'sig1vec'), ov.sig1vec=options.sig1vec; else ov.sig1vec=zeros(1,ngm); ov.sig1vec(ivec_trait1)=exp_amd(iv(cnti-1+[1:sum(ivec_trait1)]),1); cnti=cnti+sum(ivec_trait1); end
  if isfield(options,'sig2vec'), ov.sig2vec=options.sig2vec; else ov.sig2vec=zeros(1,ngm); ov.sig2vec(ivec_trait2)=exp_amd(iv(cnti-1+[1:sum(ivec_trait2)]),1); cnti=cnti+sum(ivec_trait2); end
  if isfield(options,'rhovec'), ov.rhovec=options.rhovec; else ov.rhovec=zeros(1,ngm); ov.rhovec(ivec_trait12)=erf_amd(iv(cnti-1+[1:sum(ivec_trait12)]),1); cnti=cnti+sum(ivec_trait12); end
end

% ToDo
%   Allow for specification per Gaussian mixture component which traits it includes (1, 2, or both)
