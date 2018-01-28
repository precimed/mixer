function [zmat delmat_total delmat_indep delmat_pleio] = BGMG_generate_random_stats(params,vals_r2,count_r2,nobs,Neff,H)

% If Neff or H not specified, assume sigb is in units of del (already scaled)
if ~exist('Neff','var') | isempty(Neff), Neff = 1; end 
if ~exist('H','var') | isempty(H), H = 1; end
if ~exist('dispflag','var') | isempty(dispflag), dispflag = false; end

nsnp = sum(count_r2);

randmat = rand(nobs,nsnp);
betamat1_indep = zeros(nobs,nsnp); betamat2_indep = zeros(nobs,nsnp);
betamat1_pleio = zeros(nobs,nsnp); betamat2_pleio = zeros(nobs,nsnp);

tmp1 = randmat<params.pi_1;
tmp2 = randmat<params.pi_1+params.pi_2;
tmp3 = randmat<params.pi_1+params.pi_2+params.pi_3;
ind1 = find(tmp1);
ind2 = find(tmp2&~tmp1);
ind3 = find(tmp3&~tmp2);

C0 = [params.sig0_1^2 params.sig0_1*params.sig0_2*params.rho0; params.sig0_1*params.sig0_2*params.rho0 params.sig0_2^2];
 
betamat1_indep(ind1) = params.sigb_1*randn(1,length(ind1));
betamat2_indep(ind2) = params.sigb_2*randn(1,length(ind2));
if length(ind3)>0
  C = [params.sigb_1^2 params.sigb_1*params.sigb_2*params.rhob; params.sigb_1*params.sigb_2*params.rhob params.sigb_2^2];
  tmp = mvnrnd([0 0],C,length(ind3));
  betamat1_pleio(ind3) = tmp(:,1);
  betamat2_pleio(ind3) = tmp(:,2);
end

r2vec = NaN(1,nsnp);
i0 = 0;
for bini = 1:length(vals_r2)
  r2vec(i0+[1:count_r2(bini)]) = vals_r2(bini);
  i0 = i0+count_r2(bini);
end

wmat = repmat(sqrt(Neff*H*r2vec),[nobs 1]);
delmat_indep = [sum(betamat1_indep.*wmat,2) sum(betamat2_indep.*wmat,2)];
delmat_pleio = [sum(betamat1_pleio.*wmat,2) sum(betamat2_pleio.*wmat,2)];
delmat_total = delmat_indep + delmat_pleio;
zmat = delmat_total + mvnrnd([0 0],C0,size(delmat_total,1));

