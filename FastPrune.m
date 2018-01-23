function logpvec_pruned = FastPrune(logpvec, LDmat, ldthresh);

% FastPrune     Prune SNPs based on 1KG LD structure.
% logpvec,      -log10(P) vector 
% LDmat,        1KG LD matrix (in general sparse format)
% ldthresh(optional), LD thrshold(default, 0.2)
%
% reture pruned version of logpvec with prune SNPs set to nan.
%
% Note: 
% If length of logpvec not equal row number of LDmat you get WRONG result.
% The SNPs with largest -log10(P) will be kept in each block.
%
% logpvec_pruned = FastPrune(logpvec,LDmat,ldthresh);
%
if ~exist('ldthresh', 'var'), ldthresh = 0.2; end

[sv si] = sort(abs(logpvec),'descend'); 
si = si(isfinite(sv)); 
sv = sv(isfinite(sv));
prunevec = false(size(logpvec));
if islogical(LDmat)     % .....and persumably LDmat has been set with the desired threshold! E.g., r^2>0.8 all true, the rest false.
  for i = 1:length(si)
    if ~prunevec(si(i))
       prunevec(find(LDmat(:, si(i)))) = 1; % Get rid of ("prune") these SNPs.
       prunevec(si(i)) = 0;                 % Keep this SNP.
    end
  end
else
  for i = 1:length(si)
    if ~prunevec(si(i))
       prunevec(find(LDmat(:, si(i)) >= ldthresh)) = 1; % Get rid of these SNPs.
       prunevec(si(i)) = 0;                             % Keep this SNP.
    end
  end
end
logpvec_pruned = logpvec;
logpvec_pruned(prunevec) = NaN;
