function logpvec_pruned = GenStats_FastPrune(logpvec,LDmat);

%[sv si] = sort(logpvec,'descend'); si = si(isfinite(sv)); sv = sv(isfinite(sv));
[sv si] = sort(abs(logpvec),'descend'); si = si(isfinite(sv)); sv = sv(isfinite(sv));
prunevec = false(size(logpvec));
for i = 1:length(si)
%  if mod(i,10000)==0, fprintf(1,'%d of %d\n',i,length(si)); end
  if ~prunevec(si(i))
     prunevec(find(LDmat(:,si(i)))) = 1; 
     prunevec(si(i)) = 0;
  end
end
logpvec_pruned = logpvec;
logpvec_pruned(prunevec) = NaN;
