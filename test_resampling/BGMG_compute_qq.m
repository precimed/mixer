function [qqvec1 qqvec2 logpvals_qq zvals_qq] = BGMG_compute_qq(pdfmat_z,delvals)

zstep = delvals(2)-delvals(1);
pdfvec1 = rowvec(sum(pdfmat_z,2))*zstep;
pdfvec2 = rowvec(sum(pdfmat_z,1))*zstep;
[mv mi] = min(abs(delvals));
if mv~=0
  keyboard % Not implemented
else
  zvals_qq = delvals(mi:end);
  qqvec1 = flip(cumsum(flip([pdfvec1(mi) 2*pdfvec1(mi+1:end)]))*zstep); qqvec1 = qqvec1/qqvec1(1);
  qqvec2 = flip(cumsum(flip([pdfvec2(mi) 2*pdfvec2(mi+1:end)]))*zstep); qqvec2 = qqvec2/qqvec2(1);
end

logpvals_qq = -log10(2*normcdf(-abs(zvals_qq)));

%keyboard

