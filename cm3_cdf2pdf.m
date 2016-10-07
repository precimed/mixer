function pdfmat = cm3_cdf2pdf(cdfmat)

pdfmat = zeros(size(cdfmat,1)-1,size(cdfmat,2)-1);
for i = 1:size(pdfmat,1)
  for j = 1:size(pdfmat,2)
    pdfmat(i,j) = (cdfmat(i+1,j+1)-cdfmat(i,j))-(cdfmat(i+1,j)-cdfmat(i,j))-(cdfmat(i,j+1)-cdfmat(i,j));
  end
end
