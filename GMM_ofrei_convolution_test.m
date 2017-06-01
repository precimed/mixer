if 0  % 1D convolution
limit = 20; step=0.1; x = -limit:step:limit;

conv_n = @(x, n)op_power(x, @(a,b)conv(a, b, 'same'), n);

norm1pdf  = normpdf(x, 0, 1);
norm1hist = norm1pdf / sum(norm1pdf);

power = 25;

normNpdf  = normpdf(x, 0, sqrt(power));
normNhist = normNpdf / sum(normNpdf);

plot(x, conv_n(norm1hist, power), x, normNhist);



normNhist = norm1hist;
normNhist_exp = norm1hist;
sigma     = 1;

clf;

plot(x,normNhist, x, normNhist_exp);

normNhist     = conv(normNhist, norm1hist, 'same'); normNhist = normNhist / sum(normNhist);
sigma         = sigma + 1; 
normNhist_exp = normpdf(x, 0, sqrt(sigma)); normNhist_exp = normNhist_exp / sum(normNhist_exp);
end

if 1 % 2D convolution (yep it works :))
limit = 10; step=0.5; p = -limit:step:limit;
[x,y] = meshgrid(p, p);
A = [3 -1.5; -1.5 4]; pA = reshape(mvnpdf([x(:), y(:)], 0, A), size(x)); p = p ./ sum(p(:));
B = [2 1; 1 3];       pB = reshape(mvnpdf([x(:), y(:)], 0, B), size(x)); p = p ./ sum(p(:));
C = A+B;              pC = reshape(mvnpdf([x(:), y(:)], 0, C), size(x)); p = p ./ sum(p(:));
                      pC_conv = conv2(pA, pB, 'same');
D = [0 0; 0 1];       pD = reshape(mvnpdf([x(:), y(:)], 0, D), size(x)); p = p ./ sum(p(:));
imagesc(pC);

imagesc(conv2(pC, pC, 'same'))


end