% Evaluate gain:
clear all; clc
%mex mexBVNcdf.c  % Uncomment to re-compile mexBVNcdf.c. Else, using mexBVNcdf.mexw64
sigma1 = 2;
sigma2 = 3;
rho = .9;

% Vektor version
X = (0:.1:1)';     
x = [X,-X];
for rho=.1:.01:.99
tic
nr1=BVNcdf(x,[],[sigma1,rho;rho,sigma2]);       % Matlab implementaion
slut1=toc;
tic
nr2=mvncdf(x,[],[sigma1,rho;rho,sigma2]);       % Matlab build in
slut2=toc;
tic
nr3 = mexBVNcdf(x,[0 0],[sigma1,rho;rho,sigma2]);   %mex-implementation. Acuracy is not that good on this version of the mex. This is due to the rather simple implementation of the univariate normal cdf
slut3=toc;
fprintf('FINAL (rho=%1.2f): %2.2f times faster than Matlab, Max(abs(error))=%2.5g\n',rho,slut2/slut3,max(abs(nr3-nr2)));
end  
