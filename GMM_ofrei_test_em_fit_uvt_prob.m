zvec = -10:0.1:10;
pi1t  = 0.01;
sig0t = 1.2;
sig1t = 1.8;
phi  = @(sig)normpdf(zvec, 0, sig);
pvec = (1-pi1t)*phi(sig0t) + pi1t*phi(sig1t);

clf; hold on;
qqlim=10;
plot(-log10(cumsum(pvec/sum(pvec))),  -log10(normcdf(zvec,0,1)),'r');
plot([0 qqlim],[0 qqlim], 'k');
xlim([0 10]); ylim([0 15]);

KL=@(p,q)sum(p .* (log(p) - log(q)));
sKL=@(p,q)(KL(p,q)+KL(q,p));

clf; hold on;
sig0h = 1; sig1h = 2; pi1h = 0.01;
for i=1:10
    p0 = (1-pi1h) * phi(sig0h);
    p1 = pi1h * phi(sig1h);
    gamma = p1 ./ (p0 + p1);
    sig0h = sum(pvec .* (1-gamma) .* zvec.^2) / sum(pvec .* (1-gamma));
    sig1h = sum(pvec .* gamma .* zvec.^2) ./ sum(pvec .* gamma);
    pi1h = sum(pvec .* gamma)/sum(pvec);
    plot(i, pi1h, '*');
    %fprintf('%.3e\n', KL(pvec, (1-pi1h)*phi(sig0h) + pi1h*phi(sig1h)));
end



pi=1e-3; sig0t=1.00; sig1t=1.3; pvec = (1-pi)*phi(sig0t)+pi*phi(sig1t);
cost = @(x)sKL((1-x(1))*phi(x(2))+x(1)*phi(x(3)), pvec); % x = [pi1, sig0, sig1]
x = fminsearch(cost, [0.1 1 2], struct('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8))

%0.1303    1.1998    2.2984
%0.1297    1.2002    2.3016

%sigvec = 1.0 : 0.01 : 3.0;
%pivec  = 0 : 0.01 : 1.0;
%[sig, pi] = meshgrid(sigvec, pivec);
%kl = nan(size(sig));
%for i=1:numel(sig),
%    kl(i)= KL(pvec, (1-pi(i))*phi(1)+pi(i)*phi(sig(i)));
%end
%imagesc(kl);
%[~, i] = min(kl(:)); 
%fprintf('sig1=%.3e, pi1=%.3e\n', sig(i), pi(i));
%[~, id] = min(kl); sigvec(id)
%plot(kl)
