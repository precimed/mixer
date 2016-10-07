function y = mnpdf(x,p)
%MNPDFLN Log of multinomial probability density function (log(pdf)).
%   Y = MNPDFLN(X,PROB) returns log of the pdf for the multinomial distribution with
%   probabilities PROB, evaluated at each row of X. X and PROB are M-by-K
%   matrices or 1-by-K vectors, where K is the number of multinomial bins or
%   categories.  Each row of PROB must sum to one, and the sample sizes for
%   each observation (rows of X) are given by the row sums SUM(X,2).  Y is
%   a M-by-1 vector, and MNPDF computes each row of Y using the corresponding
%   rows of the inputs, or replicates them if needed.
%
%   Example:
%    Generate a random vector with sample size 20 and probabilities P and
%    compute the multinomial pdf of X with probabilities P
%    P=[0.3,0.7];
%    X=mnrnd(20,P);
%    Y=mnpdf(X,P);
%
%   See also MNRFIT, MNRVAL, MNRND

%   Copyright 2006-2011 The MathWorks, Inc. (modified by AMD)
%   $Revision: 1.1.8.3 $  $Date: 2011/05/09 01:25:54 $

narginchk(2,2);

% If p is a column that can be interpreted as a single vector of MN
% probabilities (non-empty and sums to one), transpose it to a row.
% Otherwise, treat as a matrix with one category.
if size(p,2)==1 && size(p,1)>1 && abs(sum(p,1)-1)<=size(p,1)*eps(class(p))
    p = p';
    if size(x,2)==1 && size(x,1)>1  % transpose x if it is a non-empty column vector
        x = x';
    end
end

[m,k] = size(x);
if k < 1
    error(message('stats:mnpdf:NoCategories'));
end
n = sum(x,2);

[mm,kk] = size(p);
if kk ~= k
    error(message('stats:mnpdf:ColSizeMismatch'));
elseif mm == 1 % when m > 1
    p = repmat(p,m,1);
elseif m == 1
    m = mm;
    x = repmat(x,m,1);
    n = repmat(n,m,1);
elseif mm ~= m
    error(message('stats:mnpdf:RowSizeMismatch'));
end

outClass = superiorfloat(n,p);

xBad = any(x < 0 | x ~= round(x), 2);
pBad = any(p < 0 | 1 < p, 2) | abs(sum(p,2)-1) > size(p,2)*eps(class(p));
nBad = n < 0 | round(n) ~= n;

xPos = (x > 0); % avoid 0*log(0), but let (pi==0) & (y>0) happen
xlogp = zeros(m,k,outClass);
xlogp(xPos) = x(xPos).*log(p(xPos)); % may be bad if p was, those won't be used
xlogp = sum(xlogp, 2);

% Initialize to return zero for noninteger or negative x
y = zeros(m,1,outClass); % y(xBad) = 0;

t = ~(xBad | pBad | nBad);
y(t) = gammaln(n(t,:) + 1) - sum(gammaln(x(t,:) + 1), 2) + xlogp(t,:);

% Return NaN for invalid parameter values
y(pBad | nBad) = NaN;

