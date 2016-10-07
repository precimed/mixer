function p = BVNcdf(X,mu,omega)
%             X:     Nx2 matrix
%             mu:    1(or N)x2 vector (matrix) of means
%             omega: 2x2 covariance matrix.
if isempty(mu),     mu = [0,0]; end
if length(omega)<2, omega(1,2) = omega; omega(2,1) = omega(1,2); omega(1,1) = 1; omega(2,2) = 1; end
%if sum(diag( omega )' > sum( abs( omega - diag(diag(omega)))))~=2 ||  all(all(omega==omega.'))==0 % Test for Positive definite & symmetric covariance matrix
try
  if sum(eig(omega)<0)>0 % Test for Positive definite & symmetric covariance matrix (Edited by AMD)
      disp('The covariance matrix is not positive definite and/or symmetric! Check the inputs carefully!');
      keyboard
%      error('The covariance matrix is not positive definit and/or symmetric! Check the inputs carefully!');
  end
catch
  keyboard
end
dh=(X(:,1)-mu(:,1))/sqrt(omega(1,1));dk=(X(:,2)-mu(:,2))/sqrt(omega(2,2));
r = omega(1,2)/sqrt(omega(1,1)*omega(2,2));
p = NaN(size(dh));
p(dh == Inf & dk==Inf)   = 1;
p(dk==Inf)               = phid(dh(dk==Inf));
p(dh == Inf)             = phid(dk(dh == Inf));
p(dh == -Inf | dk==-Inf) = 0;
ind = (dh>-Inf&dh<Inf&dk>-Inf&dk<Inf);
if sum(ind)>0 %--> find P(x1<dh,x2<dk,r)
    p(ind)               = BVNcdfsub(-dh(ind),-dk(ind),r);
end
end
function p=BVNcdfsub(dh,dk,r)
if abs(r) < 0.3, lg = 3;
    %       Gauss Legendre points and weights, n =  6
    w = [0.1713244923791705 0.3607615730481384 0.4679139345726904];
    x = [0.9324695142031522 0.6612093864662647 0.2386191860831970];
elseif abs(r) < 0.75,  lg = 6;
    %       Gauss Legendre points and weights, n = 12
    w = [.04717533638651177 0.1069393259953183 0.1600783285433464...
        0.2031674267230659 0.2334925365383547 0.2491470458134029];
    x = [0.9815606342467191 0.9041172563704750 0.7699026741943050...
        0.5873179542866171 0.3678314989981802 0.1252334085114692];
else,  lg = 10;
    %       Gauss Legendre points and weights, n = 20
    w = [.01761400713915212 .04060142980038694 .06267204833410906...
        .08327674157670475 0.1019301198172404 0.1181945319615184...
        0.1316886384491766 0.1420961093183821 0.1491729864726037 0.1527533871307259];
    x = [0.9931285991850949 0.9639719272779138 0.9122344282513259...
        0.8391169718222188 0.7463319064601508 0.6360536807265150...
        0.5108670019508271 0.3737060887154196 0.2277858511416451 0.07652652113349733];
end
dim1=ones(size(dh,1),1); dim2=ones(1,lg);
hk = dh.*dk; bvn = dim1*0;
if abs(r) < 0.925, hs = ( dh.*dh + dk.*dk )/2; asr = asin(r);
    sn1=sin( asr*( 1 - x )/2 ); sn2= sin( asr*( 1 + x )/2 );
    bvn = sum((dim1*w).*exp( ( (dim1*sn1).*(hk*dim2) - hs*dim2 )./( 1 - dim1*(sn1.^2) ) ) ...
        + (dim1*w).*exp( ( (dim1*sn2).*(hk*dim2) - hs*dim2 )./( 1 - dim1*(sn2.^2) ) ),2)...
        *asr/( 4*pi ) + phid(-dh).*phid(-dk);
else, twopi = 2*pi; if r < 0, dk = -dk; hk = -hk; end
    if abs(r) < 1
        as = ( 1 - r )*( 1 + r ); a = sqrt(as); bs = ( dh - dk ).^2;
        c = ( 4 - hk )/8 ; d = ( 12 - hk )/16; asr = -( bs./as + hk )/2;
        ind = asr > -100;
        bvn(ind) = a*exp(asr(ind)).*( 1 - (c(ind).*(bs(ind)-as)).*(1-d(ind).*bs(ind)/5)/3 + (c(ind).*d(ind)).*as.^2/5 );
        ind = hk > -100; b = sqrt(bs); sp = sqrt(twopi)*phid(-b/a);
        bvn(ind) = bvn(ind) - (exp(-hk(ind)/2).*sp(ind)).*b(ind).*( 1 - c(ind).*bs(ind).*( 1 - d(ind).*bs(ind)/5 )/3 );
        a = a/2;
        for is = -1 : 2 : 1, xs = ( a + a*is*x ).^2;
            rs = sqrt( 1 - xs ); asr1 = -( (bs*dim2)./(dim1*xs) + hk*dim2 )/2;
            ind1 = (asr1 > -100) ;
            sp1 = ( 1 + (c*dim2).*(dim1*xs).*( 1 + (d*dim2).*(dim1*xs) ) );
            ep1 = exp( -(hk*dim2).*( 1 - dim1*rs )./( 2*( 1 + dim1*rs ) ) )./(dim1*rs);
            bvn = bvn + sum(a.*(dim1*w).*exp(asr1.*ind1).*( ep1.*ind1- sp1.*ind1 ),2);
        end
        bvn = -bvn/twopi;
    end
    if r > 0, bvn =  bvn + phid( -max( dh, dk ) );
    elseif r < 0, bvn = -bvn + max( 0, phid(-dh)-phid(-dk) ); end
end, p = max( 0, min( 1, bvn ) );
end

function p = phid(z)
p = erfc( -z/sqrt(2) )/2; % Normal cdf
end

