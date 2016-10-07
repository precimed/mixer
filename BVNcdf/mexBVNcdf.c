#include "stdio.h"
#include "mex.h"
#include "matrix.h"
#include "math.h"
#define Inf     1e+20
#define pi      3.14159265358979323846
#define twopi   2*pi
#define MaxQuad 10

// double phid(double x){
//     // Fejl er ca. 1e-7... Ret meget, men 300 gange hurtigere end matlab. Kan være man lige skulle forbedre denne del lidt!
//     // constants
//     double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741,a4 = -1.453152027,a5 = 1.061405429,p = 0.3275911;
//     double t,y;
//     // Save the sign of x
//     int sign = 1;
//     if (x < 0)
//         sign = -1;
//     x = fabs(x) / sqrt(2.0);
//     
//     // A&S formula 7.1.26
//     t = 1.0 / (1.0 + p*x);
//     y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * exp(-x*x);
//     return 0.5 * (1.0 + sign*y);
// }

double phid(double x){
    // Based on "BETTER APPROXIMATIONS TO CUMULATIVE NORMAL FUNCTIONS" by Craeme West, 2004--> double precision
    double XAbs = fabs(x),build,Cumnorm;
    if (XAbs > 37){
        return 0;
    } else {
        if (XAbs < 7.07106781186547){
            build = 3.52624965998911E-02 * XAbs + 0.700383064443688;
            build = build * XAbs + 6.37396220353165;
            build = build * XAbs + 33.912866078383;
            build = build * XAbs + 112.079291497871;
            build = build * XAbs + 221.213596169931;
            build = build * XAbs + 220.206867912376;
            Cumnorm = exp(-XAbs*XAbs / 2) * build;
            build = 8.83883476483184E-02 * XAbs + 1.75566716318264;
            build = build * XAbs + 16.064177579207;
            build = build * XAbs + 86.7807322029461;
            build = build * XAbs + 296.564248779674;
            build = build * XAbs + 637.333633378831;
            build = build * XAbs + 793.826512519948;
            build = build * XAbs + 440.413735824752;
            Cumnorm = Cumnorm / build;
        } else {
            build = XAbs + 0.65;
            build = XAbs + 4 / build;
            build = XAbs + 3 / build;
            build = XAbs + 2 / build;
            build = XAbs + 1 / build;
            Cumnorm = exp(-XAbs*XAbs / 2) / build / 2.506628274631;
        }
    }
    if (x > 0){
        Cumnorm = 1 - Cumnorm;}
    return Cumnorm;
}

double BVNcdf(double dh,double dk,double r ){
    // Based on Matlab code from Alan Getz
    double w[MaxQuad],x[MaxQuad],hk,bvn,hs,asr,sn,as,a,bs,c,d,b,sp,xs,rs,ep;
    int lg,i,j,is;

    if        (dh >=  Inf && dk >=  Inf)
        return 1;
    else if   (dh <=  -Inf || dk <=  -Inf)
        return 0;
    else if   (dk >= Inf)
        return phid(dh);
    else if (dh >= Inf)
        return phid(dk);
    else {
        dh=-dh; dk=-dk;
        if (fabs(r) < 0.3){             lg = 3;
        //       Gauss Legendre points and weights, n =  6
        w[0] = 0.1713244923791705; w[1] = 0.3607615730481384; w[2] = 0.4679139345726904;
        x[0] = 0.9324695142031522; x[1] = 0.6612093864662647; x[2] = 0.2386191860831970;
        } else if (fabs(r) < 0.75){     lg = 6;
        //       Gauss Legendre points and weights, n = 12
        w[0] = .04717533638651177; w[1] = 0.1069393259953183; w[2] = 0.1600783285433464;
        w[3] = 0.2031674267230659; w[4] = 0.2334925365383547; w[5] = 0.2491470458134029;
        x[0] = 0.9815606342467191; x[1] = 0.9041172563704750; x[2] = 0.7699026741943050;
        x[3] = 0.5873179542866171; x[4] = 0.3678314989981802; x[5] = 0.1252334085114692;
        } else {                        lg = 10;
        //       Gauss Legendre points and weights, n = 20
        w[0] = .01761400713915212; w[1] = .04060142980038694; w[2] = .06267204833410906;
        w[3] = .08327674157670475; w[4] = 0.1019301198172404; w[5] = 0.1181945319615184;
        w[6] = 0.1316886384491766; w[7] = 0.1420961093183821; w[8] = 0.1491729864726037; w[9] = 0.1527533871307259;
        x[0] = 0.9931285991850949; x[1] = 0.9639719272779138; x[2] = 0.9122344282513259;
        x[3] = 0.8391169718222188; x[4] = 0.7463319064601508; x[5] = 0.6360536807265150;
        x[6] = 0.5108670019508271; x[7] = 0.3737060887154196; x[8] = 0.2277858511416451; x[9] = 0.07652652113349733;
        }
        hk = dh*dk; bvn = 0;
        if (fabs(r) < 0.925){
            hs = ( dh*dh + dk*dk )/2; asr = asin(r); // findes asin i math.h? Det ser det ud til
            for(i=0;i<lg;i++){
                sn = sin( asr*( 1 - x[i] )/2 );
                bvn = bvn + w[i]*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
                sn = sin( asr*( 1 + x[i] )/2 );
                bvn = bvn + w[i]*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
            }
            bvn = bvn*asr/( 4*pi );
            bvn = bvn + phid(-dh)*phid(-dk);
        } else {
            if (r < 0){ dk = -dk; hk = -hk; }
            if (fabs(r) < 1){
                as = ( 1.0 - r )*( 1.0 + r ); a = sqrt(as); bs = ( dh - dk )*( dh - dk );
                c = ( 4.0 - hk )/8.0 ; d = ( 12.0 - hk )/16.0; asr = -( bs/as + hk )/2.0;
                if (asr > -100){ bvn = a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 ); }
                if (hk > -100){
                    b = sqrt(bs); sp = sqrt(twopi)*phid(-b/a);
                    bvn = bvn - exp(-hk/2.0)*sp*b*( 1.0 - c*bs*( 1.0 - d*bs/5.0 )/3.0 );
                }
                a = a/2.0;
                for(i=0;i<lg;i++){
                    for(j=0;j<=1;j++){
                        is = -1*(j==0)+1*(j==1);
                        xs = ( a + a*is*x[i] )*( a + a*is*x[i] );
                        rs = sqrt( 1.0 - xs ); asr = -( bs/xs + hk )/2.0;
                        if (asr > -100){
                            sp = ( 1.0 + c*xs*( 1.0 + d*xs ) );
                            ep = exp( -hk*( 1.0 - rs )/( 2.0*( 1.0 + rs ) ) )/rs;
                            bvn = bvn + a*w[i]*exp(asr)*( ep - sp );
                        }
                    }
                }
                bvn = -bvn/twopi;
            }
            if (r > 0){ bvn =  bvn + phid( -max( dh, dk ) ); }
            else if (r < 0){ bvn = -bvn + max( 0, phid(-dh)-phid(-dk) ); }
        }
        return max( 0, min( 1, bvn ) );
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *X, *mu,*omega,*P,dh,dk,r;
    mwSize dim1;

    int i;
    //input (Parameters): 
    X = mxGetPr(prhs[0]); // Pointer matrix of points to evaluate
    mu = mxGetPr(prhs[1]); // Pointer matrix of points to evaluate
    omega = mxGetPr(prhs[2]); // Pointer matrix of points to evaluate
    // Dimensions of X
    dim1 = mxGetM(prhs[0]);

    // Allocate memory for the pointers, used in the calculations
    plhs[0] = mxCreateDoubleMatrix(dim1, 1, mxREAL);
    P = mxGetPr(plhs[0]);
    // Calculate the commulative probability
    for(i=0;i<dim1;i++){
        dh = (X[i]-mu[0])/sqrt(omega[0]);
        dk = (X[i+dim1]-mu[1])/sqrt(omega[3]);
        r  =  omega[1]/sqrt(omega[0]*omega[3]);
        P[i] = BVNcdf(dh,dk,r);
    }
}