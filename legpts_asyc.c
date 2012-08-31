#include <mex.h>
#include <math.h>


#define PI 3.141592653589793
#define EPSILON 2.2204460492503131e-16
#define SQRTEPS 1.4901161193847661e-08
#define M 30
#define ITER_MAX 10

int asy1(double *nodes, double*weights, unsigned long int n)
{
    int s = n % 2;                       /* True if n is odd */
    unsigned long int n3, n4;
    int k, iter;
    double x, w, t, dt, C, theta, sint, f, fp, S, dS, dn = (double)n, p2, fn, fn5, nk, n5k, pn, jk;
    
    /* The constant out the front of the sum in evaluation (compute only once)
    // C = sqrt(pi/4)*gamma(n+1)/gamma(n+3/2)
    //  
    // Naive expansion    
    //    C = (4.0/PI);
    //    for ( k = 1; k <= n; k++ ) {
    //        C *= k/(k + 0.5);
    //    }
    //
    // Use Stirling's series */
    dS = -.125/dn;
    S = dS; 
    k = 1; 
    while ( fabs(dS/S) > EPSILON/100 ) {
        k += 1;
        dS *= -.5 * (k-1) / (k+1) / dn;
        S += dS;
    }
    double stirling[10] = {1., 1./12., 1./288., -139./51840., -571./2488320., 163879./209018880., 5246819./75246796800., -534703531./902961561600., -4483131259./86684309913600., 432261921612371./514904800886784000.};
    fn = 1; fn5 = 1; nk = dn; n5k = dn+.5;
    for ( k = 1 ; k < 10 ; k++ ) {
        fn += stirling[k]/nk;
        fn5 += stirling[k]/n5k;
        nk *= dn;
        n5k *= dn + .5;
    }
    C = exp(S)*sqrt(4.0/(dn+.5)/PI) * fn / fn5;  
    
    /* Approximations of roots near boundaries */
    double besselroots[32] = {2.404825557695773,
                            5.520078110286311,
                            8.653727912911013,
                            11.79153443901428,
                            14.93091770848779,
                            18.071063967910920,
                            21.211636629879269,
                            24.352471530749284,
                            27.493479132040250,
                            30.634606468431990,
                            33.775820213573560,
                            36.917098353664045,
                            40.058425764628240,
                            43.199791713176737,
                            46.341188371661815,
                            49.482609897397822,
                            52.624051841114984,
                            55.765510755019974,
                            58.906983926080954,
                            62.048469190227159,
                            65.189964800206866,
                            68.331469329856787,
                            71.472981603593752,
                            74.614500643701817,
                            77.756025630388066,
                            80.897555871137627,
                            84.039090776938195,
                            87.180629843641128,
                            90.322172637210500,
                            93.463718781944763,
                            96.605267950996279,
                            99.746819858680624};
    pn = dn*dn+dn+1./3.;
    
    /* Loop over each root */
    for ( k = (n+s)/2; k > 0; k-- ) {
        /* Asymptotic approximation of roots (in theta-space) */
        theta = PI*(4.0*k-1)/(4.0*n+2.0);
        n3 = n*n*n; n4 = n3*n; sint = sin(theta);
        x = (1.0 - (n-1.0)/(8*n3)-(39.0-28.0/(sint*sint))/(384.0*n4))*cos(theta);  
        t = acos(x);
        
        /* Initialise */
        dt = 1.0; 
        iter = 1;
        
        /* Skip x near the ends */
        if (x>1-4/(log(n)*n)) { 
            if (k < 30) { 
                jk = besselroots[k-1];
                t = jk/sqrt(pn)*(1.-(jk*jk-1.)/(pn*pn)/360.);
            }
            dt = 0.;
        }
        
        /* Newton iteration for roots */
        while (fabs(dt) > EPSILON) {
            feval(n, t, C, &f, &fp);        /* Evaluate at current approximation */
            dt = f/fp;                      /* Newton update */
            t += dt;                        /* Update t */
            if (iter++ > ITER_MAX) { break; } 
        }
        if (dt == 0) {
            feval(n, t, C, &f, &fp);
        }
        
        /* Convert back to x-space */
        x = cos(t);
        w = 2./(fp*fp);        
        
        /* Store nodes and weights */
        nodes[n-k] = x;
        nodes[k-1] = -x;
        weights[n-k] = w;
        weights[k-1] = w;
    }
    
    /* Enforce x = 0 for odd n */ 
    if (s == 1) {
        nodes[(n+s)/2-1] = 0;
    }

    return 0;
}

int feval(int n, double theta, double C, double* f, double* fp)
{
    double sinT = sin(theta), cotT = 1/tan(theta);
    double alpha, cosA, sinA, denom, df, dfp;
    int m;
    
    /* m = 0 */
    alpha = (n + 0.5)*theta - 0.25*PI;
    denom = sqrtl(2.0*sinT);
    cosA = cos(alpha);
    sinA = sin(alpha);
    *f = C * cosA/denom;
    *fp = C * ( 0.5*(cosA*cotT+sinA) + n*sinA) / denom;
    
    /* Loop over m until coefficient is sufficiently small */
    for ( m = 1; m <= 30; m++ ) {
        C *= (2.0*m-1.0)/(2.0*m)*(m-0.5)/(n+m+0.5);
        alpha += theta - 0.5*PI;
        denom *= 2.0*sinT;
        cosA = cos(alpha);
        sinA = sin(alpha);
        
        df = C * cosA/denom;
        dfp = C * ( (m+0.5)*(cosA*cotT+sinA) + n*sinA) / denom;
        
        *f += df;
        *fp += dfp;
        
        if (fabs(df) + fabs(dfp) < EPSILON){ break; }
    }
    
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *result_roots,*result_ders;
  int n;
  /*  check for proper number of arguments */
  if(nrhs!=1) 
    mexErrMsgTxt("One input required.");
  
  /*  create a pointer to the input matrix x  */
  n = mxGetScalar(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,n, mxREAL);
  
  /*  create a C++ pointer to a copy of the output matrix */
  result_roots = mxGetPr(plhs[0]);
  result_ders = mxGetPr(plhs[1]);
  
  /*  call the C++ subroutine */
  asy1(result_roots,result_ders,n); 
  
  return;
}