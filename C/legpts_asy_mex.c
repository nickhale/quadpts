/*******************************************************************************
 * This file is part of QUADPTS.
 * Copyright (c) 2012 Nick Hale and Alex Townsend
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

#include <mex.h>
#include <math.h>

#define PI 3.141592653589793
#define PI_HI  3.141592741012573242187
#define PI_LO -0.00000008742278000372485127293903
#define EPSILON 2.2204460492503131e-16
#define SQRTEPS 1.4901161193847661e-08
#define M 30
#define ITER_MAX 10

int asy1(double *nodes, double*weights, unsigned long int n);
int feval(int n, double theta, double C, double* f, double* fp, int k);
int feval_asy1(int n, double theta, double C, double* f, double* fp, int k);
int feval_asy2(int n, double t, double* f, double* fp);
int mycosA(int n, double theta, double* cosA, double* sinA, int k);
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static double bessj0( double x );
static double bessj1( double x );
double clenshaw( double* c, int n, double x);
double chebbessj0( double x );
double chebbessj1( double x );

int asy(double *nodes, double*weights, unsigned long int n)
{
    int s = n % 2;                       /* True if n is odd */
    unsigned long int n3, n4;
    int k, iter;
    double x, w, t, dt, C, theta, sint, f, fp, S, dS, dn = (double)n;
    double fn, fn5, nk, n5k, pn, jk, num, den, p , pp;
    
    /* The constant out the front of the sum in evaluation (compute only once)
    /* C = sqrt(pi/4)*gamma(n+1)/gamma(n+3/2) */
    /* Use Stirling's series */
    dS = -.125/dn;
    S = dS; 
    k = 1; 
    while ( fabs(dS/S) > EPSILON/100 ) {
        k += 1;
        dS *= -.5 * (k-1) / (k+1) / dn;
        S += dS;
    }
    double stirling[10] = {1., 1./12., 1./288., -139./51840., -571./2488320., 
                           163879./209018880., 5246819./75246796800., 
                          -534703531./902961561600., -4483131259./86684309913600., 
                           432261921612371./514904800886784000.};
    fn = 1; fn5 = 1; nk = dn; n5k = dn+.5;
    for ( k = 1 ; k < 10 ; k++ ) {
        fn += stirling[k]/nk;
        fn5 += stirling[k]/n5k;
        nk *= dn;
        n5k *= dn + .5;
    }
    C = exp(S)*sqrt(4.0/(dn+.5)/PI) * fn / fn5;  
    
    /* Approximations of roots near boundaries */
    double besselroots[32] = {2.404825557695773,  5.520078110286311,
                              8.653727912911013, 11.791534439014280,
                             14.930917708487790, 18.071063967910920,
                             21.211636629879269, 24.352471530749284,
                             27.493479132040250, 30.634606468431990,
                             33.775820213573560, 36.917098353664045,
                             40.058425764628240, 43.199791713176737,
                             46.341188371661815, 49.482609897397822,
                             52.624051841114984, 55.765510755019974,
                             58.906983926080954, 62.048469190227159,
                             65.189964800206866, 68.331469329856787,
                             71.472981603593752, 74.614500643701817,
                             77.756025630388066, 80.897555871137627,
                             84.039090776938195, 87.180629843641128,
                             90.322172637210500, 93.463718781944763,
                             96.605267950996279, 99.746819858680624};
    
    /* Useful constants */                         
    pn = dn*dn+dn+1./3.; 
    n3 = n*n*n; 
    n4 = n3*n; 
    
    /* Loop over each root */
    for ( k = (n+s)/2; k > 0; k-- ) {
        
        /* Initialise */
        dt = 1.0; 
        iter = 1;
        
        /* Asymptotic approximation of roots (Tricomi, 1950) */
        theta = PI*(4.0*k-1)/(4.0*n+2.0);
        sint = sin(theta);
        x = (1.0 - (n-1.0)/(8.0*n3)-(39.0-28.0/(sint*sint))/(384.0*n4))*cos(theta);  
        t = acos(x);
                
        /* Use different guesses near the ends */
        if (x > 0.5) { 
            
            /* Compute some bessel functions */
            if (k < 30) {                   /* These Bessel roots are hard coded */
                jk = besselroots[k-1];          
            }
            else {                          /* Compute more (Branders, JCP 1981) */
                p = (k-0.25)*PI;
                pp = p*p;
                num = 0.0682894897349453 + pp*(0.131420807470708 + pp*(0.0245988241803681 + pp*0.000813005721543268));
                den = p*(1.0 + pp*(1.16837242570470 + pp*(0.200991122197811 + pp*(0.00650404577261471))));
                jk = p + num/den;
            }
            
            /* Evaluate asymptotic approximations */
            if ( k <= 5 ) {                 /* Extreme boundary (Gatteschi, 1967) */
                t = jk/sqrt(pn)*(1.-(jk*jk-1.)/(pn*pn)/360.); 
                dt = 0.0;
            }
            else {                          /* Boundary (Olver, 1974) */
                p = jk/(n+.5);
                t = p + (p/tan(p)-1.)/(8.*p*(n+.5)*(n+.5));
            }
        }
        
        /* Newton iteration for roots */
        while (fabs(dt) > EPSILON) {
            feval(n, t, C, &f, &fp, k);     /* Evaluate at current approximation */
            dt = f/fp;                      /* Newton update */
            t += dt;                        /* Update t */
            if (iter++ > ITER_MAX) { break; } 
        }
        feval(n, t, C, &f, &fp, k);            /* Once more for luck */
        if ( dt!=0 ) {
            t += f/fp;
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

int feval(int n, double theta, double C, double* f, double* fp, int k)
{
    if (k > 10) {
        feval_asy1(n, theta, C, f, fp, k);
    }
    else {
        feval_asy2(n, theta, f, fp);
    }
    return 0;
}

double chebbessj0( double x )
{
   double coeffs[50] = {0.113905930122844,  -0.248029977895747,   
   0.216875874902150,  -0.210251640406686,   0.177883635103513,
  -0.121762914878623,   0.098895497651315,   0.012919886971102,
  -0.011705113604369,   0.112816231442338,  -0.085181686388300,
   0.039441066005212,  -0.025994587781630,  -0.106154546065109,
   0.071909852136308,   0.064563961364873,  -0.042499100646195,
  -0.022480172610839,   0.014332051364239,   0.005385391858862,
  -0.003324534307660,  -0.000964842744347,   0.000576985825155,
   0.000135778873812,  -0.000078709230886,  -0.000015512082673,
   0.000008723136824,   0.000001473672352,  -0.000000804536430,
  -0.000000118586276,   0.000000062900184,   0.000000008202725,
  -0.000000004230297,  -0.000000000493641,   0.000000000247703,
   0.000000000026109,  -0.000000000012756,  -0.000000000001224,
   0.000000000000583,   0.000000000000051,  -0.000000000000024,
  -0.000000000000002,   0.000000000000001,   0.000000000000000,
                   0,   0.000000000000000,  -0.000000000000000,
   0.000000000000000,   0.000000000000000,  -0.000000000000000};
   
   return clenshaw(coeffs, 50, x);
}

double chebbessj1( double x )
{
   double coeffs[50] = {0.042376943206206,  -0.130389123712321,
   0.052750018296831,  -0.074421155995637,  -0.028637713473499,
   0.017389752444885,  -0.107194432750030,   0.093954008691064,
  -0.095524857421292,   0.081871310776877,   0.035487540382714,
  -0.028040542627382,   0.091468408261080,  -0.068290226934422,
  -0.086597281912651,   0.061611441440843,   0.038365223954845,
  -0.026128637312591,  -0.010946122417318,   0.007158707791449,
   0.002256773752795,  -0.001420735583157,  -0.000357638844790,
   0.000217159662443,   0.000045317812976,  -0.000026585052558,
  -0.000004721163388,   0.000002679664528,   0.000000412920937,
  -0.000000227047734,  -0.000000030821259,   0.000000016436849,
   0.000000001989641,  -0.000000001030182,  -0.000000000112314,
   0.000000000056515,   0.000000000005596,  -0.000000000002739,
  -0.000000000000248,   0.000000000000118,   0.000000000000010,
  -0.000000000000005,  -0.000000000000000,   0.000000000000000,
   0.000000000000000,  -0.000000000000000,  -0.000000000000000,
  -0.000000000000000,                   0,                   0};
   return clenshaw(coeffs, 50, x);
}

double clenshaw( double* c, int n, double x) {
    double bk = 0., bk1 = 0., bk2 = 0.;
    int k;

    x = (4./31.)*x - 2.;
    for (k = n-1; k > 0; k-- ) { 
        bk = c[k] + x*bk1 - bk2;
        bk2 = bk1; 
        bk1 = bk;
    }
    return c[0] + .5*x*bk1 - bk2;
}


int feval_asy2(int n, double t, double* f, double* fp) {   
    double dn = (double)n, rho = dn + .5, rho2 = dn - .5;
    double sint, csct, cost, cott, tinv;
    double tB1t = 0., A2t = 0.;
    double Ja, Jab, Jb, Jbb, gt, gtdx, vals, vals2, ders, A1, tB0, denom, tmp;

    /* Evaluate bessel functions */
    Ja = chebbessj0( rho*t );
    Jb = chebbessj1( rho*t );
    Jab = chebbessj0( rho2*t );
    Jbb = chebbessj1( rho2*t );
    
    /* Evaluate functions for recurrsive definition of coefficients. */
    sint = sin(t);
    csct = 1/sint;
    cost = cos(t);
    cott = csct*cost;
    tinv = 1/t;
    gt = 0.5*(cott-tinv);
    gtdx = 0.5*(-csct*csct + tinv*tinv);
    tB0 = 0.25*gt;
    A1 = 0.125*(gtdx-gt*tinv-0.25*gt*gt);   
    
    /* first term: */
    vals = Ja;
    vals2 = Jab;
    
    /* second term: */
    tB0 = .25*gt;
    vals = vals + Jb*tB0/rho;
    vals2 = vals2 + Jbb*tB0/rho2;

    /* third term: */
    vals = vals + Ja*A1/(rho*rho); 
    vals2 = vals2 + Jab*A1/(rho2*rho2);
    
    /* higher terms: */
    
    /* combine */
    denom = sqrt(t*csct);
    *fp = n*(-cost*vals + vals2)*csct*denom;
    *f = vals*denom;
            
    return 0;
}


    
int feval_asy1(int n, double theta, double C, double* f, double* fp, int k)
{    
    double sinT = sin(theta), cosT = cos(theta), cotT = cosT/sinT;
    double alpha, cosA, sinA, denom, df, dfp, tmp;
    int m;
    
    /* m = 0 */
    denom = sqrtl(2.0*sinT);
    
    /*alpha = (n + 0.5)*theta - 0.25*PI;
    cosA = cos(alpha);
    sinA = sin(alpha);*/
    mycosA(n, theta, &cosA, &sinA, k);
    
    *f = C * cosA/denom;
    *fp = C * ( 0.5*(cosA*cotT+sinA) + n*sinA) / denom;
    
    /* Loop over m until term is sufficiently small */
    for ( m = 1; m <= M; m++ ) {
        C *= (1.0-0.5/m)*(m-0.5)/(n+m+0.5);
        denom *= 2.0*sinT;
        
        /*alpha += theta - 0.5*PI;
        cosA = cos(alpha);
        sinA = sin(alpha);*/
        
        tmp = cosA*sinT + sinA*cosT;
        sinA = sinA*sinT - cosA*cosT;
        cosA = tmp;
        
        df = C * cosA/denom;
        dfp = C * ( (m+0.5)*(cosA*cotT+sinA) + n*sinA ) / denom;
        
        *f += df;
        *fp += dfp;
        
        if (fabs(df) + fabs(dfp) < EPSILON/100){ break; }
    }
    
}

int mycosA(int n, double theta, double* cosA, double* sinA, int k) {
        int j;
        double dh, tmp = 0.0, DH, dh2, lo, hi = theta, fixsgn = (double)(1-2*(k%2));
        double k025 = (k-.25), rho = n+.5, sgn = 1.0, fact = 1.0, hi2;         
        
        /* bit shift to get a hi-lo version of theta */
        hi = (double)((float)hi);
        lo = theta - hi;  
        dh = (hi*rho-k025*PI_HI) + lo*rho - k025*PI_LO;
        /* easy way: dh = (n+0.5)*theta-(k-.25)*PI; */
        
        DH = dh; dh2 = dh*dh;
        for ( j = 0; j <= 5; j++ ) {
            tmp += sgn*DH/fact;
            sgn = -sgn;
            fact = fact*(double)((2*j+3)*(2*j+2));
            DH *= dh2;
        }
        *cosA = tmp*fixsgn;
        
        tmp = 0.0; sgn = 1.0; fact = 1.0; DH = 1.0;
        for ( j = 0; j <= 5; j++ ) {
            tmp += sgn*DH/fact;
            sgn = -sgn;
            fact = fact*(double)((2*j+2)*(2*j+1));
            DH *= dh2;
        }
        *sinA = -tmp*fixsgn;
        
    return 0;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *result_roots,*result_ders;
  int n;
  /*  Check for proper number of arguments */
  if(nrhs!=1) 
    mexErrMsgTxt("One input required.");
  
  /*  Create a pointer to the input matrix x  */
  n = mxGetScalar(prhs[0]);
  
  /*  Set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,n, mxREAL);
  
  /*  Create a C++ pointer to a copy of the output matrix */
  result_roots = mxGetPr(plhs[0]);
  result_ders = mxGetPr(plhs[1]);
  
  /*  Call the C++ subroutine */
  asy(result_roots,result_ders,n);
  
  return;
}