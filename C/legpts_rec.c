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
#define EPSILON 2.2204460492503131e-16
#define ITER_MAX 10

int rec(double *nodes, double*weights, unsigned long int n)
{
    int s = n % 2;                       /* True if n is odd */
    unsigned long int n3, n4;
    int k, iter;
    double x, w, t, dx, C, theta, sint, f, fp, S, dS, dn = (double)n;
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
        dx = 1.0; 
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
            }
            else {                          /* Boundary (Olver, 1974) */
                p = jk/(n+.5);
                t = p + (p/tan(p)-1.)/(8.*p*(n+.5)*(n+.5));
            }
            x = cos(t);
        }
       
        /* Newton iteration for roots */
        while (fabs(dx) > EPSILON) {
            feval(n, x, &f, &fp);        /* Evaluate at current approximation */
            dx = f/fp;                      /* Newton update */
            x += dx;                        /* Update t */
            if (iter++ > ITER_MAX) { break; } 
        }
        feval(n, x, &f, &fp);            /* Once more for luck */
        if ( dx!=0 ) {
            x += f/fp;
        }
        
        /* Compute weights */
        w = 2.0/(fp*fp)/(1-x*x);        
        
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

int feval(int n, double x, double* f, double* fp)
    {
    double P, Pm1 = x, Pm2 = 1.0;
    /* double PP, PPm1 = 1.0;, PPm2 = 0.0; */
    int k;
    
    /* Loop over m until term is sufficiently small */
    for ( k = 1; k < n; k++ ) {
        P = ((double)(2*k+1)*Pm1*x-(double)k*Pm2)/(double)(k+1);
        Pm2 = Pm1;
        Pm1 = P;
        /*PP = ((double)(2*k+1)*(Pm2+x*PPm1)-(double)k*PPm2)/(double)(k+1);  
        PPm2 = PPm1; 
        PPm1 = PP; */ 
    }
    *f = P;
    /**fp = PP; */
    *fp = n*(x*P-Pm2)/(1-x*x);  

}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
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
  rec(result_roots,result_ders,n);
  
  return;
}