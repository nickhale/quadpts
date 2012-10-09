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
#define SQRTEPS 1.4901161193847661e-08
#define M 30
#define ITER_MAX 10

int asy1(double *nodes, double*weights, unsigned long int n, double a, double b)
{
    int k, iter;
    double x, w, t, dt, C, C2, K, tK, tmp, f, fp, dn = (double)n, pn, jk, rho, phik;
    double dSa, dSb, dSab, dS, S, nak, nbk, nabk, fna, fnb, fnab, Cw, phi;
    
    /* The constant out the front of the sum in evaluation (compute only once) */
    /* C = 2^(2n+a+b+1)/pi*B(n+a+1,n+b+1) */
    /* Use Stirling's series */
    dSa = .5*a*a/n; dSb = .5*b*b/n; dSab = .25*(a+b)*(a+b)/n;
    dS = dSa + dSb - dSab; S = dS; k = 1;
    while ( fabs(dS/S) > EPSILON/10 ) {
        k += 1;
        tmp = -(k-1.)/(k+1.)/n;
        dSa *= tmp*a;
        dSb *= tmp*b;
        dSab *= .5*tmp*(a+b);
        dS= dSa + dSb - dSab;
        S += dS;
    }
    double stirling[10] = {1., 1./12., 1./288., -139./51840., -571./2488320., 
        163879./209018880., 5246819./75246796800., 
        -534703531./902961561600., -4483131259./86684309913600., 
        432261921612371./514904800886784000.};
    fna = 1.; fnb = 1.; fnab = 1.; nak = n+a; nbk = n+b; nabk = 2*n+a+b;
    for ( k = 1 ; k < 10 ; k++ ) {
        fna += stirling[k]/nak;
        fnb += stirling[k]/nbk;
        fnab += stirling[k]/nabk;        
        nak *= n+a;
        nbk *= n+b;
        nabk *= 2*n+a+b;
    }
    C = 2./PI * exp(S) * sqrtl(2.0*PI*(n+a)*(n+b)/(2*n+a+b))/(2*n+a+b+1) * fna*fnb/fnab; 
    C2 = .25 * C * (2*n+a+b)*(2*n+a+b+1) / ((a+n)*(b+n));
    
    /* Constant for the weights */
    Cw = 1; phi = -a*b/(double)n;
    for (k = 1; k <= 20; k++) {
        Cw = Cw + phi;
        phi = -(k+a)*(k+b)/(k+1.)/(n-k)*phi;
        if (fabs(phi) < EPSILON/10) break;
        if (k == n-1) break;
    }
    Cw *= pow(2,a+b+1);
    
    /* Tabulated Bessel roots (for approximations of roots near boundaries) */
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
    rho = dn + a + b + .5;
    
    /* Loop over each root (first half, theta <= PI/2) */
    for ( k = 1; k <= n; k++ ) {
        
        /* Asymptotic approximation of roots (in theta-space) */
        K = PI*(4.0*k-1+2*a)/(4.0*n+2.0*(a+b+1));
        tK = tan(.5*K); tmp = 2*n+a+b+1;
        t = K + ((.25-a*a)/tK-(.25-b*b)*tK)/(tmp*tmp);
        x = cos(t);
        
        if ( t > PI/2 ) break;
        
        /* Initialise */
        dt = 1.0; 
        iter = 1;
        
        /* Skip x near the ends */
        if (x>1-4/(log(n)*n)) { 
            if (k < 30) { 
                jk = besselroots[k-1];
                phik = jk/rho;
                t = phik + (.5*(a*a-.25)*(1-phik/tan(phik))/phik-.25*(a*a-b*b)*tan(.5*phik))/(rho*rho);
            } 
/*            dt = 0.; */
        }
        
        /* Newton iteration for roots */
        while (fabs(dt) > EPSILON) {
            feval(n, t, C, C2, &f, &fp, a, b);  /* Evaluate at current approximation */
            dt = f/fp;                          /* Newton update */
            t += dt;                            /* Update t */
            if (iter++ > ITER_MAX) break;
        }
        if (dt == 0) feval(n, t, C, C2, &f, &fp, a, b);
        
        /* Convert back to x-space */
        x = cos(t);
        w = Cw / (fp*fp);        
        
        /* Store nodes and weights */
        nodes[n-k] = x;
        weights[n-k] = w;
    }
    
    /* Use symmetry and flip for theta > PI/2 */
    tmp = a;
    a = b;
    b = tmp;
    
    /* Loop over each root (second half, theta > PI/2) */
    for ( k = k; k <= n; k++ ) {
        
        /* Asymptotic approximation of roots (in theta-space) */
        K = PI*(4.0*k-1+2*a)/(4.0*n+2.0*(a+b+1));
        tK = tan(.5*K); tmp = 2*n+a+b+1;
        t = K + ((.25-a*a)/tK-(.25-b*b)*tK)/(tmp*tmp);
        t = PI - t;
        x = cos(t);
        
        /* Initialise */
        dt = 1.0; 
        iter = 1;
        
        /* Skip x near the ends */
        if (x>1-4/(log(n)*n)) { 
            if (k < 30) { 
                jk = besselroots[k-1];
                phik = jk/rho;
                t = phik + (.5*(a*a-.25)*(1-phik/tan(phik))/phik-.25*(a*a-b*b)*tan(.5*phik))/(rho*rho);
            }
/*            dt = 0.; */
        }
        
        /* Newton iteration for roots */
        while (fabs(dt) > EPSILON) {
            feval(n, t, C, C2, &f, &fp, a, b);  /* Evaluate at current approximation */
            dt = f/fp;                          /* Newton update */
            t += dt;                            /* Update t */
            if (iter++ > ITER_MAX) break;
        }
        if (dt == 0) feval(n, t, C, C2, &f, &fp, a, b);
        
        /* Convert back to x-space */
        x = -cos(t);
        w = Cw/(fp*fp);        
        
        /* Store nodes and weights */
        nodes[n-k] = x;
        weights[n-k] = w;
    }
    
    return 0;
}

int feval(int n, double t, double C, double C2, double* f, double* fp, double a, double b)
{
    double cost, sint, cos5t, sin5t, cot5t, alpha, cosA, sinA, cosA2, sinA2;
    double df, dfp, sgn, tmp;
    double phi_m, phi2_m, denom_m, phi_l, phi2_l, denom_l;
    int l, m;   
    
    
    /* These are used often */
    cost = cos(t); 
    sint = sin(t);
    cos5t = cos(.5*t); 
    sin5t = sin(.5*t); 
    cot5t = cos5t/sin5t;
    
    /* m = 0 */
    alpha = .5*(2*n+a+b+1)*t - .5*(a+.5)*PI;
    sinA = sin(alpha);
    cosA = cos(alpha); 
    cosA2 = cosA*cost + sinA*sint;
    *f = cosA; 
    *fp = cosA2;
    
    /* Initialise  */
    phi_m = 1; 
    phi2_m = 1; 
    denom_m = 1;
    
    /* Loop over m */
    for ( m = 1; m <= M-1; m++ ) {
        df = 0; dfp = 0; 
        alpha = alpha + .5*t;
        
        /* sin and cos alpha_m */
        sinA = sin(alpha);
        cosA = cos(alpha); 
        sinA2 = sinA*cost - cosA*sint;
        cosA2 = cosA*cost + sinA*sint;
        
        /* Update m constants */
        phi_m *= (.5+b+m-1.)*(.5-b+m-1.)/(m)/(2.*n+a+b+m+1.);
        phi2_m *= (.5+b+m-1.)*(.5-b+m-1.)/(m)/(2.*n+a+b+m-1.);
        denom_m /= 2.*cos5t;
        
        /* Initialise l constants */
        phi_l = phi_m; 
        phi2_l = phi2_m; 
        denom_l = denom_m; 
        sgn = 1.;
        
        /* Loop over l */
        for ( l = 0; l <= m; l+=2 ) {
            
            /* Even l term */
            df = df + sgn*phi_l*denom_l*cosA;
            dfp = dfp + sgn*phi2_l*denom_l*cosA2;
            
/*            if ( l == m ) continue; /* Doesn't matter ? */
            
            /* Update l constants */
            tmp = (double)(m-l)/(l+1.)*(.5+a+l)*(.5-a+l)/(.5+b+m-l-1.)/(.5-b+m-l-1.);
            phi_l *= tmp;
            phi2_l *= tmp;
            denom_l *= cot5t;
            
            /* Odd l term */
            df += sgn*phi_l*denom_l*sinA;
            dfp += sgn*phi2_l*denom_l*sinA2;
            
            /* Update l constants */
            tmp = (m-l-1.)/(l+2.)*(.5+a+l+1.)*(.5-a+l+1.)/(.5+b+m-l-2.)/(.5-b+m-l-2.);
            phi_l *= tmp;
            phi2_l *= tmp;
            denom_l *= cot5t;
            
            /* Switch sign */
            sgn = -sgn;
        }
        
        /* Update for m */
        *f = *f + df;
        *fp = *fp + dfp;
        
        /* Quit */
        if ( fabs(df) + fabs(dfp) < EPSILON/10) break;
        
    }
    
    /* Constant out the front */
    *f *= C;
    *fp *= C2;
    
    /* Recurrence relation for derivative */
    *fp = (n*(a-b-(2.*n+a+b)*cost)**f + 2.*(n+a)*(n+b)**fp)/(2.*n+a+b)/sint;
    
    /* The factor out the front */
    tmp = pow(sin5t,a+.5)*pow(cos5t,b+.5);
    *f /= tmp;
    *fp /= tmp;
    
    return(0);
    
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    double *result_roots, *result_ders, alpha, beta;
    int n;
    
    /*  check for proper number of arguments */
    if(nrhs!=3) 
        mexErrMsgTxt("Three inputs required.");
    
    /*  create a pointer to the inputs  */
    n = mxGetScalar(prhs[0]);
    alpha = mxGetScalar(prhs[1]);
    beta = mxGetScalar(prhs[2]);    
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,n, mxREAL);
    
    /*  create a C++ pointer to a copy of the output matrix */
    result_roots = mxGetPr(plhs[0]);
    result_ders = mxGetPr(plhs[1]);
    
    /*  call the C++ subroutine */
    asy1(result_roots,result_ders,n,alpha,beta); 
    
    return;
}