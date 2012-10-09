#include <mex.h>
#include <math.h>

#define PI 3.141592653589793
#define EPSILON 2.2204460492503131e-16
#define SQRTEPS 1.4901161193847661e-08
#define M 30
#define ITER_MAX 10
#define NN 30
#define CC .5*PI

double mysin(double n, double x);
double mycos(double n, double x);
double asymbesselj(double nu, double n, double x);

int asy1(double *nodes, double*weights, unsigned long int n, double a, double b)
{
    int j, k, iter, oldk;
    double x, w, t, dt, C, C2, K, tK, tmp, f, fp, dn = (double)n, pn, jk, phik, phi;
    double dS, S, nak, nk, fna, fn, fnab, ff, Cw;
    double rho = dn + .5*(a+b+1.), rho2 = dn + .5*(a+b-1.);
    double tB1[NN], A2[NN], xcheb[NN];
    
    // The constant out the front of the sum in evaluation (compute only once)
    dS = .5*a*a/dn;
    S = dS; k = 1;
    while (fabs(dS/S) > EPSILON/10) {
        k += 1;
        dS *= -(k-1.)/(k+1.)/dn*a;
        S += dS;
    }
    double stirling[10] = {1., 1./12., 1./288., -139./51840., -571./2488320., 163879./209018880,
        5246819./75246796800., -534703531./902961561600.,
        -4483131259./86684309913600., 432261921612371./514904800886784000.};
    fn = 1.; nk = n; fna = 1.; nak = n+a;
    for ( k = 1 ; k < 10 ; k++ ) {
        fn += stirling[k]/nk;
        fna += stirling[k]/nak;
        nk *= n;
        nak *= n+a;
    }
    ff = fna/fn;
    C = exp(S)*sqrtl(.5*(dn+a)/dn)*pow(dn/rho,a)*ff;
    C2 = exp(S)*sqrtl(.5*dn/(dn+a))*pow(dn/rho2,a)*ff;
//    C2 = C*n/(n+1.0)*pow(rho/rho2,a);
    
    // Constant for the weights
    Cw = 1; phi = -a*b/(double)n;
    for (k = 1; k <= 20; k++) {
        Cw = Cw + phi;
        phi = -(k+a)*(k+b)/(k+1.)/(n-k)*phi;
        if (fabs(phi) < EPSILON/10) break;
        if (k == n-1) break;
    }
    Cw *= pow(2,a+b+1);
    
    // Chebyshev points
    for ( k=0; k<NN; k++ ) {
        xcheb[NN-k-1] = .5*CC*(1+cos(PI*k/(NN-1.)));
    }
    
    // Loop over each root (first half, theta <= PI/2)
    for ( k = 1; k <= n; k++ ) {
        
        // Asymptotic approximation of roots (in theta-space)
        K = PI*(4.0*k-1.+2.*a)/(4.0*n+2.0*(a+b+1));
        tK = tan(.5*K); tmp = 2*n+a+b+1;
        t = K + ((.25-a*a)/tK-(.25-b*b)*tK)/(tmp*tmp);
        x = cos(t);
        
        if ( t > PI/2 ) break;
        
        // Initialise
        dt = 1.0; 
        iter = 1;
        
        // Newton iteration for roots
        while (fabs(dt) > EPSILON) {
            feval(n, t, C, C2, &f, &fp, a, b, tB1, A2, xcheb);  // Evaluate at current approximation
            dt = f/fp;                          // Newton update
            t += dt;                            // Update t
            if (iter++ > ITER_MAX) break;
        }
        
        // Convert back to x-space
        x = cos(t);
        w = Cw / (fp*fp);        
        
        // Store nodes and weights
        nodes[n-k] = x;
        weights[n-k] = w;
    }
    
    // Use symmetry and flip for theta > PI/2
    tmp = a;
    a = b;
    b = tmp;
    
    // The constant out the front of the sum in evaluation (compute only once)
    dS = .5*a*a/dn;
    S = dS; j = 1;
    while (fabs(dS/S) > EPSILON/10) {
        j += 1;
        dS *= -(j-1.)/(j+1.)/dn*a;
        S += dS;
    }
    fn = 1.; nk = n; fna = 1.; nak = n+a;
    for ( j = 1 ; j < 10 ; j++ ) {
        fn += stirling[j]/nk;
        fna += stirling[j]/nak;
        nk *= n;
        nak *= n+a;
    }
    ff = fna/fn;
    C = exp(S)*sqrtl(.5*(dn+a)/dn)*pow(dn/rho,a)*ff;
    C2 = exp(S)*sqrtl(.5*dn/(dn+a))*pow(dn/rho2,a)*ff;
    
    // Loop over each root (second half, theta > PI/2)
    for ( k = k-1; k < n; k++ ) {
        
        // Asymptotic approximation of roots (in theta-space)
        K = PI*(4.0*(n-k)-1.+2.*a)/(4.0*n+2.0*(a+b+1));
        tK = tan(.5*K); tmp = 2*n+a+b+1;
        t = K + ((.25-a*a)/tK-(.25-b*b)*tK)/(tmp*tmp);
        x = cos(t);
        
        // Initialise
        dt = 1.0; 
        iter = 1;
        
        // Newton iteration for roots
        while (fabs(dt) > EPSILON) {
            feval(n, t, C, C2, &f, &fp, a, b, tB1, A2, xcheb);  // Evaluate at current approximation
            dt = f/fp;                          // Newton update
            t += dt;                            // Update t
            if (iter++ > ITER_MAX) break;
        }
        
        // Convert back to x-space
        x = -cos(t);
        w = Cw/(fp*fp);        
        
        // Store nodes and weights
        nodes[n-1-k] = x;
        weights[n-1-k] = w;
    }
    
    return 0;
}

int feval(int n, double t, double C, double C2, double* f, double* fp, double a, double b, double tB1[NN], double A2[NN], double xcheb[NN])
{
    double A = (.25-a*a), B = (.25-b*b), dn = (double)n;
    double rho = dn + .5*(a+b+1.), rho2 = dn + .5*(a+b-1.);
    double rho_3 = rho*rho*rho, rho_4 = rho_3*rho, rho2_3 = rho2*rho2*rho2, rho2_4 = rho2_3*rho2;
    double sqrtt = sqrtl(t), tB1t = 0., A2t = 0.;
    double Ja, Jab, Jb, Jbb, gt, sin5t, cos5t, gtdx, vals, vals2, valstmp, ders, A10, A1, tB0, denom;
    
    Ja = asymbesselj(a,rho,t); 
    Jab = asymbesselj(a,rho2,t);
    Jb = asymbesselj(a+1.,rho,t);
    Jbb = asymbesselj(a+1.,rho2,t);
    
    // Evaluate functions for recurrsive definition of coefficients. 
    cos5t = cos(.5*t); sin5t = sin(.5*t);
    gt = A*(cos5t/sin5t-(2./t)) - B*sin5t/cos5t;
    gtdx = A*(2./(t*t)-.5/(sin5t*sin5t)) - .5*B/(cos5t*cos5t);
    
    // first term:
    vals = Ja;
    vals2 = Jab;
    
    // second term: 
    tB0 = .25*gt;
    vals = vals + Jb*tB0/rho;
    vals2 = vals2 + Jbb*tB0/rho2;

    // third term:
    A10 = a*(A+3.*B)/24.;
    A1 = .125*gtdx - .125*(1.+2.*a)*gt/t - gt*gt/32. - A10;
    vals = vals + Ja*A1/(rho*rho); 
    vals2 = vals2 + Jab*A1/(rho2*rho2);
    // higher terms:
//    bary(tB1, xcheb, t, &tB1t);
//    bary(A2, xcheb, t, &A2t);    
    vals = vals + Jb*tB1t/rho_3 + Ja*A2t/rho_4;
    //    vals2 = vals2 + Jbb*tB1/rho2_3 + Jab*A2/rho2_4;

    valstmp = C*vals;
    denom = pow(sin(.5*t),a+.5)*pow(cos(.5*t),b+.5);
    vals = sqrtt * valstmp / denom;
    *f = vals;
    
    // Relation for derivative
    ders = ( n*(a-b-(2.*n+a+b)*cos(t))*valstmp + 2.*C2*(n+a)*(n+b)*vals2 ) / (2.*n+a+b);
    *fp = ders * sqrtt / (denom*sin(t));
    
    return(0);
    
}

double asymbesselj(double nu, double n, double x)
{
    double pp, fournu2 = 4.*nu*nu, cosx, sinx, cosp, sinp, cosW, sinW, man, wife, a, sgn, tmp, vals;
    int j, MM = 30;
    
    // avg frequency.
    pp = (.5*nu+.25)*PI;
    
    // coefficients
    a = 1.; 
    sgn = 1.0;
    man = 1.; 
    wife = 0.;
    
    for (j = 1; j <= MM/2; j+=2) { 
        
        tmp = (2.*j-1.);
        a *= (fournu2-tmp*tmp)/(8.*j)/n/x;
        wife += a;
        
        sgn = -sgn;  
        
        tmp = (2.*j+1.);
        a *= sgn*(fournu2-tmp*tmp)/(8.*(j+1))/n/x;
        man += a;
        
        sgn = -sgn;  
    }
    
    // factor out the average frequency
    cosx = mycos(n,x); sinx = mysin(n,x); 
    cosp = cos(pp); sinp = sin(pp);
    cosW = cosx*cosp + sinx*sinp;
    sinW = sinx*cosp - cosx*sinp;
    
    // besselj is the marriage of two parts.
    man *= cosW;
    wife *= sinW; 
    
    // combine
    vals = sqrt(2./PI)*(1./sqrt(n))*(1./sqrt(x))*(man - wife);
    
    return(vals);
    
}

double mycos(double n, double x)
{
    int k;
    double fl, rem, val, flx = fl*x, s;
    
    if ( n > 2 ) {
        k = floor(log2(n));
        fl = pow(2,k);
        rem = n - fl;
        flx = fl*x;
        val = cos(flx)*mycos(rem,x) - sin(flx)*mysin(rem,x);
    }
    else val = cos(n*x);
    
    return(val);
    
}

double mysin(double n, double x)
{
    int k;
    double fl, rem, val, flx = fl*x;;
    
    if ( n > 2 ) {
        k = floor(log2(n));
        fl = pow(2,k);
        rem = n - fl;
        flx = fl*x;
        val = cos(flx)*mysin(rem,x) + sin(flx)*mycos(rem,x);
    }
    else val = sin(n*x);
    
    return(val);
}

int bary ( double vals[NN] , double points[NN], double x, double *f) {
    
    int k;
    double w, u = 0.0, v = 0.0;
    
    /* Do the barycentric form in place. */
    /* Do the first node separately due to the half-weight. */
    //    if ( x == points[0] )
    //        return vals[0]
    w = 0.5 / ( x - points[0] );
    u = vals[0] * w;
    v = w;
    
    /* Do the interior nodes. */
    for ( k = 1 ; k < NN-1 ; k++ ) {
        //        if ( x == points[k] )
        //            return vals[k]
        w = (double)( 1 - 2 * ( k & 1 ) ) / ( x - points[k] );
        u += vals[k] * w;
        v += w;
    }
    
    /* Do the last node separately due to the half-weight. */
    //    if ( x == points[k] ) 
    //        return vals[k]
    w = ( 0.5 - ( k & 1 ) ) / ( x - points[k] );
    u += vals[k] * w;
    v += w;
    
    /* Return the fraction. */
    *f = u / v;
    return(0);
}

int asy2_higherterms(double a, double b, double* tB1, double* A2) {
    double A = (.25-a*a), B = (.25-b*b);
    double t[NN], v[NN], g[NN], gp[NN], gpp[NN], f[NN];
    double B0[NN], B0p[NN], B1[NN], A1[NN], A1p[NN], A1p_t[NN];
    double I[NN], J[NN], K[NN], K2[NN], tmpvec[NN];
    double sgn, tmp, tk, sin5t, cos5t, tan5t, sin5t2, sin5t4, sint;
    int j, k;
    
    // Chebyshev points and barycentric weights
    sgn = 1.;
    for ( k=0; k<NN; k++ ) {
        t[NN-k-1] = .5*CC*(1+cos(PI*k/(NN-1.)));
        v[NN-k-1] = sgn;
        sgn = -sgn;
    }
    v[NN] = .5*v[NN];
    v[0] = .5*v[0];
    
    // Basic vectors
    g[0] = 0.; gp[0] = -A/6.-.5*B; gpp[0] = 0.; f[0] = -A/12.-.25*B;    
    for ( k=1; k<NN; k++ ) {
        tk = t[k]; 
        cos5t = cos(.5*tk);
        sin5t = sin(.5*tk);
        tan5t = tan(.5*tk);
        sin5t2 = sin5t*sin5t;
        sin5t4 = sin5t2*sin5t2;
        g[k] = A*(1./tan5t-2./tk)-B*tan5t;
        gp[k] = A*(2./(tk*tk)-.5/(sin5t*sin5t))-.5*B/(cos5t*cos5t);
        sint = sin(tk);
        gpp[k] = A*(-4./(tk*tk*tk)+.25*sint/sin5t4) - 4.*B*sin5t4/(sint*sint*sint);
        f[k] = A*(1./(tk*tk) - .25/sin5t2) - .25*B/(cos5t*cos5t);
    }
    
    // B0 and A1
    tmp = a*(A+3.*B)/24.;
    B0[0] = .25*(-A/6.-.5*B);
    B0p[0] = 0.;
    A1[0] = 0.;
    A1p[0] = 0.;
    A1p_t[0] = -A/720.-A*A/576.-A*B/96.-B*B/64.-B/48.+a*(A/720.+B/48.);
    for ( k=1; k<NN; k++ ) {
        B0[k] = .25*g[k]/t[k];
        B0p[k] = .25*(gp[k]-g[k]/t[k])/t[k];
        A1[k] = .125*gp[k] - .5*(1.+2.*a)*B0[k] - g[k]*g[k]/32. - tmp;
        A1p[k] = .125*gpp[k] - .5*(1.+2.*a)*B0p[k] - gp[k]*g[k]/16.;
        A1p_t[k] = A1p[k]/t[k];
    }
    
    // Integrals
    cumsum(A1p_t, I);
    for ( k=0; k<NN; k++ ) tmpvec[k] = f[k]*A1[k];
    cumsum(tmpvec, J); 
    
    // B1
    tB1[0] = 0.;
    B1[0] = A/720.+A*A/576.+A*B/96.+B*B/64.+B/48.+a*(A*A/576.+B*B/64.+A*B/96.)-a*a*(A/720.+B/48.);
    for ( k=1; k<NN; k++ ) {
        tB1[k] =  -.5*A1p[k] - (.5+a)*I[k] + .5*J[k];
        B1[k] = tB1[k]/t[k];
    }
    
    // Integrals
    for ( k=0; k<NN; k++ ) tmpvec[k] = f[k]*tB1[k];
    cumsum( tmpvec, K );
    diff( tB1, K2 );
    
    // A2
    A2[0] = 0.;
    tmp = .5*K2[0] - (.5+a)*B1[0] - .5*K[0];
    for ( k=1; k<NN; k++ ) A2[k] = .5*K2[k] - (.5+a)*B1[k] - .5*K[k] - tmp;
    
    for ( k=0; k<NN; k++ ) {
        tB1[k] = tB1[k];
        A2[k] = A2[k];
    }
    return(0);
    
}

/* indefinite integral */
int cumsum( double f[NN], double g[NN]) {
    int j, k;
    double t[NN], B[NN][NN] = {0}, T[NN][NN], Ti[NN][NN];
    double c = 1./(NN-1.0), s, sgn;
    
    // vals --> coeffs
    sgn = 1.;
    for ( j=0; j<NN; j++ ) {
        Ti[j][0] = sgn*c;
        for ( k=1; k<NN-1; k++ ) {
            Ti[j][k] = c*(cos(j*(NN-k-1.)*PI/(NN-1.))+cos(j*(NN+k-1.)*PI/(NN-1.)));
        }
        Ti[j][NN-1] = c;
        sgn = -sgn;
    }
    for ( k=0; k<NN; k++ ) {
        Ti[0][k] *= .5;
        Ti[NN-1][k] *= .5;
    }
    
    // coeff cumsum matrix
    for ( j=2; j<NN; j++ ) {
        B[j-1][j] = -1./(2.*(j-1.));
        B[j][j-1] = 1./(2.*j);
    }  
    sgn = -1.;
    for ( j=1; j<NN; j++ ) {
        B[0][j] = sgn*(B[j-1][j] + B[j+1][j]);
        sgn = -sgn;
    }
    B[0][0] = 1.;
    B[1][0] = 1.;
    
    // coeffs --> vals
    for ( j=0; j<NN; j++ ) {
        for ( k=0; k<NN; k++ ) {
            T[j][k] = cos(k*(NN-j-1.)*PI/(NN-1.));
        }
    }
    
    // Do the multiplications
    multiply( Ti, f, g);
    multiply( B, g, t);
    multiply( T, t, g);
    
    // Constant factor for domain [0 C]
    for ( j=0; j<NN; j++ ) g[j] *= .5*CC;
    
    return(0);
}

/* differentiate */
int diff( double f[NN], double g[NN] ) {    
    int j,k;
    double t[NN], v[NN], s, sgn, D[NN][NN];
    // Chebyshev points and barycentric weights
    sgn = 1.;
    for ( k=0; k<NN; k++ ) {
        t[k] = -cos(PI*k/(NN-1.));
        v[k] = sgn;
        sgn = -sgn;
    }
    v[NN-1] = .5*v[NN-1];
    v[0] = .5*v[0];
    
    // make the diffmat
    for ( j=0; j<NN; j++ ) {
        for ( k=0; k<NN; k++ ) {
            D[j][k] = v[k]/v[j]/(t[j] - t[k]);
        }
        D[j][j] = 0.0;
        s = 0.;
        for ( k=0; k<NN; k++ ) {
            s += D[j][k];
        }
        D[j][j] = -s;
    }
    
    // Do the multiply
    multiply( D, f, g);
    
    // Constant factor for domain [0 C]
    for ( j=0; j<NN; j++ ) g[j] *= 4. / PI;
    
    return(0);
}

// Simple multiply function
int multiply( double A[NN][NN], double f[NN], double c[NN] ) {
    int j, k;
    double s;
    
    // Loop
    for ( j=0; j<NN; j++) {
        s = 0.;
        for ( k=0; k<NN; k++ ) {
            s += A[j][k]*f[k];        
        }
        c[j] = s;
    }
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