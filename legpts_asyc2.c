#include <mex.h>
#include <math.h>


/* #define PI 3.141592653589793 */
#define EPSILON 2.2204460492503131e-40

#define PI 3.141592653589793238462643l


#define M 60
#define ITER_MAX 10

/* USAGE 
 
 ./legpts_asy1 N - compute nodes and weights of degree N. Display on screen.
 ./legpts_asy1 N x.txt w.txt - compute node and weights, and print to file.
 
*/

int asy1(double *nodes, double*weights, unsigned long int n)
{
    int s = n % 2;
    unsigned long int n3, n4;
    int k, iter;
    long double theta, x, winv, C, sint, f, fp, dx;
    
    C = (4.0/PI);
    for ( k = 1; k <= n; k++ ) {
        C *= k/(k + 0.5);
    }

    for ( k = (n+s)/2; k > 0; k-- ) {
        theta = PI*(4.0*k-1)/(4.0*n+2.0);
        n3 = n*n*n; n4 = n3*n; sint = sin(theta);
        x = (1.0 - (n-1.0)/(8*n3)-(39.0-28.0/(sint*sint))/(384.0*n4))*cos(theta);  
        dx = 1.0;
        iter = 1;

        while (fabs(dx) > 2.0*EPSILON) {
            feval(n, x, C, &f, &fp, &winv); 
            dx = -f/fp;                     
            x += dx;          
            if (iter++ > ITER_MAX) { break; } 
        }

        nodes[n-k] = (double)x;
        nodes[k-1] = -(double)x;
        winv = 1.0/winv;
        weights[n-k] = (double)winv;
        weights[k-1] = (double)winv;
    }

    if (s == 1) {
        nodes[(n+s)/2-1] = 0.0;
    }

    return 0;
}

int feval(unsigned long int n, long double x, long double C, long double* f, long double* fp, long double* winv)
{
    long double theta = acos(x), sinT = sin(theta), cotT = x/sinT;
    long double alpha, cosA, sinA, denom;
    int m;
    
    alpha = (n + 0.5)*theta - 0.25*PI;
    denom = sqrtl(2.0*sinT);
    cosA = cosl(alpha);
    sinA = sinl(alpha);
    *f = C * cosA/denom;
    *fp = C * ( 0.5*(cosA*cotT+sinA) + n*sinA) / denom;

    for ( m = 1; m <= M; m++ ) {
        C *= (2.0*m-1.0)/(2.0*m)*(m-0.5)/(n+m+0.5);
        alpha += theta - 0.5*PI;
        denom *= 2.0*sinT;
        cosA = cos(alpha);
        sinA = sin(alpha);
        *f += C * cosA/denom;
        *fp += C * ( (m+0.5)*(cosA*cotT+sinA) + n*sinA) / denom;
        
        if (fabs(C) < EPSILON*EPSILON){ break; }
    }

    *winv = .5*(*fp)*(*fp);    
    *fp /= sinT;
    
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
  plhs[1] = mxCreateDoubleMatrix(n,1, mxREAL);
  
  /*  create a C++ pointer to a copy of the output matrix */
  result_roots = mxGetPr(plhs[0]);
  result_ders = mxGetPr(plhs[1]);
  
  /*  call the C++ subroutine */
  asy1(result_roots,result_ders,n); 
  
  return;
}