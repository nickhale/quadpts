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
#define EXP1 2.718281828459046
#define M 30
#define ITER_MAX 10

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
double clenshaw( double* c, int n, double x);
double chebbessj0( double x );
double chebbessj1( double x );
int besselj_taylor(double mu, double x, double z, double *j);
double besselj(double mu, double x);

double besselj(double mu, double x)
{
    if (mu < 0) {
        return -besselj(-mu, x);
    }
    else if (mu == 0.) {
        return chebbessj0( x );
    }
    else if ( mu == 1.) {
        return chebbessj1( x );
    }
    else if ( pow((EXP1*x)/(2*mu),mu) < EPSILON ) {
        return pow((EXP1*x)/(2*mu),mu);
    }
    else if ( mu > 10. ) {
        return (2.*(mu+1.)/x)*besselj( mu+1., x) - besselj( mu+2., x);
    }
    else {
        return (2.*(mu-1.)/x)*besselj( mu-1., x) - besselj( mu-2., x);
    }
}

int besselj_taylor(double mu, double t, double z, double *J) {
    int j,k,l,lmax = 10;
    int lchoosej, sgn;
    double c, s, s2; 
    
    *J = besselj(mu,z);
    return 0;
    
    c = 1.;
    s = chebbessj0( z );    
    for ( l = 1; l <= lmax; l++ ) {
        c *= -t/(double)(2*l);
        s2 = 0.;
        sgn = 1;
        lchoosej = 1;
        for ( j = 0; j <= l; j++ ) {
           s2 += (double)(sgn*lchoosej)*besselj((double)(-l+2*j),z);
           
           sgn *= -1;
           lchoosej *= (l-j)/(j+1.);
        }
        mexPrintf("%16.16f\n",c);
        s += c*s2;
    }
    
    *J = s;
    return 0;
}
               


double chebbessj0_old( double x )
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
   
   return clenshaw(coeffs, 50, 2./31.*x-1.);
}

double chebbessj0( double x )
{
   double coeffs[30] = {0.21675296296201045174313830389014,
                             -0.73197997451675389391923209392453,
                              0.22797763497031857388689019863932,
                              0.15284420981300738125325845173708,
                            -0.035295616033358979845090857076833,
                           -0.0099799992651306256792765875736422,
                            0.0018142609366514010419902709448759,
                            0.0003235891046425889772220595605353,
                          -0.00004841898485851272568247586094686,
                        -0.0000062916348373359738400767136988387,
                        0.00000079942476664286739130267762585159,
                       0.000000081592939856988078921076285919465,
                     -0.0000000090060168465171482923300334739558,
                    -0.00000000075626888288033084632885102808322,
                    0.000000000073775551261977522865908802683463,
                   0.0000000000052601044177443783363134989761516,
                 -0.00000000000045967798937305825795095961441697,
                -0.000000000000028468224107434042976887575897462,
                0.0000000000000022529894502335234532480058641737,
               0.00000000000000012330257543067426162304448507358,
            -0.0000000000000000089162571254206527010196106770777,
            -0.0000000000000000004370789114627443101549370861892,
           0.000000000000000000029094922022045915240836975369384,
          0.0000000000000000000012914245410749442185935885510156,
       -0.000000000000000000000079637265917628544002457121314317,
      -0.0000000000000000000000032293263813559717959165648165844,
      0.00000000000000000000000018548302717794476652598170488439,
    0.0000000000000000000000000069227609604047599589831775451484,
  -0.00000000000000000000000000037144957948853240739662579526224,
 -0.000000000000000000000000000012863369539271733425846532217877};
 return clenshaw(coeffs, 30, 2./5.*x-1.);
}

double chebbessj1( double x )
{
   double coeffs[30] = {0.1284551818460639559996165057131,
                             -0.26022312539939923941511375346869,
                             -0.32867361592127520313615266371342,
                              0.10454109055311047880391056435423,
                             0.038152487629942511871667620455564,
                           -0.0084048807536382567003801782916386,
                           -0.0017675094305799908454387298390052,
                           0.00030357174228846830117312224376546,
                          0.000044589555418507427004803699992487,
                         -0.000006309760806013143194723266294448,
                        -0.0000007102154103115846437486386391493,
                       0.000000085637327129795935698154712385195,
                      0.0000000078024604299104507568326769727303,
                    -0.00000000082043459676868790821360894454719,
                   -0.000000000062735952044990044987373716210529,
                   0.0000000000058515773654603478845696247289584,
                  0.00000000000038530096794249504838824312898447,
                -0.000000000000032300898514797817202651808191342,
               -0.0000000000000018668799186079360973921773757747,
               0.00000000000000014214956856492052415261667849562,
             0.0000000000000000073192279383126793066919948260834,
           -0.00000000000000000051054544180991906063316071531991,
          -0.000000000000000000023697774261425138825268836475081,
          0.0000000000000000000015251857780890269937688432598129,
        0.000000000000000000000064437294353830434183774921444486,
      -0.0000000000000000000000038497275295236734422588484518447,
     -0.00000000000000000000000014923327190038186991290092391468,
    0.0000000000000000000000000083194820176135439148972849450031,
   0.00000000000000000000000000029739429087342291692752577774426,
 -0.000000000000000000000000000015567380220870463958716645537429};
   return clenshaw(coeffs, 30, 2./5.*x-1.);
}


double clenshaw( double* c, int n, double x) {
    double bk = 0., bk1 = 0., bk2 = 0.;
    int k;

    x *= 2.;
    for (k = n-1; k > 0; k-- ) { 
        bk = c[k] + x*bk1 - bk2;
        bk2 = bk1; 
        bk1 = bk;
    }
    return c[0] + .5*x*bk1 - bk2;
}



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *j, mu, x, z;
  /*  Check for proper number of arguments */
  if(nrhs!=3) 
    mexErrMsgTxt("Two inputs required.");
  
  /*  Create a pointer to the input matrix x  */
  mu = mxGetScalar(prhs[0]);
  x = mxGetScalar(prhs[1]);
  z = mxGetScalar(prhs[2]);
  
  /*  Set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
  
  /*  Create a C++ pointer to a copy of the output matrix */
  j = mxGetPr(plhs[0]);
  
  /*  Call the C++ subroutine 
  *j = besselj(mu,x,z); */
  besselj_taylor(mu,x,z,j);  
  
  return;
}