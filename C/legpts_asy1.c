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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.141592653589793 // Pi (double)
#define EPSILON 2.2204460492503131e-16 // double

#define M 30
#define ITER_MAX 10

/* USAGE 
 
 ./legpts_asy1 N - compute nodes and weights of degree N. Display on screen.
 ./legpts_asy1 N x.txt w.txt - compute node and weights, and print to file.
 
*/

int main(int argc, char *argv[])
{
    unsigned long int n = atoi(argv[1]); // Convert command line input from str to int
    int s = n % 2;                       // True if n is odd
    unsigned long int n3, n4;
    int k, iter;
    double x, w, t, dt, C, theta, sint, f, fp, S, dS, dn = (double)n, p2, fn, fn5, nk, n5k;
    // Allocate memory for nodes and weights
    double *nodes = (double *) malloc(n * sizeof (double));
    double *weights = (double *) malloc(n * sizeof (double));
    FILE *fid, *fid2; 
    
    // The constant out the front of the sum in evaluation (compute only once)
    // C = sqrt(pi/4)*gamma(n+1)/gamma(n+3/2)
    //  
    // Naive expansion    
//    C = (4.0/PI);
//    for ( k = 1; k <= n; k++ ) {
//        C *= k/(k + 0.5);
//    }
    //
    // Use Stirling's series
    dS = -.125/dn;
    S = dS; 
    k = 1; 
    while ( fabs(dS/S) > EPSILON/100 ) {
        k += 1;
        dS *= -.5 * (k-1) / (k+1) / dn;
        S += dS;
    }
    double stirling[10] = {1., 1./12., 1./288., -139./51840., 163879./209018880., 5246819./75246796800., -534703531./902961561600., -4483131259./86684309913600., 432261921612371./514904800886784000.};
    fn = 1; fn5 = 1; nk = dn; n5k = dn+.5;
    for ( k = 1 ; k < 10 ; k++ ) {
        fn += stirling[k]/nk;
        fn5 += stirling[k]/n5k;
        nk *= dn;
        n5k *= dn + .5;
    }
    C = exp(S)*sqrt(4.0/(dn+.5)/PI) * fn / fn5;  
    
    // Loop over each root
    for ( k = (n+s)/2; k > 0; k-- ) {
        // Asymptotic approximation of roots (in theta-space)
        theta = PI*(4.0*k-1)/(4.0*n+2.0);
        n3 = n*n*n; n4 = n3*n; sint = sin(theta);
        x = (1.0 - (n-1.0)/(8*n3)-(39.0-28.0/(sint*sint))/(384.0*n4))*cos(theta);  
        t = acos(x);
        
        // Initialise
        dt = 1.0; 
        iter = 1;
        
        // Newton iteration for roots
        while (fabs(dt) > 10*EPSILON) {
            feval(n, t, C, &f, &fp);        // Evaluate at current approximation
            dt = f/fp;                      // Newton update
            t += dt;                        // Update t
            if (iter++ > ITER_MAX) { break; } 
        }
        
        // Convert back to x-space
        x = cos(t);
        w = 2./(fp*fp);        
        
        // Store nodes and weights
        nodes[n-k] = x;
        nodes[k-1] = -x;
        weights[n-k] = w;
        weights[k-1] = w;
    }
    
    // Enforce x = 0 for odd n
    if (s == 1) {
        nodes[(n+s)/2-1] = 0;
    }
    
    // Print
    if (argc < 3){             // Print the output to screen
        printf( "x = \n");
        for ( k = 0; k < n; k++ ) {
            printf( "%16.16f\n", nodes[k]);
        }
        printf( "w = \n");
        for ( k = 0; k < n; k++ ) {
            printf( "%16.16f\n", weights[k]);
        }
    }
    if (argc > 2) {           // Print x to file
        fid = fopen(argv[2],"w");
        for ( k = 0; k < n; k++ ) {
            fprintf( fid, "%16.16f\n", nodes[k]);
        } 
        fclose(fid);
    }
    if (argc > 3) {           // Print w to file
        fid2 = fopen(argv[3],"w");
        for ( k = 0; k < n; k++ ) {
            fprintf( fid2, "%16.16f\n", weights[k]);
        } 
        fclose(fid2);
    }
    
    // Free memory
    free(nodes);   
    nodes = NULL;
    free(weights); 
    weights = NULL;

    return 0;
}

int feval(int n, double theta, double C, double* f, double* fp, double* winv)
{
    double sinT = sin(theta), cotT = cos(theta);
    double alpha, cosA, sinA, denom;
    int m;
    
    // m = 1
    alpha = (n + 0.5)*theta - 0.25*PI;
    denom = sqrtl(2.0*sinT);
    cosA = cos(alpha);
    sinA = sin(alpha);
    *f = C * cosA/denom;
    *fp = C * ( 0.5*(cosA*cotT+sinA) + n*sinA) / denom;

    // Loop over m until coefficient is sufficiently small
    for ( m = 1; m <= M; m++ ) {
        C *= (2.0*m-1.0)/(2.0*m)*(m-0.5)/(n+m+0.5);
        alpha += theta - 0.5*PI;
        denom *= 2.0*sinT;
        cosA = cos(alpha);
        sinA = sin(alpha);
        *f += C * cosA/denom;
        *fp += C * ( (m+0.5)*(cosA*cotT+sinA) + n*sinA) / denom;
        
        if (fabs(C) < EPSILON){ break; }
    }
        
}
