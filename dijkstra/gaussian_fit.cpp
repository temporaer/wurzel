////////////////////////////////////////////////////////////////////////////////////
//  Example program that shows how to use levmar in order to fit the three-
//  parameter exponential model x_i = p[0]*exp(-p[1]*i) + p[2] to a set of
//  data measurements; example is based on a similar one from GSL.
//
//  Copyright (C) 2008  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////////
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "gaussian_fit.hpp"
#include "../../../plain/levmar-2.5/levmar.h"
 
#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif
 
 
/* the following macros concern the initialization of a random number generator for adding noise */
#undef REPEATABLE_RANDOM
#define DBL_RAND_MAX (double)(RAND_MAX)
 
#ifdef _MSC_VER // MSVC
#include <process.h>
#define GETPID  _getpid
#elif defined(__GNUC__) // GCC
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#else
#warning Do not know the name of the function returning the process id for your OS/compiler combination
#define GETPID  0
#endif /* _MSC_VER */
 
#ifdef REPEATABLE_RANDOM
#define INIT_RANDOM(seed) srandom(seed)
#else
#define INIT_RANDOM(seed) srandom((int)GETPID()) // seed unused
#endif
 
/* Gaussian noise with mean m and variance s, uses the Box-Muller transformation */
double gNoise(double m, double s)
{
double r1, r2, val;
 
  r1=((double)random())/DBL_RAND_MAX;
  r2=((double)random())/DBL_RAND_MAX;
 
  val=sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
 
  val=s*val+m;
 
  return val;
}
 
/* model to be fitted to measurements: x_i = p[0]*exp(-p[1]*i) + p[2], i=0...n-1 */
void expfunc(double *p, double *x, int m, int n, void *data)
{
register int i;
 
  for(i=0; i<n; ++i){
    x[i]=p[0]*exp(-p[1]*i) + p[2];
  }
}
void expfunc_coor(double *p, double *x, int m, int n, void *data)
{
register int i;
  double* d = (double*)data;
  //double startval = *d;
  //d++;
 
  for(i=0; i<n; i+=1){

    //x[i]=p[0]*exp(-p[1]*d[i]) + p[2];
      //x[i]=p[0]*exp(-p[1]*d[i]);
      x[i]=p[0]*exp(-p[1]*d[i]);

  }
}
/* Jacobian of expfunc_coor() */
void jacexpfunc_coor(double *p, double *jac, int m, int n, void *data)
{   
register int i, j;
  double* d = (double*)data;
 
  /* fill Jacobian row by row */
  for(i=j=0; i<n; ++i){
    jac[j++]=exp(-p[1]*d[i]);
    jac[j++]=-p[0]*d[i]*exp(-p[1]*d[i]);
    //jac[j++]=-d[i]*exp(-p[0]*d[i]);
    //jac[j++]=1.0;
  }
}
 
/* Jacobian of expfunc() */
void jacexpfunc(double *p, double *jac, int m, int n, void *data)
{   
register int i, j;
 
  /* fill Jacobian row by row */
  for(i=j=0; i<n; ++i){
    jac[j++]=exp(-p[1]*i);
    jac[j++]=-p[0]*i*exp(-p[1]*i);
    jac[j++]=1.0;
  }
}
void fit_gauss_curve(double* p, double* x, int n, double* coor){
	const int m = 2;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

	/* initial parameters estimate: (1.0, 0.0, 0.0) */
	//p[0]=1.0; p[1]=0.0; p[2]=0.0;

	/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

	/* invoke the optimization function */
	int ret=dlevmar_der(expfunc_coor, jacexpfunc_coor, p, x, m, n, 1000, opts, info, NULL, NULL, coor); // with analytic Jacobian
	//int ret=dlevmar_dif(expfunc_coor, p, x, m, n, 1000, opts, info, NULL, NULL, coor); // without Jacobian
	//printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
	//printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], p[2]);
}
 
//int main()
//{
//const int n=40, m=3; // 40 measurements, 3 parameters
//double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ], coor[2*n];
//register int i;
//int ret;
 
  /* generate some measurement using the exponential model with
   * parameters (5.0, 0.1, 1.0), corrupted with zero-mean
   * Gaussian noise of s=0.1
   */
  //INIT_RANDOM(0);
  //for(i=0; i<n; ++i){
    //x[i]=(5.0*exp(-0.1*i) + 1.0) + gNoise(0.0, 0.1);
    //coor[2*i]   = i;
    //coor[2*i+1] = 2;
  //}
 
  //fit_gauss_curve(p, x, n, coor);
 
  //exit(0);
//}
