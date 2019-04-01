/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

 #ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

/**
 * \file lens_distortion.cpp
 * \brief Functions for lens distortion model basic operations
 * \author Luis Alvarez \n \n
 */


#include <stdio.h>
#include <math.h>
#include "ami_lens_distortion.h"
#include "ami_utilities.h"
#include "ami_pol.h"

/*******************************************************************************
  FUNCTION TO COMPUTE THE REAL ROOTS OF A CUBIC POLYNOMIAL. IT RETURNS
  THE NUMBER OF ROOTS FOUND SORTED BY MAGNITUD
*******************************************************************************/
AMI_DLL_CPP int ami_RootCubicPolynomial(
  double *a, /* POLINOMIAL COEFICIENTS a[0]+a[1]x+a[2]x^2 +... */
  int N, /* DEGREE OF POLINOMIAL (IT HAS TO BE 3) */
  double *x) /* POLINOMIAL ROOTS */
{
  double a1,a2,a3,Q,R,S,T,D,A;
  if(N!=3 || a[3]==0) return(-100000);
  a1=a[2]/a[3];
  a2=a[1]/a[3];
  a3=a[0]/a[3];
  Q=(3*a2-a1*a1)/9.;
  R=(9*a1*a2-27*a3-2*a1*a1*a1)/54.;
  D=Q*Q*Q+R*R;
  if(D>0){
    S=R+sqrt(D);
    T=R-sqrt(D);
    if(S>0) S=pow(S,(double)1./3.);
    else S=-pow(-S,(double)1./3.);
    if(T>0) T=pow(T,(double)1./3.);
    else T=-pow(-T,(double)1./3.);
    x[0]=S+T-a1/3.;
    return(1);
  }
  else{
    double PI2=acos(-1.);
    if(Q!=0) A=acos(R/sqrt(-Q*Q*Q));
    else A=0;
    Q=2.*sqrt(-Q);
    x[0]=Q*cos(A/3.)-a1/3.;
    x[1]=Q*cos(A/3.+2.*PI2/3.)-a1/3.;
    x[2]=Q*cos(A/3+4.*PI2/3.)-a1/3.;

    if(fabs(x[0])>fabs(x[1])){ Q=x[1]; x[1]=x[0]; x[0]=Q; }
    if(fabs(x[0])>fabs(x[2])){ Q=x[2]; x[2]=x[0]; x[0]=Q; }
    if(fabs(x[1])>fabs(x[2])){ Q=x[2]; x[2]=x[1]; x[1]=Q; }

    return(3);
  }
}

/*******************************************************************************
  EVALUATION OF A POLYNOM USING HORNER ALGORITHM
*******************************************************************************/
AMI_DLL_CPP double ami_polynomial_evaluation(
  double *a, /* POLYNOM COEFICIENT */
  int Na, /* POLYNOM DEGREE */
  double x) /* POINT WHERE THE POLYNOM IS EVALUATED */
{
  double sol=a[Na];
  int i;
  for(i=Na-1;i>-1;i--) sol=sol*x+a[i];
  return(sol);
}

/*******************************************************************************
 COMPUTE THE LENS DISTORTION MODEL IN A POINT
*******************************************************************************/
AMI_DLL_CPP void ami_lens_distortion_model_evaluation(
  double *a, // INPUT POLINOMIAL DISTORTION MODEL
  int Na, // INPUT DEGREE OF POLINOMIAL DISTORTION MODEL
  double xc,double yc,  // INPUT CENTER OF DISTORTION
  double x_input,double y_input,  // INPUT POINT
  double *x_output,double *y_output  // OUTPUT UNDISTORTED POINT
  )
{
  double norm=sqrt((x_input-xc)*(x_input-xc)+(y_input-yc)*(y_input-yc));
  double A=ami_polynomial_evaluation(a,Na,norm);
  *x_output=xc+(x_input-xc)*A;
  *y_output=yc+(y_input-yc)*A;
}

/*******************************************************************************
FUNCION AUXILIAR ami_inverse_lens_distortion()
*******************************************************************************/
AMI_DLL_CPP double ami_inverse_lens_distortion_newton_raphson(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na) /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
{
  int i;
  double paso,d,*b,*b2,root;
  if(a[Na]==0.) return(-1);
  /* WE ALLOCATE MEMORY */
  b=(double*)malloc( sizeof(double)*(Na+2) );
  b2=(double*)malloc( sizeof(double)*(Na+2) );
  d=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  b[0]=-1; b[1]=a[0];
  paso=d;
  for(i=2;i<(Na+2);i++){
    b[i]=a[i-1]*paso;
    paso*=d;
  }
  for(i=0;i<(Na+2);i++) b2[i]=b[Na+1-i];

  // WE IMPROVE SOLUTION USING NEWTON RAPHSON
  double norm2;
  double xn=x,yn=y,xp,yp,xn2=0.,yn2=0.,xn3,yn3,error,tol=1e-8;
  ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn,yn,&xp,&yp);
  norm2=(xp-x)*(xp-x)+(yp-y)*(yp-y);
  error=tol+1;
  int iter=0;
  root=1;
  double lambda=1;
  while(error>tol && ++iter<100){
    xn2=x0+(xn-x0)*(root+1e-6);
    yn2=y0+(yn-y0)*(root+1e-6);
    xn3=x0+(xn-x0)*(root-1e-6);
    yn3=y0+(yn-y0)*(root-1e-6);
    ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn2,yn2,&xp,&yp);
    double norm3=(xp-x)*(xp-x)+(yp-y)*(yp-y);
    ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn3,yn3,&xp,&yp);
    double norm4=(xp-x)*(xp-x)+(yp-y)*(yp-y);
    double derivative=(norm3-norm4)/2e-6;
    double derivative2=(norm3+norm4-2.*norm2)/1e-12;
    if(derivative2==0) break;
    double root2=root-lambda*derivative/derivative2;
    error=fabs(root2-root);
    xn2=x0+(xn-x0)*root2;
    yn2=y0+(yn-y0)*root2;
    ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn2,yn2,&xp,&yp);
    norm3=(xp-x)*(xp-x)+(yp-y)*(yp-y);
    if(norm3<norm2){
      root=root2;
      norm2=norm3;
      if(lambda<1.) lambda*=2.;
    }
    else {
      lambda/=2.;
    }
  }
  *xt=xn2; *yt=yn2;

  free(b); free(b2);
  return(norm2);
}

/*******************************************************************************
  FUNCTION TO INVERSE THE LENS DISTORTION TRANSFORMATION
*******************************************************************************/
AMI_DLL_CPP int ami_inverse_lens_distortion(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na) /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
{
  int i,Nr;
  double paso,d,*b,*b2,*rx,*ry,root;
  if(a[Na]==0.) return(-1);
  /* WE ALLOCATE MEMORY */
  b=(double*)malloc( sizeof(double)*(Na+2) );
  b2=(double*)malloc( sizeof(double)*(Na+2) );
  rx=(double*)malloc( sizeof(double)*(Na+2) );
  ry=(double*)malloc( sizeof(double)*(Na+2) );
  /* WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS */
  d=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  b[0]=-1; b[1]=a[0];
  paso=d;
  for(i=2;i<(Na+2);i++){
    b[i]=a[i-1]*paso;
    paso*=d;
  }

  if(Na==2){
    Nr=ami_RootCubicPolynomial(b,3,rx);
    for(i=0;i<Na;i++) ry[i]=0;
  }
  else{
    for(i=0;i<(Na+2);i++) b2[i]=b[Na+1-i];
    Nr=ami_polynomial_root(b2,Na+1,rx,ry);
  }
  /* WE SELECT THE REAL ROOT NEAR TO 1 */
  root=10e5;
  for(i=0;i<Nr;i++){
    if(fabs(ry[i])<0.00000000001 && fabs(root-1)>fabs(rx[i]-1)) root=rx[i];
  }
  if(Nr==0){
    root=1.;
  }

  /* WE TRANSFORM THE POINT COORDINATES */
  *xt=x0+(x-x0)*root;
  *yt=y0+(y-y0)*root;

  double x2,y2;
  double error=ami_inverse_lens_distortion_newton_raphson(x,y,x0,y0,&x2,&y2,a,Na);

  if(error<2.){
    if( ((x-x2)*(x-x2)+(y-y2)*(y-y2)) < ((x-(*xt))*(x-(*xt))+(y-(*yt))*(y-(*yt))) ){
      *xt=x2; *yt=y2;
    }
  }

  free(b); free(rx); free(ry); free(b2);
  return(0);
}

/*******************************************************************************
  FUNCTION TO INVERSE THE LENS DISTORTION TRANSFORMATION
*******************************************************************************/
AMI_DLL_CPP int ami_inverse_lens_distortion_fast(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na,/* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
  double dl1r /*COEFICIENT INTERPOLATED FROM VECTOR OF MAX DISTORSION DISTANCES*/)
{
  /* WE TRANSFORM THE POINT COORDINATES */
  *xt=x0+(x-x0)*dl1r;
  *yt=y0+(y-y0)*dl1r;

  double x2,y2;
  double error=ami_inverse_lens_distortion_newton_raphson(x,y,x0,y0,&x2,&y2,a,Na);

  if(error<2.){
    if(((x-x2)*(x-x2)+(y-y2)*(y-y2)) < ((x-(*xt))*(x-(*xt))+(y-(*yt))*(y-(*yt))))
    {
      *xt=x2; *yt=y2;
    }
  }

  return(0);
}
