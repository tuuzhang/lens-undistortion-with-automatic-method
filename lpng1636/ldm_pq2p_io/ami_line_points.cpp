/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

/**
 * \file line_points.cpp
 * \brief line_points class AMI_DLL_H basic methods
 * \author Luis Alvarez \n \n
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include "ami_line_points.h"
using namespace ami;
using namespace std;



/**
 * \fn double line_points::points_to_equation()
 * \brief LINE ESTIMATION FROM POINT RETURN -1 IF IT DOES NOT WORK. OTHERWISE
          RETURN THE AVERAGE OF THE SQUARED DISTANCE OF THE POINTS TO THE LINE
 * \author Luis Alvarez
 */

/****************************************************************/
/* LINE ESTIMATION FROM POINT RETURN -1 IF IT DOES NOT WORK. OTHERWISE
RETURN THE AVERAGE OF THE SQUARED DISTANCE OF THE POINTS TO THE LINE */
/****************************************************************/
AMI_DLL_CPP double line_points::points_to_equation()
{
  int i,j,k;
  long double suu,suv,svv,um,vm,h,r[4][3],min,paso,norma;
  long double cero=10e-100;
  int N=points.size();

  if(N<2){
    printf("Numero de puntos para el Calculo de la recta 2D menor que 2\n");
    return(-1.);
  }

  suu=0; suv=0; svv=0; um=0; vm=0;
  for(i=0;i<N;i++){
    um+=points[i].x;
    vm+=points[i].y;
  }
  um/=N; vm/=N;
  for(i=0;i<N;i++){
    suu+=(points[i].x-um)*(points[i].x-um);
    svv+=(points[i].y-vm)*(points[i].y-vm);
    suv+=(points[i].x-um)*(points[i].y-vm);
  }
  suu/=N; svv/=N; suv/=N;
  if(fabs(suv)<= cero){
    if(suu<svv && svv>cero){
      //a=1; b=0; c=-um;
      rect.set_a((double) 1.);
      rect.set_b((double) 0.);
      rect.set_c((double) -um);
      return(0);
    }
    if(svv<suu && suu>cero){
      //a=0; b=1; c=-vm;
      rect.set_a((double) 0.);
      rect.set_b((double) 1.);
      rect.set_c((double) -vm);
      return(0);
    }
    printf("No se pudo calcular la recta 2D\n");
    return(-1);
  }

  r[2][1]=r[3][1]=r[0][0]=r[1][0]=1.;
  h=0.5*(suu-svv)/suv;
  if(h>0){
    r[0][1]=-h-sqrt(1.+h*h);
    r[0][2]=-(um+r[0][1]*vm);
    r[1][1]=-1./r[0][1];
    r[1][2]=-(um+r[1][1]*vm);
    r[2][0]=h+sqrt(1.+h*h);
    r[2][2]=-(r[2][0]*um+vm);
    r[3][0]=-1./r[2][0];
    r[3][2]=-(r[3][0]*um+vm);
  }
  else{
    r[0][1]=-h+sqrt(1+h*h);
    r[0][2]=-(um+r[0][1]*vm);
    r[1][1]=-1./r[0][1];
    r[1][2]=-(um+r[1][1]*vm);
    r[2][0]=h-sqrt(1+h*h);
    r[2][2]=-(r[2][0]*um+vm);
    r[3][0]=-1./r[2][0];
    r[3][2]=-(r[3][0]*um+vm);
  }

  for(j=0;j<4;j++){
    norma=sqrt(r[j][0]*r[j][0]+r[j][1]*r[j][1]);
    for(i=0;i<3;i++)
      r[j][i]/=norma;
  }

  min=0.; k=0;
  for(i=0;i<N;i++){
    paso=r[0][0]*points[i].x+r[0][1]*points[i].y+r[0][2];
    min+=paso*paso;
  }
  for(j=1;j<4;j++){
    h=0;
    for(i=0;i<N;i++){
      paso=r[j][0]*points[i].x+r[j][1]*points[i].y+r[j][2];
      h+=paso*paso;
    }
    if(h<min){
      k=j;
      min=h;
    }
  }

  rect.set_a((double) r[k][0]);
  rect.set_b((double) r[k][1]);
  rect.set_c((double) r[k][2]);
  return(min/N);
}
 
  /**
 * \fn double line_points::distance(point2d<double> &p )
 * \brief Return the min distance of a point to the collection of points of the line points
 * \author Luis Alvarez
 */
AMI_DLL_CPP double line_points::distance(point2d<double> &p /** point2d */){
  double min=1e20;
  for(int m=0;m<(int)points.size();m++){
    double temp=(points[m].x-p.x)*(points[m].x-p.x)+(points[m].y-p.y)*(points[m].y-p.y);
    if(temp<min) min=temp;
  }
  if(min>0) return(sqrt(min));
  else return(0.);
}