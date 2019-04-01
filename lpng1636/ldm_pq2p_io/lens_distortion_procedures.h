/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef LENS_DISTORTION_PROCEDURE_H
#define LENS_DISTORTION_PROCEDURE_H

#include "ami_point2d.h"
#include "ami_lens_distortion_model.h"
#include "ami_line_points.h"
#include "ami_image.h"
#include "ami_utilities.h"
#include <vector>
using namespace std;

//------------------------------------------------------------------------------
/**
 * \fn double model_center_estimation_2p(vector< line_points > &lines,
                                  lens_distortion_model &d, int w, int h,
                                  vector<bool> v)
 * \brief Function to optimize the lens distortion model parameters and its center
 * \param [in] lines Primitives detected in the image
 * \param [in,out] d Lens distortion model
 * \param [in] w Image width
 * \param [in] h Image height
 * \param [in] v Bool vector which indicates the parameters to estimate
 * \return Returns the distance between all the points and their associated lines
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
double model_center_estimation_2p(vector< line_points > &lines,
                                  lens_distortion_model &d, int w, int h,
                                  vector<bool> v);

//------------------------------------------------------------------------------

int build_l1r_vector(std::vector<double> &l1r, 
                     double max_distance_corner,int Na, double *a);

//------------------------------------------------------------------------------

int build_l1r_quotient_vector(std::vector<double> &l1r,
                              double max_distance_corner, int Na, double *a);

//------------------------------------------------------------------------------

enum ImageAmplification{FIT_ALL, FIT_WIDTH, FIT_HEIGHT};

//------------------------------------------------------------------------------

//Computing l1r for division model using Newton-Raphson
static void compute_division_l1r(vector<double> &l1r, 
                                 const vector<double> &d1, 
                                 double max_distance_corner)
{
  l1r.resize((int)(max_distance_corner+1.5));
    l1r[0]=1;
    double r0=1;
    for(int m=1;m<(int)l1r.size();m++){
      //WE COMPUTE THE INVERSE USING NEWTON-RAPHSON
      double h=1e-5;
      for(int i=0;i<1000;i++){
        //EVALUATION OF LENS DISTORTION MODEL AND DERIVATIVE
        double r2=r0*r0;
        double sum=d1[0];
        for(int k=1;k<(int)d1.size();k++){
          sum+=d1[k]*r2;
          r2*=r0*r0;
        }
        double f_r0=r0/sum;
        //DERIVATIVE
        r2=(r0+h)*(r0+h);
        sum=d1[0];
        for(int k=1;k<(int)d1.size();k++){
          sum+=d1[k]*r2;
          r2*=(r0+h)*(r0+h);
        }
        double f_r0_h=(r0+h)/sum;
        double derivative=(f_r0_h-f_r0)/h;
        // WE COMPUTE THE NEW ROOT
        double r1=r0-(f_r0-m)/derivative;
        if(fabs(r0-r1)<fabs(r0)*1e-5){
          r0=r1;
          break;
        }
        r0=r1;
      }
      l1r[m]=r0/m;
    }
}

//------------------------------------------------------------------------------

//Computing the scale factor according to the selected value of 
//ImageAmplification
static double compute_scale(const ImageAmplification amp, int width, int height, 
                            const lens_distortion_model d)
{
  point2d<double> ldm_center = d.get_distortion_center();
  double scale = 0., scale2 = 0.;
  switch(amp)
  {
    case FIT_ALL: //All image is included in corrected one
    {
      ami::point2d<double> temp2 = 
        d.evaluation(ami::point2d<double>((double) width,(double) height));
      scale=(temp2-ldm_center ).norm()/ldm_center.norm();

      temp2=d.evaluation( ami::point2d<double>((double)width,(double)0.));
      scale2=(temp2-ldm_center ).norm()/ldm_center.norm();
      if(scale2<scale) scale=scale2;

      temp2=d.evaluation( ami::point2d<double>((double) 0.,(double)height));
      scale2=(temp2-ldm_center ).norm()/ldm_center.norm();
      if(scale2<scale) scale=scale2;

      temp2=d.evaluation( ami::point2d<double>((double) 0.,(double) 0.));
      scale2=(temp2-ldm_center ).norm()/ldm_center.norm();
      break;
    }
    case FIT_WIDTH: //Image is fitted to keep width
    {
      ami::point2d<double> temp2 = 
        d.evaluation( ami::point2d<double>((double) width,ldm_center.y));
      scale=(temp2-ldm_center ).norm()/(ldm_center.x);

      temp2=d.evaluation( ami::point2d<double>((double) 0.,ldm_center.y));
      scale2=(temp2-ldm_center ).norm()/(ldm_center.x);
      break;
    }
    case FIT_HEIGHT: //Image is fitted to keep height
    {
      ami::point2d<double> temp2 = 
        d.evaluation( ami::point2d<double>(ldm_center.x,(double) height));
      scale=(temp2-ldm_center ).norm()/(ldm_center.y);

      temp2=d.evaluation( ami::point2d<double>(ldm_center.x,(double) 0.));
      scale2=(temp2-ldm_center ).norm()/(ldm_center.y);
      break;
    }
  }
  if(scale2<scale) scale=scale2;
  return scale;
}

//------------------------------------------------------------------------------
/**
 * \fn  template <class  U>
        ami::image<U> lens_distortion_procedures::undistort_image_inverse(
        ami::image<U> input_image,
        lens_distortion_model &d,
        const double &image_amplification_factor)
 * \brief ESTIMATE AN UNDISTORTED IMAGE USING DISTORTION MODEL 
          (inverse method)
 * \author Luis Alvarez
 */

template <class  U>
ami::image<U> undistort_image_inverse(ami::image<U> input_image,
                                      const lens_distortion_model &d,
                                      ImageAmplification amp
                                     )
{
  int width0=input_image.width();
  int height0=input_image.height();
  int size =width0*height0;
  int width=width0,height=height0;


  // WE CREATE OUTPUT IMAGE
  ami::image<U> output_image(width,height,0,0,0);

  //CALCULATE MAXIMUM DISTANCE FROM CENTRE TO A CORNER
  point2d<double> ldm_center=d.get_distortion_center();
  double max_distance_corner = sqrt(update_rsqmax(ldm_center,width,height));

  //BUILD INTERMEDIATE VECTOR
  vector<double> l1r;
  double *a;
  int Na;
  vector<double> d1=d.get_d();
  Na=2*(d1.size()-1);
  if(d1.size()<2) {output_image=input_image; return(output_image);}
  a=(double*)malloc(sizeof(double)*(Na+1));
  a[0]=d1[0];
  for(int i=1;i<(int)d1.size();i++){a[2*i-1]=0.; a[2*i]=d1[i]; }
  while(Na>0 && a[Na]==0) --Na;

  //We update the max_distance_corner according to lens distorsion max
  //displacement and the type of distortion model
  if(d.get_type()==POLYNOMIAL)
  {
    double step =0;
    double power=max_distance_corner;
    for(int k=0;k<=Na;k++){
      step+=power*a[k];
      power*=max_distance_corner;
    }
    if(step>max_distance_corner) max_distance_corner=step;
  }
  else
  {
    if(d1.size()==2){
    max_distance_corner=max_distance_corner/(d1[0]+d1[1]*max_distance_corner*
                                                           max_distance_corner);
    }
    else{
      max_distance_corner=max_distance_corner/(d1[0]+d1[1]*max_distance_corner*
                                                max_distance_corner+d1[2]*
                                                max_distance_corner*
                                                max_distance_corner*
                                                max_distance_corner*
                                                max_distance_corner);
    }
  }
  
  // WE BUILD THE LENS DISTORTION INVERSE VECTOR
  if(d.get_type()==POLYNOMIAL)
  {
    if(Na<2) {
      output_image=input_image;
      return(output_image);
    }
    if(build_l1r_vector(l1r,max_distance_corner,Na,a)==-1) {
      output_image=input_image;
      return(output_image);
    }
  }
  else
  {
    if(d1.size()==2){
      if(build_l1r_quotient_vector(l1r,max_distance_corner,Na,a)<0){
        output_image=input_image;
        return(output_image);
      }
    }
    else{
      compute_division_l1r(l1r, d1, max_distance_corner);
    }
  }
  
  // WE FIT IMAGE SCALING
  double scale = compute_scale(amp, width, height, d);
  ami::point2d<double> t = ldm_center*scale - ldm_center;

  int nc,n2,i,j;
  #ifdef _OPENMP
  #pragma omp parallel for \
   shared(width,height,width0,height0,output_image,input_image,size)\
   private(nc,i,j,n2)
  #endif
  for (nc=0;nc<3;nc++)
  {
    n2=nc*size;
    for(i=0;i<height;i++){
      for(j=0;j<width;j++){
        ami::point2d<double> temp(j*scale-t.x, i*scale-t.y);
        double distance_centre= (ldm_center-temp).norm();

        //INTERPOLATION
        int ind=(int)distance_centre;
        if(ind+1>=(int)l1r.size()) continue;
        double dl1r=l1r[ind]+(distance_centre-ind)*(l1r[ind+1]-l1r[ind]);
        ami::point2d<double> p;

        p.x=ldm_center.x+(temp.x-ldm_center.x)*dl1r;
        p.y=ldm_center.y+(temp.y-ldm_center.y)*dl1r;

        int m = (int)p.y;
        int n = (int)p.x;
        if(0<=m && m<height0 && 0<=n && n<width0)
        {
          //COLOUR INTERPOLATION
          double di=p.y-m;
          double dj=p.x-n;
          unsigned int k=i*width+j;
          unsigned int k0=m*width0+n;
          double accum=0;
          double w_accum=0;
          double w=((1.-di)*(1.-dj));

          accum+=(double)w*input_image[k0+n2];
          w_accum+=w;


          if( (di*(1.-dj))>0. && (m+1)<height0)
          {
            k0=(m+1)*width0+n;
            w=(di*(1.-dj));
            accum+=(double)w*input_image[k0+n2];
            w_accum+=w;
          }
          if( ((1-di)*(dj))>0. && (n+1)<width0)
          {
            k0=(m)*width0+n+1;
            w=(1.-di)*(dj);
            accum+=(double)w*input_image[k0+n2];
            w_accum+=w;
          }
          if( ((di)*(dj))>0. && (n+1)<width0 && (m+1)<height0)
          {
            k0=(m+1)*width0+n+1;
            w=(di)*(dj);
            accum+=(double)w*input_image[k0+n2];
            w_accum+=w;
          }
          if(w_accum>0.) output_image[k+n2]=(U) (accum/w_accum);
        }
      }
    }
  }
  free(a);
  return(output_image);
}

#endif
