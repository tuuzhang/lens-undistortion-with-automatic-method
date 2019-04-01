/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file image_contours.h
 * \brief Subpixel image contour class AMI_DLL_H definition
 * \author Luis Alvarez \n \n
*/
#ifndef SUBPIXEL_IMAGE_CONTOURS_H
#define SUBPIXEL_IMAGE_CONTOURS_H

#include <stdlib.h>
//#include "../auxiliary/ami_image.h"
#include "ami_point2d.h"


#ifdef __AMIDEBUG__
  #include "wxAmiDebugLog.h"
#endif
using namespace std;

namespace ami
{

template<class T >
class image;

// class AMI_DLL_H OF SUBPIXEL PRECISION CONTOURS
/**
 * \class  image_contours
 * \brief class to store subpixel contours extracted from soccer stadium images.
 * \author Luis Alvarez
 */
class AMI_DLL_H   image_contours{
  bool *c /** Image contour c[i]=1 in contours points */;
  float *x /** Image with horizontal subpixel precision contour value */;
  float *y /** Image with vertical subpixel precision contour value*/;
  float *d /** Distance of the center line point to the white line border
          contour*/;
  float *cosine /** x orientation of the center lines*/;
  float *sine /** y orientation of the center lines*/;
  int width/** Image width*/;
  int height /** Image height */;
  vector<int> index; /** Index vector of the edges position */

public:
 
  /**
  * \fn image_contours()
  * \brief Constructor without taking memory
  * \author Luis Alvarez
  */
  image_contours(){
    width = 0;
    height = 0;
    c = NULL;
    x = NULL;
    y = NULL;
    d = NULL;
    cosine = NULL;
    sine = NULL;
  }
 
  /** 
   * \fn image_contours(int width_c,int height_c)
   * \brief Constructor taking memory
   * \author Luis Alvarez
   */
  image_contours(int width_c,int height_c);
  
  /** \fn ~image_contours()
   * \brief Destructor to free memory
   * \author Luis Alvarez
   */
  ~image_contours();
  
  /** 
   * \fn image_contours & operator=(const image_contours &contours)
   * \brief Assignment operator
   * \author Daniel Santana-Cedrés
   */
  image_contours & operator=(const image_contours &contours);

  /**
  * \fn bool *get_c()
  * \brief Return array c to identity contour points
  * \author Luis Alvarez
  */
  bool *get_c(){return c;}

  const bool *get_c() const {return c;}

  /**
  * \fn float *get_x()
  * \brief Return array x of subpixel x coordinate location
  * \author Luis Alvarez
  */
  float *get_x(){return x;}

  const float *get_x() const {return x;}

  /**
  * \fn float *get_y()
  * \brief Return array y of subpixel y coordinate location
  * \author Luis Alvarez
  */
  float *get_y(){return y;}

  const float *get_y() const {return y;}

  /**
  * \fn float *get_d()
  * \brief Return array d of distance to a contour pixel to the boundary of
  *        contour pixel area
  * \author Luis Alvarez
  */
  const float *get_d() const {return d;}

  /**
  * \fn const float *get_cosine()
  * \brief Return array cosine of x coordinate contour point orientation
  * \author Luis Alvarez
  */
  float *get_cosine(){return cosine;}

  const float *get_cosine() const {return cosine;}

  /**
  * \fn float *get_sine()
  * \brief Return array sine of y coordinate contour point orientation
  * \author Luis Alvarez
  */
  float *get_sine(){return sine;}

  const float *get_sine() const {return sine;}

  /**
  * \fn const int get_width() const
  * \brief Return image width
  * \author Luis Alvarez
  */
  int get_width() const {return width;}

  /**
  * \fn const int get_height() const
  * \brief Return image height
  * \author Luis Alvarez
  */
  int get_height() const {return height;}

  /**
  * \fn vector<int> get_index()
  * \brief This function returns the index vector of the edges
  * \author Luis Alvarez  
  */
  const vector<int>& get_index()const {return index;}

  /**
  * \fn void set_index(vector<int> &index2)
  * \brief This method assigns an index vector with the position of the edges
  * \author Luis Alvarez  
  */
  void set_index(vector<int> &index2){index=index2;}

  /**
  * \fn void clean(const int neighborhood_radius,const int min_neighbor_points,
            const double min_orientation_value,const int min_distance_point)
  * \brief  The method cleans the contours according to their curvature and the
            stability of their orientation
  * \author Luis Alvarez  
  */
  void clean(const int neighborhood_radius,const int min_neighbor_points,
             const double min_orientation_value,const int min_distance_point);

};


}

#endif
