/*
 * Copyright (c) 2010-2011, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */


/**
 * \file image.h
 * \brief Class to store multichannel image
 * \author Luis Alvarez, Pedro Henríquez \n \n
*/
#ifndef image_H
#define image_H

#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <typeinfo>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_png.h"

namespace ami {

// ITERATOR ALIAS

typedef std::vector<unsigned char>::iterator usIt;
typedef std::vector<int>::iterator iIt;
typedef std::vector<short>::iterator sIt;
typedef std::vector<float>::iterator fIt;
typedef std::vector<double>::iterator dIt;

typedef std::vector<unsigned char>::const_iterator uscIt;
typedef std::vector<int>::const_iterator icIt;
typedef std::vector<short>::const_iterator scIt;
typedef std::vector<float>::const_iterator fcIt;
typedef std::vector<double>::const_iterator dcIt;

}


namespace ami {

////////////////////////////////////////////////////////////////////////////////
/**
 * \class image
 * \brief Class to store multichannel images and basic methods
 * \author Luis Alvarez
 */
template <class T>
class image{
  std::vector <T> image_          /** std vector to allocate image "channel 1 + 
                                  channel 2 + ..." */;
  int width_                      /** Image width */;
  int height_                     /** Image height */;
  int nChannels_                  /** Number of image channels */;
  std::vector <int> roi_          /** ROI : Region Of Interest =(x0,x1,y0,y1,c0,c1)
                                  subwindow corners used (if defined) to operate only
                                  in a subwindow area */;

  std::vector <float> origin_     /** World coordinate origin */;
  std::vector <float> pixel_size_ /** Pixel size */;

  public:

  // CONSTRUCTOR - DESTRUCTOR METHODS
  image();
  image(const int width,const int height);
  image(const T *vector,const int width,const int height,const int channels);
  image(const int width,const int height, const int nChannels, const T &a);
  image(const int width,const int height, const T &a, const T &b, const T &c);
  image(const int width,const int height, const int nChannels);
  image(std::string name /** Input file name */ );


  image(const image<bool> &img_in, bool isbool);
  image(const image<T> &img_in);
  ~image();

  // CLASS ELEMENT ACCESS METHODS
  inline std::vector <T> *get_image(){return &image_;}
  inline int width() const {return width_;}
  inline int height() const {return height_;}
  inline int nChannels() const {return nChannels_;}
  inline int size() const {return image_.size();}
  inline T& operator[](const int &i);
  inline T& operator()(const int &x,const int &y,const int &channel);
  inline const T& operator[](const int &i) const;

  // BASIC OPERATION METHODS
  int write(std::string name);
  int read(std::string name);
  image & operator=(const image &image2);
  void clear();
  void init(const int width,const int height, const int nChannels);
  void set_nchannels(int nchannel1);
  void set_size(int width,int height);

  // BASIC METHODS TO MANIPULATE SUBWINDOWS (ROI : REGION OF INTEREST )
  std::vector <int> get_roi() const{return roi_;}
  template <class U> void  get_roi_image(image<U> &image2) const;
  image<T>  const get_roi_image(const std::vector <int> &roi2) const;

  // BASIC METHODS TO MANIPULATE IMAGE WITH ITERATORS
  typedef typename std::vector <T>::iterator iterator;
  inline iterator begin() {return(image_.begin());}
  inline iterator end() {return(image_.end());}


};

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> void image<T>::init(const int width,const int height,
                                              const int nChannels)
 * \brief Initializes the image taking memory
 * \author Carlos Falcon
*/
template <class T>
void image<T>::init(const int width,const int height, const int nChannels)
{
  width_     = width;
  height_    = height;
  nChannels_ = nChannels;
  int size   = width_*height_*nChannels_;
  image_.resize(size);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(const T *vector,const int width,
                                          const int height,const int channels)
 * \brief Image from T vector
 * \author Carlos Falcon
*/

template <class T>
image<T>::image(const T *vector,const int width,const int height,
                const int channels)
{
  width_=width ;
  height_=height;
  nChannels_=channels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
  size=width_*height_;
  #ifdef _OPENMP
  #pragma omp parallel \
    shared(size, vector)
  #endif
  for(int n=0;n<nChannels();n++){
    int n2_=n*size;
    #ifdef _OPENMP
    #pragma omp for nowait
    #endif
    for(int i=0;i<size;i++)
    {
      image_[i+n2_]=vector[i+n2_];
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(const image<T> &img_in)
 * \brief Image copy constructor
 * \author Carlos Falcon
*/

template <class T>
image<T>::image(const image<T> &img_in)
{
  width_=img_in.width() ;
  height_=img_in.height();
  nChannels_=img_in.nChannels();
  int size=width_*height_*nChannels_;
  image_.resize(size);
  size=width_*height_;
  #ifdef _OPENMP
  #pragma omp parallel \
    shared(size)
  #endif
  for(int n=0;n<img_in.nChannels();n++){
    int n2_=n*size;
    #ifdef _OPENMP
    #pragma omp for nowait
    #endif
    for(int i=0;i<size;i++)
    {
      image_[i+n2_]=img_in[i+n2_];
    }
  }

  if(img_in.get_roi().size() == 6){
    roi_.resize(6);
    roi_.at(0)=img_in.get_roi().at(0);
    roi_.at(1)=img_in.get_roi().at(1);
    roi_.at(2)=img_in.get_roi().at(2);
    roi_.at(3)=img_in.get_roi().at(3);
    roi_.at(4)=img_in.get_roi().at(4);
    roi_.at(5)=img_in.get_roi().at(5);
  }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>  const image<T>::
        get_roi_image(const std::vector <int> &roi2) const
 * \brief Function to get a subwindow image
 * \author Luis Alvarez
*/
template <class T>
image<T>  const image<T>::get_roi_image(const std::vector <int> &roi2) const
{
  int x0,x1,y0,y1,c0,c1;

  // IF THERE IS NOT ROI DEFINED WE RETURN AN EMPTY IMAGE
  if(roi2.size()<6){
    return(image<T>());
  }
  else{
    x0=roi2.at(0);
    x1=roi2.at(1);
    y0=roi2.at(2);
    y1=roi2.at(3);
    c0=roi2.at(4);
    c1=roi2.at(5);
  }

  int width,height;
  if( x1<= width_ &&  y1<= height_ && x0<x1 && y0<y1 &&  c1<= nChannels_ &&
      c0<c1){
    width=x1-x0;
    height=y1-y0;
  }
  else{
    return(image<T>());
  }

  image<T> image2(width,height,c1-c0);

  int size_=width_*height_;
  int size=width*height;
  int i,n,j,i2,i2_,n2,n2_;

  #ifdef _OPENMP
  #pragma omp parallel \
  shared(image2,x0,x1,y0,y1,c0,c1,width,size,size_) \
  private(n,i,j,n2,n2_,i2,i2_)
  #endif
  for(n=c0;n<c1;n++){
    n2=(n-c0)*size;
    n2_=n*size_;
    #ifdef _OPENMP
    #pragma omp for nowait
    #endif
    for(i=y0;i<y1;++i){
      i2_=n2_+i*width_;
      i2=n2+(i-y0)*width-x0;
      for(j=x0;j<x1;++j){
        image2[i2+j]= image_[i2_+j] ;
      }
    }
  }
  return(image2);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> template <class U> void  image<T>::
        get_roi_image(image<U> &image2) const
 * \brief Function to fill a subwindow image
 * \author Luis Alvarez
*/

template <class T> template <class U> void  image<T>::
  get_roi_image(image<U> &image2) const
{
  int x0,x1,y0,y1,c0,c1;

  // IF THERE IS NOT ROI DEFINED WE TAKE THE WHOLE IMAGE
  if(roi_.size()<6){
    x0=0;
    x1=width_;
    y0=0;
    y1=height_;
    c0=0;
    c1=nChannels_;
  }
  else{
    x0=roi_.at(0);
    x1=roi_.at(1);
    y0=roi_.at(2);
    y1=roi_.at(3);
    c0=roi_.at(4);
    c1=roi_.at(5);
  }

  if( x1<= width_ &&  y1<= height_ && x0<=x1 && y0<=y1 &&  c1<= nChannels_ &&
      c0<=c1){
    int width=x1-x0;
    int height=y1-y0;

    // WE ALLOCATE MEMORY IF NEEDED
    if(image2.width()!=width || image2.height()!=height ||
       image2.nChannels()!=nChannels_){
      image2=image<U>(width,height,c1-c0);
    }

    int size_=width_*height_;
    int size=width*height;
    int i,n,j,i2,i2_,n2,n2_;

    #ifdef _OPENMP
    #pragma omp parallel \
    shared(x0,x1,y0,y1,c0,c1,width,size,size_) \
    private(n,i,j,n2,n2_,i2,i2_)
    #endif
    for(n=c0;n<c1;n++){
      n2=(n-c0)*size;
      n2_=n*size_;
      #ifdef _OPENMP
      #pragma omp for nowait
      #endif
      for(i=y0;i<y1;++i){
        i2_=n2_+i*width_;
        i2=n2+(i-y0)*width-x0;
        for(j=x0;j<x1;++j){
          image2[i2+j]=(U) image_[i2_+j] ;
        }
      }
    }
    return;
  }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T> & image<T>::operator=(const image &image2)
 * \brief Operator = (equality of images of different size is not allowed )
 * \author Luis Alvarez
*/
template <class T>
image<T> & image<T>::operator=(const image &image2)
{
  if( width_==image2.width() && height_==image2.height() &&
      nChannels_==3 && image2.nChannels()==1){
    if(image_.size()!=(unsigned int) (3*width_*height_) ) image_.resize(3*width_*height_);
    int size_=width_*height_;
    int size_2=2*size_;
    for(int m=0;m<size_;m++){
      image_[m]         = image2[m];
      image_[m+size_]   = image2[m];
      image_[m+size_2]  = image2[m];
    }
    return *this;
  }

  if( width_!=image2.width() || height_!=image2.height() ||
      nChannels_!=image2.nChannels()){
    width_=image2.width();
    height_=image2.height();
    nChannels_=image2.nChannels();
    image_.resize(image2.size());
  }

  int k,k_end=image2.size();
  #ifdef _OPENMP
  #pragma omp parallel for shared(k_end) private(k)
  #endif
  for(k=0;k<k_end;++k)
    image_[k]=  image2[k];

  return *this;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> int image<T>::write( std::string name )
 * \brief Function to write an image to disk.
 * \param [in] name Filename of the image.
 * \author Pedro Henriquez and Daniel Santana-Cedrés
*/
template <class T>
int image<T>::write(std::string name /** Image file name */)
{
  int pos=name.find_last_of('.');
  int size_ = width_*height_;
  unsigned char *red = new unsigned char[size_],
                *green = new unsigned char[size_],
                *blue = new unsigned char[size_];

  if(pos == (int)std::string::npos) return -1;

  for(int i=0; i<size_; i++)
  {
    red[i]   = image_[i];
    green[i] = image_[i+size_];
    blue[i]  = image_[i+size_*2];
  }

  std::string extension=name.substr(pos+1);
  if( (extension == std::string("png")) || (extension == std::string("PNG")))
  {
    int output_value = ami_write_png(strdup(name.c_str()),red,green,blue,width_,
                                     height_);
    delete []red;
    delete []green;
    delete []blue;
    return output_value;
  }
  delete []red;
  delete []green;
  delete []blue;
  printf("WRITE::Unrecognized image format\n");
  return -1;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> T& image<T>::operator()(const int &x,const int &y,
                                                  const int &channel)
 * \brief Operator () to acces image value
 * \author Javier Martin
*/
template <class T>
T& image<T>::operator()(const int &x,const int &y,const int &channel){
  #ifdef IMAGE_DEBUG
    if((x+y*width_+channel*width_*height_)>=(image_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
                                        image_.size(),x+y*width_+
                                        channel*width_*height_);
      int j; scanf("%d",&j);
      exit(0);
    }
  #endif
  return image_[x+y*width_+channel*width_*height_];
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> const T& image<T>::operator[](const int &i) const
 * \brief Operator [] to acces image value
 * \author Luis Alvarez
*/
template <class T>
const T& image<T>::operator[](const int &i) const
{
  #ifdef IMAGE_DEBUG
    if(i>=(int)(image_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
             (int)image_.size(),i);
      int j; scanf("%d",&j);
      exit(0);
    }
  #endif
  return image_[i];
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> T& image<T>::operator[](const int &i)
 * \brief Operator [] to acces image value
 * \author Luis Alvarez
*/
template <class T>
T& image<T>::operator[](const int &i)
{
  #ifdef IMAGE_DEBUG
    if(i>=(int)(image_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
             (int)image_.size(),i);
      int j; scanf("%d",&j);
      exit(0);
    }
  #endif
  return  (image_.at(i));
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(int width,int height)
 * \brief Image constructor taking memory
 * \author Luis Alvarez
*/
template <class T>
image<T>::image(int width,int height)
{
  width_=width ;
  height_=height;
  nChannels_=1;
  image_.resize(width*height);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::image(const int width,const int height, const T &a, const T &b,
                       const T &c)
 * \brief Image constructor taking memory and initialiting the channels: R with
          a value, G with b value and B with c value
 * \author Carlos Falcon
*/
template <class T>
image<T>::image(const int width,const int height, const T &a, const T &b,
                const T &c)
{
  width_=width ;
  height_=height;
  nChannels_=3;
  int size=width_*height_;
  image_.resize(size*3);
  int i;
  #ifdef _OPENMP
  #pragma omp parallel for shared(size) private(i)
  #endif
  for(i=0;i<size;i++)
    image_.at(i)=a;
  #ifdef _OPENMP
  #pragma omp parallel for shared(size) private(i)
  #endif
  for(i=size;i<2*size;i++)
    image_.at(i)=b;
  #ifdef _OPENMP
  #pragma omp parallel for shared(size) private(i)
  #endif
  for(i=2*size;i<3*size;i++)
    image_.at(i)=c;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::image(const int width,const int height, const int nChannels,
                       const T &a)
 * \brief Image constructor taking memory and initialiting with a value
 * \author Luis Alvarez
*/
template <class T>
image<T>::image(const int width,const int height,const int nChannels,const T &a)
{
  width_=width ;
  height_=height;
  nChannels_=nChannels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
  int i;
  #ifdef _OPENMP
  #pragma omp parallel for shared(size) private(i)
  #endif
  for(i=0;i<size;i++) image_.at(i)=a;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(const int width,const int height,
                                          const int nChannels)
 * \brief Image constructor taking memory and initialiting with a value
 * \author Luis Alvarez
*/
template <class T>
image<T>::image(const int width,const int height, const int nChannels)
{
  width_=width;
  height_=height;
  nChannels_=nChannels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::~image()
 * \brief Basic constructor
 * \author Luis Alvarez
*/
template <class T>
image<T>::image() : width_(0), height_(0), nChannels_(0) {}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::image()
 * \brief Basic destructor
 * \author Luis Alvarez
*/
template <class T>
image<T>::~image(){}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn  template <class T> image<T>::image(std::string name)
 * \brief Constructor from an image file
 * \param [in] name Input file name.
 * \author Pedro Henriquez
*/
template <class T>
image<T>::image(std::string name /** INPUT FILE NAME */ )
{
  read(name);
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn void image<T>::clear()
* \brief Function to put the image to 0 and clear its vectors.
* \author Pedro Henríquez
*/
template <class T>
void image<T>::clear()
{
  width_=0;
  height_=0;
  nChannels_=0;
  // roi_clear();
  image_.clear();
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn void image<T>::set_nchannels(int nchannel1)
* \brief Change number of channels
* \param [in] nchannel1 New number of channels.
* \author Pedro Henríquez
*/
template <class T>
void image<T>::set_nchannels(int nchannel1)
{
  nChannels_=nchannel1;
  image_.resize(nChannels_*width_*height_);
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn void image<T>::set_size(int width,int height)
* \brief Change the image dimensions
* \param [in] width,height New image dimensions.
* \author Pedro Henríquez
*/
template <class T>
void image<T>::set_size(int width,int height)
{
  width_=width;
  height_=height;
  image_.resize(nChannels_*width_*height_);
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn int image<T>::read (std::string name)
* \brief Read an image selecting the library depend on the image format, returns
          0 when it can't load the image.
* \author Pedro Henríquez and Daniel Santana-Cedrés
*/
template <class T>
int image<T>::read (std::string name)
{
  unsigned char *red, *green, *blue;
  int width, height;
  // READ THE IMAGE WITH IO_PNG
  int output_value = ami_read_png(strdup(name.c_str()),&red,&green,&blue,
                                  &width,&height);
  if(output_value!=0)
  {
    return output_value;
  }
  // INITIALIZE ATTRIBUTES AND COUNTERS
  width_    = width;
  height_   = height;
  nChannels_= 3;
  int size_ = width_ * height_;
  image_.resize(size_*nChannels_);
  // FILL THE IMAGE VECTOR WITH THE INFORMATION FROM IO_PNG
  for(int i=0; i<size_; i++)
  {
    image_[i]         = (T)red[i];
    image_[i+size_]   = (T)green[i];
    image_[i+size_*2] = (T)blue[i];
  }
  // FREE MEMORY
  free(red);
  free(green);
  free(blue);
  return 0;
}

}

#endif
