/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


 #ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file lens_distortion.h
 * \brief Functions for lens distortion model basic operations
 * \author Luis Alvarez \n \n
*/
#ifndef LENS_DISTORTION_H
#define LENS_DISTORTION_H

#ifdef __AMIDEBUG__
  #include "wxAmiDebugLog.h"
#endif

/**
 * \fn void ami_lens_distortion_model_evaluation(double *a,int Na, double xc,
        double yc,double x_input,double y_input,double *x_output,
        double *y_output)
 *  \brief  Compute the lens distortion model in a point
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] a Input polynomial distortion model
 *  \param [in] Na Input degree of polynomial distortion model
 *  \param [in] xc x coordinate of the distortion center
 *  \param [in] yc y coordinate of the distortion center
 *  \param [in] x_input x coordinate of the input point
 *  \param [in] y_input y coordinate of the input point
 *  \param [out] x_output x coordinate of the output undistorted point
 *  \param [out] y_output y coordinate of the output undistorted point
 *  \return Returns the coordinates of the undistorted point
 * \author Luis Alvarez
 */
AMI_DLL_H void ami_lens_distortion_model_evaluation(double *a,int Na, double xc,
                                                    double yc,double x_input,
                                                    double y_input,
                                                    double *x_output,
                                                    double *y_output);

//------------------------------------------------------------------------------

/**
 * \fn int ami_RootCubicPolynomial(double *a,int N,double *x)
 *  \brief  Function to compute the real roots of a cubic polynomial.
 *  \pre Any parameter can be null.
 *  \pre N has to be 3.
 *  \post It returns the number of roots found sorted by magnitude
 *  \param [in]  a Polynomial coefficients a[0]+a[1]x+a[2]x^2 +...
 *  \param [in]  N Degree of polinomial (it has to be 3)
 *  \param [out] x Polynomial roots
 *  \return Returns the number of roots found sorted by magnitude
 * \author Luis Alvarez
 */
AMI_DLL_H int ami_RootCubicPolynomial(double *a,int N,double *x);

//------------------------------------------------------------------------------

/**
 * \fn double ami_polynomial_evaluation(double *a,int Na,double x)
 *  \brief  Evaluation of a polynom using horner algorithm.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \post It returns the number of roots found sorted by magnitud.
 *  \param[in] a Polinomial coeficients a[0]+a[1]x+a[2]x^2 +...
 *  \param[in] Na Polynom degree
 *  \param[in] x Point where the polynom is evaluated
 *  \return Returns the evaluation of x in polynomial a
 * \author Luis Alvarez
 */
AMI_DLL_H double ami_polynomial_evaluation(double *a,int Na,double x);

//------------------------------------------------------------------------------

/**
 * \fn int ami_inverse_lens_distortion(double x,double y,double x0,double y0,
        double *xt,double *yt,double *a,int Na)
 *  \brief  Function to inverse the lens distortion transformation.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] x x coordinate of the point to inverse
 *  \param [in] y y coordinate of the point to inverse
 *  \param [in] x0 x coordinate of the image center
 *  \param [in] y0 y coordinate of the image center
 *  \param [out] xt x coordinate of the inverse point transformed
 *  \param [out] yt y coordinate of the inverse point transformed
 *  \param [in] a Lens distortion model polynom
 *  \param [in] Na Degree of the lens distortion model polynom
 *  \return Returns 0 if the function finishes properly
 * \author Luis Alvarez
 */
AMI_DLL_H int ami_inverse_lens_distortion(double x,double y,double x0,double y0,
                                          double *xt,double *yt,double *a,
                                          int Na);

//------------------------------------------------------------------------------

/**
 * \fn int ami_inverse_lens_distortion(double x,double y,double x0,double y0,
        double *xt,double *yt,double *a,int Na, double dl1r)
 *  \brief  Function to inverse the lens distortion transformation.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] x x coordinate of the point to inverse
 *  \param [in] y y coordinate of the point to inverse
 *  \param [in] x0 x coordinate of the image center
 *  \param [in] y0 y coordinate of the image center
 *  \param [out] xt x coordinate of the inverse point transformed
 *  \param [out] yt y coordinate of the inverse point transformed
 *  \param [in] a Lens distortion model polynom
 *  \param [in] Na Degree of the lens distortion model polynom
 *  \param [in] dl1r Coefficient interpolated from vector of max distorsion 
 *              distances
 *  \return 0
 * \author Luis Alvarez and Pedro Henriquez
 */
int ami_inverse_lens_distortion_fast(double x,double y,double x0,double y0,
                                     double *xt,double *yt, double *a, int Na,
                                     double dl1r);


//------------------------------------------------------------------------------

/**
 * \fn int ami_inverse_lens_distortion(double x,double y,double x0,double y0,
        double *xt,double *yt,double *a,int Na, double dl1r)
 *  \brief  Function to inverse the lens distortion using Newton-Raphson algorithm.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] x x coordinate of the point to inverse
 *  \param [in] y y coordinate of the point to inverse
 *  \param [in] x0 x coordinate of the image center
 *  \param [in] y0 y coordinate of the image center
 *  \param [out] xt x coordinate of the inverse point transformed
 *  \param [out] yt y coordinate of the inverse point transformed
 *  \param [in] a Lens distortion model polynom
 *  \param [in] Na Degree of the lens distortion model polynom
 *  \return Returns the interpolated point or -1 otherwise
 * \author Luis Alvarez and Pedro Henriquez
 */
double ami_inverse_lens_distortion_newton_raphson(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na); /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
#endif
