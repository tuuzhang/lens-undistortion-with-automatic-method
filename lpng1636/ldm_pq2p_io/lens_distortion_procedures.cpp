/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#include "lens_distortion_procedures.h"
#include "ami_lens_distortion.h"
#include "ami_utilities.h"
#include <iostream>
#include <algorithm>
using namespace std;

//------------------------------------------------------------------------------
/** \class E_accumulators 
 * \brief Class to store the accumulators and methods for computing the gradient
 *        and the hessian used in the minimization process
 * \author Daniel Santana-Cedrés
 */
class E_accumulators{
  /**
   * Attributes for the different accumulators
   */
  double E_p1_p2_xc_yc, E_p1plush_p2_xc_yc, E_p1_p2plush_xc_yc,
         E_p1_p2_xcplush_yc, E_p1_p2_xc_ycplush, E_p1minush_p2_xc_yc,
         E_p1_p2minush_xc_yc, E_p1_p2_xcminush_yc, E_p1_p2_xc_ycminush,
         E_p1plush_p2plush_xc_yc, E_p1plush_p2_xcplush_yc,
         E_p1plush_p2_xc_ycplush, E_p1_p2plush_xcplush_yc,
         E_p1_p2plush_xc_ycplush, E_p1_p2_xcplush_ycplush, E1_p1_p2_xc_yc;
  public:
    /**
     * \brief Resetting the accumulators.
     * \return An E_accumulators object with the attributes setted to 0.
     */
    void reset();
    
    /**
     * \brief Method for computing the gradient of E.
     * \param[out] gradE Vector to store the result of computing the gradient.
     * \param[in] v Boolean vector which indicates the parameters to store in 
     *              gradE.
     * \param[in] h1 Value for computing the gradient of the distortion 
     *               parameters.
     * \param[in] h2 Value for computing the gradient of the distortion center.
     * \return The computed gradient inside gradE.
     */
    void compute_gradient(double *gradE, const vector<bool> v, const double h1,
                          const double h2);
    
    /**
     * \brief Procedure for computing the Hessian matrix
     * \param[out] hessianE Matrix to store the result of computing the hessian.
     * \param[in] gamma Value for convergence.
     * \param[in] v Boolean vector which indicates the parameters to store in 
     *              hessianE.
     * \param[in] h1 Value for computing the hessian of the distortion 
     *               parameters.
     * \param[in] h2 Value for computing the hessian of the distortion center.
     * \return The computed hessian inside hessianE
     */
    void compute_hessian(double **hessianE, const double gamma, 
                         const vector<bool> v,const double h1, const double h2);
    
    /**
     * \brief This method computes the accumulators according to the distance 
     *        between the points and their associated lines, by modifying the 
     *        lens distortion model.
     * \param[in,out] lines The set of detected primitives.
     * \param[in] ldm The estimated lens distortion model.
     * \param[in] dc0 Lens distortion center.
     * \param[in] h1 Value for the variation of the parameters of the lens 
     *               distortion model.
     * \param[in] h2 Value for the variation of the coordinates of the lens 
     *               distortion center.
     * \param[in] w Image width.
     * \param[in] h Image height.
     * \return The updated values of the accumulators.
     */
    void compute_accumulators(vector< line_points > &lines,
                              lens_distortion_model ldm, point2d<double> dc0, 
                              const double h1, const double h2, double prevp1, 
                              double prevp2, int w, int h);
    
    /**
     * \brief Get the value of the current energy.
     * \return The value of the current energy.
     */
    double get_E_p1_p2_xc_yc();
    
    /**
     * \brief Set the value of the new energy obtained by the optimization 
     *        process.
     * \param[in] val The new value for the energy.
     */
    void set_E1_p1_p2_xc_yc(double val);
    
    /**
     * \brief Get the value of the new energy.
     * \return The value of the new energy.
     */
    double get_E1_p1_p2_xc_yc();
};
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/**
 * \fn double distortion_points_to_line_equation(const lens_distortion_model &d,
                                                 line_points &l)
 * \brief Line estimation from point after applying a lens distortion model to 
          the points. return -1 if it does not work. Otherwise return the 
          squared distance of the points to the line after applying a lens 
          distortion model
 * \param[in] d Input distortion model
 * \param[in,out] l Input/output set of primitives
 * \author Luis Alvarez
 */
//Based on the book: O. Faugeras, Three-dimensional computer vision, MIT Press, 
//ISBN: 0-262-06158-9, 1993.
double distortion_points_to_line_equation(const lens_distortion_model &d,
                                          line_points &l)
{
  int i,j,k;
  long double suu,suv,svv,um,vm,h,r[4][3],min,step,norm;
  long double zero=10e-100;
  int N=l.get_points().size();
  double a,b,c;

  if(N<2){
    cout << "Number of points for computing the 2D line less than 2" << endl;
    return(-1.);
  }

  // WE CREATE AN AUXILIAR VECTOR OF POINTS
  vector< point2d<double> > points(l.get_points().size());
  for(i=0;i<((int)points.size());i++) 
    points[i]=d.evaluation(l.get_points()[i]);

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
  if(fabs(suv)<= zero){
    if(suu<svv && svv>zero){
      a=1.;
      b=0.;
      c=-um;
      l.set_abc(a,b,c);
      return(0.);
    }
    if(svv<suu && suu>zero){
      a=0.;
      b=1.;
      c=-vm;
      l.set_abc(a,b,c);
      return(0.);
    }
    cout << "It was not possible to compute the 2D straight line" << endl;
    return(-1);
  }

  r[2][1]=r[3][1]=r[0][0]=r[1][0]=1.;
  h=0.5*(suu-svv)/suv;
  long double root = sqrt(1.+h*h);
  r[0][1] = (h>0) ? -h-root : -h+root;
  r[2][0] = -r[0][1];
  r[0][2] = -(um+r[0][1]*vm);
  r[1][1] = -1./r[0][1];
  r[1][2] = -(um+r[1][1]*vm);
  r[2][2] = -(r[2][0]*um+vm);
  r[3][0] = -1./r[2][0];
  r[3][2] = -(r[3][0]*um+vm);
  
  for(j=0;j<4;j++){
    norm=sqrt(r[j][0]*r[j][0]+r[j][1]*r[j][1]);
    for(i=0;i<3;i++)
     r[j][i]/=norm;
  }

  min=0.; k=0;
  for(i=0;i<N;i++){
    step=r[0][0]*points[i].x+r[0][1]*points[i].y+r[0][2];
    min+=step*step;
  }
  for(j=1;j<4;j++){
    h=0;
    for(i=0;i<N;i++){
      step=r[j][0]*points[i].x+r[j][1]*points[i].y+r[j][2];
      h+=step*step;
    }
    if(h<min){
      k=j;
      min=h;
    }
  }

  l.set_abc((double) r[k][0],(double)r[k][1],(double)r[k][2]);

  return min;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Updating the lens distortion model and incrementing the value of the 
//accumulator
void update_ldm_and_increment_accumulator(lens_distortion_model &ldm, 
                                          const double n1, const double dk1, 
                                          const double n2, const double dk2,
                                          double &accum,line_points &lp)
{
  ldm.get_d()[1] = n1/dk1;
  ldm.get_d()[2] = n2/dk2;
  accum         += distortion_points_to_line_equation(ldm,lp);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void update_den_ks(point2d<double> dc, int w, int h, double &denk1, double &denk2)
{
  double dmi = update_rsqmax(dc,w,h);
  double r2 = (sqrt(dmi))/2; double r2_2 = r2*r2; double r2_4 = r2*r2*r2*r2;
  denk1 = -12*r2_2;
  denk2 = -12*r2_4;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void E_accumulators::reset()
{
  E_p1_p2_xc_yc            = 
  E_p1plush_p2_xc_yc       = 
  E_p1_p2plush_xc_yc       = 
  E_p1_p2_xcplush_yc       = 
  E_p1_p2_xc_ycplush       = 
  E_p1minush_p2_xc_yc      = 
  E_p1_p2minush_xc_yc      = 
  E_p1_p2_xcminush_yc      = 
  E_p1_p2_xc_ycminush      = 
  E_p1plush_p2plush_xc_yc  = 
  E_p1plush_p2_xcplush_yc  = 
  E_p1plush_p2_xc_ycplush  = 
  E_p1_p2plush_xcplush_yc  = 
  E_p1_p2plush_xc_ycplush  = 
  E_p1_p2_xcplush_ycplush  = 
  E1_p1_p2_xc_yc           = 0.;
}
//------------------------------------------------------------------------------
    
//------------------------------------------------------------------------------
//Method for computing the gradient of E
void E_accumulators::compute_gradient(double *gradE, const vector<bool> v,
                                      const double h1, const double h2)
{
  //Polynomial
  //Gradient vector
  /*-(E(p1+h1,p2,xc,yc) - E(p1,p2,xc,yc))/h1
    -(E(p1,p2+h1,xc,yc) - E(p1,p2,xc,yc))/h1
    -(E(p1,p2,xc+h2,yc) - E(p1,p2,xc,yc))/h2
    -(E(p1,p2,xc,yc+h2) - E(p1,p2,xc,yc))/h2*/
  gradE[0] = (v[0]) ? -((E_p1plush_p2_xc_yc - E_p1_p2_xc_yc)/h1) : 0.;
  gradE[1] = (v[1]) ? -((E_p1_p2plush_xc_yc - E_p1_p2_xc_yc)/h1) : 0.;
  gradE[2] = (v[2]) ? -((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2) : 0.;
  gradE[3] = (v[3]) ? -((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2) : 0.;
}
//------------------------------------------------------------------------------
    
//------------------------------------------------------------------------------
//Method for computing the hessian of E
void E_accumulators::compute_hessian(double **hessianE, double gamma, 
                                     vector<bool> v, const double h1, 
                                     const double h2)
{
  //Polynomial
  //Hessian matrix
  hessianE[0][0] = ((E_p1plush_p2_xc_yc + E_p1minush_p2_xc_yc -
                      2*E_p1_p2_xc_yc)/(h1*h1))+gamma;

  hessianE[0][1] = hessianE[1][0] = (1/h1)*
                  (((E_p1plush_p2plush_xc_yc - E_p1plush_p2_xc_yc)/h1)-
                    ((E_p1_p2plush_xc_yc - E_p1_p2_xc_yc)/h1));

  hessianE[0][2] = hessianE[2][0] = (1/h1)*
                  (((E_p1plush_p2_xcplush_yc - E_p1plush_p2_xc_yc)/h2)-
                    ((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2));

  hessianE[0][3] = hessianE[3][0] = (1/h1)*
                  (((E_p1plush_p2_xc_ycplush - E_p1plush_p2_xc_yc)/h2)-
                    ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));

  hessianE[1][1] = ((E_p1_p2plush_xc_yc + E_p1_p2minush_xc_yc -
                      2*E_p1_p2_xc_yc)/(h1*h1))+gamma;

  hessianE[1][2] = hessianE[2][1] = (1/h1)*
                  (((E_p1_p2plush_xcplush_yc - E_p1_p2plush_xc_yc)/h2)-
                    ((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2));

  hessianE[1][3] = hessianE[3][1] = (1/h1)*
                  (((E_p1_p2plush_xc_ycplush - E_p1_p2plush_xc_yc)/h2)-
                    ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));

  hessianE[2][2] = ((E_p1_p2_xcplush_yc + E_p1_p2_xcminush_yc -
                      2*E_p1_p2_xc_yc)/(h2*h2))+gamma;

  hessianE[2][3] = hessianE[3][2] = (1/h2)*
                  (((E_p1_p2_xcplush_ycplush - E_p1_p2_xcplush_yc)/h2)-
                    ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));

  hessianE[3][3] = ((E_p1_p2_xc_ycplush + E_p1_p2_xc_ycminush -
                      2*E_p1_p2_xc_yc)/(h2*h2))+gamma;
  if(!v[0])
  {
    hessianE[0][0]=1.+gamma;
    hessianE[0][1]=hessianE[0][2]=hessianE[0][3]=
    hessianE[1][0]=hessianE[2][0]=hessianE[3][0]=0.;
  }
  if(!v[1])
  {
    hessianE[1][1]=1.+gamma;
    hessianE[1][0]=hessianE[1][2]=hessianE[1][3]=
    hessianE[0][1]=hessianE[2][1]=hessianE[3][1]=0.;
  }
  if(!v[2])
  {
    hessianE[2][2]=1.+gamma;
    hessianE[2][0]=hessianE[2][1]=hessianE[2][3]=
    hessianE[0][2]=hessianE[1][2]=hessianE[3][2]=0.;
  }
  if(!v[3])
  {
    hessianE[3][3]=1.+gamma;
    hessianE[3][0]=hessianE[3][1]=hessianE[3][2]=
    hessianE[0][3]=hessianE[1][3]=hessianE[2][3]=0.;
  }
}
//------------------------------------------------------------------------------    
    
//------------------------------------------------------------------------------
void E_accumulators::compute_accumulators(vector< line_points > &lines,
                          lens_distortion_model ldm, point2d<double> dc0,
                          const double h1, const double h2, double prevp1,
                          double prevp2, int w, int h)
{
  point2d<double> dc;
  double denk1, denk2;
  //POLYNOMIAL
  if(ldm.get_type()==POLYNOMIAL)
  {
    for(int i=0; i<(int)lines.size(); i++)
    {
      dc = dc0;
      //E(p1,p2,xc,yc)
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*prevp2), denk1, 
        (4*prevp2-prevp1), denk2, E_p1_p2_xc_yc, lines[i]); 

      //E(p1+h1,p2,xc,yc)
      update_ldm_and_increment_accumulator(ldm, ((prevp1+h1)-16*prevp2),
        denk1, (4*prevp2-(prevp1+h1)), denk2, E_p1plush_p2_xc_yc, lines[i]);

      //E(p1,p2+h1,xc,yc)
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*(prevp2+h1)), 
        denk1, (4*(prevp2+h1)-prevp1), denk2, E_p1_p2plush_xc_yc, lines[i]);

      //E(p1,p2,xc+h2,yc)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*prevp2), denk1, 
        (4*prevp2-prevp1), denk2, E_p1_p2_xcplush_yc, lines[i]);

      //E(p1,p2,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*prevp2), denk1, 
        (4*prevp2-prevp1), denk2, E_p1_p2_xc_ycplush, lines[i]);

      //E(p1-h1,p2,xc,yc)
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, ((prevp1-h1)-16*prevp2), 
        denk1, (4*prevp2-(prevp1-h1)), denk2, E_p1minush_p2_xc_yc,lines[i]);

      //E(p1,p2-h1,xc,yc)
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*(prevp2-h1)), 
        denk1, (4*(prevp2-h1)-prevp1), denk2, E_p1_p2minush_xc_yc,lines[i]);

      //E(p1,p2,xc-h2,yc)
      dc.x = dc0.x-h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*prevp2), denk1, 
        (4*prevp2-prevp1), denk2, E_p1_p2_xcminush_yc, lines[i]);

      //E(p1,p2,xc,yc-h2)
      dc.x = dc0.x;
      dc.y = dc0.y-h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*prevp2), denk1, 
        (4*prevp2-prevp1), denk2, E_p1_p2_xc_ycminush, lines[i]);

      //E(p1+h1,p2+h1,xc,yc)
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm,((prevp1+h1)-16*(prevp2+h1)),
        denk1, (4*(prevp2+h1)-(prevp1+h1)), denk2, E_p1plush_p2plush_xc_yc, 
        lines[i]);

      //E(p1+h1,p2,xc+h2,yc)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, ((prevp1+h1)-16*prevp2), 
        denk1, (4*prevp2-(prevp1+h1)), denk2, E_p1plush_p2_xcplush_yc, 
        lines[i]);

      //E(p1+h1,p2,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, ((prevp1+h1)-16*prevp2), 
        denk1, (4*prevp2-(prevp1+h1)), denk2, E_p1plush_p2_xc_ycplush, 
        lines[i]);

      //E(p1,p2+h1,xc+h2,yc)
      dc.x = dc0.x+h2;
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*(prevp2+h1)), 
        denk1, (4*(prevp2+h1)-prevp1), denk2, E_p1_p2plush_xcplush_yc, 
        lines[i]);

      //E(p1,p2+h1,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*(prevp2+h1)), 
        denk1, (4*(prevp2+h1)-prevp1), denk2, E_p1_p2plush_xc_ycplush, 
        lines[i]);

      //E(p1,p2,xc+h2,yc+h2)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, (prevp1-16*prevp2), denk1, 
        (4*prevp2-prevp1), denk2, E_p1_p2_xcplush_ycplush, lines[i]);
    }
  }
  else //DIVISION
  {
    for(int i=0; i<(int)lines.size(); i++)
    {
      dc = dc0;
      //E(p1,p2,xc,yc)
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2))), denk1, 
        (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1))), denk2,E_p1_p2_xc_yc,
        lines[i]);

      //E(p1+h1,p2,xc,yc)
      update_ldm_and_increment_accumulator(ldm, 
        (((-(prevp1+h1))/(1+prevp1+h1))+((16*prevp2)/(1+prevp2))), denk1, 
        (((-4*prevp2)/(1+prevp2))+((prevp1+h1)/(1+prevp1+h1))), denk2, 
        E_p1plush_p2_xc_yc, lines[i]);

      //E(p1,p2+h1,xc,yc)
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*(prevp2+h1))/(1+prevp2+h1))), denk1, 
        (((-4*(prevp2+h1))/(1+prevp2+h1))+(prevp1/(1+prevp1))), denk2,
        E_p1_p2plush_xc_yc, lines[i]);

      //E(p1,p2,xc+h2,yc)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2))), denk1, 
        (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1))), denk2, 
        E_p1_p2_xcplush_yc, lines[i]);

      //E(p1,p2,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);;
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2))), denk1,
        (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1))), denk2, 
        E_p1_p2_xc_ycplush, lines[i]);

      //E(p1-h1,p2,xc,yc)
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-(prevp1-h1))/(1+prevp1-h1))+((16*prevp2)/(1+prevp2))), denk1, 
        (((-4*prevp2)/(1+prevp2))+((prevp1-h1)/(1+prevp1-h1))), denk2, 
        E_p1minush_p2_xc_yc, lines[i]);

      //E(p1,p2-h1,xc,yc)
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*(prevp2-h1))/(1+prevp2-h1))), denk1,
        (((-4*(prevp2-h1))/(1+prevp2-h1))+(prevp1/(1+prevp1))), denk2, 
        E_p1_p2minush_xc_yc, lines[i]);

      //E(p1,p2,xc-h2,yc)
      dc.x = dc0.x-h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2))), denk1,
        (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1))), denk2, 
        E_p1_p2_xcminush_yc, lines[i]);

      //E(p1,p2,xc,yc-h2)
      dc.x = dc0.x;
      dc.y = dc0.y-h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2))), denk1, 
        (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1))), denk2, 
        E_p1_p2_xc_ycminush, lines[i]);

      //E(p1+h1,p2+h1,xc,yc)
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-(prevp1+h1))/(1+prevp1+h1))+((16*(prevp2+h1))/(1+prevp2+h1))),
        denk1, 
        (((-4*(prevp2+h1))/(1+prevp2+h1))+((prevp1+h1)/(1+prevp1+h1))), 
        denk2, E_p1plush_p2plush_xc_yc, lines[i]);

      //E(p1+h1,p2,xc+h2,yc)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-(prevp1+h1))/(1+prevp1+h1))+((16*prevp2)/(1+prevp2))), denk1,
        (((-4*prevp2)/(1+prevp2))+((prevp1+h1)/(1+prevp1+h1))), denk2, 
        E_p1plush_p2_xcplush_yc, lines[i]);

      //E(p1+h1,p2,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-(prevp1+h1))/(1+prevp1+h1))+((16*prevp2)/(1+prevp2))), denk1,
        (((-4*prevp2)/(1+prevp2))+((prevp1+h1)/(1+prevp1+h1))), denk2,
        E_p1plush_p2_xc_ycplush, lines[i]);

      //E(p1,p2+h1,xc+h2,yc)
      dc.x = dc0.x+h2;
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*(prevp2+h1))/(1+prevp2+h1))), denk1,
        (((-4*(prevp2+h1))/(1+prevp2+h1))+(prevp1/(1+prevp1))), denk2,
        E_p1_p2plush_xcplush_yc, lines[i]);

      //E(p1,p2+h1,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*(prevp2+h1))/(1+prevp2+h1))), denk1,
        (((-4*(prevp2+h1))/(1+prevp2+h1))+(prevp1/(1+prevp1))), denk2,
        E_p1_p2plush_xc_ycplush, lines[i]);

      //E(p1,p2,xc+h2,yc+h2)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      update_den_ks(dc,w,h,denk1,denk2);
      update_ldm_and_increment_accumulator(ldm, 
        (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2))), denk1,
        (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1))), denk2,
        E_p1_p2_xcplush_ycplush, lines[i]);
    }
  }
}
//------------------------------------------------------------------------------
double E_accumulators::get_E_p1_p2_xc_yc()
{
  return E_p1_p2_xc_yc;
}
//------------------------------------------------------------------------------
void E_accumulators::set_E1_p1_p2_xc_yc(double val)
{
  E1_p1_p2_xc_yc = val;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
double E_accumulators::get_E1_p1_p2_xc_yc()
{
  return E1_p1_p2_xc_yc;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/**
 * \fn double model_center_estimation_2p(vector< line_points > &lines,
                                  lens_distortion_model &d, int w, int h,
                                  vector<bool> v)
 * \brief Function to optimize the lens distortion model parameters and its 
 *        center
 * \param [in] lines Primitives detected in the image
 * \param [in,out] d Lens distortion model
 * \param [in] w Image width
 * \param [in] h Image height
 * \param [in] v Bool vector which indicates the parameters to estimate
 * \return Returns the distance between all the points and their associated 
 *         lines
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
double model_center_estimation_2p(vector< line_points > &lines,
                                  lens_distortion_model &d, int w, int h,
                                  vector<bool> v)
{
  lens_distortion_model ldm;
  point2d<double> distortion_center = d.get_distortion_center();
  ldm.set_type(d.get_type());
  const double h1 = 1e-4, h2 = 1e-2;
  E_accumulators Ec;
  //Ec.reset();

  double gamma = 10.0;

  double prevk1 = d.get_d()[1], prevk2 = d.get_d()[2];
  double nextk1 = 1.0,          nextk2 = 1.0;
  double dmi = update_rsqmax(distortion_center,w,h);
  //We compute p values using the relation r1=2r2
  double r1 = sqrt(dmi), r1_2 = r1*r1, r1_4 = r1*r1*r1*r1;
  double r2 = r1/2,      r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
  double prevp1 = 0.0, prevp2 = 0.0, nextp1 = 0.0, nextp2 = 0.0;
  if(d.get_type()==POLYNOMIAL)
  {
    //p1 = k14r_2^2 + k216r_2^4
    prevp1 = prevk1*4*r2_2 + prevk2*16*r2_4;
    //p2 = k1r_2^2 + k2r_2^4
    prevp2 = prevk1*r2_2 + prevk2*r2_4;
  }
  else
  {
    //p1 = k14r_2^2 + k216r_2^4
    prevp1 = (1/(1+prevk1*r1_2+prevk2*r1_4))-1;
    //p2 = k1r_2^2 + k2r_2^4
    prevp2 = (1/(1+prevk1*r2_2+prevk2*r2_4))-1;
  }
  
  //double nextp1 = 0.0, nextp2 = 0.0;
  double d0[4] = {prevp1, prevp2, distortion_center.x, distortion_center.y};
  double d1[4];
  point2d<double> prevdc = distortion_center;
  point2d<double> nextdc = distortion_center;
  nextdc.x += h2;
  nextdc.y += h2;

  int convergence_it = 0;
  ldm.get_d().resize(3);
  ldm.get_d()[0] = 1.;
  double *gradE;
  ami2_malloc1d(gradE,double,4);
  double **hessianE;
  ami2_malloc2d(hessianE,double,4,4);

  const double TOL = 1e-4;
  
  while(( ((fabs(prevk1-nextk1))     > ((fabs(nextk1)+1e-30)*TOL))
        ||((fabs(prevk2-nextk2))     > ((fabs(nextk2)+1e-30)*TOL))
        ||((fabs(prevdc.x-nextdc.x)) > ((fabs(nextdc.x)+1e-30)*TOL))
        ||((fabs(prevdc.y-nextdc.y)) > ((fabs(nextdc.y)+1e-30)*TOL)))
        && convergence_it<=100)
  {
    //----1ST WE COMPUTE THE Es
    //INITIALIZATION OF THE ACCUMULATORS
    Ec.reset();
    prevdc = point2d<double>(d0[2],d0[3]);
    
    //COMPUTE THE ACCUMULATORS FOR THE COLLECTION OF LINES
    Ec.compute_accumulators(lines,ldm,prevdc,h1,h2,prevp1,prevp2,w,h);
    
    //COMPUTE THE GRADIENT AND HESSIAN
    Ec.compute_gradient(gradE, v, h1, h2);
    Ec.compute_hessian(hessianE, gamma, v, h1, h2);

    //WE SOLVE THE SYSTEM
    if(ami2_gauss(hessianE,gradE,4)!=0)
    {
      cout << "The system has not solution" << endl;
      exit(EXIT_FAILURE);
    }

    //WE COMPUTE THE ERROR WITH THE RESULT OF THE SYSTEM
    //d1 = d0 + z
    d1[0] = d0[0] + gradE[0];
    d1[1] = d0[1] + gradE[1];
    d1[2] = d0[2] + gradE[2];
    d1[3] = d0[3] + gradE[3];
    nextdc.x = d1[2];
    nextdc.y = d1[3];
    ldm.set_distortion_center(nextdc);
    
    double denk1, denk2;
    update_den_ks(prevdc,w,h,denk1,denk2);
    if(d.get_type()==POLYNOMIAL)
    {
      prevk1 = (prevp1-16*prevp2)/(denk1);
      prevk2 = (4*prevp2-prevp1)/(denk2);
    }
    else
    {
      prevk1 = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(denk1);
      prevk2 = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(denk2);
    }
    
    update_den_ks(nextdc,w,h,denk1,denk2);
    nextp1 = d1[0];
    nextp2 = d1[1];
    if(d.get_type()==POLYNOMIAL)
    {
      nextk1 = ldm.get_d()[1] = (nextp1-16*nextp2)/(denk1);
      nextk2 = ldm.get_d()[2] = (4*nextp2-nextp1)/(denk2);
    }
    else
    {
      nextk1 = ldm.get_d()[1] = (((-nextp1)/(1+nextp1))+
                                 ((16*nextp2)/(1+nextp2)))/(denk1);
      nextk2 = ldm.get_d()[2] = (((-4*nextp2)/(1+nextp2))+
                                  (nextp1/(1+nextp1)))/(denk2);
    }
    
    Ec.set_E1_p1_p2_xc_yc(0.0);
    for(int i=0; i<(int)lines.size(); i++)
      Ec.set_E1_p1_p2_xc_yc(Ec.get_E1_p1_p2_xc_yc() + 
                            distortion_points_to_line_equation(ldm,lines[i]));

    //TEMPORARY MODEL FOR CHECKING THE INVERTIBILITY
    update_den_ks(point2d<double>(d1[2],d1[3]),w,h,denk1,denk2);
    lens_distortion_model tmp_ldm;
    tmp_ldm.get_d().resize(3);
    tmp_ldm.get_d()[0] = 1.;
    if(d.get_type()==POLYNOMIAL)
    {
      tmp_ldm.get_d()[1] = (d1[0]-16*d1[1])/(denk1);
      tmp_ldm.get_d()[2] = (4*d1[1]-d1[0])/(denk2);
    }
    else
    {
      tmp_ldm.get_d()[1] = (((-d1[0])/(1+d1[0]))+
                            ((16*d1[1])/(1+d1[1])))/(denk1);
      tmp_ldm.get_d()[2] = (((-4*d1[1])/(1+d1[1]))+(d1[0]/(1+d1[0])))/(denk2);
    }
    tmp_ldm.set_distortion_center(point2d<double>(d1[2],d1[3]));
    tmp_ldm.set_type(d.get_type());
    
    if(Ec.get_E1_p1_p2_xc_yc() < Ec.get_E_p1_p2_xc_yc() && 
       check_invertibility(tmp_ldm,w,h)==true)
    {
      prevp1 = nextp1;
      prevp2 = nextp2;
      d0[0]  = d1[0];
      d0[1]  = d1[1];
      d0[2]  = d1[2];
      d0[3]  = d1[3];
      gamma /= 10;
    }
    else
    {
      gamma *= 10;
    }
    convergence_it++;
  }
  free(gradE);
  ami2_free2d(hessianE);

  //WE UPDATE THE FINAL MODEL
  double denk1, denk2;
  update_den_ks(point2d<double>(d0[2],d0[3]),w,h,denk1,denk2);
  if(d.get_type()==POLYNOMIAL)
  {
    d.get_d()[1] = (d0[0]-16*d0[1])/(denk1);
    d.get_d()[2] = (4*d0[1]-d0[0])/(denk2);
  }
  else
  {
    d.get_d()[1] = (((-d0[0])/(1+d0[0]))+((16*d0[1])/(1+d0[1])))/(denk1);
    d.get_d()[2] = (((-4*d0[1])/(1+d0[1]))+(d0[0]/(1+d0[0])))/(denk2);
  }
  
  d.set_distortion_center(point2d<double>(d0[2],d0[3]));

  //AND RETURN THE ERROR VALUE
  return std::min(Ec.get_E1_p1_p2_xc_yc(), Ec.get_E_p1_p2_xc_yc());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/**
 * \fn int build_l1r_vector(vector<double> &l1r,point2d<double> &dc, 
                            double max_distance_corner,int Na, double *a)
 * \brief Build an intermediate vector with values of L-1(r) for d = (0 to 
 *        max_distance_corner)
 * \param [in] [out] l1r vector to store the roots
 * \param [in] dc distortion center
 * \param [in] max_distance_corner Maximum distance from distortion center to a 
 *             corner
 * \param [in] a Lens distortion model polynom
 * \param [in] Na Degree of the lens distortion model polynom
 * \author Luis Alvarez, Pedro Henriquez
 */
int build_l1r_vector(std::vector<double> &l1r, 
                     double max_distance_corner,int Na, double *a)
{
  //BUILD INTERMEDIATE VECTOR WITH L-1(r) FROM CENTER TO FURTHEST CORNER
  if(a[Na]==0.) return(-1);
  l1r.resize((int)(max_distance_corner+1.5));

  // AUXILIARY VARIABLES
  double *b,*b2,root2,root=1.;

  // WE ALLOCATE MEMORY
  b=(double*)malloc( sizeof(double)*(Na+2) );
  b2=(double*)malloc( sizeof(double)*(Na+1) );

  // WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS AND THE DERIVATIVE */
  for(int i=1; i<(Na+2); i++){ b[i]=a[i-1];}
  for(int i=0; i<(Na+1); i++){ b2[i]=a[i]*(i+1);}

  // WE BUILD THE VECTOR OF INVERSE LENS DISTORTION FUNCTION
  for(int dist=1; dist<(int)l1r.size(); dist++)
  {
    // WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS AND THE DERIVATIVE */
    b[0]=-dist;

    //NEWTON-RAPHSON TO COMPUTE THE POLYNOMIAL ROOT
    for(int k=0;k<10000;k++){
      double pol_eval=ami_polynomial_evaluation(b,Na+1,root);
      double pol_der=ami_polynomial_evaluation(b2,Na,root);
      if(pol_der==0.) break;
      root2=root-pol_eval/pol_der;
      if(fabs(root-root2)<(fabs(root)*1e-12)){
        root=root2;
        break;
      }
      root=root2;
    }

    //PUSH RESULT
    l1r[dist]=root/dist;
  }
  free(b);
  free(b2);
  l1r[0]=l1r[1];

  return 0;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/**
 * \fn int build_l1r_quotient_vector(vector<double> &l1r,point2d<double> &dc,
                                     double max_distance_corner,int Na,
                                     double *a)
 * \brief Build an intermediate vector with values of L-1(r) for
          d = (0 to max_distance_corner)
 * \param [in] [out] l1r vector to store the roots
 * \param [in] dc distortion center
 * \param [in] max_distance_corner Maximum distance from distortion center to a
                corner
 * \param [in] a Lens distortion model polynom
 * \param [in] Na Degree of the lens distortion model polynom
 * \author Luis Alvarez, Pedro Henriquez
 */
int build_l1r_quotient_vector(std::vector<double> &l1r,
                              double max_distance_corner, int Na, double *a)
{
  //BUILD INTERMEDIATE VECTOR WITH L-1(r) FROM CENTER TO FURTHEST CORNER
  if(a[Na]==0. || Na<2) return(-1);
  l1r.resize((int)(max_distance_corner+1.5));

  // WE BUILD THE VECTOR OF INVERSE LENS DISTORTION FUNCTION
  double root,disc,root1,root2;
  int dist=(int)(max_distance_corner+0.5);
  disc= 1.-4*dist*dist*a[2];
  if(disc<0) return(-2);
  disc=sqrt(disc);
  root1=(1-disc)/(2*dist*a[2]);
  root2=(1+disc)/(2*dist*a[2]);
  root = (root1>0) ? root1 : root2;
  l1r[dist] = root/dist;
  
  while( (dist--)>0)
  {
    disc= 1.-4*dist*dist*a[2];
    if(disc<0) return(-2);
    disc=sqrt(disc);
    root1=(1-disc)/(2*dist*a[2]);
    root2=(1+disc)/(2*dist*a[2]);
    root = (fabs(root-root1)<fabs(root-root2)) ? root1: root2;
    l1r[dist] = root/dist;
  }

  l1r[0]=l1r[1];

  return 0;
}
//------------------------------------------------------------------------------