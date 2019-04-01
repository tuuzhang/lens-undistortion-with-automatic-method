/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

#include "line_extraction.h"
#include "lens_distortion_procedures.h"
#include "ami_utilities.h"
#include <algorithm>
//#include <cstdio>
#include<assert.h>
static const double ami_pi=acos((double) -1.);
static const int MIN_POINTS_IN_LINE=20;
static const int MAX_DISTANCE_BETWEEN_LINES=10;

//------------------------------------------------------------------------------
//This function corrects the orientation of the edge.
//Parameters: the edge position, the orientation (sin,cos) and the lens
//            distortion model
 /**
 * \fn point2d<float> orientation_update(point2d<double> p, float sine, 
 *                                      float cosine, lens_distortion_model ldm)
 * \brief Correction of the edge orientation (sin,cos) using the provided lens 
          distortion model
 * \param[in] p The edge point
 * \param[in] sine Value of the sine component
 * \param[in] cosine Value of the cosine component
 * \param[in] ldm The lens distortion model
 * \return A vector with the corrected orientation
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
point2d<float> orientation_update(point2d<double> p, float sine, float cosine,
                                 lens_distortion_model ldm)
{
  point2d<float> corrected_orientation;
  double a,b,norm;
  //Point plus the orientation
  point2d<double> p_ori(p.x + cosine, p.y + sine);
  //We apply the model to the original point
  point2d<double> p_prime(ldm.evaluation(p));
  //We apply the model to the point plus the orientation
  point2d<double> p_ori_prime(ldm.evaluation(p_ori));

  //We compute the new orientation and the norm
  a = p_ori_prime.x - p_prime.x;
  b = p_ori_prime.y - p_prime.y;
  norm = sqrt(a*a + b*b);
  corrected_orientation.x = ((a*a + b*b) <= 0) ? cosine : a/norm;
  corrected_orientation.y = ((a*a + b*b) <= 0) ? sine   : b/norm;
  return corrected_orientation;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Updating the image contours input object according to
// the initial lens distortion model
/**
 * \fn void update_contours(lens_distortion_model ini_ldm,
                     const image_contours &contours,
                     image_contours &contours_modified,
                     const vector<int> &index
                    )
 * \brief Updating the contours object with the provided
          lens distortion model. If the model is empty, just makes a copy.
 * \param[in] ini_ldm The lens distortion model
 * \param[in] contours The input contours object
 * \param[out] contours_modified The output contours object
 * \param[in] index Vector with the indexes of the edge points
 * \return The updated contours object
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void update_contours(lens_distortion_model ini_ldm,
                     const image_contours &contours,
                     image_contours &contours_modified,
                     const vector<int> &index
                    )
{
  if(ini_ldm.is_identity())
  {
    // The modified model has the default value. we only copy
    contours_modified = contours;
  }
  else
  {
    // Else, we update the contours object using the provided model
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<(int)index.size(); i++)
    {
      int j = index[i];
      // We update the contour component of contours object
      contours_modified.get_c()[j] = true;
      // CORRECT THE POSITION
      point2d<double> ori(contours.get_x()[j], 
                          contours.get_y()[j]);
      point2d<double> corrected(ini_ldm.evaluation(ori));
      contours_modified.get_x()[j] = corrected.x;
      contours_modified.get_y()[j] = corrected.y;

      // Correct the orientation
      point2d<float> corrected_orientation = orientation_update(ori,
                                      contours.get_sine()[j],
                                      contours.get_cosine()[j],
                                      ini_ldm);
      
      // We update the value of the sin and cos
      contours_modified.get_cosine()[j] = corrected_orientation.x;
      contours_modified.get_sine()[j]   = corrected_orientation.y;
      if(contours_modified.get_sine()[j]<0){
        contours_modified.get_cosine()[j] *=-1;
        contours_modified.get_sine()[j]   *=-1;
      }
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Filling the voting matrix corresponding to the k slice
/**
 * \fn void voting(const image_contours &contours_modified, 
            lens_distortion_model it_ldm,
            double max_norm,
            double p,
            double step_angle,
            int height_score,
            int width_score,
            float angle_point_orientation_max_difference,
            float angle_resolution,
            float distance_resolution,
            double *sine, double *cosine,
            vector<vector<float> > &score_k,
            const vector<int> &index
           )
 * \brief Filling the voting matrix corresponding to the k slice
 * \param[in] contours_modified A contours object modified with the 
 *                                initial lens distortion model 
 *                                (from update_contours)
 * \param[in] it_ldm The lens distortion model corresponding to the current 
 *                   slice
 * \param[in] step_angle Step for the angle dimension in the Hough voting matrix
 * \param[in] (height_score,width_score) Size of the voting matrix
 * \param[in] angle_increment Increment for the angle interval
 * \param[in] distance_resolution Resolution of the distance interval (width of 
 *                                the voting matrix)
 * \param[in] (sine,cosine) Arrays of samples of sine and cosine functions
 * \param[out] score_k The voting matrix
 * \param[in] index Vector with the indexes of the edge points
 * \return The filled voting matrix
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void voting(const image_contours &contours_modified, 
            lens_distortion_model it_ldm,
            double step_angle,
            int height_score,
            int width_score,
            float angle_increment,
            float distance_resolution,
            double *sine, double *cosine,
            vector<vector<float> > &score_k,
            const vector<int> &index
           )
{
  // We do the voting process only with the points inside the vector index
  for(int ind=0; ind<(int)index.size(); ind++)
  {
    double distance,x2,y2,x3,y3;
    int q,d,angle;
    q = index[ind];
    // We compute the distortion model according to the type
    x3 = x2 = contours_modified.get_x()[q];
    y3 = y2 = contours_modified.get_y()[q];
    if(it_ldm.get_d()[1]!=0)
    {
      point2d<double> res = it_ldm.evaluation(point2d<double>(x2,y2));
      x3 = res.x;
      y3 = res.y;
    }
    // Orientation correction
    point2d<float> corrected_orientation = orientation_update(
                                              point2d<double>(x2,y2),
                                              contours_modified.get_sine()[q],
                                              contours_modified.get_cosine()[q],
                                              it_ldm);
    // We estimate the angle interval
    if(corrected_orientation.y>=0)
    {
      angle = (int) ((ami_pi-atan2(corrected_orientation.y,
                      corrected_orientation.x))/step_angle);
    }
    else
    {
      angle = (int) ((ami_pi-atan2(-corrected_orientation.y,
                      -corrected_orientation.x))/step_angle);
    }
    if(angle==height_score) angle=height_score-1;

    int l_min=(int) (angle-angle_increment);
    int l_max=(int) (angle+angle_increment);
    int id = 2;

    for(int l=l_min; l<l_max; l++)
    {
      int angle_index = l, sine_sign = 1; 
      if(l<0) angle_index = height_score+l;
      if(l>=height_score)
      {
        angle_index = l-height_score;
        sine_sign = -1;
      }
      
      distance=-cosine[angle_index]*y3-sine[angle_index]*x3;
      for(int nd=-id;nd<=id;nd++)
      {
        d= width_score/2+((int)(distance/distance_resolution+0.5))+nd;
        if(d>=0 && d<width_score)
        {
          double rect_distance = fabs((distance+cosine[angle_index]*y3+
                            sine_sign*sine[angle_index]*x3)+
                            nd*distance_resolution);
          score_k[d][angle_index]+=1./
                            (1.+rect_distance);
        }
      }
    }
  } // End index
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Updating the value of the distance according to the higher average sum
void update_distance(const vector<vector<float> > &score_k, int &d0, int ang, int Ri, 
                    int width_score)
{
  float sum=0., sum_l=0., sum_r=0.;
  float mean=0., mean_l=0., mean_r=0.;
  int n_ele=2*Ri+1;
  // Sums
  // Left (centered on d0-1)
  if(d0>1)
  {
    sum_l = score_k[d0-2][ang] + score_k[d0-1][ang] + score_k[d0][ang];
    if (n_ele!=0) mean_l = sum_l/n_ele;
  }
  // Right (centered on d0+1)
  if(d0<width_score-2)
  {
    sum_r = score_k[d0][ang] + score_k[d0+1][ang] + score_k[d0+2][ang];
    if(n_ele!=0) mean_r = sum_r/n_ele;
  }
  // Center (centered on d0)
  if(d0>0 && d0<width_score-1)
  {
    sum = score_k[d0-1][ang] + score_k[d0][ang] + score_k[d0+1][ang];
    if(n_ele!=0) mean = sum/n_ele;
  }
  // We update d0 according to the highest value
  float val_l = sum_l*mean_l, val = sum*mean, val_r = sum_r*mean_r;
  float max_val = max(val_l,max(val,val_r));
  if(max_val == val_l)
    d0 = d0-1;
  if(max_val == val_r)
    d0 = d0+1;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Marking a region with first_max_score inside the hough space 
// given by a distance interval [d_min,dmax] and the angle ang
void region_marking(vector<vector<float> > &score_k, int d_min, int d_max, 
                    int ang, float first_max_score, int width_score)
{
  //We chech the limits of the minimum and maximum distances
  if(d_min<0) d_min = 0;
  if(d_max>=width_score) d_max=width_score-1;
  for(int dist=d_min; dist<=d_max; dist++)
    score_k[dist][ang] = first_max_score;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Updating the value of the marking radius
void update_marking_radius(int kang, float angle_resolution, 
                           float distance_point_line_max, int &Ri)
{
  int f=5;
  if((kang*angle_resolution*f)<(distance_point_line_max+2))
    Ri = distance_point_line_max+2;
  else
    Ri = kang*angle_resolution*f;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Maximum selection inside the voting matrix
/**
 * \fn void maximum_selection(vector<vector<float> > &score_k,
                              const vector<int> &i_pos,
                              const vector<int> &j_pos,
                              double &max_score,
                              double &votation_score,
                              int width_score,
                              int height_score,
                              int nlines_plus,
                              float distance_resolution,
                              float distance_point_line_max,
                              float angle_point_orientation_max_difference,
                              float angle_resolution,
                              double *sine, double *cosine,
                              vector<line_points> &lines
                             )
 * \brief Selection the n most voted lines inside
 * \param[in] score_k The voting matrix
 * \param[in] (i_pos,j_pos) Indexes of the position of the lines with a 
 *                          significant amount of votes
 * \param[out] max_score Maximum score of each selected line
 * \param[out] votation_score Total amount of votes of the slice
 * \param[in] (width_score,height_score) Size of the voting matrix
 * \param[in] nlines_plus Number of lines to detect
 * \param[in] distance_resolution Resolution of the distance of the Hough voting
 *                                matrix
 * \param[in] distance_point_line_max Maximum allowed distance between the edge 
 *                                    point and the associated line
 * \param[in] angle_point_orientation_max_difference Maximum angle difference 
 *     between the orientation of the edge point and the orientation of the line
 * \param[in] angle_resolution Resolution of the angle of the Hough voting 
 *                             matrix
 * \param[in] (sine,cosine) Orientation of the edge points
 * \param[out] lines line_points vector in which we store the detected lines
 * \return The n most voted lines
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void maximum_selection(vector<vector<float> > &score_k,
                       const vector<int> &i_pos,
                       const vector<int> &j_pos,
                       double &max_score,
                       double &votation_score,
                       int width_score,
                       int height_score,
                       int nlines_plus,
                       float distance_resolution,
                       float distance_point_line_max,
                       float angle_point_orientation_max_difference,
                       float angle_resolution,
                       double *sine, double *cosine,
                       vector<line_points> &lines
                      )
{
  float first_max_score = 0;
  vector<int> m(nlines_plus);
  vector<int> n(nlines_plus);
  for(int l=0;l<nlines_plus;l++)
  {
    max_score=0;
    //We select a local maximum inside the slice (with a neigborhood of one)
    for(int ind_pos=0; ind_pos<(int)i_pos.size(); ind_pos++)
    {
      int i = i_pos[ind_pos];
      int j = j_pos[ind_pos];
      int iu,id,jl,jr;
      iu = id = jl = jr = 0;
      iu = (i==0)              ? height_score-1 : i-1;
      id = (i==height_score-1) ? 0              : i+1;
      jl = (j==0)              ? width_score-1  : j-1;
      jr = (j==width_score-1)  ? 0              : j+1;
      if(score_k[j][i]>=max_score && score_k[j][i]!=first_max_score &&
         (score_k[j][i]>=score_k[jl][iu] || score_k[jl][iu]==first_max_score) &&
         (score_k[j][i]>=score_k[j][iu]  || score_k[j][iu]==first_max_score)  &&
         (score_k[j][i]>=score_k[jr][iu] || score_k[jr][iu]==first_max_score) &&
         (score_k[j][i]>=score_k[jl][i]  || score_k[jl][i]==first_max_score)  &&
         (score_k[j][i]>=score_k[jr][i]  || score_k[jr][i]==first_max_score)  &&
         (score_k[j][i]>=score_k[jl][id] || score_k[jl][id]==first_max_score) &&
         (score_k[j][i]>=score_k[j][id]  || score_k[j][id]==first_max_score)  &&
          (score_k[j][i]>=score_k[jr][id] || score_k[jr][id]==first_max_score))
      {
        max_score = score_k[j][i];
        m[l] = i; n[l] = j;
      }
    } // End ind for
    //We save the line parameters
    lines[l].set_a((float)sine[m[l]]);
    lines[l].set_b((float)cosine[m[l]]);
    lines[l].set_c((float)(n[l]-width_score/2)*distance_resolution);
    //We increment the votation score of the slice
    votation_score+=max_score;
    //**************************************************************************
    //We take the first maximum value+0.01
    if(l==0)
      first_max_score = max_score+0.01;
    //We set to first_max_score a neighborhood of the maximum score point
    int line_angle = m[l];
    int d0_top, d0_bottom;
    d0_top = d0_bottom = n[l];
    int alpha_min=(int) (line_angle-angle_point_orientation_max_difference/
                                                            angle_resolution);
    int alpha_max=(int) (line_angle+angle_point_orientation_max_difference/
                                                            angle_resolution);
    int Ri_top,Ri_bottom;
    Ri_top = Ri_bottom = distance_point_line_max+2;
    int d_min_top, d_min_bottom;
    d_min_top = d_min_bottom = d0_top-Ri_top;
    int d_max_top, d_max_bottom;
    d_max_top = d_max_bottom = d0_top+Ri_top;
    bool first_pos=true, first_neg=true;
    int kang_top, kang_bottom;
    kang_top = kang_bottom = 0;
    //int f = 5;
    score_k[n[l]][m[l]] = first_max_score;
    
    for(int ang_top=line_angle, ang_bottom=line_angle; 
        ang_top<alpha_max || ang_bottom>=alpha_min; 
        ang_top++, ang_bottom--)
    {
      int alpha_top = ang_top;
      int alpha_bottom = ang_bottom;
      if(ang_top<alpha_max)
      {
        //Checking the angle top limit
        if(ang_top>=height_score)
        {
          alpha_top = ang_top - height_score;
          if(first_pos)
          {
            d0_top = width_score/2 - (d0_top - width_score/2);
            first_pos = false;
          }
        }
        //We update top distance center according to the maximum average 
        //sum
        update_distance(score_k, d0_top, alpha_top, Ri_top, width_score);
        //We update the minimum and maximum distances using the distance updated 
        //and the radius
        d_min_top    = d0_top - Ri_top; 
        d_max_top    = d0_top + Ri_top;
        //We mark the region in the score matrix given by the distances and the 
        //angle
        region_marking(score_k, d_min_top, d_max_top, alpha_top, 
                       first_max_score, width_score);
         //We update the value of marking radius for the next iteration
        update_marking_radius(kang_top, angle_resolution, 
                            distance_point_line_max, Ri_top);
        kang_top++;
      }
      
      if(ang_bottom>=alpha_min)
      {
        //Checking the angle bottom limit
        if(ang_bottom<0)
        {
          alpha_bottom = height_score+ang_bottom;
          if(first_neg)
          {
            d0_bottom = width_score/2 + (width_score/2 - d0_bottom);
            first_neg = false;
          }
        }
        //We update bottom distance center according to the maximum average 
        //sum
        update_distance(score_k, d0_bottom, alpha_bottom, Ri_bottom, 
                        width_score);
        //We update the minimum and maximum distances using the distance updated 
        //and the radius
        d_min_bottom = d0_bottom - Ri_bottom;
        d_max_bottom = d0_bottom + Ri_bottom;
        //We mark the region in the score matrix given by the distances and the 
        //angle
        region_marking(score_k, d_min_bottom, d_max_bottom, alpha_bottom, 
                      first_max_score, width_score);
        //We update the value of marking radius for the next iteration
        update_marking_radius(kang_bottom, angle_resolution, 
                              distance_point_line_max, Ri_bottom);
        kang_bottom++;
      }
    }
  } // End lines loop
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Associate points to lines according to thier distance
/**
 * \fn void capturing_points(const image_contours &contours_modified,
                      const image_contours &contours,
                      float angle_point_orientation_max_difference,
                      image_primitives &image_primitive,
                      image_primitives &image_primitive_corrected,
                      image_primitives &image_primitive_original,
                      int nlines_plus,
                      float max_distance,
                      const vector<int> &index
                     )
                     
 * \brief Associate points to lines according to their distance
 * \param[in] contours_modified contours object modified with the 
 *                               input lens distortion model 
 *                               (from update_contours)
 * \param[in] contours Original contours object
 * \param[in] angle_point_orientation_max_difference Maximum angle difference 
 *     between the orientation of the edge point and the orientation of the line
 * \param[in,out] image_primitive Set of detected primitives
 * \param[out] image_primitive_corrected Set of corrected primitives
 * \param[out] image_primitive_original Set of original primitives
 * \param[in] nlines_plus Total amount of lines to detect
 * \param[in] max_distance Maximum distance between a point an its associated 
 *                         line
 * \param[in] index Indexes of the edge points
 * \return The primitives with the captured points
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void capturing_points(const image_contours &contours_modified,
                      const image_contours &contours,
                      float angle_point_orientation_max_difference,
                      image_primitives &image_primitive,
                      image_primitives &image_primitive_corrected,
                      image_primitives &image_primitive_original,
                      int nlines_plus,
                      float max_distance,
                      const vector<int> &index
                     )
{
  double x2,y2;
  double dot_product_min=cos(ami_pi*angle_point_orientation_max_difference/
                             180.);
  lens_distortion_model ld = image_primitive.get_distortion();
  for(int k=0;k<nlines_plus;k++) 
    image_primitive.get_lines()[k].get_points().clear();
  for(int ind=0; ind<(int)index.size(); ind++)
  {
    int q=index[ind];
    x2=contours_modified.get_x()[q];
    y2=contours_modified.get_y()[q];
    double x2ori = contours.get_x()[q];
    double y2ori = contours.get_y()[q];
    point2d<double> p2(x2,y2);
    point2d<double> p2ori(x2ori,y2ori);
    point2d<double> p2d=ld.evaluation(p2);
    point2d<float> corrected_orientation = orientation_update(p2,
                                            contours_modified.get_sine()[q],
                                            contours_modified.get_cosine()[q],
                                            image_primitive.get_distortion());
    for(int k=0;k<nlines_plus;k++)
    {
      if(fabs(image_primitive.get_lines()[k].get_b()*corrected_orientation.x
          -image_primitive.get_lines()[k].get_a()*corrected_orientation.y )<
          dot_product_min) continue;
      if(fabs(image_primitive.get_lines()[k].evaluation(p2d)) < max_distance)
      {
         image_primitive.get_lines()[k].get_points().push_back(p2);
         image_primitive_corrected.get_lines()[k].get_points().push_back(p2d);
         image_primitive_original.get_lines()[k].get_points().push_back(p2ori);
         break;
      }
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Erasing a point ip in the line il in the image primitives objects
void erase_point(image_primitives &image_primitive_original,
                 image_primitives &image_primitive,
                 image_primitives &image_primitive_corrected,
                 int il, int &ip
                )
{
  image_primitive_original.get_lines()[il].get_points().erase(
             image_primitive_original.get_lines()[il].get_points().begin()+ip);
  image_primitive.get_lines()[il].get_points().erase(
             image_primitive.get_lines()[il].get_points().begin()+ip);
  image_primitive_corrected.get_lines()[il].get_points().erase(
             image_primitive_corrected.get_lines()[il].get_points().begin()+ip);
  ip--;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Computing the orientation sign of a point
double compute_orientation_sign(const line_points l,
                                int width, int ip, 
                                const image_contours &contours)
{
  point2d<double> current_point = l.get_points()[ip];
  int pos  = width*floor(current_point.y+0.5) + floor(current_point.x+0.5);
  double a = l.get_a();
  double b = l.get_b();
  return (contours.get_cosine()[pos]*b + contours.get_sine()[pos]*-a);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Ensuring the consistence of lines orientation
/**
 * \fn void ensure_consistent_line_orientations(
 *                               image_primitives &image_primitive_original,
                                 image_primitives &image_primitive,
                                 image_primitives &image_primitive_corrected,
                                 const image_contours &contours,
                                 int width
                                )
                     
 * \brief Ensuring the consistence of lines orientation, according to the  
 *        predominant orientation of their points
 * \param[out] image_primitive Set of detected primitives
 * \param[out] image_primitive_corrected Set of corrected primitives
 * \param[in,out] image_primitive_original Set of original primitives
 * \param[in] contours Original contours object
 * \param[in] width Image width
 * \return The debugged primitives
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void ensure_consistent_line_orientations(image_primitives &image_primitive,
                          image_primitives &image_primitive_corrected,
                          image_primitives &image_primitive_original,
                          const image_contours &contours,
                          int width
                         )
{
  for(int il=0; il<(int)image_primitive_original.get_lines().size(); il++)
  {
    int pos_count=0, neg_count=0;
    for(int ip=0; 
        ip<(int)image_primitive_original.get_lines()[il].get_points().size(); 
        ip++)
    {
      double orientation_sign = 
        compute_orientation_sign(image_primitive_original.get_lines()[il], 
                                 width, ip, contours);
      if(orientation_sign>0)
        pos_count++;
      if(orientation_sign<0)
        neg_count++;
    }
    if(pos_count!=0 && neg_count!=0)
    {
      double s = (pos_count>neg_count) ? 1 : -1;
      for(int ip=0; 
          ip<(int)image_primitive_original.get_lines()[il].get_points().size(); 
          ip++)
      {
        double orientation_sign = 
          compute_orientation_sign(image_primitive_original.get_lines()[il], 
                                   width, ip, contours);
        if(s*orientation_sign<0)
          erase_point(image_primitive_original,image_primitive,
                      image_primitive_corrected,il,ip);
      }
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Removing lines with less than 20 points and
// joining very close lines
/**
 * \fn void remove_and_join(image_primitives &image_primitive,
                            image_primitives &image_primitive_corrected,
                            image_primitives &image_primitive_original,
                            float angle_point_orientation_max_difference,
                            float distance_point_line_max,
                            int width,
                            const image_contours &contours
                           )
                     
 * \brief Removing short lines and to join lines which are very close
 * \param[in,out] image_primitive Set of detected primitives
 * \param[in,out] image_primitive_corrected Set of corrected primitives
 * \param[in,out] image_primitive_original Set of original primitives
 * \param[in] angle_point_orientation_max_difference Maximum angle difference 
 *     between the orientation of the edge point and the orientation of the line
 * \param[in] distance_point_line_max Maximum allowed distance between the edge 
 *                                    point and the associated line
 * \param[in] width Image width
 * \param[in] contours Original contours object
 * \return The set of primitives with the short lines removed and with the close
 *         lines joined
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void remove_and_join(image_primitives &image_primitive,
                     image_primitives &image_primitive_corrected,
                     image_primitives &image_primitive_original,
                     float angle_point_orientation_max_difference,
                     float distance_point_line_max,
                     int width,
                     const image_contours &contours
                    )
{
  double dot_product_min=cos(2.*ami_pi*angle_point_orientation_max_difference/
                             180.);
  if(dot_product_min<0.95) dot_product_min=0.95;
  
  vector<line_points> &vlip  = image_primitive.get_lines(),
                      &vlipc = image_primitive_corrected.get_lines(),
                      &vlipo = image_primitive_original.get_lines();

  float min_line_points=0.05*vlip[0].get_points().size();
  if(min_line_points<20) min_line_points=20;
  float distance_line_line_min=2.*distance_point_line_max;
  if(distance_line_line_min>10) distance_line_line_min=10.;

  for(int k=0;k<(int)vlip.size();k++)
  {
    line_points &linek  = vlip[k],
                &linekc = vlipc[k],
                &lineko = vlipo[k];
    if(linek.get_points().size()<min_line_points)
    {
      vlip.erase(vlip.begin()+k);
      vlipc.erase(vlipc.begin()+k);
      vlipo.erase(vlipo.begin()+k);
      k--;
      continue;
    }
    for(int l=k+1;l<(int)vlip.size();l++)
    {
      line_points &linel  = vlip[l];
      line_points &linelc = vlipc[l];
      line_points &linelo = vlipo[l];
      // We check the number of points
      if(!linel.get_points().size()>0)
        continue;
      double aux = linek.get_a()*linel.get_a() + linek.get_b()*linel.get_b();
      if(!(fabs(aux)>dot_product_min))
        continue;
      // We check the orientation of the first point of each primitive
      double a1 = linekc.get_a(), a2 = linelc.get_a(), b1 = linekc.get_b(), 
             b2 = linelc.get_b();
      point2d<double> point = lineko.get_points()[0];
      int pos = width*point.y + point.x;
      double sign1 = b1*contours.get_cosine()[pos]-a1*contours.get_sine()[pos];
      point = linelo.get_points()[0];
      pos = width*point.y + point.x;
      double sign2 = b2*contours.get_cosine()[pos]-a2*contours.get_sine()[pos];
      if(sign1*sign2<=0)
          continue;
      // We check the average distance of the points to the line
      aux=0;
      int m;
      for(m=0;m<(int)linel.get_points().size();m++)
      {
        point2d<double> p2d = linelc.get_points()[m];
        double dist=fabs(linek.evaluation(p2d));
        aux+=dist;
        if(dist>(5.*distance_point_line_max))
          break;
      }
      if(m < (int)linel.get_points().size()) 
        continue;
      if((aux/(int)linel.get_points().size())<MIN_POINTS_IN_LINE)
      {
        //We check the distance of the points of a line to the points of the
        //other line. This distance should be big to avoid joining parallel 
        //lines
        aux=0;
        for(m=0;m<(int)linel.get_points().size();m++)
        {
          point2d<double> p2d = linelc.get_points()[m];
          double a,b,c;
          linek.get_abc(a,b,c);
          double d = a*p2d.x + b*p2d.y + c;
          point2d<double> np2d(p2d.x-d*a, p2d.y-d*b);
          aux+=linekc.distance(np2d);
        }
        if((aux/linelc.get_points().size())>MAX_DISTANCE_BETWEEN_LINES)
        {
          // We add the points of the line to the line point structure
          linek.get_points().insert(linek.get_points().end(),
                                    linel.get_points().begin(),
                                    linel.get_points().end());

          linekc.get_points().insert(linekc.get_points().end(),
                                     linelc.get_points().begin(),
                                     linelc.get_points().end());

          lineko.get_points().insert(lineko.get_points().end(),
                                     linelo.get_points().begin(),
                                     linelo.get_points().end());
          // We remove the line points structure
          vlip.erase(vlip.begin()+l);
          vlipc.erase(vlipc.begin()+l);
          vlipo.erase(vlipo.begin()+l);
          l--;
        }
      }
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Computing the line equation using its points
/**
 * \fn void recompute_line_equations(image_primitives &image_primitive,
                                    image_primitives &image_primitive_corrected,
                                    image_primitives &image_primitive_original,
                                    bool lens_distortion_estimation
                                    )
                     
 * \brief Recomputing the equations of the lines using their points
 * \param[in,out] image_primitive Set of detected primitives
 * \param[in,out] image_primitive_corrected Set of corrected primitives
 * \param[out] image_primitive_original Set of original primitives
 * \param[in] lens_distortion_estimation Flag for considering the lens 
 *                                       distortion
 * \return The set of primitives with the recomputed equations
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
void recompute_line_equations(image_primitives &image_primitive,
                              image_primitives &image_primitive_corrected,
                              image_primitives &image_primitive_original,
                              bool lens_distortion_estimation
                             )
{
  double a,b,c;
  if (lens_distortion_estimation==true)
  {
    for(int i=0;i<(int)image_primitive.get_lines().size();i++)
    {
      if(image_primitive.get_lines()[i].get_points().size()>2)
      {
        image_primitive_corrected.get_lines()[i].points_to_equation();
        image_primitive_corrected.get_lines()[i].get_abc(a,b,c);
        image_primitive.get_lines()[i].set_abc(a,b,c);
        image_primitive_original.get_lines()[i].set_abc(a,b,c);
      }
    }
  }
  else
  {
    for(int i=0;i<(int)image_primitive.get_lines().size();i++)
      if(image_primitive.get_lines()[i].get_points().size()>2)
      {
        image_primitive.get_lines()[i].points_to_equation();
        image_primitive.get_lines()[i].get_abc(a,b,c);
        image_primitive_corrected.get_lines()[i].set_abc(a,b,c);
        image_primitive_original.get_lines()[i].set_abc(a,b,c);
      }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
 /**
 * \fn double line_equation_distortion_extraction_improved_hough(
                              const image_contours &contours,
                              image_primitives &image_primitive,
                              float distance_point_line_max, int nlines=100,
                              float angle_resolution=0.1,
                              float distance_resolution=1.,
                              float initial_distortion_parameter=0.0,
                              float final_distortion_parameter=1.0,
                              float distortion_parameter_resolution=0.1,
                              float angle_point_orientation_max_difference=2.,
                              bool lens_distortion_estimation=true,
                              lens_distortion_model ini_ldm = 
                              lens_distortion_model())
 * \brief Computation of lines using an improved version of Hough which includes
          1 parameter lens distortion division model
 * \param[in] contours Contour information
 * \param[out] image_primitive Image primitives where lines and distortion model
 *                             is defined
 * \param[in] distance_point_line_max Maximum distance allowed between points
 *                                    and associated lines
 * \param[in] nlines Number of lines to return
 * \param[in] angle_resolution Angle resolution
 * \param[in] distance_resolution Distance resolution
 * \param[in] initial_distortion_parameter Initial value of normalized
 *                                         distortion parameter
 * \param[in] final_distortion_parameter Final value normalized distortion
 *                                       parameter
 * \param[in] distortion_parameter_resolution Distortion parameter
 *                                            discretization step
 * \param[in] angle_point_orientation_max_difference Maximum difference (in
 *                    degrees) of the point orientation angle and the line angle
 * \param[in] lens_distortion_estimation Boolean to control if we estimate the
 *                                       lens distortion model
 * \param[in] ini_ldm Initial distortion model
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
double line_equation_distortion_extraction_improved_hough(
                                  const image_contours &contours,
                                  image_primitives &image_primitive,
                                  float distance_point_line_max,
                                  int nlines,
                                  float angle_resolution,
                                  float distance_resolution ,
                                  float initial_distortion_parameter,
                                  float final_distortion_parameter,
                                  float distortion_parameter_resolution,
                                  float angle_point_orientation_max_difference,
                                  bool lens_distortion_estimation,
                                  lens_distortion_model ini_ldm)
{
  int l;
  int width_score,height_score,depth_score;
  int width=contours.get_width();
  int height=contours.get_height();
  double best_distortion_parameter=0;
  double *sine,*cosine;
  double max_norm=(double) width*width+height*height;
  point2d<double> dc((double) width/2.,(double) height/2.);
  double dmi=0.;

  // We check discretization parameter values: distance and angle resolution
  assert(distance_resolution>0 && angle_resolution>0 && angle_resolution<=180);
  // We define score volume size
  width_score=(int) (1.2*(2.*sqrt(max_norm)/distance_resolution+2));
  height_score=(int) (180./angle_resolution);
  depth_score=distortion_parameter_resolution>0 ?
    (int)(1.+(final_distortion_parameter-initial_distortion_parameter)/
          distortion_parameter_resolution) : 1;

  // We define orientation discretization vector
  sine   =(double*)malloc(sizeof(double)*height_score);
  cosine =(double*)malloc(sizeof(double)*height_score);
  double step_angle=angle_resolution*ami_pi/180.;
  for(l=0;l<height_score;l++) {
      cosine[l]=cos((double) l*step_angle);
      sine[l]=sin((double) l*step_angle);
  }
  
  //****************************************************************************
  // Update the contours object with the initial lens distortion model
  image_contours contours_modified(width, height);
  const vector<int> &index = contours.get_index();
 
  update_contours(ini_ldm,contours,contours_modified,index);
  
  //We update dmi
  if(ini_ldm.get_type()==DIVISION)
    dmi = ini_ldm.is_identity() ?
              dc.norm()*dc.norm() :
              update_rsqmax(ini_ldm.get_distortion_center(), width, height);

  //****************************************************************************
  // We fill score volume
  double max_votation_score = -1;
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int k=0; k<depth_score; k++)
  {
    // Vector of the actual slice
    vector<vector<float> > score_k(width_score, 
                                   std::vector<float>(height_score,0.));

    // We compute the lens distortion parameter value for each iteration
    double kdistortion = 0.0;
    double p = initial_distortion_parameter+(k*distortion_parameter_resolution);
    if(ini_ldm.get_type()==DIVISION)
    {
      kdistortion = -p / (dmi + dmi*p);
    }
    else
    {
      kdistortion = p/max_norm;
    }

    // We build the model for this iteration
    lens_distortion_model it_ldm;
    it_ldm.set_distortion_center(dc);
    it_ldm.get_d().resize(2);
    it_ldm.get_d()[0] = 1.;
    it_ldm.get_d()[1] = kdistortion;
    it_ldm.set_type(ini_ldm.get_type());
    
    voting(contours_modified,it_ldm,step_angle,height_score,width_score,
           angle_point_orientation_max_difference/angle_resolution,
           distance_resolution,sine, cosine,score_k,index);

    // We create two vectors with the positions inside the score matrix with a
    // score higher than 10
    vector<int> i_pos(width_score*height_score,0);
    vector<int> j_pos(width_score*height_score,0);
    int pind=0;
    for(int i=0;i<height_score;i++)
    {
      for(int j=0;j<width_score;j++)
      {
        if(score_k[j][i]>10)
        {
          i_pos[pind] = i;
          j_pos[pind] = j;
          pind++;
        }
      }
    }
    i_pos.resize(pind);
    j_pos.resize(pind);
    //**************************************************************************
    // We select the maximum of the slice
    double votation_score=0., max_score;
    vector<line_points> lines(nlines);
    
    maximum_selection(score_k,i_pos,j_pos,max_score,votation_score,
                      width_score,height_score,nlines,
                      distance_resolution,distance_point_line_max,
                      angle_point_orientation_max_difference,
                      angle_resolution,sine,cosine,lines);
    //**************************************************************************
    // We select the maximun distortion level
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    if(votation_score>max_votation_score)
    {
      max_votation_score=votation_score;
      image_primitive.set_lines(lines);
      best_distortion_parameter = kdistortion;
    }
  } // End k loop
  //****************************************************************************
  // Image_primitives objects for corrected points and original points
  image_primitives image_primitive_corrected, image_primitive_original;
  image_primitive_corrected.set_lines(image_primitive.get_lines());
  image_primitive_original.set_lines(image_primitive.get_lines());

  // We fill the image primitive lens distortion model
  lens_distortion_model ld;
  ld.set_distortion_center(dc);
  ld.get_d().resize(2);
  ld.get_d()[0]=1.;
  ld.get_d()[1]=best_distortion_parameter;
  ld.set_type(ini_ldm.get_type());
  image_primitive.set_distortion(ld);
  image_primitive_original.set_distortion(ld);
  //****************************************************************************
  // We compute the points of the line points
  // following the distance of the points to the line
  float max_distance = (distance_point_line_max+0.5*distance_resolution);
  capturing_points(contours_modified,contours,
                   angle_point_orientation_max_difference,
                   image_primitive,image_primitive_corrected,
                   image_primitive_original,nlines,
                   max_distance,index);

  //****************************************************************************
  // Debugging of the primitives through the orientation of the points
  ensure_consistent_line_orientations(image_primitive,image_primitive_corrected,
                                      image_primitive_original,contours,width);

  //****************************************************************************
  // We recompute the line equations
  recompute_line_equations(image_primitive,image_primitive_corrected,
                           image_primitive_original,lens_distortion_estimation);

  //****************************************************************************
  // We remove lines with a small number of points and we merge lines
  // which are too close
  remove_and_join(image_primitive,image_primitive_corrected,
                  image_primitive_original,
                  angle_point_orientation_max_difference,
                  distance_point_line_max,width,contours);
  //****************************************************************************
  // After removing the lines, we recompute the equations
  recompute_line_equations(image_primitive,image_primitive_corrected,
                           image_primitive_original,lens_distortion_estimation);
  
  //****************************************************************************
  free(sine); free(cosine);
  image_primitive = image_primitive_original;
  return max_votation_score;
}
//------------------------------------------------------------------------------