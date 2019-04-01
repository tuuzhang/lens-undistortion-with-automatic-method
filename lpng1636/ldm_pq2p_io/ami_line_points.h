/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file line_points.h
 * \brief line_points class AMI_DLL_H to store line equation and associated
          points
 * \author Luis Alvarez \n \n
*/

#ifndef LINE_POINTS_H
#define LINE_POINTS_H
#include "ami_line.h"
#include "ami_point2d.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace ami;

namespace ami
{

/**
 * \class  line_points
 * \author Luis Alvarez
 * \brief Class  to store together a collection of aligned points and line
            equation and basic method
*/
class AMI_DLL_H  line_points{
  ami::line rect/** Line defined by the equation (ax+by+c=0)*/;
  std::vector< ami::point2d<double> > points/** Line points */;


public:
  /**
  * \fn line_points()
  * \brief Constructor
  * \author Luis Alvarez
  */
  line_points(){ rect.set_abc((double) 0.,(double) 0,(double) 0);}

  /**
  * \fn line_points()
  * \brief Constructor from a line_points object
  * \author Luis Alvarez
  */
  line_points(const line_points &line){ rect = line.rect; points = line.points;}

  /**
  * \fn line_points(double a,double b,double c)
  * \author Luis Alvarez
  * \param a,b,c
  */
  line_points(double a/** Coefficient 'a' of the line */,
              double b/** Coefficient 'b' of the line */,
              double c/** Coefficient 'c' of the line */){ rect.set_abc(a,b,c);}

  /**
  * \fn ~line_points()
  * \author Luis Alvarez
  */
  ~line_points(){};

  /**
  * \fn std::vector<point2d> get_points()
  * \author Luis Alvarez
  */
  std::vector<point2d<double> > &get_points() {return points;}
 
  /**
  * \fn const std::vector<point2d> &get_points()
  * \author Luis Alvarez
  */
  const std::vector<point2d<double> > &get_points() const {
    return points;} /** RETURN p*/

  /**
  * \fn  const ami::line &get_rect()  const
  * \author Luis Alvarez
  */
  const  ami::line &get_rect()  const {return rect;}
  
  /**
  * \fn const double &get_a()
  * \author Luis Alvarez
  */
  const double &get_a() const {return rect.get_a();}
  
  /**
  * \fn const double &get_b()
  * \author Luis Alvarez
  */
  const double &get_b() const {return rect.get_b();}
  
  /**
  * \fn const double &get_c()
  * \author Luis Alvarez
  */
  const double &get_c() const {return rect.get_c();}

 /**
  * \fn get_abc(double &a2,double &b2,double &c2)
  * \author Luis Alvarez
  */
  void get_abc(double &a2/** a coefficient of the line */,
               double &b2/** b coefficient of the line */,
               double &c2/** c coefficient of the line */)
  {a2=rect.get_a(); b2=rect.get_b(); c2=rect.get_c();}

  /**
  * \fn void set_a(double a2)
  * \author Luis Alvarez
  */
  void set_a(double a2/** a coefficient of the line */){rect.set_a(a2);}
 
  /**
  * \fn void set_b(double b2)
  * \author Luis Alvarez
  */
  void set_b(double b2/** b coefficient of the line */){rect.set_b(b2);}
 
  /**
  * \fn void set_c(double c2)
  * \author Luis Alvarez
  */
  void set_c(double c2/** c coefficient of the line */){rect.set_c(c2);}

  /**
  * \fn void set_abc(double a,double b,double c){rect.set_abc(a,b,c);}
  * \author Luis Alvarez
  */
  void set_abc(double a/** a coefficient of the line */,
               double b/** b coefficient of the line */,
               double c/** c coefficient of the line */){rect.set_abc(a,b,c);}
  
  /**
  * \fn void set_points(std::vector<point2d<double> > p2)
  * \author Luis Alvarez
  */
  void set_points(std::vector<point2d<double> > p2/** Vector of points */){
    points=p2;};

  double points_to_equation()/** COMPUTES RECT FROM x,y*/;
  
  /**
  * \fn double evaluation(const point2d<double> &p) const
  * \author Luis Alvarez
  */
  double evaluation(const point2d<double> &p /** 2D point */) const {
    return (rect.get_a()*p.x+rect.get_b()*p.y+rect.get_c());
  }

  line_points & operator= (const line_points &line)
  {
    rect = line.rect;
    points = line.points;
    return *this;
  }
  
  double distance(point2d<double> &p /** point2d */);

};

}

#endif
