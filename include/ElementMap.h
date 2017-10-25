#ifndef ELEMENTMAP_H
#define ELEMENTMAP_H

//const double _MASS	=	(1.672621777e-27);	//kg electron
//const double _Q		=	(1.602176565e-19);	//C
const double _MASS	=	(1.672621777e-27);	//kg proton
const double _Q		=	(1.602176565e-19);	//C
const double C_light	=	299792458.0;		//m/s
const double Epsilon	=	8.854187817e-12;	//F/m

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace Eigen;
using namespace std;
class ElementMap
{
public:
  virtual ~ElementMap();
  ElementMap();
  ElementMap(double);
  ElementMap(double,double);
  virtual void Map();
  virtual MatrixXd getmap() const {return R;}
  virtual MatrixXd getmap()       {return R;}
  virtual const vector<double> getfield(const vector<double> &xyz);
  void setlength(double);
  void setbeta(double);
  void setgamma(double);
  virtual void print();
protected:
  MatrixXd R=MatrixXd::Identity(6,6);
  double beta,gamma;
  double gradient,field,k,kx,ky,h,Bro,aperture;
  double length,location;
};
#endif
