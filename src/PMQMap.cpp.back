#include "PMQMap.h"

PMQMap::~PMQMap(){}
PMQMap::PMQMap(){}
PMQMap::PMQMap(double length,double radius,double aperture,double gamma)//mm,m,mm
{
  this->length=length;
  this->radius=radius*1000;
  this->aperture=aperture;
  setgamma(gamma);
  Map();
}
void BendingMap::Map()
{
  angle=length/radius;
  h=angle/abs(radius)/abs(angle);
  kx=abs(h);
  ky=0;
  length=abs(radius)*angle;
//xx
  R(0,0)=cos(kx*length);
  R(0,1)=sin(kx*length)/(kx*1000);
  R(1,0)=-kx*1000*sin(kx*length);
  R(1,1)=cos(kx*length);
//yy
  R(2,2)=1;
  R(2,3)=length/1000;
  R(3,2)=0;
  R(3,3)=1;
//zz
  R(4,4)=1;
  R(4,5)=(-h*h*(kx*length*beta*beta-sin(kx*length))/(kx*kx*kx)+length/pow(gamma,2)*(1-h+h/kx/kx))/1000;
  R(5,4)=0;
  R(5,5)=1;
//zx
  R(4,0)=-h*sin(kx*length)/kx;
  R(4,1)=-h*(1-cos(kx*length))/kx/kx/1000;
  R(5,0)=0;
  R(5,1)=0;
//xz
  R(0,4)=0;
  R(0,5)=h*(1-cos(kx*length))/1000;
  R(1,4)=0;
  R(1,5)=h*sin(kx*length)/kx;
//
}
void BendingMap::setangle(double angle)
{
  this->angle=angle/180*M_PI;
}
void BendingMap::setradius(double radius)
{
  this->radius=radius*1000;//m transfer to mm
}
void BendingMap::setaperture(double aperture)
{
  this->aperture = aperture;
}
void BendingMap::print()
{
  cout<<"The transfer matrix of bending(radius = "<<radius/1000<<" m,angle = "<<angle/M_PI*180<<" deg) is\n"<<R<<endl;
}
