#include "DriftMap.h"
DriftMap::~DriftMap(){}
DriftMap::DriftMap(){}
DriftMap::DriftMap(double length,double gamma):ElementMap(length)
{
  setgamma(gamma);
  Map();
}
void DriftMap::Map()
{
  //xx
  R(0,0)=1;
  R(0,1)=length;
  R(1,0)=0;
  R(1,1)=1;
  //yy
  R(2,2)=1;
  R(2,3)=length;
  R(3,2)=0;
  R(3,3)=1;
  //zz
  R(4,4)=1;
  R(4,5)=length/gamma/gamma;
  R(5,4)=0;
  R(5,5)=1;
  //zx
  //xz
}
const vector<double> DriftMap::getfield(const vector<double> &xyz)
{
  return vector<double>{0.0,			//Ex
    0.0,			//Ey
    0.0,			//Ez
    0.0,			//Bx
    0.0,			//By
    0.0};			//Bz
}
void DriftMap::print()
{
  cout<<"The transfer matrix of drift(length = "<<length<<" m) is\n"<<R<<endl;
}
