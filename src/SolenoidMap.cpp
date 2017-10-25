#include "SolenoidMap.h"


SolenoidMap::SolenoidMap(double length,double field,double gamma)
{
  setlength(length);
  setfield(field);
  setgamma(gamma);
  Map();
}
void SolenoidMap::Map()
{
/*
  Bro=_MASS/_Q*gamma*beta*C_light;//T m
  k=field/(2*Bro);
  //cout<<"GAMMA"<<_MASS/_Q*gamma*0.001/field<<"   "<<beta<<endl;
  //xx
  R(0,0)=cos(k*length)*cos(k*length);		R(0,1)=sin(k*length)*cos(k*length)/k;
  R(1,0)=-k*sin(k*length)*cos(k*length);	R(1,1)=cos(k*length)*cos(k*length);
  //yy
  R(2,2)=cos(k*length)*cos(k*length);		R(2,3)=sin(k*length)*cos(k*length)/k;
  R(3,2)=-k*sin(k*length)*cos(k*length);	R(3,3)=cos(k*length)*cos(k*length);
  //zz
  R(4,4)=1;					R(4,5)=length/gamma/gamma;
  R(5,4)=0;					R(5,5)=1;
  //xy
  R(0,2)=sin(k*length)*cos(k*length);		R(0,3)=sin(k*length)*sin(k*length)/k;
  R(1,2)=-k*sin(k*length)*sin(k*length);	R(1,3)=sin(k*length)*cos(k*length);
  //yx
  R(2,0)=-sin(k*length)*cos(k*length);		R(2,1)=-sin(k*length)*sin(k*length)/k;
  R(3,0)=k*sin(k*length)*sin(k*length);		R(3,1)=-sin(k*length)*cos(k*length);
*/

//ONLY for the motion throught a longitudinal field. NOT consider the frige field
  Bro=gamma*_MASS*beta*C_light/_Q;//T m
  k=field/(Bro);
  //xx
  R(0,0)=1;	R(0,1)=sin(k*length)/k;
  R(1,0)=0;	R(1,1)=cos(k*length);
  //yy
  R(2,2)=1;	R(2,3)=sin(k*length)/k;
  R(3,2)=0;	R(3,3)=cos(k*length);
  //zz
  R(4,4)=1;					R(4,5)=length/gamma/gamma;
  R(5,4)=0;					R(5,5)=1;
  //xy
  R(0,2)=0;	R(0,3)=(cos(k*length)-1)/k;
  R(1,2)=0;	R(1,3)=-sin(k*length);
  //yx
  R(2,0)=0;	R(2,1)=-(cos(k*length)-1)/k;
  R(3,0)=0;	R(3,1)=sin(k*length);

  
  //zx
  R(4,0)=0;		R(4,1)=0;
  R(5,0)=0;		R(5,1)=0;
  //xz
  R(0,4)=0;		R(0,5)=0;
  R(1,4)=0;		R(1,5)=0;
}
const vector<double> SolenoidMap::getfield(const vector<double> &xyz)
{
  return vector<double>{0.0,0.0,0.0,0.0,0.0,field};
}

void SolenoidMap::setfield(double f)
{
  field = f;
}

void SolenoidMap::print()
{
  cout<<"The transfer matrix of Solenoid(length = "<<length<<" mm,Field = "<<field<<" T) is\n"<<R<<endl;
  cout<<R.block<4,4>(0,0).determinant()<<endl;
  cout<<"R.block<2,2>(0,0) "<<R.block<2,2>(0,0).determinant()<<endl;
  cout<<"R.block<2,2>(0,2) "<<R.block<2,2>(0,2).determinant()<<endl;
  cout<<"R.block<2,2>(2,0) "<<R.block<2,2>(2,0).determinant()<<endl;
  cout<<"R.block<2,2>(2,2) "<<R.block<2,2>(2,2).determinant()<<endl;
}
