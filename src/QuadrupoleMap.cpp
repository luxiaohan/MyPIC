#include "QuadrupoleMap.h"

QuadrupoleMap::~QuadrupoleMap()
{
}

QuadrupoleMap::QuadrupoleMap(){}

QuadrupoleMap::QuadrupoleMap(double length,double gradient,double gamma)//(m,T/m)
{
  setlength(length);
  setgradient(gradient);
  setgamma(gamma);
  Map();
}
void QuadrupoleMap::Map()
{
  Bro=gamma*_MASS*C_light*beta/_Q;//T m
  //cout<<"BRO  "<<Bro<<"   "<<gamma<<"  "<<beta<<endl;
  //cout<<"Foucus length  "<<1/(gradient * length /Bro)<<endl;
  k=sqrt(abs(gradient/Bro));
  //cout<<k*k<<endl;
  if(gradient>0)
  {
    //xx
    R(0,0)=cos(k*length);	R(0,1)=sin(k*length)/k;
    R(1,0)=-k*sin(k*length);	R(1,1)=cos(k*length);
    //yy
    R(2,2)=cosh(k*length);	R(2,3)=sinh(k*length)/k;
    R(3,2)=k*sinh(k*length);	R(3,3)=cosh(k*length);
    //zz
    R(4,4)=1;			R(4,5)=length/gamma/gamma;
    R(5,4)=0;			R(5,5)=1;
  }
  else if(gradient<=0)
  {
    //xx
    R(0,0)=cosh(k*length);	R(0,1)=sinh(k*length)/k;
    R(1,0)=k*sinh(k*length);	R(1,1)=cosh(k*length);
    //yy
    R(2,2)=cos(k*length);	R(2,3)=sin(k*length)/k;
    R(3,2)=-k*sin(k*length);	R(3,3)=cos(k*length);
    //zz
    R(4,4)=1;			R(4,5)=length/gamma/gamma;
    R(5,4)=0;			R(5,5)=1;
  }
  //
}

const vector<double> QuadrupoleMap::getfield(const vector<double> &xyz)
{
  return vector<double>{0.0,			//Ex
    			0.0,
			0.0,
			gradient * xyz.at(2), 
			gradient * xyz.at(0),
			0.0};

}


void QuadrupoleMap::setgradient(double gra)
{
  gradient=gra;//T/m
}


void QuadrupoleMap::print(){
  cout<<"The transfer matrix of quadrupole(length = "<<length<<" mm,gradient = "<<gradient*1000<<" T/m) is\n"<<R<<endl;
}




