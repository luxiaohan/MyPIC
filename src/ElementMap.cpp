#include "ElementMap.h"

ElementMap::~ElementMap(){}
ElementMap::ElementMap(){}
ElementMap::ElementMap(double length)
{//(m,T/m)
  this->length=length;
};


void ElementMap::setlength(double length)
{
  this->length=length;
}

void ElementMap::setbeta(double a)
{
  if(a<=1&&a>=0)
  {
    beta=a;
    gamma=1/sqrt(1-beta*beta);
  }else{}
}

void ElementMap::setgamma(double a)
{
  if(a>1)
  {
    gamma=a;
    beta=sqrt(1-1/gamma/gamma);
  }else{}
}

void ElementMap::Map(){}
const vector<double> ElementMap::getfield(const vector<double> &xyz){}
void ElementMap::print()
{
  cout<<"Element cout"<<endl;
}
