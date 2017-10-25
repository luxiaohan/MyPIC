#ifndef THINLENSMAP_H
#define THINLENSMAP_H

#include "ElementMap.h"
class ThinlensMap:public ElementMap
{
public:
  ThinlensMap(double focallengthx,double focallengthy);//mm
  void Map();
  void print();
protected:
double fx,fy;
MatrixXd R2= MatrixXd::Zero(3,3);
};
#endif
