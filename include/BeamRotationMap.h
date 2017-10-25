#ifndef BEAMROTATIONMAP_H
#define BEAMROTATIONMAP_H

#include "ElementMap.h"
class BeamRotationMap:public ElementMap
{
public:
  virtual ~BeamRotationMap();
  BeamRotationMap();
  BeamRotationMap(double,double,double);//angle//°
  void Map();
  void print();
protected:
  double angleXY,angleXZ,angleZY;
  MatrixXd Rxy=MatrixXd::Zero(6,6);
  MatrixXd Rxz=MatrixXd::Zero(6,6);
  MatrixXd Rzy=MatrixXd::Zero(6,6);
};
#endif
