#ifndef GAPMAP_H
#define GAPMAP_H

#include "ElementMap.h"
#include "Field.h"
#include "MyFunc.h"
#include "DriftMap.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
class GapMap: public ElementMap
{
public:
  virtual ~GapMap();
  GapMap();
  GapMap(double length,double E0T,double phase,double frequency,double gamma);//(mm,MV,deg,mm)
  GapMap(double length,const string &p,double phase,double frequency,double gamma);
  GapMap(double tubelength,double gaplength,double radius,double phase_at_entry,double frequency,double gamma,int nharm);
  void Map();
  double getEf();
  void FieldMethod(double &x,double &px,double &y,double &py,double &z,double &pz);
  void print();
  MatrixXd getmap();
  MatrixXd getmap(double x,double px,double y,double py,double z,double dp);
  
  void MakeMap(double tubelength,double gaplength,double radius,double frequency,int nharm);
  void ULSCF(double lamada,double tube,double gap,double radius,vector<double> poten,int nharm,vector<double> &coeff,vector<double> &coeff1);

protected:
  MatrixXd R1=MatrixXd::Identity(6,6);

  Field F1;
  double Energygain,wavelength,frequency,phase,E0T;
  int Slice=0,SliceMatrix=0,phaseSYNflag=1,E0TFlag=0,FieldFlag=0;

  
};

#endif
