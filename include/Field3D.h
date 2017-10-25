#ifndef FIELD3D_H 
#define FIELD3D_H

#include "ElementMap.h"
#include <string>
#include <vector>
using namespace std;
using vd	=vector<double>;
using v2d	=vector<vector<double> > ;
class Field3D:public ElementMap
{
public:
  virtual ~Field3D();
  Field3D();
  Field3D(const string &path,double,double,double,double,double,int);//(mm,T/m)
  const vector<double> getfield(const vector<double> &xyz,double time);
  void synchro(double,double,double);
  void SetOffset(double);
  void ScanOffset(double);
  virtual void Map();
  virtual void print();
protected:
  v2d cor;
private:
  string filename;
  int type=1;
  double ke,kb;
  double locationEnd;
  double frq;
  vd zero=vd{0,0,0,0,0,0};
  //type0 get field
  double d1,d2;
  v2d gridLocation;
  v2d eField;
  v2d mField;
  vd fieldofParticle=vd{0,0,0,0,0,0};
  //type1
  int nx,ny,nz;
  double locationMax;
  double xmin,xmax,ymin,ymax,zmax;
  v2d Field;
  double gridNum;
  double synPhase,offsetPhase;
  double spacex,spacey,spacez;
  //type2
  int nr;	//z had been defined before
  double rmax;
  double spacer;

};

#endif
