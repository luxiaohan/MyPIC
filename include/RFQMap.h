#ifndef RFQMAP_H
#define RFQMAP_H


#include "ElementMap.h"

class RFQMap : public ElementMap
{
  public:
    virtual ~RFQMap();
    RFQMap();
    RFQMap(double length,double voltage,double phase,double radius,double gamma,double frequency=499.8e6,double AAcc=10000000.0,double AFocus=1.0,char particletype='e');
    void initial(double length,double voltage,double phase,double radius,double gamma,double frequency=499800000,double AAcc=10000.0,double AFocus=1,char particletype='e');
    virtual void Map();
    virtual void print();
    void settype(int,int,int);
    void clear();
  protected:
    vector<MatrixXd> dR;
    MatrixXd Rt1=MatrixXd::Identity(2,2);
    MatrixXd Rt2=MatrixXd::Identity(2,2);
    MatrixXd Rt3=MatrixXd::Identity(2,2);
    double step=0.1;
    double gammain,gammaout,gammasyn,betain,betaout,betasyn,tsyn,zsyn;
    int N,S,type,pretype,nexttype;
    double x,dx,y,dy,dz,dp;
    double kx1,kx2,ky1,ky2,kz1,kz2,C1,C2,C3,Q,Mass,omiga;
    double voltage,phase,AAcc,AFocus,radius;
    char particletype;
    inline int sign(double x)
    {
      if(x < 0) return -1;
      else if(x==0) return 0;
      else return 1;
    }
};
#endif


