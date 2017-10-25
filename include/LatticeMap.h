#ifndef LATTICEMAP_H
#define LATTICEMAP_H

#include "AllElement.h"
#include <vector>
#include <string>

using namespace std;
using namespace Eigen;

struct Element
{
  bool operator== (const Element& e)
  {
    if(name==e.name&&p.size()==e.p.size())
    {
      for(int i=0;i<p.size();i++)
      {
        if(p[i]!=e.p[i])
          return false;
      }
      return true;
    }
    return false;
  }
  bool operator!= (const Element& e)
  {
    if(name==e.name&&p.size()==e.p.size())
    {
      for(int i=0;i<p.size();i++)
      {
        if(p[i]!=e.p[i])
          return true;
      }
      return false;
    }
    return true;
  }
  string name;
  string path;
  vector<double> p;
  double position=0,length=0;
  int ordernumber=0;
};
class LatticeMap
{
public:
  virtual ~LatticeMap();
  LatticeMap();
  LatticeMap(const LatticeMap& );
  LatticeMap& operator= (const LatticeMap&);
  bool operator== (const LatticeMap&);
  void callelement(Element &v,ifstream &infile);
  void add(const int &num,const Element &v);
  void modify(const int &num,const Element &v);
  double read(const char *p,double EGamma,int numberdivided=1,int numberrepeat=1);

  void del(int Nelement);
  void deleteall();
  void print();
  void prelattice(const char *in,const char *out,int n);
  void prelattice1(const char *in,const char *out);
  void addEGamma(const double num);

  MatrixXd getMat(int L,double &x,double &dx,double &y,double &dy,double &z,double &dp);
  const vector<double> getFieldBase(const vector<double> &xyz,double tolerrance,double time);
  const vector<double> getFieldPre(const vector<double> &xyz,double time);
  const vector<double> getField(const vector<double> &xyz,double time);
  const int getEleNum(double z);
  double    getRadius(double z);
  vector<MatrixXd>		Mat;
  vector<Element>		ELE;

  vector<double> EBeta;
  vector<double> EGamma;
  double fitness;
  double latticeLength;
//protected:
  //Element
  Element elemtemp;

  vector<QuadrupoleMap>		FQ;
  vector<BendingMap>		BM;
  vector<DriftMap>		DL;
  vector<SolenoidMap>		SN;
  vector<RFQMap>		RM;
  vector<RFQ>			RQ;
  vector<Field3D>		FIELD3D;

  
  //Common
  int NDivide;

  int Mfail=0,particle_count=0,totalparticalnumber;
  int i=0,j,k;
  
  //MPI
  int myid,numprocs;
  
  //RFQ_old
  int positionpre,positionnow,positionnext;
  double rfqtypepre,rfqtypenow,rfqtypenext;
  
  //Field
  int bendmorder=0,quadmorder=0,driftorder=0,solenorder=0,rfqorder=0,fieldorder=0;
private:
  template <typename T1,typename T2>
  void pushback(const T1 &ele, T2 &vecEle,const double Egamma,int &order,Element &);
};

#endif
