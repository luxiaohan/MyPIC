#ifndef MAP_H
#define MAP_H

#ifdef WIN32
#include <direct.h>  
#include <io.h>  
#else
#include <stdarg.h>  
#include <sys/stat.h>  
#include <sys/types.h>
#endif  

#include <cstdlib>
#include <cstdio>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <string>
#include <vector>
#include <ctime>
#include <algorithm>
#include "AllElement.h"
#include "LatticeMap.h"
#include "BeamMap.h"
#include "MyFunc.h"
#include "Distribution.h"
#include <cmath>
#include <Eigen/Dense>
#include <mpich/mpi.h>
using std::vector;
using Eigen::MatrixXd;

class Mapmain
{
  public:
    virtual ~Mapmain();
    Mapmain();
    void latticeinput(char *p,int numberdivided=1);
    void prelattice(const char *,const char *,int);
    void caculate();
    void twissscan(int ,double,double,double,double);

    void twisscaculate();
    void twisstrack(double btmp, double atmp,char a,int enumber=0);
    void twisstrack(double bxtmp, double axtmp,double bytmp, double aytmp,int enumber=0);
    //MatrixXd twissmapx(const MatrixXd &R);//beta alpha gamma,in x direction test
    //MatrixXd twissmapy(const MatrixXd &R);//beta alpha gamma,in y direction test

    double getalphax();
    double getalphay();
    double getbetax();
    double getbetay();
    double getgammax();
    double getgammay();
    double GetPhaseAdvance_X(const MatrixXd &R);
    double GetPhaseAdvance_Y(const MatrixXd &R);
    double GetPhaseAdvance_Z(const MatrixXd &R);
    
    //TRACK
    void preparticle(const string &p);
    void postparticle();
    void particleread(const string &p);
    void ptrackonce();
    void ptrackall(int spacechargeflag = 1,char *p ="particles.dat");
    
    //SPACE CHARGE
    void spacechargepromote(double);
    void spacechargepre();
    double scfunx(double);
    double scfuny(double);
    double scfunz(double);
    double testfun(double);
    double Romberg(int);
    
    //GA TWISS
    void Gtwissscan(int period,double betax1,double betax2,double alphax1,double alphax2,double betay1,double betay2,double alphy1,double alphay2);
    void Gmutate(vector<double> &chromo);
    void Gselection();
    bool GSortbyfitness(const vector<double> &v1,const vector<double> &v2);
    void Gfitness(vector<double> &v1);
    void Gfitness2(vector<double> &v1);
    
    //GA LATTICE
    int Glatticescan(LatticeMap e1,LatticeMap e2,const vector<double> &v1,const vector<double> &v2,int FIT=2);
    int Glatticescan(LatticeMap e1,LatticeMap e2,const vector<double> &v1,const double &PX,const double &PY);
    int Glatticescan(LatticeMap e1,LatticeMap e2,const double &P1,const double &P2);
    void GLfitness(LatticeMap &v);
    void GLfitness1(LatticeMap &v);
    void GLfitness2(LatticeMap &v);
    void GLfitness3(LatticeMap &v);
    void GLfitness4(LatticeMap &v);
    
  protected:
    int BEAMLOSE=0;
    string line,pline;
    MatrixXd R1,R2;
    int DEBUG=0;
    
    //LATTICE
    LatticeMap lattice1;
    vector<double> location;
    int i=0,j,k,m;
    int _ELECOUNT;
    int Mcount=0,Mfail=0,particle_count=0,totalparticlenumber;
    double EGammaIni=1.0;
    double EBetaIni=sqrt(1.0-1.0/EGammaIni/EGammaIni);
    
    //BEAM
    BeamMap beam;
    double alpha,beta,gamma,TAlphax,TBetax,TGammax,TAlphay,TBetay,TGammay;
    double PhaseAdvance_X=0,PhaseNow_X=0,PhaseAdvance_Y=0,PhaseNow_Y=0;
    Vector3d Twiss_X,Twiss_Y,twiss;
    
    //SPACE CHARGE
    double X,Y,Z,gam,muspace[3];
    
    //MPI
    int myid,numprocs,namelen,NDivide;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    //GA Lattice
    int FitType;
    vector<double> twiss_In;
    vector<double> twiss_Out;
    double emitGLx,emitGLy;
    double PhaseExpect_X,PhaseExpect_Y;
    
    //GA TWISS
    int simpnum;
    double Param[200][4];
    int ParamFlag=0;
    double limit1[4],limit2[4];
    
    //Constant
    const double Dielectric_const=8.8541878176e-12;//F/m
};

#endif
