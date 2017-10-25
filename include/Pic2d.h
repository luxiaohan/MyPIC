#ifndef PIC2D_H
#define PIC2D_H
#include <vector>
#include <cmath>
#include <fftw3.h>
#include "LatticeMap.h"
#include "BeamMap.h"
using std::vector;
using v2= vector<vector<double> >;
class PIC2D
{
  public:
    /*
       2015-10-28@LiuZhicong
       */
    PIC2D();
    PIC2D(double parnumber,double gridnumx,double gridnumy,double time,double Egamma,double macroNum,int sli);
    ~PIC2D();
    int Initial();

    //PIC Solver
    int PIC2D_ICEFROG();
    int PIC2D_FIREFROG();
    int PIC2D_Weight();
    int PIC2D_FFT();
    int PIC2D_InternalField();
    int PIC2D_ExternalField(double time);
    int PIC2D_ExternalFieldPre(double time);
    int PIC2D_SpaceCharge() 
    {
      PIC2D_Weight(); 
      PIC2D_FFT(); 
      PIC2D_InternalField();
    }
    int PIC2D_PUSH(double time);




    //import
    int ReadParticle(const vector<double> &particleParameter);
    double ReadLattice(const char *p);
    int CalEmit();
    int CalMatchR();

    //get var
    vector<vector<double> > &xyz() 	 {return particle_xy;}
    vector<vector<double> > &phi() 	 {return phi_grid;}
    vector<vector<double> > &mf() 	 {return mField;}
    vector<vector<double> > &ef() 	 {return eField;}
    vector<vector<double> > &ex() 	 {return ex_grid;}
    vector<vector<double> > &ey() 	 {return ey_grid;}
    vector<vector<double> > &rho() 	 {return rho_grid;}
    vector<vector<double> > &tww() 	 {return ParaTransfer;}
  protected:
    //some number
    const double Huge		=1;
    const double Macrocharge	=1.602176565e-19*Huge;
    const double Macromass	=1.672621777e-27*Huge;	//unit:C,kg
    double Egamma=1.0001,Ebeta=sqrt(1.0-1/Egamma/Egamma);
    double phaseX=0,phaseY=0;

    double location=0;

    //Property of grid
    const int particlenumber=10000;
    double timeStep;//unit: second
    const double timeStepFix=1e-9;
    const int slice=10;
    const int numofgridx=256;
    const int numofgridy=256;
    double stepx=0.001,stepy=0.001,vt=stepx*stepy;//unit:m

    vector<int> countx	=vector<int>(particlenumber);
    vector<int> county	=vector<int>(particlenumber);
    vector<double> gridx	=vector<double>(numofgridx);
    vector<double> gridy	=vector<double>(numofgridy);

    vector<double> v1	=vector<double>(particlenumber); //temporary weight 
    vector<double> v2	=vector<double>(particlenumber); //the squence of vari:
    vector<double> v3	=vector<double>(particlenumber); //   :particlenumber !!
    vector<double> v4	=vector<double>(particlenumber);

    vector<vector<double> > particle_xy	=vector<vector<double>  >(particlenumber,   vector<double>(6));
    vector<vector<double>  >   rho_grid	=vector<vector<double>  >(numofgridx,   vector<double>(numofgridy));
    vector<vector<double>  >   phi_grid	=vector<vector<double>  >(numofgridx,   vector<double>(numofgridy));
    vector<vector<double>  >   ex_grid	=vector<vector<double>  >(numofgridx,   vector<double>(numofgridy));
    vector<vector<double>  >   ey_grid	=vector<vector<double>  >(numofgridx,   vector<double>(numofgridy));
    vector<vector<double>  >   eFieldInner	=vector<vector<double>  >(particlenumber,   vector<double>(3));
    vector<vector<double>  >   mFieldInner	=vector<vector<double>  >(particlenumber,   vector<double>(3));
    vector<vector<double>  >   eField	=vector<vector<double>  >(particlenumber,   vector<double>(3));
    vector<vector<double>  >   mField	=vector<vector<double>  >(particlenumber,   vector<double>(3));

    //FFT
    double *in,*out,*outf;
    fftw_plan g2f,f2g;

    //other
    BeamMap beam;
    vector<vector<double>  >   ParaTransfer;
    LatticeMap MyLattice;
  private:
    //PIC_Solver
    int PIC2D_PushVelocityFirstHalf();
    int PIC2D_PushVelocityLastHalf();
    int PIC2D_PushVelocity();
    int PIC2D_PushPosition();

    //PIC2D_Push
    int coun=0;
    int eleNumTemp=-1;
};
#endif
