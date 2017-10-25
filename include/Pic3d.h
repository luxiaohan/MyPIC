#ifndef PIC3D_H
#define PIC3D_H
#include <vector>
#include <cmath>
#include <fftw3.h>
#include "LatticeMap.h"
#include "BeamMap.h"
using std::vector;
using v2d= vector<vector<double> > ;
using v1d= vector<double> ;
using v3d= vector<vector<vector<double> > > ;
class PIC3D
{
  public:
    /*
       2015-12-30@LiuZhicong
       */
    PIC3D();
    PIC3D(double parnumber,
	double gridnumx,
	double gridnumy,
	double gridnumz,
	double steptime,
	double macroNum,
	int frq);
    ~PIC3D();
    int Initial();

    //PIC Solver
    int PIC3D_ICEFROG();
    int PIC3D_FIREFROG();
    int PIC3D_Weight();
    int PIC3D_FFT();
    int PIC3D_InternalField();
    int PIC3D_ExternalField(double time);
    int PIC3D_ExternalFieldPre(double time);
    int PIC3D_SpaceCharge() 
    {
      PIC3D_Weight(); 
      PIC3D_FFT(); 
      PIC3D_InternalField();
    }
    int PIC3D_PUSH(double time,int post=1);




    //import
    int ReadParticle(const vector<double> &particleParameter,double gam,double dw,double length_z=0);
    int ReadParticle(const vector<double> &particleParameter,double gam);
    double ReadLattice(const char *p);	//return the length of the lattice
    int CalEmit();
    int CalMatchR();

    //get var
    v2d &xyz() 	 {return particle_xy;}
    BeamMap 	&be() 	 {return beam;}
    LatticeMap 	&Lat() 	 {return MyLattice;}
    vector<int> &pLS()	 {return particleLost;}
    int particleLiveNumber;
    int particleFixedNumber;
    int particleInMeshNumber;
    v3d &phi() 	 {return phi_grid;}
    v2d &mf() 	 {return mField;}
    v2d &ef() 	 {return eField;}
    v3d &ex() 	 {return ex_grid;}
    v3d &ey() 	 {return ey_grid;}
    v3d &rho() 	 {return rho_grid;}
    v2d &tww() 	 {return ParaTransfer;}
    v2d &tw2() 	 {return ParaTransfer2;}
  protected:
    //some number
    double Huge		=1;
    double BaseCharge	=1.602176565e-19;
    double BaseMass	=1.672621777e-27;
    double BaseEnergyInMeV	=938.27204;
    double Macrocharge		=BaseCharge*Huge;
    double Macromass		=BaseMass*Huge;	//unit:C,kg
    double Egamma=1.0001,Ebeta=sqrt(1.0-1/Egamma/Egamma);
    double dw=0;
    double frq=0;

    double location=0;

    //Property of grid
    int particlenumber=10000;
    double timeStep;//unit: second
    double timeStepFix=1e-9;
    int slice=10;
    int numofgridx=256;
    int numofgridy=256;
    int numofgridz=10;
    double stepx=0.001,stepy=0.001,stepz=0.001,vt=stepx*stepy*stepz;//unit:m

    vector<int> countx	=vector<int>(particlenumber);
    vector<int> county	=vector<int>(particlenumber); 
    vector<int> countz	=vector<int>(particlenumber);
    vector<double> gridx	=vector<double>(numofgridx);
    vector<double> gridy	=vector<double>(numofgridy);
    vector<double> gridz	=vector<double>(numofgridz);

    v2d vs 		=v2d(particlenumber, v1d(8));
    v2d particle_xy	=v2d(particlenumber, v1d(6));
    vector<int> particleLost	 =  vector<int>(particlenumber,0);
    vector<int> particleOutOfMesh=  vector<int>(particlenumber,0);
    vector<int> particleFixed	 =  vector<int>(particlenumber,0);

    v3d rho_grid	=v3d(numofgridx,   v2d(numofgridy,v1d(numofgridz)));
    v3d phi_grid	=v3d(numofgridx,   v2d(numofgridy,v1d(numofgridz)));
    v3d ex_grid		=v3d(numofgridx,   v2d(numofgridy,v1d(numofgridz)));
    v3d ey_grid		=v3d(numofgridx,   v2d(numofgridy,v1d(numofgridz)));
    v3d ez_grid		=v3d(numofgridx,   v2d(numofgridy,v1d(numofgridz)));

    v2d  eFieldInner	=v2d(particlenumber,   v1d(3));
    v2d  mFieldInner	=v2d(particlenumber,   v1d(3));
    v2d  eField	=v2d(particlenumber,  v1d(3));
    v2d  mField	=v2d(particlenumber,  v1d(3));

    //FFT
    double *out;
    double *in,*outf;
    double *in_xy,*out_xy;
    fftw_plan g2f_xy_many,f2g_xy_many;
    fftw_complex *in_z,*out_z,*outf_z;
    fftw_plan g2f_z_many,f2g_z_many;
    fftw_r2r_kind kind_xy_forward[2] = {FFTW_RODFT00,FFTW_RODFT00};
    fftw_r2r_kind kind_xy_backward[2]= {FFTW_RODFT00,FFTW_RODFT00};
    //many

    //other
    BeamMap beam;
    LatticeMap MyLattice;
    v2d  ParaTransfer;
    v2d  ParaTransfer2;
  private:
    //PIC_Solver
    int PIC3D_PushVelocityFirstHalf();
    int PIC3D_PushVelocityLastHalf();
    int PIC3D_PushVelocity();
    int PIC3D_PushPosition(int criParLost,int criPost);	//particle lost ON/OFF, post process ON/OFF

    //PIC3D_Push
    double timeline;
    int coun=0;
    int eleNumTemp=-1;

    //particle lost
    bool RFQLongtitudeLostFlag=0;
    bool RFQLongtitudeLostSwitch=1;

    //field synchronized
    vector<int> fieldPassed;
    double ScanPhase(int order);
    int scanFlag=0,tempOrder;
    int LostFlag=1;
    double ph;
    
    //getfield
    double getfield(v3d &,double,double,double);

    //benchmark
    void benchmark_flux();
    void benchmark_phi_theory_x();
    void benchmark_phi_theory_z();
    void benchmark_electri_theory_x();
    void benchmark_electri_theory_z();
};
#endif
