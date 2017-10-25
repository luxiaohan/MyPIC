#include "Pic3d.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "MyFunc.h"
#include "Distribution.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <algorithm>
#include "FFT.h"

using namespace std;
//GS UNIT
/*
   2015-10-28@LiuZhicong
   */

PIC3D::PIC3D()
{ 
  Initial();
}

PIC3D::PIC3D(double parnumber,double gridnumx,double gridnumy,double gridnumz,double time,double macroNum,int frqe):
  particlenumber(parnumber),
  numofgridx(gridnumx),
  numofgridy(gridnumy),
  numofgridz(gridnumz),
  timeStep(time),
  timeStepFix(time),
  Huge(macroNum),
  frq(frqe)
{ 
  Initial();
}

PIC3D::~PIC3D()
{
  fftw_free(in);
  fftw_free(out);
  fftw_free(outf);
  fftw_free(in_xy);
  fftw_free(in_z);
  fftw_free(out_xy);
  fftw_free(out_z);
  fftw_free(outf_z);
  fftw_destroy_plan(g2f_xy_many);
  fftw_destroy_plan(f2g_xy_many);
  fftw_destroy_plan(g2f_z_many);
  fftw_destroy_plan(f2g_z_many);
}

int PIC3D::Initial()
{
  cout<<"Initialing..."<<endl;
  in  	= (double*)             fftw_malloc(sizeof(double)*numofgridx*numofgridy*numofgridz);
  out 	= (double*)             fftw_malloc(sizeof(double)*numofgridx*numofgridy*numofgridz);
  outf 	= (double*)             fftw_malloc(sizeof(double)*numofgridx*numofgridy*numofgridz);
  in_xy = (double*)             fftw_malloc(sizeof(double)*numofgridx*numofgridy*numofgridz);
  out_xy= (double*)             fftw_malloc(sizeof(double)*numofgridx*numofgridy*numofgridz);
  in_z  = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)*numofgridx*numofgridy*numofgridz);
  out_z = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)*numofgridx*numofgridy*numofgridz);
  outf_z= (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)*numofgridx*numofgridy*numofgridz);
  for(int i=0;i<numofgridx;++i)
  {
    for(int j=0;j<numofgridy;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {
	in_z[i*numofgridy*numofgridz+j*numofgridz+k][1]=0.0;
      }
    }
  }
  //1
  fftw_plan_with_nthreads(omp_get_max_threads());
  int n_xy[]={numofgridx,numofgridy};
  g2f_xy_many=fftw_plan_many_r2r(2, n_xy, numofgridz,
      in_xy  ,n_xy,
      numofgridz, 1,
      out_xy,n_xy ,
      numofgridz, 1,
      kind_xy_forward, FFTW_MEASURE);
  //2
  fftw_plan_with_nthreads(omp_get_max_threads());
  int n_z[]={numofgridz};
  g2f_z_many=fftw_plan_many_dft(1, n_z, numofgridx*numofgridy,
      in_z , n_z ,
      1,numofgridz,
      out_z, n_z ,
      1,numofgridz,
      FFTW_FORWARD, FFTW_MEASURE);
  //3
  fftw_plan_with_nthreads(omp_get_max_threads());
  f2g_z_many=fftw_plan_many_dft(1, n_z, numofgridx*numofgridy,
      out_z  , n_z ,
      1,numofgridz,
      outf_z, n_z ,
      1,numofgridz,
      FFTW_BACKWARD, FFTW_MEASURE);
  //4
  fftw_plan_with_nthreads(omp_get_max_threads());
  f2g_xy_many=fftw_plan_many_r2r(2, n_xy, numofgridz,
      out_xy,n_xy,
      numofgridz, 1,
      in_xy,n_xy ,
      numofgridz, 1,
      kind_xy_backward, FFTW_MEASURE);


  /*
     fftw_plan_with_nthreads(omp_get_max_threads());
     f2g=fftw_plan_r2r_3d(
     numofgridx,numofgridy,numofgridz,
     in,	out,
     FFTW_RODFT00,
     FFTW_RODFT00,
     FFTW_R2HC,
     FFTW_MEASURE);
     fftw_plan_with_nthreads(omp_get_max_threads());
     f2g=fftw_plan_r2r_3d(
     numofgridx,numofgridy,numofgridz,
     out,	outf,
     FFTW_RODFT00,
     FFTW_RODFT00,
     FFTW_HC2R,
     FFTW_MEASURE);
     */
  beam.initial('p',Egamma,particle_xy.size());
  remove("lostparticle.dat");
  return 1;
}

int PIC3D::ReadParticle(const vector<double> &particleParameter,double gam,double dw,double length_z)
{
  if(particleParameter.size()!=6) 
  {
    cout<<"particle read in error!"<<endl;
    return -1;
  }
  Egamma=gam;
  Ebeta =sqrt(1.0-1.0/Egamma/Egamma);
  cout<<"Gamma   =  "<<Egamma<<"\nBeta    =  "<<Ebeta<<"\nSpeed   =  "<<(Ebeta*C_light)<<endl;
  KV2D_T(particleParameter,particle_xy);
  //KV2D_A(particleParameter,particle_xy);
  double gammatemp,betatemp;
  default_random_engine generator;
  uniform_real_distribution<double> distribution(-length_z,0);
  uniform_real_distribution<double> distribution2(1-dw,1+dw);
  for(int i=0;i<particlenumber;i++) 
  {
    particle_xy[i][1]*=Ebeta*C_light;
    particle_xy[i][3]*=Ebeta*C_light;
    particle_xy[i][4]=distribution(generator);
    //gammatemp=1+(Egamma-1)*distribution2(generator);
    //betatemp =sqrt(1.0-1.0/gammatemp/gammatemp);

    particle_xy[i][5]=distribution2(generator)*Ebeta*C_light;
  }


  particleLiveNumber=particlenumber; 
  if(particlenumber==1)
  {
    for(int i=0;i<6;i++)
    {
      if(i%2==0)
	particle_xy[0][i]=0.0/2;
      else
	particle_xy[0][i]=0;
    }
    particle_xy[0][3]=0;
    particle_xy[0][4]=-length_z/2;
    particle_xy[0][5]=Ebeta*C_light;
  }
  cout<<"Distribution   : KV"<<endl;
  cout<<"Particle number: "<<particle_xy.size()<<endl;
  ofstream pout("data/initialKV.dat");
  printvector(pout,particle_xy);
  pout.close();
  beam.caculate_emittance(particle_xy,particleLost,frq);
  return 0;
}

int PIC3D::ReadParticle(const vector<double> &particleParameter,double gam)
{
  if(particleParameter.size()<9)
  {
    cout<<"particleParameter.size()<9  in read particle"<<endl;
    return -1;
  }
  Egamma=gam;
  Ebeta =sqrt(1.0-1.0/Egamma/Egamma);
  double speed=Ebeta*C_light;
  cout<<"Gamma   =  "<<Egamma<<"\nBeta    =  "<<Ebeta<<"\nSpeed   =  "<<(Ebeta*C_light)<<endl;
  KV3D_T(particleParameter,particle_xy);
  ofstream pout("data/initialKV.dat");

  for(int i=0;i<particlenumber;i++) 
  {
    particle_xy[i][1]*=speed;
    particle_xy[i][3]*=speed;
    particle_xy[i][5]*=speed;
    particle_xy[i][5]+=speed;
  }
    printvector(pout,particle_xy);
  particleLiveNumber=particlenumber; 
  if(particlenumber==1)
  {
    for(int i=0;i<6;i++)
    {
      if(i%2==0)
	particle_xy[0][i]=0.0/2;
      else
	particle_xy[0][i]=0;
    }
    particle_xy[0][3]=0;
    particle_xy[0][5]=Ebeta*C_light;
  }
  cout<<"Distribution   : KV"<<endl;
  cout<<"Particle number: "<<particle_xy.size()<<endl;
  //ofstream pout("data/initialKV.dat");
  //printvector(pout,particle_xy);
  pout.close();
  beam.caculate_emittance(particle_xy,particleLost,frq);
  return 0;
}

double PIC3D::ReadLattice(const char *p)
{
  MyLattice.read(p,Egamma,1);
  return MyLattice.latticeLength;
}

int PIC3D::PIC3D_ICEFROG()
{
  //before cycle, we push x in halfstep.		FOR ICEFROG!!!
  timeStep=-timeStep;
  PIC3D_ExternalField(0);
  PIC3D_PushVelocityLastHalf(); 
  timeStep=timeStepFix;
}

int PIC3D::PIC3D_FIREFROG()
{
  //afte cycle, we push dx in halfstep			FOR ICEFROG!!!
  PIC3D_PushVelocityFirstHalf();
}


int PIC3D::PIC3D_Weight()
{
  if(particlenumber==0) return -1;
#pragma omp single
  {
    for(int i=0;i<numofgridx;++i)
      for(int j=0;j<numofgridy;++j)
	for(int k=0;k<numofgridz;++k)
	{
	  rho_grid[i][j][k]=0;
	}
    //caculate the step
    vector<double> total(3),ave(6),sigma(3); 
    for(int i=0;i<particle_xy.size();i++)
    {
      for(int j=0;j<3;++j)
      {
	total[j]+=particle_xy[i][2*j];
      }
    }
    for(int i=0;i<3;++i)
    {
      ave[i]=total[i]/particle_xy.size();
    }

    for(int i=0;i<particle_xy.size();i++)
    {
      for(int j=0;j<3;++j)
      {
	sigma[j]		+=(particle_xy[i][2*j]-ave[j])		*	(particle_xy[i][2*j]-ave[j]);
      }
    }
    for(int i=0;i<3;++i)
    {
      sigma[i]/=particle_xy.size();
      //sigma[i]=beam.get
      sigma[i]=sqrt(sigma[i]);
    }
    beam.caculate_emittance(particle_xy,particleLost,frq);
    sigma[0]=beam.getsigmax();
    sigma[1]=beam.getsigmay();
    sigma[2]=beam.getsigmaz();
    ave[0]=beam.particle_ave[0];
    ave[1]=beam.particle_ave[2];
    ave[2]=beam.particle_ave[4];

    if(particlenumber==1)
    {
      sigma[0]=0.001; 
      sigma[1]=0.001;
      sigma[2]=0.001;
    }
    stepx=sqrt(sigma[0])*2*3.5/(numofgridx);
    stepy=sqrt(sigma[1])*2*3.5/(numofgridy);
    if(MyLattice.ELE[MyLattice.getEleNum(ave[2])].name=="rfq")
    {
      stepz=  beam.particle_ave[5]/ MyLattice.ELE[MyLattice.getEleNum(ave[2])].p[4]/numofgridz;
    }
    else
    {
      stepz=sqrt(sigma[2])*2*3.5/(numofgridz);
    }
    
    if(particlenumber==1)
    {
      stepx=1; 
      stepy=1;
      stepz=1;
    }
    
    /*
       ave[0]=  particle_xy[0][0]-stepx/2-numofgridx/4.0*stepx; 
       ave[1]=  particle_xy[0][2]-stepy/2-numofgridy/4.0*stepy;
       ave[2]=  particle_xy[0][4]-stepz/2-numofgridz/4.0*stepz;
       cout<<ave[2]<<"  "<<particle_xy[0][4]<<endl;
       */
    //cout<<"step  : "<<stepx<<"  "<<stepy<<"  "<<stepz<<endl;
    //stepz=0.001; 

    for(int i=0;i<gridx.size();i++)
    {
      //gridx[i]=ave[0]+stepx*(i+0.5-numofgridx/2);
      gridx[i]=ave[0]+stepx*(i-numofgridx/2);
    }
    for(int i=0;i<gridy.size();++i)
    {
      //gridy[i]=ave[1]+stepy*(i+0.5-numofgridy/2);
      gridy[i]=ave[1]+stepy*(i-numofgridy/2);
    }
    for(int i=0;i<gridz.size();++i)
    {
      //gridz[i]=ave[2]+stepz*(i+0.5-numofgridz/2);
      gridz[i]=ave[2]+stepz*(i-numofgridz/2);
    }
    //get which gird the particle locate
    vt=stepx*stepy*stepz;
  }

  int idx,idy,idz;
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
    if(particleLost[i]) continue;
    countx[i]=floor((particle_xy[i][0]-gridx[0])/stepx);
    county[i]=floor((particle_xy[i][2]-gridy[0])/stepy);
    countz[i]=floor((particle_xy[i][4]-gridz[0])/stepz);
    //cout<<countx[i]<<county[i]<<" "<<countz[i]<<endl;
    if(   (countx[i]>=numofgridx-1)
	||(county[i]>=numofgridy-1)
	||(countz[i]>=numofgridz)
	||countx[i]<0||county[i]<0||countz[i]<0)
    {
      particleOutOfMesh[i]=1;
      for(int j=0;j<8;++j)
      {
	vs[i][j]=  0;
      }
    }
    else
    {
      particleOutOfMesh[i]=0;
      /*
	 if(countx[i]>=numofgridx) countx[i]=numofgridx-1;
	 if(county[i]>=numofgridy) county[i]=numofgridy-1;  
	 if(countz[i]>=numofgridz) countz[i]=numofgridz-1;
	 */
      for(int j=0;j<8;++j)
      {
	/*  cout<<countx[i]+ ((j%2==0)? 1 : 0 )<<endl;
	    cout<<county[i]+ ((j%4<2)? 1 : 0 )<<endl;
	    cout<<countz[i]+ ((j<0)? 1 : 0 )<<endl;*/
	if(countz[i]==numofgridz-1)
	{
	  idz=(j  <=3) 	? -countz[i] : 0 ;
	}
	else
	{
	  idz=(j  <=3) 	? 1 : 0;
	}

	vs[i][j]=  	
	  abs(particle_xy[i][0]-gridx[countx[i]+ ((j%2==0)	? 1 : 0) ]) * 
	  abs(particle_xy[i][2]-gridy[county[i]+ ((j%4< 2) 	? 1 : 0) ]) *
	  abs(particle_xy[i][4]-gridz[countz[i]+ idz ])
	  /vt;

	//cout<<"vs "<<i<<" "<<j<<" "<<vs[i][j]<<endl;
      }
    }
    /*
       vs[i][3]=   (particle_xy[i][0]-gridx[countx[i]])  *   (particle_xy[i][2]-gridy[county[i]])   /vt;
       vs[i][2]=abs(particle_xy[i][0]-gridx[countx[i]+1])*   (particle_xy[i][2]-gridy[county[i]])   /vt;
       vs[i][1]=   (particle_xy[i][0]-gridx[countx[i]])  *abs(particle_xy[i][2]-gridy[county[i]+1]) /vt;
       vs[i][0]=abs(particle_xy[i][0]-gridx[countx[i]+1])*abs(particle_xy[i][2]-gridy[county[i]+1]) /vt;
       */
  }
#pragma omp single
  {
    particleInMeshNumber=0;
    for(int i=0;i<particlenumber;++i)
    {
      if(particleLost[i] or particleOutOfMesh[i]) continue;
      ++particleInMeshNumber;
      rho_grid[countx[i]  ][county[i]  ][countz[i]]  		+=	vs[i][0]*Macrocharge;
      rho_grid[countx[i]+1][county[i]  ][countz[i]]		+=	vs[i][1]*Macrocharge;
      rho_grid[countx[i]  ][county[i]+1][countz[i]]		+=	vs[i][2]*Macrocharge;
      rho_grid[countx[i]+1][county[i]+1][countz[i]]		+=	vs[i][3]*Macrocharge;
      if(countz[i]==numofgridz-1)
      {
	rho_grid[countx[i]  ][county[i]  ][0]		+=	vs[i][4]*Macrocharge;
	rho_grid[countx[i]+1][county[i]  ][0]		+=	vs[i][5]*Macrocharge;
	rho_grid[countx[i]  ][county[i]+1][0]		+=	vs[i][6]*Macrocharge;
	rho_grid[countx[i]+1][county[i]+1][0]		+=	vs[i][7]*Macrocharge;
      }
      else
      {
	rho_grid[countx[i]  ][county[i]  ][countz[i]+1]		+=	vs[i][4]*Macrocharge;
	rho_grid[countx[i]+1][county[i]  ][countz[i]+1]		+=	vs[i][5]*Macrocharge;
	rho_grid[countx[i]  ][county[i]+1][countz[i]+1]		+=	vs[i][6]*Macrocharge;
	rho_grid[countx[i]+1][county[i]+1][countz[i]+1]		+=	vs[i][7]*Macrocharge;
      }
    }
//    
//       double sum =0;
//       for(int i=0;i<numofgridx;++i)
//       {
//       for(int j=0;j<numofgridy;++j)
//       {
//       for(int k=0;k<numofgridz;++k)
//       {
//       sum+=rho_grid[i][j][k];
//       }
//       }
//       }
//       cout<<"total charge  "<<sum<<"  number= "<<sum/Macrocharge<<endl;
  }

  return 1;
}

int PIC3D::PIC3D_FFT()
{
  //now we start fft
  /*FFT
    2015-10-29@LiuZhicong
    */
  int type=1;
  if(type==1)//DST in x and y ,DFT in z
  {
    //forward Fourier
#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  in_xy[i*numofgridy*numofgridz+j*numofgridz+k]=rho_grid[i][j][k]/stepx/stepy/stepz;
	}
      }
    }

#pragma omp barrier
#pragma omp single 
    fftw_execute(g2f_xy_many);
#pragma omp barrier
    int coun365;
#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  coun365=i*numofgridy*numofgridz+j*numofgridz+k;
	  in_z[coun365][0]=out_xy[coun365];
	}
      }
    }

#pragma omp barrier
#pragma omp single 
    fftw_execute(g2f_z_many);
#pragma omp barrier
    /*
       for(int k=0;k<numofgridz;++k)
       {
       for(int i=0;i<numofgridx;++i)
       {
       for(int j=0;j<numofgridy;++j)
       {
       in_xy[i*numofgridy+j]=in[i*numofgridy*numofgridz+j*numofgridz+k];
       }
       }
       fftw_execute(g2f_xy);
       for(int i=0;i<numofgridx;++i)
       {
       for(int j=0;j<numofgridy;++j)
       {
       in[i*numofgridy*numofgridz+j*numofgridz+k]=out_xy[i*numofgridy+j];
       }
       }
       }

       for(int i=0;i<numofgridx;++i)
       {
       for(int j=0;j<numofgridy;++j)
       {
       for(int k=0;k<numofgridz;++k)
       {
       in_z[i*numofgridy+j][0]=in[i*numofgridy*numofgridz+j*numofgridz+k];
       in_z[i*numofgridy+j][1]=0;
       }

       fftw_execute(g2f_z);
       for(int k=0;k<numofgridz;++k)
       {
       out[i*numofgridy*numofgridz+j*numofgridz+k][0]=out_z[i*numofgridy+j][0];
       out[i*numofgridy*numofgridz+j*numofgridz+k][1]=out_z[i*numofgridy+j][1];
       }
       } 
       }
       */
    /* if want to use this part, change the type of "out" from fftw_complex to double
#pragma omp barrier
#pragma omp single 
fftw_execute(g2f);
#pragma omp barrier
*/
    //devide by K
    double K_rho2phi=0,knx=0,kny=0,knz=0;
    int coun445;
#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  if(1)
	  {
	    knx=M_PI/2*(i+1)/(numofgridx+1)/stepx;
	    kny=M_PI/2*(j+1)/(numofgridy+1)/stepy;
	    knz=M_PI  *(k  )/(numofgridz)/stepz;
	  }
	  else if(0)// back up
	  {
	    knx=M_PI/2*(i+1)/(numofgridx+1)/stepx;
	    kny=M_PI/2*(j+1)/(numofgridy+1)/stepy;
	    knz=M_PI/1*(k  )/(numofgridz  )/stepz;
	  }
	  else if(1)
	  {
	    if(i==0&&j==0&&k==0) continue;
	    knx=M_PI/2*(i)/(numofgridx)/stepx;
	    kny=M_PI/2*(j)/(numofgridy)/stepy;
	    knz=M_PI  *(k)/(numofgridz)/stepz;
	  }
	  //knz=M_PI*(k+1)/(numofgridz+1)/stepz;
	  /*
	     if(k<=numofgridz/2)
	     {
	     knz=M_PI*(k+1)/(numofgridz+1)/stepz;
	     }
	     else if(k>=(numofgridz+1)/2-1)
	     {
	     int temp= numofgridz-k;
	     knz=M_PI*(temp+1)/(numofgridz+1)/stepz;
	     }
	     else
	     {
	     cout<<"error in fft z"<<endl;
	     }
	     */
	  //r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
	  //Sk=	(2.0 * sin( pi*(i-1)/ (2 * Grid_nx)  )/ Grid_dx)**2 + &
	  //  	(2.0 * sin( pi*(j-1)/ (2 * Grid_ny)  )/ Grid_dy)**2 + &
	  //	(2.0 * sin( pi*(k-1)/      Grid_nz   )/ Grid_dz)**2	
	  if(1)
	  {
	    K_rho2phi=	pow(2*sin(knx*stepx)/stepx,2)+
	                pow(2*sin(kny*stepy)/stepy,2)+
	                pow(2*sin(knz*stepz)/stepz,2);
	  }
	  else if(0)
	  {
	    K_rho2phi=4*(knx*knx+kny*kny+knz*knz);
	  }
	  //this one is very like the one above ,and according to my caclulate ,this one is the analysic solution.but the charge is not normalize to 1 but to 1.00004. Confusing!
	  coun445=i*numofgridy*numofgridz+j*numofgridz+k;
	  out_z[coun445][0]/=K_rho2phi*Epsilon;
	  out_z[coun445][1]/=K_rho2phi*Epsilon;
	}
      }
    }
#pragma omp single
    if(0)
    {
      out_z[0][0]=0;
      out_z[0][1]=0;
    }    
    //backward Fourier
    /*
#pragma omp barrier
#pragma omp single 
fftw_execute(f2g);
#pragma omp barrier
*/
#pragma omp barrier
#pragma omp single 
    fftw_execute(f2g_z_many);
#pragma omp barrier

#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  coun365=i*numofgridy*numofgridz+j*numofgridz+k;
	  out_xy[coun365]=outf_z[coun365][0];
	}
      }
    }

#pragma omp barrier
#pragma omp single 
    fftw_execute(f2g_xy_many);
#pragma omp barrier


#pragma omp for
    for(int i=0;i<numofgridx;i++)
    {
      for(int j=0;j<numofgridy;j++)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  phi_grid[i][j][k]=in_xy[i*numofgridy*numofgridz+j*numofgridz+k]/(4*(numofgridx+1)*(numofgridy+1)*(numofgridz));
	  //The fftw library is not normalized
	  //FFTW_RODFT00 computes an RODFT00 transform, i.e. a DST-I. (Logical N=2*(n+1), inverse is FFTW_RODFT00.)
	}
      }
    }
  }
  else if(type==2) // in stand
  {
    //forward Fourier
#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  in_xy[i*numofgridy*numofgridz+j*numofgridz+k]=rho_grid[i][j][k]/stepx/stepy/stepz;
	}
      }
    }

#pragma omp barrier
#pragma omp single 
    fftw_execute(g2f_xy_many);
#pragma omp barrier
    int coun365;
#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  coun365=i*numofgridy*numofgridz+j*numofgridz+k;
	  in_z[coun365][0]=out_xy[coun365];
	}
      }
    }

#pragma omp barrier
#pragma omp single 
    fftw_execute(g2f_z_many);
#pragma omp barrier
    /*
       for(int k=0;k<numofgridz;++k)
       {
       for(int i=0;i<numofgridx;++i)
       {
       for(int j=0;j<numofgridy;++j)
       {
       in_xy[i*numofgridy+j]=in[i*numofgridy*numofgridz+j*numofgridz+k];
       }
       }
       fftw_execute(g2f_xy);
       for(int i=0;i<numofgridx;++i)
       {
       for(int j=0;j<numofgridy;++j)
       {
       in[i*numofgridy*numofgridz+j*numofgridz+k]=out_xy[i*numofgridy+j];
       }
       }
       }

       for(int i=0;i<numofgridx;++i)
       {
       for(int j=0;j<numofgridy;++j)
       {
       for(int k=0;k<numofgridz;++k)
       {
       in_z[i*numofgridy+j][0]=in[i*numofgridy*numofgridz+j*numofgridz+k];
       in_z[i*numofgridy+j][1]=0;
       }

       fftw_execute(g2f_z);
       for(int k=0;k<numofgridz;++k)
       {
       out[i*numofgridy*numofgridz+j*numofgridz+k][0]=out_z[i*numofgridy+j][0];
       out[i*numofgridy*numofgridz+j*numofgridz+k][1]=out_z[i*numofgridy+j][1];
       }
       } 
       }
       */
    /* if want to use this part, change the type of "out" from fftw_complex to double
#pragma omp barrier
#pragma omp single 
fftw_execute(g2f);
#pragma omp barrier
*/
    //devide by K
    double K_rho2phi=0,knx=0,kny=0,knz=0;
    int coun445;
#pragma omp for
    for(int i=0;i<numofgridx;++i)
    {
      for(int j=0;j<numofgridy;++j)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  if(0)
	  {
	    //if(i==0&&j==0&&k==0) continue;
	    knx=M_PI/2*(i+1)/(numofgridx+1)/stepx;
	    kny=M_PI/2*(j+1)/(numofgridy+1)/stepy;
	    knz=M_PI/2*(k+1)/(numofgridz+1)/stepz;
	  }
	  else if(1)
	  {
	    if(i==0&&j==0&&k==0) continue;
	    knx=M_PI/2*(i)/(numofgridx)/stepx;
	    kny=M_PI/2*(j)/(numofgridy)/stepy;
	    knz=M_PI/2*(k)/(numofgridz)/stepz;
	  }
	  /*
	     if(k<=numofgridz/2)
	     {
	     knz=M_PI*(k+1)/(numofgridz+1)/stepz;
	     }
	     else if(k>=(numofgridz+1)/2-1)
	     {
	     int temp= numofgridz-k;
	     knz=M_PI*(temp+1)/(numofgridz+1)/stepz;
	     }
	     else
	     {
	     cout<<"error in fft z"<<endl;
	     }
	     */
	  //r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
	  //Sk=	(2.0 * sin( pi*(i-1)/ (2 * Grid_nx)  )/ Grid_dx)**2 + &
	  //  	(2.0 * sin( pi*(j-1)/ (2 * Grid_ny)  )/ Grid_dy)**2 + &
	  //		(2.0 * sin( pi*(k-1)/      Grid_nz   )/ Grid_dz)**2	
	  K_rho2phi=	pow(2*sin(knx*stepx)/stepx,2)+
	    pow(2*sin(kny*stepy)/stepy,2)+
	    pow(2*sin(knz*stepz)/stepz,2);
	  //K_rho2phi=4*(knx*knx+kny*kny+knz*knz);
	  //this one is very like the one above ,and according to my caclulate ,this one is the analysic solution.but the charge is not normalize to 1 but to 1.00004. Confusing!
	  coun445=i*numofgridy*numofgridz+j*numofgridz+k;
	  out_z[coun445][0]/=K_rho2phi*Epsilon;
	  out_z[coun445][1]/=K_rho2phi*Epsilon;
	}
      }
    }
    /*
       out_z[0][0]=0;
       out_z[0][1]=0;
       */

    //backward Fourier
    /*
#pragma omp barrier
#pragma omp single 
fftw_execute(f2g);
#pragma omp barrier
*/



#pragma omp for
    for(int i=0;i<numofgridx;i++)
    {
      for(int j=0;j<numofgridy;j++)
      {
	for(int k=0;k<numofgridz;++k)
	{
	  phi_grid[i][j][k]=outf[i*numofgridy*numofgridz+j*numofgridz+k]/(4*(numofgridx+1)*(numofgridy+1)*(numofgridz));
	  //The fftw library is not normalizedi
	  //FFTW_RODFT00 computes an RODFT00 transform, i.e. a DST-I. (Logical N=2*(n+1), inverse is FFTW_RODFT00.)
	  //FFTW_HC2R computes the reverse of FFTW_R2HC, above. (Logical N=n, inverse is FFTW_R2HC.)
	}
      }
    }
  }

  //Diferentiate to get the e field
#pragma omp for
  for(int i=1;i<numofgridx-1;++i)
  {
    for(int j=1;j<numofgridy-1;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {
	ex_grid[i][j][k]=-(phi_grid[i+1][j][k]-phi_grid[i-1][j][k])/stepx/2;//unit:V/m;
	ey_grid[i][j][k]=-(phi_grid[i][j+1][k]-phi_grid[i][j-1][k])/stepy/2;//divide gird step;
	if(k==0)
	  ez_grid[i][j][k]=-(phi_grid[i][j][k+1]-phi_grid[i][j][numofgridz-1])/stepz/2;
	else if(k==numofgridz-1)
	  ez_grid[i][j][k]=-(phi_grid[i][j][0]-phi_grid[i][j][k-1])/stepz/2;
	else
	  ez_grid[i][j][k]=-(phi_grid[i][j][k+1]-phi_grid[i][j][k-1])/stepz/2;
	//should I caculate the eField at half grid?
      }
    }
  }
  /*
     ofstream exout("ex_grid.dat");
     ofstream eyout("ey_grid.dat");
     for(int j=0;j<numofgridy;j++)
     {
     for(int i=0;i<numofgridx;i++)
     {
     exout<<ex_grid[i][j]<<"   ";
     eyout<<ey_grid[i][j]<<"   ";
     }
     exout<<endl;
     eyout<<endl;
     }
     exout.close();
     eyout.close();
     */
  return 0;
}

int PIC3D::PIC3D_InternalField()
{
  double e1=0,e2=0,e3=0,b1=0,b2=0,b3=0;
  double gamma,beta;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    if(particleLost[i] || particleOutOfMesh[i]) continue;
    if(!(countx[i]>=numofgridx||county[i]>=numofgridy||countz[i]>=numofgridz
	  ||countx[i]<0||county[i]<0||countz[i]<0))
    {
      e1  =   ex_grid[countx[i]][county[i]][countz[i]]      	* vs[i][0]
	+     ex_grid[countx[i]+1][county[i]][countz[i]]    	* vs[i][1]
	+     ex_grid[countx[i]][county[i]+1][countz[i]]    	* vs[i][2]
	+     ex_grid[countx[i]+1][county[i]+1][countz[i]]  	* vs[i][3]
	+     ex_grid[countx[i]][county[i]][countz[i]+1]      	* vs[i][4]
	+     ex_grid[countx[i]+1][county[i]][countz[i]+1]    	* vs[i][5]
	+     ex_grid[countx[i]][county[i]+1][countz[i]+1]    	* vs[i][6]
	+     ex_grid[countx[i]+1][county[i]+1][countz[i]+1]  	* vs[i][7];
      e2  =   ey_grid[countx[i]][county[i]][countz[i]]      	* vs[i][0]
	+     ey_grid[countx[i]+1][county[i]][countz[i]]    	* vs[i][1]
	+     ey_grid[countx[i]][county[i]+1][countz[i]]    	* vs[i][2]
	+     ey_grid[countx[i]+1][county[i]+1][countz[i]]  	* vs[i][3]
	+     ey_grid[countx[i]][county[i]][countz[i]+1]      	* vs[i][4]
	+     ey_grid[countx[i]+1][county[i]][countz[i]+1]    	* vs[i][5]
	+     ey_grid[countx[i]][county[i]+1][countz[i]+1]    	* vs[i][6]
	+     ey_grid[countx[i]+1][county[i]+1][countz[i]+1]  	* vs[i][7];
      e3  =   ez_grid[countx[i]][county[i]][countz[i]]      	* vs[i][0]
	+     ez_grid[countx[i]+1][county[i]][countz[i]]    	* vs[i][1]
	+     ez_grid[countx[i]][county[i]+1][countz[i]]    	* vs[i][2]
	+     ez_grid[countx[i]+1][county[i]+1][countz[i]]  	* vs[i][3]
	+     ez_grid[countx[i]][county[i]][countz[i]+1]      	* vs[i][4]
	+     ez_grid[countx[i]+1][county[i]][countz[i]+1]    	* vs[i][5]
	+     ez_grid[countx[i]][county[i]+1][countz[i]+1]    	* vs[i][6]
	+     ez_grid[countx[i]+1][county[i]+1][countz[i]+1]  	* vs[i][7];
      if(particlenumber==1)
      {
	//benchmark_flux();
	benchmark_phi_theory_x();
	benchmark_phi_theory_z();
	benchmark_electri_theory_x();
	benchmark_electri_theory_z();
	ofstream ff("fieldtempz");
	for(int j=0;j<gridz.size();++j)
	{
	  //if(pow(gridx[j]-particle_xy[i][0],2)<1e-7) continue;
	  ff<<j
	    <<"  "
	    //<<sqrt(pow((ex_grid[j][county[i]][countz[i]]+ex_grid[j][county[i]][countz[i]+1]+ex_grid[j][county[i]+1][countz[i]+1]+ex_grid[j][county[i]+1][countz[i]])/2,2))
	    <<ez_grid[countx[i]][county[i]][j]
	    <<"  "
	    <<Macrocharge/Epsilon/(4*M_PI*(pow(gridz[j]-particle_xy[i][4],2)))
	    <<endl;
	}
	ff.close();
	ff.open("fieldtempx");
	int numt=100;
	double lengthy=gridx.back()-gridx[0];
	for(int j=1;j<numt;++j)
	{
	  double x=gridx[0]+j*lengthy/numt,y=particle_xy[0][2],z=particle_xy[0][4];
	  //if(pow(gridx[j]-particle_xy[i][0],2)<1e-7) continue;
	  ff<<j
	    <<"  "
	    //<<sqrt(pow((ex_grid[j][county[i]][countz[i]]+ex_grid[j][county[i]][countz[i]+1]+ex_grid[j][county[i]+1][countz[i]+1]+ex_grid[j][county[i]+1][countz[i]])/2,2))
	    <<getfield(ex_grid,x,y,z)
	    <<"  "
	    <<Macrocharge/Epsilon/(4*M_PI*(pow(x-particle_xy[i][0],2)))
	    <<endl;
	}
	ff.close();
	ff.open("fieldphi");
	for(int j=0;j<gridx.size();++j)
	{
	  double x=gridx[j],y=particle_xy[0][2],z=particle_xy[0][4];
	  //if(pow(gridx[j]-particle_xy[i][0],2)<1e-7) continue;
	  ff<<setw(13)
	    <<gridx[j]
	    <<"  "
	    <<setw(13)
	    <<phi_grid[j][32][16]
	    <<"  "
	    <<setw(13)
	    <<sqrt(pow(getfield(phi_grid,x,y,z),2))
	    <<"  "
	    <<setw(13)
	    <<Macrocharge/Epsilon/(4*M_PI*(abs(gridx[j]-particle_xy[i][0])))
	    <<endl;
	}
	ff.close();
	cout<<"e1  "<<ex_grid[countx[i]][county[i]][countz[i]]<<"  "<<ey_grid[countx[i]][county[i]][countz[i]]<<"  "<<ey_grid[countx[i]][county[i]][countz[i]];
	cout<<" ecal "<<9.0e9*Macrocharge/(stepx*stepx+stepy*stepy+stepz*stepz)*4/sqrt(3)<<endl;
	cout<<countx[i]<<"  "<<county[i]<<"  "<<countz[i]<<endl;
	cout<<"e2  "<<ex_grid[countx[i]-1][county[i]-1][countz[i]-1]<<"  "<<ey_grid[countx[i]-1][county[i]-1][countz[i]-1]<<"  "<<ey_grid[countx[i]-1][county[i]-1][countz[i]-1];
	cout<<" ecal "<<9.0e9*Macrocharge/(stepx*stepx+stepy*stepy+stepz*stepz)*4/9/sqrt(3)<<endl;
      }

      e1=e1/C_light;
      e2=e2/C_light;
      e3=e3/C_light;
      beta=particle_xy[i][4]/C_light;
      gamma=1.0/sqrt(1-beta*beta);
      //Lorentz Transform
      eFieldInner[i][0]=C_light	*gamma*(e1+beta*b2);
      eFieldInner[i][1]=C_light	*gamma*(e2-beta*b1);
      eFieldInner[i][2]=C_light	*e3;
      mFieldInner[i][0]=	gamma*(b1-beta*e2);
      mFieldInner[i][1]=	gamma*(b2+beta*e1);
      mFieldInner[i][2]=	b3;

    }
    else
    {
      eFieldInner[i][0]=0;
      eFieldInner[i][1]=0;
      eFieldInner[i][2]=0;
      mFieldInner[i][0]=0;
      mFieldInner[i][1]=0;
      mFieldInner[i][2]=0;
    }	
  }
  //#pragma omp single
  //cout<<eFieldInner[2][0]<<"  "<<eFieldInner[2][1]<<"  "<<eFieldInner[2][2]<<endl;
}

int PIC3D::PIC3D_ExternalFieldPre(double time)
{
  vector<double> exterFieldTemp(6);
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
    if(particleLost[i]) continue;
    exterFieldTemp=MyLattice.getFieldPre(particle_xy[i],time);		//Ex,Ey,Ez,Bx,By,Bz
    //	Bz is positive when point at direction s
    //cout<<exterFieldTemp[5]<<" x "<<endl;
    for(int j=0;j<3;j++)
    {
      eField[i][j]  =  (eFieldInner[i][j]+exterFieldTemp[j])	*Macrocharge/Macromass/Egamma;//the field had time Q and divided by M already!!
      mField[i][j]  =  (mFieldInner[i][j]+exterFieldTemp[j+3])	*Macrocharge/Macromass/Egamma;//the field had time Q and divided by M already!!
    }
    mField[i][2]=-mField[i][2];//for  (Bz is positive when point at direction s)
  }
}

int PIC3D::PIC3D_ExternalField(double time)
{
  double factorTemp=0;
  vector<double> exterFieldTemp(6);
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
    if(particleLost[i]) continue;
    exterFieldTemp=MyLattice.getField(particle_xy[i],time);		//Ex,Ey,Ez,Bx,By,Bz
    //	Bz is positive when point at direction s
    factorTemp=Macrocharge/Macromass*sqrt(1-particle_xy[i][5]*particle_xy[i][5]/C_light/C_light);
    //cout<<"factor= "<<factorTemp<<endl;
    for(int j=0;j<3;j++)
    {
      eField[i][j]  =  (eFieldInner[i][j]+exterFieldTemp[j])	* factorTemp;//the field had time Q and divided by M and gamma already!!
      mField[i][j]  =  (mFieldInner[i][j]+exterFieldTemp[j+3])	* factorTemp;//the field had time Q and divided by M and gamma already!!
    }
    /*
       if(i==2) {
       cout<<"Inner  "<<eFieldInner[i][0]<<", Exter  "<<exterFieldTemp[0]<<" cal "<<9.0e9*Macrocharge/stepx/stepx<<endl;
       }
       */
    /*
       mField[i][0]=-mField[i][0];//for  (Bz is positive when point at direction s)
       mField[i][1]=-mField[i][1];//for  (Bz is positive when point at direction s)
       mField[i][2]=-mField[i][2];//for  (Bz is positive when point at direction s)
       */
  }
}

int PIC3D::PIC3D_PUSH(double time,int post)
{
  //Adjust timeStep according to Element Length
#pragma omp single
  {
    timeline=time;
    if(0&&MyLattice.getEleNum(particle_xy[0][4])!=eleNumTemp)
    {
      eleNumTemp=MyLattice.getEleNum(particle_xy[0][4]);
      int Slice=MyLattice.ELE[eleNumTemp].length/(timeStepFix*particle_xy[0][5]);
      timeStep=MyLattice.ELE[eleNumTemp].length/(Slice+1)/particle_xy[0][5];
      //timeStep=timeStepFix;
      cout<<Slice+1<<"  Sle  "<<timeStep<<"  "<<eleNumTemp<<" ele "<<"   ";
      coun=0;
    }
    else
    {
      //coun++;
    }
  }
#pragma omp barrier
  /*
     int zhengze=0;
     if(zhengze==1)
     {
     for(int i=0;i<particlenumber;++i)
     {
     particle_xy[i][1]-=-mField.at(i).back()*particle_xy[i][2]/(Egamma)/2;
     particle_xy[i][3]-= mField.at(i).back()*particle_xy[i][0]/(Egamma)/2;
  //cout<<"Macro    "<<Ax<<"  "<<Macrocharge/(Egamma*Macromass)<<endl;
  }
  }*/

  PIC3D_ExternalField(time);
  //PIC3D_PushVelocityFirstHalf();
#pragma omp single
  {
    coun++;
    if(post==1)
      CalEmit();
    else if(post==0)
      beam.caculate_emittance(particle_xy,particleLost,frq);
  }
  /*
     if(coun==0)
     {
     PIC3D_ExternalField();
     }*/
  PIC3D_PushVelocity();
  //PIC3D_PushVelocityLastHalf();
  PIC3D_PushPosition(LostFlag,post);
  /*
     if(zhengze==1)
     {
     for(int i=0;i<particlenumber;i++)
     {
     particle_xy[i][1]+=-mField.at(i).back()*particle_xy[i][2]/(Egamma)/2;
     particle_xy[i][3]+= mField.at(i).back()*particle_xy[i][0]/(Egamma)/2;
     }
     }*/
}
int PIC3D::PIC3D_PushVelocity()
{
  double particle_xTemp=0,particle_yTemp=0,particle_zTemp=0,cc=0,mFieldTotal=0,FFactor=0;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    if(particleLost[i]||particleFixed[i]) continue;

    //1:half-step accleration in the electric field
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=eField[i][j]*timeStep/2;//*FFactor//unit:m/s
    }

    //2:particle velocity rotate in magnetic field
    mFieldTotal=sqrt(pow(mField[i][0],2) + pow(mField[i][1],2) + pow(mField[i][2],2));
    if(mFieldTotal>1e-10)
    {
      //FFactor=tan(mFieldTotal*timeStep/2.0)  /  (mFieldTotal*timeStep/2.0);  
      FFactor=1;
      cc=FFactor*1.0/(1+pow(timeStep/2,2)*pow(mFieldTotal,2));
      //cout<<cc*(  1+pow(timeStep/2,2)*(pow(mField[i][2],2) - pow(mField[i][0],2) - pow(mField[i][1],2)  )  )*particle_xy[i][5]-particle_xy[i][5]<<endl;
      particle_xTemp=particle_xy[i][1];
      particle_yTemp=particle_xy[i][3];
      particle_zTemp=particle_xy[i][5];
      particle_xy[i][1]=	cc*(
	  (  1+pow(timeStep/2,2)*(  pow(mField[i][0],2) - pow(mField[i][1],2) - pow(mField[i][2],2)  )  )	*	particle_xTemp	+
	  timeStep*(	mField[i][2]+mField[i][0]*mField[i][1]*timeStep/(2)  )			*	particle_yTemp	+
	  timeStep*((-1)*	mField[i][1]+mField[i][0]*mField[i][2]*timeStep/(2)  )			*	particle_zTemp
	  );

      particle_xy[i][3]=	cc*(
	  timeStep*((-1)*	mField[i][2]+mField[i][0]*mField[i][1]*timeStep/(2)  )			*	particle_xTemp	+
	  (  1+pow(timeStep/2,2)*(0-pow(mField[i][0],2) + pow(mField[i][1],2) - pow(mField[i][2],2)  )  )	*	particle_yTemp	+
	  timeStep*(	mField[i][0]+mField[i][1]*mField[i][2]*timeStep/(2)  )			*	particle_zTemp
	  );

      particle_xy[i][5]=	cc*(
	  timeStep*(		mField[i][1]+mField[i][0]*mField[i][2]*timeStep/(2)  )			*	particle_xTemp	+
	  timeStep*((-1)*	mField[i][0]+mField[i][1]*mField[i][2]*timeStep/(2)  )			*	particle_yTemp	+
	  (  1+pow(timeStep/2,2)*(pow(mField[i][2],2) - pow(mField[i][0],2) - pow(mField[i][1],2)  )  )	*	particle_zTemp
	  );
      //cout<<particle_xy[i][5]-particle_zTemp<<endl;
    }
    //3:curvilinear accleration
    //This is a linac disign,so the R is infinite,so this step won't make any change


    //4:half-step accleration in the electric field again
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=eField[i][j]*timeStep/2;//FFactor*//unit:m/s
    }
  }
}

int PIC3D::PIC3D_PushVelocityFirstHalf()
{
  double particle_xTemp=0,particle_yTemp=0,particle_zTemp=0,cc=0,mFieldTotal=0,FFactor=0,mFieldTimeStep=timeStep/2;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    if(particleLost[i]) continue;
    //1:half-step accleration in the electric field
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=eField[i][j]*timeStep/2;//unit:m/s
    }

    //2:particle velocity rotate in magnetic field
    mFieldTotal=sqrt(pow(mField[i][0],2) + pow(mField[i][1],2) + pow(mField[i][2],2));
    if(mFieldTotal>1e-10)
    {
      FFactor=tan(mFieldTotal*mFieldTimeStep/2.0)  /  (mFieldTotal*mFieldTimeStep/2.0);   
      FFactor=1;
      cc=FFactor*1.0/(1+pow(mFieldTimeStep/2,2)*pow(mFieldTotal,2));
      particle_xTemp=particle_xy[i][1];
      particle_yTemp=particle_xy[i][3];
      particle_zTemp=particle_xy[i][5];
      particle_xy[i][1]=	cc*(
	  (  1+pow(mFieldTimeStep/2,2)*(  pow(mField[i][0],2) - pow(mField[i][1],2) - pow(mField[i][2],2)  )  )	*	particle_xTemp	+
	  mFieldTimeStep*(	mField[i][2]+mField[i][0]*mField[i][1]*mFieldTimeStep/(2)  )			*	particle_yTemp	+
	  mFieldTimeStep*((-1)*	mField[i][1]+mField[i][0]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_zTemp
	  );

      particle_xy[i][3]=	cc*(
	  mFieldTimeStep*((-1)*	mField[i][2]+mField[i][0]*mField[i][1]*mFieldTimeStep/(2)  )			*	particle_xTemp	+
	  (  1+pow(mFieldTimeStep/2,2)*(0-pow(mField[i][0],2) + pow(mField[i][1],2) - pow(mField[i][2],2)  )  )	*	particle_yTemp	+
	  mFieldTimeStep*(	mField[i][0]+mField[i][1]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_zTemp
	  );

      particle_xy[i][5]=	cc*(
	  mFieldTimeStep*(	mField[i][1]+mField[i][0]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_xTemp	+
	  mFieldTimeStep*((-1)*	mField[i][0]+mField[i][1]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_yTemp	+
	  (  1+pow(mFieldTimeStep/2,2)*(-pow(mField[i][0],2) - pow(mField[i][1],2) + pow(mField[i][2],2)  )  )	*	particle_zTemp
	  );
    }
  }
}

int PIC3D::PIC3D_PushVelocityLastHalf()
{
  double particle_xTemp=0,particle_yTemp=0,particle_zTemp=0,cc=0,mFieldTotal=0,FFactor=0,mFieldTimeStep=timeStep/2;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    if(particleLost[i]) continue;
    mFieldTotal=sqrt(pow(mField[i][0],2) + pow(mField[i][1],2) + pow(mField[i][2],2));
    if(mFieldTotal>1e-10)
    {
      FFactor=tan(mFieldTotal*mFieldTimeStep/2.0)  /  (mFieldTotal*mFieldTimeStep/2.0);   
      FFactor=1;
      //2:particle velocity rotate in magnetic field
      cc=FFactor*1.0/(1+pow(mFieldTimeStep/2,2)*pow(mFieldTotal,2));
      particle_xTemp=particle_xy[i][1];
      particle_yTemp=particle_xy[i][3];
      particle_zTemp=particle_xy[i][5];
      particle_xy[i][1]=	cc*(
	  (  1+pow(mFieldTimeStep/2,2)*(  pow(mField[i][0],2) - pow(mField[i][1],2) - pow(mField[i][2],2)  )  )	*	particle_xTemp	+
	  mFieldTimeStep*(	mField[i][2]+mField[i][0]*mField[i][1]*mFieldTimeStep/(2)  )			*	particle_yTemp	+
	  mFieldTimeStep*((-1)*	mField[i][1]+mField[i][0]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_zTemp
	  );


      particle_xy[i][3]=	cc*(
	  mFieldTimeStep*((-1)*	mField[i][2]+mField[i][0]*mField[i][1]*mFieldTimeStep/(2)  )			*	particle_xTemp	+
	  (  1+pow(mFieldTimeStep/2,2)*(0-pow(mField[i][0],2) + pow(mField[i][1],2) - pow(mField[i][2],2)  )  )	*	particle_yTemp	+
	  mFieldTimeStep*(	mField[i][0]+mField[i][1]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_zTemp
	  );


      particle_xy[i][5]=	cc*(
	  mFieldTimeStep*(	mField[i][1]+mField[i][0]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_xTemp	+
	  mFieldTimeStep*((-1)*	mField[i][0]+mField[i][1]*mField[i][2]*mFieldTimeStep/(2)  )			*	particle_yTemp	+
	  (  1+pow(mFieldTimeStep/2,2)*(-pow(mField[i][0],2) - pow(mField[i][1],2) + pow(mField[i][2],2)  )  )	*	particle_zTemp
	  );
    }
    //3:curvilinear accleration
    //This is a linac disign,so the R is infinite,so this step won't make any change


    //4:half-step accleration in the electric field again
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=eField[i][j]*timeStep/2;//unit:m/s
    }
    //if(eField[i][2]>1)
    //cout<<i<<"  "<<eField[i][2]<<"  "<<particle_xy[i][4]<<"  "<<particle_xy[i][5]<<endl;
  }
}

int PIC3D::PIC3D_PushPosition(int criParticleLost,int criPost)
{
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    if(particleLost[i]||particleFixed[i]) continue;
    //5:advance the position
    if(particle_xy[i][4]<0)
    {
      particle_xy[i][4]+=particle_xy[i][5]/1*timeStep;
    }else
    {
      particle_xy[i][0]+=particle_xy[i][1]/1*timeStep;//unit:m
      particle_xy[i][2]+=particle_xy[i][3]/1*timeStep;
      particle_xy[i][4]+=particle_xy[i][5]/1*timeStep;
    }
    if(0)
    {
      if(abs(particle_xy[i][0])>=stepx*(numofgridx/2-0.5))
      {
	particle_xy[i][0]=(particle_xy[i][0])/abs(particle_xy[i][0])*stepx*(numofgridx/2-1);
      }
      if(abs(particle_xy[i][2])>=stepy*(numofgridy/2-0.5))
      {
	particle_xy[i][2]=(particle_xy[i][2])/abs(particle_xy[i][2])*stepy*(numofgridy/2-1);
      }
    }
    else if(criParticleLost)//boundary Lost
    {
      double beta = beam.particle_ave[5] / C_light;
      double r    = sqrt(particle_xy[i][0]*particle_xy[i][0]+particle_xy[i][2]*particle_xy[i][2]);
      double ape  = MyLattice.getRadius(particle_xy[i][4]);
      //if(ape<1e88)
      //  cout<<r<<"  "<<ape<<"  "<<(r<ape ? 1:0)<<endl;
      if(r>ape)
      {
	//cout<<"  "<<beam.particle_ave[5]/325.0e6<<"  ";
	ofstream lost("lostparticle.dat",ios::app);
	for(int j=0;j<6;++j)
	{
	  lost<<particle_xy[i][j]<<"   ";
	}
	lost<<"Transvers  "<<MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name<<endl;
	lost.close();
	//cout<<"Particle Lost here for Transvers : "<<MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name<<endl;
	particleLost[i]=1;
#pragma omp atomic
	--particleLiveNumber;
      }

      if(RFQLongtitudeLostSwitch&&MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name=="rfq")
      {
	RFQLongtitudeLostFlag	= true;
	RFQLongtitudeLostSwitch	= false;
      }
      if(RFQLongtitudeLostFlag&&MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name=="rfq")
      {
	if(abs(particle_xy[i][4]-beam.particle_ave[4])>beam.particle_ave[5]/MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].p[4]*1.5)
	{
	  if(1)
	  {
	    //cout<<"  "<<beam.particle_ave[5]/325.0e6<<"  ";
	    cout<<"Particle Lost here for Longtitude: "<<particle_xy[i][4]<<"  "<<MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].p[3]<<endl;
	    particleLost[i]=1;
#pragma omp atomic
	    --particleLiveNumber;
	  }
	  else if(0)
	  {
	    if(particle_xy[i][4]>beam.particle_ave[4])
	    {
	      if(particle_xy[i][5]>beam.particle_ave[5])
	      {
		particle_xy[i][5]=2*beam.particle_ave[5]-particle_xy[i][5];
	      }
	    }
	    else if(particle_xy[i][4]<beam.particle_ave[4])
	    {	  
	      if(particle_xy[i][5]<beam.particle_ave[5])
	      {
		particle_xy[i][5]=2*beam.particle_ave[5]-particle_xy[i][5];
	      }
	    }
	  }
	}
      }
    else if(MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name=="rfq")
    {
      //stop

    }
    if(0&&particle_xy[i][4]>0.5)
    {
      particleFixed[i]=1;
      ++particleFixedNumber;
    }
  }
}

#pragma omp barrier
#pragma omp single
if(criPost)
{
  for(int i=0;i<particlenumber;++i)
  {
    tempOrder=MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].ordernumber;
    if(MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name=="field"
	&&int(MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].p[1])!=2
	&&find(fieldPassed.begin(),fieldPassed.end(),tempOrder)==fieldPassed.end()
      )
    {
      fieldPassed.push_back(tempOrder);
      scanFlag=1;
      ifstream offsetin;
      offsetin.open("fieldoffset.dat");
      string str;
      double phasein;
      int    orderin;
      if(offsetin.is_open())
      {
	while(offsetin>>orderin>>phasein)
	{
	  if(orderin==tempOrder)
	  {
	    ph=phasein;
	    scanFlag=0;
	    break;
	  }
	}
      }
      offsetin.close();
      if(scanFlag==1)
      {
#pragma omp parallel num_threads(1)
	{
	  ph=ScanPhase(tempOrder);
	}
	ofstream offsetout("fieldoffset.dat",ios::app);
	offsetout<<tempOrder<<"  "<<ph<<endl;
	offsetout.close();
      }
      cout<<"ph offset = "<<ph<<endl;
      int a;
      if(tempOrder<2)
      {

      }
      //cin>>a;
      MyLattice.FIELD3D[tempOrder].SetOffset(ph);
      scanFlag=1;
      break;
    }
  }
}
}





int PIC3D::CalEmit()
{
  beam.caculate_emittance(particle_xy,particleLost,frq);
  //cout<<"emittancex=   "<<beam.getemittancex()<<endl;
  if(1)
  {
    double beta = beam.eBetaAverage;
    double gamma= beam.eGammaAverage;
    ParaTransfer.push_back({
	beam.particle_ave[4],    	//ave of position z
	(gamma-1)*BaseEnergyInMeV,
	sqrt(beam.getsigmax()),
	beam.getsigmaxdx()/sqrt(beam.getsigmax()),
	sqrt(beam.getsigmay()),
	beam.getsigmaydy()/sqrt(beam.getsigmay()),
	sqrt(beam.z_sigma),	//sigma_phi
	beam.getsigmazdz()/sqrt(beam.getsigmaz()),
	sqrt(beam.sigmaz),	//sigma_dw
	beam.emittancex*beta*gamma,
	beam.emittancey*beta*gamma,
	beam.emittancez/360/frq*1e9*BaseEnergyInMeV*1e3,//	keV,ns
	beam.TBetax,
	beam.TAlphax,
	beam.TBetay,
	beam.TAlphay,
	beam.TBetaz,
	beam.TAlphaz,
	particleLiveNumber,
	double(particleLiveNumber)/particlenumber,
	beam.dppmax,
	beam.dwwmax,
	beam.dbbmax

    });
  }
  else if(0)
    ParaTransfer.push_back({particle_xy[0][0],particle_xy[0][1],particle_xy[0][2],particle_xy[0][3],particle_xy[0][4],particle_xy[0][5]});
  else
    ParaTransfer.push_back({particleLiveNumber});
}

int PIC3D::CalMatchR()
{
  double K=	2*pow(Macrocharge,2)*particlenumber
    /
    (  pow(Egamma*Ebeta,3) * Macromass * pow(C_light,2)  );

  double kz0=Macrocharge  *  MyLattice.getField( vector<double>(6,0),0).back()
    /
    (2*Egamma * Macromass);
  double kz1=Macrocharge  *  MyLattice.getField( vector<double>(6,0),0).back()
    /
    (2*M_PI*Egamma * Macromass);

  double u=	K/(2*kz0);
  double rb=	sqrt(  sqrt(pow(u,2)+1/kz0)  +  u  );
  double r=Egamma * Macromass*sqrt(particle_xy[0][1]*particle_xy[0][1]+particle_xy[0][3]*particle_xy[0][3])/(Macrocharge  *  MyLattice.getField( vector<double>(6,0),0).back());
  cout<<r<<"   "<<kz1<<endl;
}


double PIC3D::ScanPhase(int order)
{
  LostFlag=0;//boundary Lost
  cout<<"ScanPhase , field NO."<<order<<endl;
  const v2d 	particleT=particle_xy;
  const vector<int>	particleL=particleLost;
  const double	particleN=particleLiveNumber;
  const double 	timeStart=timeline;
  v1d energyScanPhase;
  v1d offsetScanPhase;
  double phase=0;
  //PIC3D *accScanPhase=new PIC3D;
  //*accScanPhase=*this;
  while(phase<360)
  {
    double timetemp=timeStart;
    particle_xy=particleT;
    particleLost=particleL;
    particleLiveNumber=particleN;
    beam.caculate_emittance(particle_xy,particleLost,frq);
    while(1)
    {
      timetemp+=timeStep;
      MyLattice.FIELD3D[order].ScanOffset(phase);
      if(0)
      {
	PIC3D_SpaceCharge();
      }
      {
	PIC3D_PUSH(timetemp,0);
      }
      int count=0;
      for(int i=0;i<particlenumber;++i)
      {
	if(particleLost[i]||particleFixed[i]) continue;
	//particle loss in scan need to considered!!!
	if(   MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].ordernumber	!=order
	    ||MyLattice.ELE[MyLattice.getEleNum(particle_xy[i][4])].name	!="field")
	  count++;
      }
      if(count==particleLiveNumber) break;
    }
    beam.caculate_emittance(particle_xy,particleLost,frq);
    int size=energyScanPhase.size();
    if(size>2)
    {
      if(beam.eGammaAverage<energyScanPhase[size-1]&&energyScanPhase[size-1]>energyScanPhase[size-2])
      {
	timeline=timeStart;
	particle_xy=particleT;
	particleLost=particleL;
	particleLiveNumber=particleN;
	LostFlag=1;//boundary Lost
	beam.caculate_emittance(particle_xy,particleLost,frq);
	return offsetScanPhase[size-2];
      }
      if(beam.eGammaAverage>energyScanPhase[size-1]&&energyScanPhase[size-1]<energyScanPhase[size-2])
      {
	timeline=timeStart;
	particle_xy=particleT;
	particleLost=particleL;
	particleLiveNumber=particleN;
	beam.caculate_emittance(particle_xy,particleLost,frq);
	return offsetScanPhase[size-2]+180;
      }
    }
    energyScanPhase.	push_back(beam.eGammaAverage);
    offsetScanPhase.	push_back(phase);
    //cout<<"phase :"<<phase<<"  "<<(energyScanPhase.back()-1)*BaseEnergyInMeV<<endl;
    phase+=1;
  }
  auto it = max_element(energyScanPhase.begin(),energyScanPhase.end());
  int lo= it- energyScanPhase.begin();
  cout<<offsetScanPhase[lo]<<"  "<<(energyScanPhase[lo]-1)*BaseEnergyInMeV<<"  "<<beam.particle_ave[4]<<endl;
  cin>>phase;

  //delete accScanPhase;
  timeline=timeStart;
  particle_xy=particleT;
  particleLost=particleL;
  particleLiveNumber=particleN;
  beam.caculate_emittance(particle_xy,particleLost,frq);
  return offsetScanPhase[lo];
}

double PIC3D::getfield(v3d &field,double x,double y,double z)
{
  double xn,yn,zn,sum=0;
  int xGrid,yGrid,zGrid,indexGrid;
  vd v(8);
  xn=(x-gridx[0])/stepx;
  yn=(y-gridy[0])/stepy;
  zn=(z-gridz[0])/stepz;
  int  cx=floor((x-gridx[0])/stepx);
  int  cy=floor((y-gridy[0])/stepy);
  int  cz=floor((z-gridz[0])/stepz);
  if(   (cx>=numofgridx-1)
      ||(cy>=numofgridy-1)
      ||(cz>=numofgridz)
      ||cx<0||cy<0||cz<0)
    /*
       if(	x>gridx.back()	||y>gridy.back()	||z>gridz.back()||
       x<gridx[0]	||y<gridy[0]		||z<gridz[0])
       */
  {
    cout<<cx<<"  "<<cy<<"  "<<cz<<"getfield error in pic3d"<<endl;
    for(int j=0;j<8;++j)
    {
      v[j]=  0;
    }
    return 0;
  }
       else
       {
	 double sumv=0,a,b;
	 bool c;
	 for(int i=0;i<8;++i)
	 {
	   a = ceil(xn);
	   b = floor(xn);
	   c = i%2<1;
	   if(a!=b)
	   {
	     xGrid=(c?a:b);
	   }
	   else
	   {
	     xGrid=(c?a:a+1);
	   }

	   a = ceil(yn);
	   b = floor(yn);
	   c = i%4<2;
	   if(a!=b)
	   {
	     yGrid=(c?a:b);
	   }
	   else
	   {
	     yGrid=(c?a:a+1);
	   }

	   a = ceil(zn);
	   b = floor(zn);
	   c = i<4;
	   if(a!=b)
	   {
	     zGrid=(c?a:b);
	   }
	   else
	   {
	     zGrid=(c?a:a+1);
	   }

	   v[i]=abs(xn-xGrid) * abs(yn-yGrid) * abs(zn-zGrid);

	   a = ceil(xn);
	   b = floor(xn);
	   c = !(i%2<1);
	   if(a!=b)
	   {
	     xGrid=(c?a:b);
	   }
	   else
	   {
	     xGrid=(c?a:a+1);
	   }

	   a = ceil(yn);
	   b = floor(yn);
	   c = !(i%4<2);
	   if(a!=b)
	   {
	     yGrid=(c?a:b);
	   }
	   else
	   {
	     yGrid=(c?a:a+1);
	   }

	   a = ceil(zn);
	   b = floor(zn);
	   c = !(i<4);
	   if(a!=b)
	   {
	     zGrid=(c?a:b);
	   }
	   else
	   {
	     zGrid=(c?a:a+1);
	   }

	   sum+=field[xGrid][yGrid][zGrid]*v[i];
	 }
	 return sum;
       }
}

void PIC3D::benchmark_flux()
{
  int div=100;
  double r=0.001;
  int rnum=100;
  double rlimit=0.9*(gridx.back()-particle_xy[0][0]);

  double rstep=(rlimit-r)/rnum;
  double phi=0,sit=0;
  double x,y,z,ex,ey,ez,e;
  ofstream fluxout("benchmark_flux.dat");
  while(r<rlimit)
  {
    r+=rstep;
    double etotal=0;
    phi=0;
    for(int i=0;i<div;++i)
    {
      phi=double(i)/div*M_PI*2;
      for(int j=0;j<div;++j)
      {
	sit=(double(j)+0.5)/div*M_PI;
	x=r*sin(sit)*cos(phi)	+particle_xy[0][0];
	y=r*sin(sit)*sin(phi)	+particle_xy[0][2];
	z=r*cos(sit)		+particle_xy[0][4];

	ex=getfield(ex_grid,x,y,z);
	ey=getfield(ey_grid,x,y,z);
	ez=getfield(ez_grid,x,y,z);

	//if(i==0&&j==div/2) cout<<ex<<"  "<<ey<<"  "<<ez<<endl;
	//if(i==0&&j==div-1) cout<<ex<<"  "<<ey<<"  "<<ez<<endl;

	/*
	   ex=Macrocharge/r/r/(4*M_PI*Epsilon)*sin(sit)*cos(phi);
	   ey=Macrocharge/r/r/(4*M_PI*Epsilon)*sin(sit)*sin(phi);
	   ex=Macrocharge/r/r/(4*M_PI*Epsilon)*cos(sit);
	   */
	e=ex*sin(sit)*cos(phi)+
	  ey*sin(sit)*sin(phi)+
	  ez*cos(sit);
	//cout<<e<<"  e"<<endl;
	//e=1;
	e*=	(r*1.0/div*M_PI)   *   (r*sin(sit-0.5/div*M_PI)*2*M_PI/div+r*sin(sit+0.5/div*M_PI)*2*M_PI/div)/2;
	etotal+=e;
      }
    }
    fluxout<<r<<"  "<<etotal<<"  "<<Macrocharge/Epsilon<<endl;
    cout<<r<<"  "<<etotal*3<<"  "<<Macrocharge/Epsilon<<endl;
  }
}

void PIC3D::benchmark_phi_theory_x()
{
  cout<<"benchmark_phi_theory X!"<<endl;
  if(particlenumber!=1) return;
  ofstream ff("field_phi_theory_x.dat");

  int dotnumber=gridx.size()*3;
  double x,y,z,xc,yc,zc,r;
  const double x0	=particle_xy[0][0];
  const double y0	=particle_xy[0][2];
  const double z0	=particle_xy[0][4];
  const double xleft 	=gridx[0];
  const double xright	=gridx.back();
  const double xspace	=stepx*(gridx.size()+1);
  const double yspace	=stepy*(gridy.size()+1);
  const double zspace	=stepz*(gridz.size());
  int nue=100;
  int factorx,factory;

  double field_theory,field_theory_x,field_theory_y,field_theory_z;

  for(int j=0;j<dotnumber;++j)
  {
    x=0.95*(xleft+j*(xright-xleft)/(dotnumber-1));
    y=y0;
    z=z0;
    field_theory=0;
    field_theory_x=0;
    field_theory_y=0;
    field_theory_z=0;
    for(int xn=-nue;xn<=nue-1;++xn)
    {
      for(int yn=-nue;yn<=nue-1;++yn)
      {
	for(int zn=0;zn<1;++zn)
	{
	  xc=x0+xn*xspace;
	  yc=y0+yn*yspace;
	  zc=z0+zn*zspace;
	  r = sqrt(pow(x-xc,2)+pow(y-yc,2)+pow(z-zc,2));
	  factorx=(xn%2==0)?1:-1;
	  factory=(yn%2==0)?1:-1;
	  field_theory	+=factorx*factory*Macrocharge/(4*M_PI*Epsilon*r);
	  if(0)
	  {
	    field_theory =factorx*factory*Macrocharge/(4*M_PI*Epsilon*r);
	    field_theory_x+=(x-xc)/r*field_theory;
	    field_theory_y+=(y-yc)/r*field_theory;
	    field_theory_z+=(z-zc)/r*field_theory;
	  }
	}
      }
    }
    ff<<setw(13)
      <<x
      <<"  "
      <<setw(13)
      <<-gridx.back()-stepx

      <<setw(13)
      <<getfield(phi_grid,x,y,z)
      <<"  "

      <<setw(13)
      <<field_theory;

    if(0)
    {
      ff<<"  "      
	<<setw(13)
	<<field_theory_x

	<<"  "      
	<<setw(13)
	<<field_theory_y

	<<"  "      
	<<setw(13)
	<<field_theory_z;
    }
    ff<<endl;
  }
  ff.close();
}

void PIC3D::benchmark_phi_theory_z()
{
  cout<<"benchmark_phi_theory Z!"<<endl;
  if(particlenumber!=1) return;
  ofstream ff("field_phi_theory_z.dat");

  int dotnumber=gridz.size()*20;
  double x,y,z,xc,yc,zc,r;
  const double x0	=particle_xy[0][0];
  const double y0	=particle_xy[0][2];
  const double z0	=particle_xy[0][4];
  const double zleft 	=gridz[0];
  const double zright	=gridz.back();
  const double xspace	=stepx*(gridx.size()+1);
  const double yspace	=stepy*(gridy.size()+1);
  const double zspace	=stepz*(gridz.size());
  int nue=100;
  //int factorx,factory;

  double field_theory,field_theory_x,field_theory_y,field_theory_z;

  for(int j=0;j<dotnumber;++j)
  {
    x=x0;
    y=y0;
    z=0.95*(zleft+j*(zright-zleft)/(dotnumber-1));
    field_theory=0;
    field_theory_x=0;
    field_theory_y=0;
    field_theory_z=0;
    for(int xn=0;xn<1;++xn)
    {
      for(int yn=0;yn<1;++yn)
      {
	for(int zn=-nue;zn<=nue;++zn)
	{
	  xc=x0+xn*xspace;
	  yc=y0+yn*yspace;
	  zc=z0+zn*zspace;
	  r = sqrt(pow(x-xc,2)+pow(y-yc,2)+pow(z-zc,2));
	  field_theory	+=Macrocharge/(4*M_PI*Epsilon*r);
	  if(0)
	  {
	    field_theory =Macrocharge/(4*M_PI*Epsilon*r);
	    field_theory_x+=(x-xc)/r*field_theory;
	    field_theory_y+=(y-yc)/r*field_theory;
	    field_theory_z+=(z-zc)/r*field_theory;
	  }
	}
      }
    }
    ff<<setw(13)
      <<z-z0
      <<"  "
      <<setw(13)
      <<-gridz.back()-stepz

      <<setw(13)
      <<getfield(phi_grid,x,y,z)
      <<"  "

      <<setw(13)
      <<field_theory;

    if(0)
    {
      ff<<"  "      
	<<setw(13)
	<<field_theory_x

	<<"  "      
	<<setw(13)
	<<field_theory_y

	<<"  "      
	<<setw(13)
	<<field_theory_z;
    }
    ff<<endl;
  }
  ff.close();
}

void PIC3D::benchmark_electri_theory_x()
{
  cout<<"benchmark_electri_theory X!"<<endl;
  if(particlenumber!=1) return;
  ofstream ff("field_electri_theory_x.dat");

  int dotnumber=gridx.size()*3;
  double x,y,z,xc,yc,zc,r;
  const double x0	=particle_xy[0][0];
  const double y0	=particle_xy[0][2];
  const double z0	=particle_xy[0][4];
  const double left 	=gridx[0];
  const double right	=gridx.back();
  const double xspace	=stepx*(gridx.size()+1);
  const double yspace	=stepy*(gridy.size()+1);
  const double zspace	=stepz*(gridz.size());
  int nue=200;
  //int factorx,factory;

  double field_theory,field_theory_x,field_theory_y,field_theory_z;

  for(int j=0;j<dotnumber;++j)
  {
    if(0)
    {
      cout<<j<<"/"<<dotnumber<<endl;
    }
    x=0.95*(left+j*(right-left)/(dotnumber-1));
    y=y0;
    z=z0;
    field_theory=0;
    field_theory_x=0;
    field_theory_y=0;
    field_theory_z=0;
    for(int xn=-nue;xn<=nue-1;++xn)
    {
      for(int yn=-nue;yn<=nue-1;++yn)
      {
	for(int zn=0;zn<1;++zn)
	{
	  xc=x0+xn*xspace;
	  yc=y0+yn*yspace;
	  zc=z0+zn*zspace;
	  r = sqrt(pow(x-xc,2)+pow(y-yc,2)+pow(z-zc,2));
	  if(1)
	  {
	    field_theory =Macrocharge/(4*M_PI*Epsilon*r*r);
	    field_theory_x+=(x-xc)/r*field_theory;
	    field_theory_y+=(y-yc)/r*field_theory;
	    field_theory_z+=(z-zc)/r*field_theory;
	  }
	}
      }
    }
    ff<<setw(13)
      <<x-x0
      <<"  "
      <<setw(13)
      <<-gridx.back()-stepx

      <<setw(13)
      <<getfield(ex_grid,x,y,z)
      <<"  "

      <<setw(13)
      <<getfield(ey_grid,x,y,z)
      <<"  "

      <<setw(13)
      <<getfield(ez_grid,x,y,z);

    if(1)
    {
      ff<<"  "      
	<<setw(13)
	<<field_theory_x

	<<"  "      
	<<setw(13)
	<<field_theory_y

	<<"  "      
	<<setw(13)
	<<field_theory_z;
    }
    ff<<endl;
  }
  ff.close();
}

void PIC3D::benchmark_electri_theory_z()
{
  cout<<"benchmark_electri_theory Z!"<<endl;
  if(particlenumber!=1) return;
  ofstream ff("field_electri_theory_z.dat");

  int dotnumber=gridz.size()*20;
  double x,y,z,xc,yc,zc,r;
  const double x0	=particle_xy[0][0];
  const double y0	=particle_xy[0][2];
  const double z0	=particle_xy[0][4];
  const double zleft 	=gridz[0];
  const double zright	=gridz.back();
  const double xspace	=stepx*(gridx.size()+1);
  const double yspace	=stepy*(gridy.size()+1);
  const double zspace	=stepz*(gridz.size());
  int nue=2000;
  //int factorx,factory;

  double field_theory,field_theory_x,field_theory_y,field_theory_z;

  for(int j=0;j<dotnumber;++j)
  {
    x=x0;
    y=y0;
    z=0.95*(zleft+j*(zright-zleft)/(dotnumber-1));
    field_theory=0;
    field_theory_x=0;
    field_theory_y=0;
    field_theory_z=0;
    for(int xn=0;xn<1;++xn)
    {
      for(int yn=0;yn<1;++yn)
      {
	for(int zn=-nue;zn<=nue;++zn)
	{
	  xc=x0+xn*xspace;
	  yc=y0+yn*yspace;
	  zc=z0+zn*zspace;
	  r = sqrt(pow(x-xc,2)+pow(y-yc,2)+pow(z-zc,2));
	  if(1)
	  {
	    field_theory =Macrocharge/(4*M_PI*Epsilon*r*r);
	    field_theory_x+=(x-xc)/r*field_theory;
	    field_theory_y+=(y-yc)/r*field_theory;
	    field_theory_z+=(z-zc)/r*field_theory;
	  }
	}
      }
    }
    ff<<setw(13)
      <<z-z0
      <<"  "
      <<setw(13)
      <<-gridz.back()-stepz

      <<setw(13)
      <<getfield(ex_grid,x,y,z)
      <<"  "

      <<setw(13)
      <<getfield(ey_grid,x,y,z)
      <<"  "

      <<setw(13)
      <<getfield(ez_grid,x,y,z);

    if(1)
    {
      ff<<"  "      
	<<setw(13)
	<<field_theory_x

	<<"  "      
	<<setw(13)
	<<field_theory_y

	<<"  "      
	<<setw(13)
	<<field_theory_z;
    }
    ff<<endl;
  }
  ff.close();
}

