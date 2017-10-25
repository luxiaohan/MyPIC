#include "Pic2d.h"

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
using namespace std;
//GS UNIT
/*
   2015-10-28@LiuZhicong
   */

PIC2D::PIC2D(){  Initial();}

PIC2D::PIC2D(double parnumber,double gridnumx,double gridnumy,double time,double gam,double macroNum,int sli):
  particlenumber(parnumber),
  numofgridx(gridnumx),
  numofgridy(gridnumy),
  timeStep(time),
  timeStepFix(time),
  Egamma(gam),
  Huge(macroNum),
  slice(sli+1)
{ 
  Initial();
  cout<<slice<<endl;
  cout<<"PicGamma=  "<<Egamma<<"\nSpeed   =  "<<(Ebeta*C_light)<<endl;
}

PIC2D::~PIC2D()
{
  fftw_free(in);
  fftw_free(out);
  fftw_free(outf);
  fftw_destroy_plan(g2f);
  fftw_destroy_plan(f2g);
}

int PIC2D::Initial()
{
  {
  in  	= (double*) 	fftw_malloc(sizeof(double)*numofgridx*numofgridy);
  out 	= (double*) 	fftw_malloc(sizeof(double)*numofgridx*numofgridy);
  outf 	= (double*) 	fftw_malloc(sizeof(double)*numofgridx*numofgridy);
  }
  fftw_plan_with_nthreads(omp_get_max_threads());
  g2f=fftw_plan_r2r_2d(numofgridx,numofgridy,in,out,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE);
  fftw_plan_with_nthreads(omp_get_max_threads());
  f2g=fftw_plan_r2r_2d(numofgridx,numofgridy,out,outf,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE);
  beam.initial('e',Egamma,particle_xy.size());
  return 1;
}

int PIC2D::ReadParticle(const vector<double> &particleParameter)
{
  KV2D_T(particleParameter,particle_xy);
  //KV2D_A(particleParameter,particle_xy);
  for(int i=0;i<4;i++)
  {
    if(i%2==0)
    particle_xy[0][i]=0.005;
    else
    particle_xy[0][i]=0.00;
  }

  for(int i=0;i<particlenumber;i++) 
  {
    particle_xy[i][5]=Ebeta*C_light;
  }

  cout<<"Distribution : KV"<<endl;
  ofstream pout("data/initialKV.dat");
  printvector(pout,particle_xy);
  pout.close();
  return 0;
}

double PIC2D::ReadLattice(const char *p)
{
  return MyLattice.read(p,Egamma,1);
}

int PIC2D::PIC2D_ICEFROG()
{
  //before cycle, we push x in halfstep.		FOR ICEFROG!!!
  timeStep=-timeStep;
  PIC2D_ExternalField(0.0);
  PIC2D_PushVelocityLastHalf(); 
  timeStep=timeStepFix;
}

int PIC2D::PIC2D_FIREFROG()
{
  //afte cycle, we push dx in halfstep			FOR ICEFROG!!!
  PIC2D_PushVelocityFirstHalf();
}


int PIC2D::PIC2D_Weight()
{
  if(particlenumber==0) return -1;
#pragma omp single
{
  for(int i=0;i<numofgridx;i++)
    for(int j=0;j<numofgridy;j++)
    {
      rho_grid[i][j]=0;
    }
  //caculate the step
  int girdfixedflag=1;
  if(girdfixedflag==1)
  {
    for(int i=0;i<gridx.size();i++)
    {
      //gird fixed
      gridx[i]=stepx*(i+0.5-numofgridx/2);
      gridy[i]=stepy*(i+0.5-numofgridy/2);
    }
  }
  else if(girdfixedflag==0)
  {
    /*nomalize the coordinate of particle to the grid*/
    double xmin=0,xmax=0,xave=0,ymin=0,ymax=0,yave=0;
    //    xmin=min(&particle_xy[i][0],particlenumber);
    //    xmax=max(&particle_xy[i][0],particlenumber);
    //    xave=ave(&particle_xy[i][0],particlenumber);
    //    ymin=min(&particle_xy[i][2],particlenumber);
    //    ymax=max(&particle_xy[i][2],particlenumber);
    //    yave=ave(&particle_xy[i][2],particlenumber);
    //get the coordinate of grid

    int cut;
    if(numofgridx<=8)
    {
      cut=1;
    }
    else if(numofgridx>8&&numofgridx<=32)
    {
      cut=3;
    }
    else if(numofgridx>32)
    {
      cut=3;
    }
    stepx=(xmax-xmin)/(numofgridx/cut-3);
    stepy=(ymax-ymin)/(numofgridy/cut-3);
    
    if(xmax-xmin<1) stepx=1e-3;
    if(ymax-ymin<1) stepy=1e-3;
    
    if(particlenumber==1)
    {
      stepx=1;
      stepy=1;
    }

    for(int i=0;i<gridx.size();i++)
    {
      //gird change with the particle
      gridx[i]=xave+stepx*(i+0.5-numofgridx/2);
      gridy[i]=yave+stepy*(i+0.5-numofgridy/2);
    }
  }
  //get which gird the particle locate
  vt=stepx*stepy;
}
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
    countx[i]=floor((particle_xy[i][0]-gridx[0])/stepx);
    county[i]=floor((particle_xy[i][2]-gridy[0])/stepy);
    if(countx[i]>=numofgridx) countx[i]=numofgridx-1;
    if(county[i]>=numofgridy) county[i]=numofgridy-1;
    v4[i]=   (particle_xy[i][0]-gridx[countx[i]])  *   (particle_xy[i][2]-gridy[county[i]])   /vt;
    v3[i]=abs(particle_xy[i][0]-gridx[countx[i]+1])*   (particle_xy[i][2]-gridy[county[i]])   /vt;
    v2[i]=   (particle_xy[i][0]-gridx[countx[i]])  *abs(particle_xy[i][2]-gridy[county[i]+1]) /vt;
    v1[i]=abs(particle_xy[i][0]-gridx[countx[i]+1])*abs(particle_xy[i][2]-gridy[county[i]+1]) /vt;
  }
#pragma omp single
{

  for(int i=0;i<particlenumber;++i)
  {
    rho_grid[countx[i]  ][county[i]  ]	+=	v1[i]*Macrocharge;
    rho_grid[countx[i]+1][county[i]  ]	+=	v2[i]*Macrocharge;
    rho_grid[countx[i]  ][county[i]+1]	+=	v3[i]*Macrocharge;
    rho_grid[countx[i]+1][county[i]+1]	+=	v4[i]*Macrocharge;
  }
}
/*
#pragma omp single
{
  double rhosum=0;
  for(int i=0;i<gridx.size();++i)
  {
    for(int j=0;j<gridy.size();++j)
    {
      rhosum+=rho_grid[i][j];
    }
  }
  cout<<"Rho total = "<<rhosum<<endl;
}
*/
    //cout<<rho_grid[countx[i]][county[i]]+rho_grid[countx[i]+1][county[i]]+rho_grid[countx[i]][county[i]+1]+rho_grid[countx[i]+1][county[i]+1]<<endl;
  
  return 1;
}

int PIC2D::PIC2D_FFT()
{
  //now we start fft
  /*FFT
    2015-10-29@LiuZhicong
    */
  int type=1;
  if(type==1)//DST
  {
    //forward Fourier
#pragma omp for
    for(int i=0;i<numofgridx;i++)
    {
      for(int j=0;j<numofgridy;j++)
      {
	in[i*numofgridy+j]=rho_grid[i][j]/stepx/stepy;
      }
    }
#pragma omp barrier
#pragma omp single 
    fftw_execute(g2f);
#pragma omp barrier
    //devide by K
    double K_rho2phi=0,knx=0,kny=0;
#pragma omp for
    for(int i=0;i<numofgridx;i++)
    {
      for(int j=0;j<numofgridy;j++)
      {
	//if(i==0&&j==0) continue;
	knx=M_PI/2*(i+1)/(numofgridx+1);
	kny=M_PI/2*(j+1)/(numofgridy+1);
	//knx=M_PI/2*(i)/(numofgridx);
	//kny=M_PI/2*(j)/(numofgridy);
	//K_rho2phi=pow(2*sin(knx*stepx)/stepx,2)+pow(2*sin(kny*stepy)/stepy,2);
	K_rho2phi=pow(2*sin(knx)/stepx,2)+pow(2*sin(kny)/stepy,2);
	//K_rho2phi=4*(knx*knx+kny*kny);
	//this one is very like the one above ,and according to my caclulate ,this one is the analysic solution.but the charge is not normalize to 1 but to 1.00004. Confusing!
	//cout<<"i = "<<i<<" ,j = "<<j<<" ,K = "<<K_rho2phi<<endl;
	out[i*numofgridy+j]/=(K_rho2phi*Epsilon);
      }
    }
    //out[0]=0;

    //backward Fourier
#pragma omp barrier
#pragma omp single 
    fftw_execute(f2g);
#pragma omp barrier
#pragma omp for
    for(int i=0;i<numofgridx;i++)
    {
      for(int j=0;j<numofgridy;j++)
      {
	phi_grid[i][j]=outf[i*numofgridy+j]/(4*(numofgridy+1)*(numofgridy+1));
      }
    }
  }

  //Diferentiate to get the e field
#pragma omp for
  for(int i=1;i<numofgridx-1;i++)
  {
    for(int j=1;j<numofgridy-1;j++)
    {
      ex_grid[i][j]=-(phi_grid[i+1][j]-phi_grid[i-1][j])/stepx/2;//unit:V/m;
      ey_grid[i][j]=-(phi_grid[i][j+1]-phi_grid[i][j-1])/stepy/2;//divide gird step;
      //should I caculate the eField at half grid?
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

int PIC2D::PIC2D_InternalField()
{
  double e1=0,e2=0,e3=0,b1=0,b2=0,b3=0;
  //printvector(cout,ex_grid);
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
    e1  =   ex_grid[countx[i]][county[i]]      *v1[i]
      +     ex_grid[countx[i]+1][county[i]]    *v2[i]
      +     ex_grid[countx[i]][county[i]+1]    *v3[i]
      +     ex_grid[countx[i]+1][county[i]+1]  *v4[i];
    e2  =   ey_grid[countx[i]][county[i]]      *v1[i]
      +     ey_grid[countx[i]+1][county[i]]    *v2[i]
      +     ey_grid[countx[i]][county[i]+1]    *v3[i]
      +     ey_grid[countx[i]+1][county[i]+1]  *v4[i];
    e3  =   0;
    /*
    ofstream ff("temp1");
    for(int j=0;j<gridy.size();++j)
    {
      ff<<j
	<<"  "
	<<sqrt(pow(ex_grid[j][county[i]],2)+pow(ey_grid[j][county[i]],2))
	//<<sqrt(pow((ex_grid[j][county[i]]+ex_grid[j][county[i]+1])/2,2)+pow((ey_grid[j][county[i]]+ey_grid[j][county[i]+1])/2,2))
	<<"  "
	<<Macrocharge/Epsilon/(2*M_PI*sqrt(pow(gridx[j]-0.005,2)))
	<<endl;
    }
    ff.close();
    cout<<"e1  "<<ex_grid[countx[i]][county[i]]<<"  "<<ey_grid[countx[i]][county[i]];
    cout<<"  ecal = "<<Macrocharge/Epsilon/(2*M_PI*sqrt(stepx*stepx+stepy*stepy))/sqrt(2)*2<<endl;
    cout<<"e2  "<<ex_grid[countx[i]-1][county[i]-1]<<"  "<<ey_grid[countx[i]][county[i]];
    cout<<"  ecal = "<<Macrocharge/Epsilon/(2*M_PI*sqrt(stepx*stepx+stepy*stepy))/sqrt(2)*2/3<<endl;
    */
    e1=e1/C_light;
    e2=e2/C_light;
    e3=e3/C_light;
    //Lorentz Transform
    
    eFieldInner[i][0]=C_light	*Egamma*(e1+Ebeta*b2);
    eFieldInner[i][1]=C_light	*Egamma*(e2-Ebeta*b1);
    eFieldInner[i][2]=C_light	*e3;
    mFieldInner[i][0]=		Egamma*(b1-Ebeta*e3);
    mFieldInner[i][1]=		Egamma*(b2+Ebeta*e3);
    mFieldInner[i][2]=		b3;
    /*
    eFieldInner[i][0]=0;
    eFieldInner[i][1]=0;
    eFieldInner[i][2]=0;
    mFieldInner[i][0]=0;
    mFieldInner[i][1]=0;
    mFieldInner[i][2]=0;
    */
  }
}

int PIC2D::PIC2D_ExternalFieldPre(double time)
{
  vector<double> exterFieldTemp(6);
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
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

int PIC2D::PIC2D_ExternalField(double time)
{
  vector<double> exterFieldTemp(6);
#pragma omp for
  for(int i=0;i<particlenumber;i++)
  {
    exterFieldTemp=MyLattice.getField(particle_xy[i],time);		//Ex,Ey,Ez,Bx,By,Bz
    //	Bz is positive when point at direction s
    for(int j=0;j<3;j++)
    {
      eField[i][j]  =  (eFieldInner[i][j]+exterFieldTemp[j])	*Macrocharge/Macromass/Egamma;//the field had time Q and divided by M already!!
      mField[i][j]  =  (mFieldInner[i][j]+exterFieldTemp[j+3])	*Macrocharge/Macromass/Egamma;//the field had time Q and divided by M already!!
    }
    mField[i][2]=-mField[i][2];//for  (Bz is positive when point at direction s)
  }
}

int PIC2D::PIC2D_PUSH(double time)
{
  //Adjust timeStep according to Element Length
#pragma omp single
  {
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

/*  if(coun==0)
  {
    PIC2D_ExternalFieldPre();
  }else*/
  {
    PIC2D_ExternalField(time);
  }
  PIC2D_PushVelocityFirstHalf();
#pragma omp single
  {
    if(0)
    {
      if(coun%(slice*4)==0)   CalEmit(); 
    }else
    {
      CalEmit();
    }
  coun++;
  }
  /*
  if(coun==0)
  {
    PIC2D_ExternalField();
  }*/
  PIC2D_PushVelocityLastHalf();
  PIC2D_PushPosition();
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
int PIC2D::PIC2D_PushVelocity()
{
  double particle_xTemp=0,particle_yTemp=0,particle_zTemp=0,cc=0,mFieldTotal=0,FFactor=0;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    mFieldTotal=sqrt(pow(mField[i][0],2) + pow(mField[i][1],2) + pow(mField[i][2],2));
    if(mFieldTotal<1e-10)
    {
      continue;
    }
    FFactor=tan(mFieldTotal*timeStep/2.0)  /  (mFieldTotal*timeStep/2.0);
    FFactor=1;
    //1:half-step accleration in the electric field
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=FFactor*eField[i][j]*timeStep/2;//unit:m/s
    }

    //2:particle velocity rotate in magnetic field
    cc=FFactor*1.0/(1+pow(timeStep/2,2)*pow(mFieldTotal,2));
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
	timeStep*(	mField[i][1]+mField[i][0]*mField[i][2]*timeStep/(2)  )			*	particle_xTemp	+
	timeStep*((-1)*	mField[i][0]+mField[i][1]*mField[i][2]*timeStep/(2)  )			*	particle_yTemp	+
	(  1+pow(timeStep/2,2)*(-pow(mField[i][0],2) - pow(mField[i][1],2) + pow(mField[i][2],2)  )  )	*	particle_zTemp
	);

    //3:curvilinear accleration
    //This is a linac disign,so the R is infinite,so this step won't make any change


    //4:half-step accleration in the electric field again
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=FFactor*eField[i][j]*timeStep/2;//unit:m/s
    }
  }
}

int PIC2D::PIC2D_PushVelocityFirstHalf()
{
  double particle_xTemp=0,particle_yTemp=0,particle_zTemp=0,cc=0,mFieldTotal=0,FFactor=0,mFieldTimeStep=timeStep/2;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    mFieldTotal=sqrt(pow(mField[i][0],2) + pow(mField[i][1],2) + pow(mField[i][2],2));
    if(mFieldTotal<1e-10)
    {
      continue;
    }
    FFactor=tan(mFieldTotal*mFieldTimeStep/2.0)  /  (mFieldTotal*mFieldTimeStep/2.0);
    FFactor=1;
    //1:half-step accleration in the electric field
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=eField[i][j]*timeStep/2;//unit:m/s
    }

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
}

int PIC2D::PIC2D_PushVelocityLastHalf()
{
  double particle_xTemp=0,particle_yTemp=0,particle_zTemp=0,cc=0,mFieldTotal=0,FFactor=0,mFieldTimeStep=timeStep/2;
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    mFieldTotal=sqrt(pow(mField[i][0],2) + pow(mField[i][1],2) + pow(mField[i][2],2));
    if(mFieldTotal<1e-10)
    {
      continue;
    }
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

    //3:curvilinear accleration
    //This is a linac disign,so the R is infinite,so this step won't make any change


    //4:half-step accleration in the electric field again
    for(int j=0;j<3;++j)
    {
      particle_xy[i][2*j+1]+=eField[i][j]*timeStep/2;//unit:m/s
    }
  }
}

int PIC2D::PIC2D_PushPosition()
{
#pragma omp for
  for(int i=0;i<particlenumber;++i)
  {
    //5:advance the position
    particle_xy[i][0]+=particle_xy[i][1]/1*timeStep;//unit:m
    particle_xy[i][2]+=particle_xy[i][3]/1*timeStep;
    particle_xy[i][4]+=particle_xy[i][5]/1*timeStep;
    if(1)
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
  }
}
int PIC2D::CalEmit()
{
  {
    for(int i=0;i<particle_xy.size();++i)
    {
      beam.setparticle(i,	
	  particle_xy[i][0],
	  particle_xy[i][1]/(Ebeta*C_light),
	  particle_xy[i][2],
	  particle_xy[i][3]/(Ebeta*C_light),
	  particle_xy[i][4],
	  particle_xy[i][5]/(Ebeta*C_light)
	  );
    }
    beam.caculate_emittance();
    //cout<<"emittancex=   "<<beam.getemittancex()<<endl;
    
    double betaX   = beam.getsigmax()/beam.getemittancex();
    phaseX        += 1.0/betaX * 0.001 /M_PI*180; //*length
    double betaY   = beam.getsigmay()/beam.getemittancey();
    phaseY        += 1.0/betaY * 0.001 /M_PI*180; //*length
    ParaTransfer.push_back({
      particle_xy[0][4],
      beam.getemittancex(),
      sqrt(beam.getsigmax()),
      beam.getsigmaxdx()/sqrt(beam.getsigmax()),
      beam.getemittancey(),
      sqrt(beam.getsigmay()),
      beam.getsigmaydy()/sqrt(beam.getsigmay()),
      phaseX,
      phaseY});
  }
}

int PIC2D::CalMatchR()
{

  double K=	2*pow(Macrocharge,2)*particlenumber
    /
    (  pow(Egamma*Ebeta,3) * Macromass * pow(C_light,2)  );

  double kz0=Macrocharge  *  MyLattice.getField( vector<double>(6,0) ,0).back()
    /
    (2*Egamma * Macromass);
  double kz1=Macrocharge  *  MyLattice.getField( vector<double>(6,0) ,0).back()
    /
    (2*M_PI*Egamma * Macromass);

  double u=	K/(2*kz0);
  double rb=	sqrt(  sqrt(pow(u,2)+1/kz0)  +  u  );
  double r=Egamma * Macromass*sqrt(particle_xy[0][1]*particle_xy[0][1]+particle_xy[0][3]*particle_xy[0][3])/(Macrocharge  *  MyLattice.getField( vector<double>(6,0),0 ).back());
  cout<<r<<"   "<<kz1<<endl;
}

