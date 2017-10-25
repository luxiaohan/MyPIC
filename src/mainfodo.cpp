#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <fftw3.h>
#include "Pic2d.h"
#include "Pic3d.h"
using namespace std;

const double epsilon=8.854187818e-12;//unit: F/m
const double macrocharge=1.602176565e-19,macromass=1.672621777e-27;

int main(int argc,char** argv)
{
  const int numofgridx=256;
  const int numofgridy=256;
  double KneticE=0.15;	//MeV
  double MSta=938.282046;
  //double Egamma=1.0005,     Ebeta=sqrt(1-1.0/Egamma/Egamma);

  double Egamma=1.0+KneticE/MSta,     Ebeta=sqrt(1-1.0/Egamma/Egamma);
  double speed=Ebeta*C_light;
  double elelength=0.2;
  double err=-1e-8;
  double slice=200*(1+err);
  double timeNow=0;
  double steptime=elelength/slice/speed;
  cout<<"TimeStep=  "<<steptime<<endl;
  cout<<"stepLength="<<steptime*speed<<endl;
  double totaltime=1e-5;
  int stepnumber=totaltime/steptime;
  
  double emitx=1e-3;
  double emity=1e-3;
  
  double rxinitial=0.0186601000307;
  double ryinitial=0.00795824136235;
  
  int particlenumber=40000;
  double current = 0.32; //A
  double qDensity = current / speed / particlenumber;
  double huge = qDensity/1.602176565e-19;
  PIC2D acc1(particlenumber,numofgridx,numofgridy,steptime,Egamma,huge,slice);
  //particle number,  grid number,  timestep,  Initial gamma,   Huge number
  srand(234);
  

  acc1.ReadParticle(vector<double>{
      emitx*speed,  pow(rxinitial,2)/emitx/speed,	0, 
      emity*speed,  pow(ryinitial,2)/emity/speed,	0});
      
  elelength = acc1.ReadLattice("latticefodo.txt");
  char addresspout[99];
  char addressmout[99];

  acc1.PIC2D_ICEFROG();
  cout<<setiosflags(ios::left);
#pragma omp parallel num_threads(4) 
  {
    for(int i=0;i<stepnumber;i++)
    {
      acc1.PIC2D_SpaceCharge();
      acc1.PIC2D_PUSH(timeNow);
#pragma omp single
      {
	timeNow+=steptime;
	cout<<i<<"/"<<stepnumber<<"   "
	<<setw(10)<<acc1.tww().back().at(0)<<"   "
	<<setw(10)<<acc1.tww().back().at(1)<<"  "
	<<setw(10)<<1.0/sqrt(1-acc1.tww().back().at(1)*acc1.tww().back().at(1)/C_light/C_light)*MSta-MSta<<"MeV  "<<endl;
	if(i%100==0)
	{
	  sprintf(addressmout,"data/Ptc%06d.dat",i);
	  ofstream pout(addressmout);
	  if(pout)  printvector(pout,acc1.xyz());
	  pout.close();
	}
      }
      if(acc1.tww().back()[0]>elelength) 
      {
	//cout<<acc1.tww().back()[0]<<endl;
	break;
      }
    }
  }
  acc1.PIC2D_FIREFROG();
  //acc1.CalMatchR();
  string addout="x2_"+to_string(stepnumber)+".dat";
  ofstream twwout(addout);
  printvector(twwout,acc1.tww());

  //
  /*
     ofstream xyout("acc1.xyz().dat");
     for(int i=0;i<particlenumber;i++)
     {
     xyout<<setw(10)<<acc1.xyz()[i][0]<<"   "<<setw(10)<<acc1.xyz()[i][1]<<"   "<<setw(10)<<acc1.xyz()[i][2]<<"   "<<setw(10)<<acc1.xyz()[i][3]<<endl;
     }
     gout.close();
     */


  
  //benchmark rho
  /*
  vector<vector<double>  >   rho_test(numofgridx,   vector<double>(numofgridy));
  ofstream ggout("rhobench.dat");
  double sumrho=0;
  for(int j=1;j<numofgridy-1;j++)
  {
    for(int i=1;i<numofgridy-1;i++)
    {
      rho_test[i][j]=-(
                        (acc1.phi()[i-1][j]-2*acc1.phi()[i][j]+acc1.phi()[i+1][j])/1
                       +(acc1.phi()[i][j-1]-2*acc1.phi()[i][j]+acc1.phi()[i][j+1])/1
                       )*epsilon;
      ggout<<setw(10)<<rho_test[i][j]<<"   ";
      sumrho+=rho_test[i][j];
    }
    ggout<<endl;
  }
  cout<<"sum rho bench=  "<<sumrho<<endl;
  ggout.close();
  */
/*
  //benchmark phi ana
  vector<vector<double>  >   phibench(numofgridx,   vector<double>(numofgridy));
  vector<vector<double>  >   exbench(numofgridx,   vector<double>(numofgridy));
  vector<vector<double>  >   eybench(numofgridx,   vector<double>(numofgridy));
  ofstream phibenchout("phibench.dat");
  double x,y,x1,y1;
  double centerlengthx=numofgridx/2-0.5;
  double centerlengthy=numofgridy/2-0.5;
  for(int j=0;j<numofgridy;j++)
  {
  for(int i=0;i<numofgridx;i++)
  {
  x=double(i)-numofgridx/2+0.5;
  y=double(j)-numofgridy/2+0.5;

  phibench[i][j]=0;
  exbench[i][j]=0;
  eybench[i][j]=0;
  for(int m=-1;m<2;m++)
  for(int n=-1;n<2;n++)
  {
  x1=x+m*(numofgridx+2);
  y1=y+n*(numofgridy+2);
  x1=x1*0.001;
  y1=y1*0.001;
  if((abs(m)+abs(n))%2==0)
  {
  phibench[i][j]+=-0.5/M_PI/epsilon*10000*macrocharge*log(sqrt(x1*x1+y1*y1));
  exbench[i][j]+=0.5/M_PI/epsilon*10000*macrocharge/(x1*x1+y1*y1)*x1;
  eybench[i][j]+=0.5/M_PI/epsilon*10000*macrocharge/(x1*x1+y1*y1)*y1;
  //cout<<x1*x1+y1*y1<<"   "<<log(sqrt(x1*x1+y1*y1))<<"   "<<phibench[i][j]<<endl;
  }
  else
  {
  phibench[i][j]-=-0.5/M_PI/epsilon*10000*macrocharge*log(sqrt(x1*x1+y1*y1));
  exbench[i][j]+=0.5/M_PI/epsilon*10000*macrocharge/(x1*x1+y1*y1)*x1;
  eybench[i][j]+=0.5/M_PI/epsilon*10000*macrocharge/(x1*x1+y1*y1)*y1;
  }
  }
  phibenchout<<setw(10)<<(phibench[i][j])<<"   ";
  }
  phibenchout<<endl;
  }
  phibenchout.close();

  //benchmark phi
  ofstream gout("phi.dat");
  for(int j=0;j<numofgridy;j++)
  {
  for(int i=0;i<numofgridx;i++)
  {
  gout<<setw(10)<<acc1.phi()[i][j]<<"   ";
  }
  gout<<endl;
}
gout.close();

ofstream conout("contrast.dat");
for(int j=numofgridy/3;j<2*numofgridy/3;j++)
{
  for(int i=numofgridx/3;i<2*numofgridx/3;i++)
  {
    conout<<setw(10)<<acc1.ex()[i][j]-exbench[i][j]<<"   ";
  }
  conout<<endl;
}
conout.close();
*/
cout<<"COMPLETE!"<<endl;
}
