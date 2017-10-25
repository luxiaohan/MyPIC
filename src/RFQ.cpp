#include "RFQ.h"
#include "MyFunc.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <gsl/gsl_sf_bessel.h>

RFQ::~RFQ()
{
}
RFQ& RFQ::operator=(const RFQ &rh)
{
//do nothing
  
}


RFQ::RFQ(const string &path,double posit,double length,double frequency)
{
  position	=	posit;
  this->length  	=	length;
  this->frequency	=	frequency;
  syn_phase 		=	0;//M_PI/2;
  ifstream rfStream;
  rfStream.open(path);
  if(!rfStream.is_open()) cout<<"RFQ Structure path ERROR!"<<endl;
  string str;
  vector<string> strVec;
  vd vtemp(8);
  for(int i = 0;rfStream.peek()!=EOF;++i)
  {
    getline(rfStream,str);
    //cout<<str<<endl;
    StringSplit(str,strVec);
    if(i==0)
      cell_position.push_back(position);
    else
      cell_position.push_back(cell_position.back()+cell_length.back());
    cell_length.push_back(stod(	strVec[1])/100);
    cell_position_end.push_back(cell_length.back()+cell_position.back());
    UL.push_back(	stod(	strVec[2])*1e3);  //kV to V
    radius.push_back(	stod(	strVec[3])/100);  //cm to m
    aperture.push_back(	stod(	strVec[4])/100);
    mod.push_back(	stod(	strVec[5]));

    for(int j=0;j<8;j++)
    {
      vtemp[j]=stod(strVec.at(j+6));
    }
   /* 
    for(int j=2;j<8;j++)
    {
      vtemp[j]=0;
    }
    */
    /*
    if(i==10) 
    {
      cout<<"rfq Structure error on LINE : "<<i+1<<"\nThe term is "<<strVec.size()<<endl;
      for(int j=0;j<vtemp.size();++j)
      {
	cout<<vtemp.at(j)<<"\t"<<j<<endl;
      }
      int a;
      cin>>a;
    }
    */
    factor.push_back(vtemp);
  }
    //ofstream a("8item.txt");
    //printvector(a,factor);
  cout<<"RFQ cell number : "<<UL.size()<<endl;
  cout<<"RFQ length      : "<<cell_position_end.back()<<endl;
}

const vector<double> RFQ::getfield(const vector<double> &xyz,double time)
{
  double Er,Ez,Et,Ex,Ey;
  int    eleNum=	getEleNum(xyz[4]);
  if(eleNum==-1) return vector<double>{0,0,0,0,0,0};
  double theta	=	atan2(xyz[2],xyz[0]);
  double radiu	=	sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2]);
  double z	=	xyz[4]-cell_position[eleNum];
  double k	=	M_PI/cell_length[eleNum];
  if(radiu<ERR)
  {
    Er=0.0;
    Et=0.0;
    Ez=-abs(UL[eleNum])/2.0*(
	factor[eleNum][1]*sin(k*z)*k+
	factor[eleNum][6]*sin(3*k*z)*3*k)*
      pow(-1.0,eleNum+1);
    //cout<<"z  = "<<z<<"  XYZ: "<<xyz[4]<<"  eleNum :"<<eleNum<<" cellposi : "<<cell_position[eleNum]<<endl;
      //cout<<"Ez = "<<Ez<<endl;
  }
  else
  {
    double b41,b43,b22,b62;
    b41=gsl_sf_bessel_In(4,  k*radiu);
    b43=gsl_sf_bessel_In(4,3*k*radiu);
    b22=gsl_sf_bessel_In(2,2*k*radiu);
    b62=gsl_sf_bessel_In(6,2*k*radiu);
    Er =   abs(UL[eleNum])/2.0* ( 
	factor[eleNum][0]*2*radiu		/pow(radius[eleNum],2)*cos(2*theta)+
      	factor[eleNum][2]*6*pow(radiu,5) 	/pow(radius[eleNum],6)*cos(6*theta)+
	factor[eleNum][1]*gsl_sf_bessel_In(1,  k*radiu)  *k*cos(  k*z)*pow(-1.0,eleNum+1)+
	factor[eleNum][6]*gsl_sf_bessel_In(1,3*k*radiu)*3*k*cos(3*k*z)*pow(-1.0,eleNum+1)+
	factor[eleNum][3]*gsl_sf_bessel_In(5,  k*radiu)  *k*cos(  k*z)*cos(4*theta)*pow(-1.0,eleNum+1)+
	factor[eleNum][7]*gsl_sf_bessel_In(5,3*k*radiu)*3*k*cos(3*k*z)*cos(4*theta)*pow(-1.0,eleNum+1)+
	factor[eleNum][4]*gsl_sf_bessel_In(3,2*k*radiu)*2*k*cos(2*k*z)*cos(2*theta)*pow(-1.0,eleNum+1)+
	factor[eleNum][5]*gsl_sf_bessel_In(7,2*k*radiu)*2*k*cos(2*k*z)*cos(6*theta)*pow(-1.0,eleNum+1));
    Et =  -abs(UL[eleNum])/2.0* (
	factor[eleNum][0]*pow(radiu/radius[eleNum],2)*sin(2*theta)*2 +
	factor[eleNum][2]*pow(radiu/radius[eleNum],6)*sin(6*theta)*6 +
	factor[eleNum][3]*b41*cos(  k*z)   *sin(4*theta)*4*pow(-1.0,eleNum+1)+
	factor[eleNum][7]*b43*cos(3*k*z)   *sin(4*theta)*4*pow(-1.0,eleNum+1)+
	factor[eleNum][4]*b22*cos(2*k*z)   *sin(2*theta)*2*pow(-1.0,eleNum+1)+
	factor[eleNum][5]*b62*cos(2*k*z)   *sin(6*theta)*6*pow(-1.0,eleNum+1))
        /radiu;
    Ez=   -abs(UL[eleNum])/2.0* (
	factor[eleNum][1]*gsl_sf_bessel_In(0,  k*radiu)*sin(  k*z)  *k                +
	factor[eleNum][6]*gsl_sf_bessel_In(0,3*k*radiu)*sin(3*k*z)*3*k                +
	factor[eleNum][3]*b41*sin(  k*z)  *k *cos(4* theta)+ 
	factor[eleNum][7]*b43*sin(3*k*z)*3*k *cos(4* theta)+
	factor[eleNum][4]*b22*sin(2*k*z)*2*k *cos(2* theta)+
	factor[eleNum][5]*b62*sin(2*k*z)*2*k *cos(6* theta))*pow(-1.0,eleNum+1);
  }
  double phaseCos=cos(2*M_PI*frequency*time+syn_phase);
  //cout<<frequency<<"   "<<time<<"   \n";
  //cout<< "  cos   "<<phaseCos<<"   "<<2*M_PI*frequency*time+syn_phase<<endl;
  Ex=Er*cos(theta)-Et*sin(theta);
  Ex*=phaseCos;
  Ey=Er*sin(theta)+Et*cos(theta);
  Ey*=phaseCos;
  Ez*=phaseCos;
  //if(abs(Ez)>1) cout<<radiu<<"  "<<k<<" "<<Ex<<"  "<<Ey<<"  "<<Ez<<endl;
  return vector<double>{Ex,Ey,Ez,0.0,0.0,0.0};
}

double RFQ::getradius(double z)
{
  int eleNum=	getEleNum(z);
  //cout<<aperture[eleNum]<<endl;
  if(eleNum==-1)  return 1e99;
  else		  return aperture[eleNum];
}

int RFQ::getEleNum(double z)
{
  int count = 0;
  while(count<cell_position.size()
      &&z > cell_position_end[count])
  {
    ++count;
  }
  if(count==cell_position.size())
  {
    //cout<<"Particle in RFQ error: exceed the RFQ length but still in RFQ"<<endl;
    count=-1;
  }
  return count;
}


