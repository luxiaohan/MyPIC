#include "Field3D.h"
#include "MyFunc.h"
#include <fstream>

Field3D::~Field3D(){}
Field3D::Field3D(const string &path,double position,double length,double phase,double bfactor,double efactor,int type)
{
  filename=path;
  this->location=	position;
  this->length	=	length;
  this->type	=	type;
  ke		=	efactor;
  kb		=	bfactor;
  locationEnd	=	position+length;
  synPhase	=	phase/180*M_PI;
  frq		=	325e6;

  ifstream fin;
  string 		str;
  vector<string> 	strVec;
  if(type==0)
  {
    fin.open(path);
    if(!fin.is_open())
    {
      cout<<"field file error : "<<path<<endl;
      return;
    }
    vector<double>	xyzTemp(3);
    vector<double>	eFieldTemp(3);
    vector<double>	mFieldTemp(3);
    for(int i = 0;i<2;++i)
    {
      //to skip several lines
      getline(fin,str);
    }
    for(int i = 0;fin.peek()!=EOF;++i)
    {
      getline(fin,str);
      //cout<<str<<endl;
      StringSplit(str,strVec);
      xyzTemp[0]	=	stod(	strVec[0])	/1000;		//mm to m
      xyzTemp[1]	=	stod(	strVec[1])	/1000;
      xyzTemp[2]	=	stod(	strVec[2])	/1000	+position;
      gridLocation.push_back(xyzTemp);

      eFieldTemp[0]=	stod(	strVec[3])*ke	/1000	;
      eFieldTemp[1]=	stod(	strVec[4])*ke	/1000;
      eFieldTemp[2]=	stod(	strVec[5])*ke	/1000;
      eField.push_back(eFieldTemp);

      mFieldTemp[0]=	stod(	strVec[6])*kb	/1000;
      mFieldTemp[1]=	stod(	strVec[7])*kb	/1000;
      mFieldTemp[2]=	stod(	strVec[8])*kb	/1000;
      mField.push_back(mFieldTemp);
    }
    fin.close();
  }
  else if(type == 1)
  {
    vector<string> filedex={".edx",".edy",".edz",".bdx",".bdy",".bdz"};
    for(int dire=0;dire<6;dire++)
    {
      string filepath=path+filedex[dire];
      fin.open(filepath);
      if(!fin.is_open())
      {
	cout<<"field type1 file error : "<<filepath<<endl;
	return;
      }
      if(dire==0)
      {
	fin>>nz>>zmax;
	fin>>nx>>xmin>>xmax;
	fin>>ny>>ymin>>ymax;
	Field.resize((nx+1)*(ny+1)*(nz+1),vector<double>(6));
	locationMax=location+zmax;
	spacex=(xmax-xmin)/nx;
	spacey=(ymax-ymin)/ny;
	spacez=(zmax	 )/nz;
	getline(fin,str);
	getline(fin,str);
      }
      else
      {
	for(int i = 0;i<4;++i)
	{
	  //to skip several lines
	  getline(fin,str);
	}
      }
      if(stoi(str)!=1)	cout<<"mark error : "<<filepath<<"  "<<str<<endl;


      double temp=0;
      int index=0;
      for(int k=0;k<=nz;++k)
      {
	for(int j=0;j<=ny;++j)
	{
	  for(int i=0;i<=nx;++i)
	  {
	    fin>>temp;

	    index=k*(ny+1)*(nx+1)+j*(nx+1)+i;
	    Field[index][dire]=temp;
	  }
	}
      }
      fin.close();
      gridNum=Field.size();

    }
    for(int i=0;i<gridNum;++i)
    {
      for(int j=0;j<3;++j)
      {
	Field[i][j]=Field[i][j]*1.0e6*ke;
      }
      for(int j=3;j<6;++j)
      {
	Field[i][j]=Field[i][j]*kb;
      }
    }
    string fout=filename+"_fieldout"+".dat";
    ofstream fieldout(fout);

    for(int k=0;k<=nz;++k)
    {
      for(int j=0;j<=ny;++j)
      {
	for(int i=0;i<=nx;++i)
	{
	  fieldout<<i*spacex+xmin<<"  ";
	  fieldout<<j*spacey+ymin<<"  ";
	  fieldout<<k*spacez<<"  ";
	  double index=k*(ny+1)*(nx+1)+j*(nx+1)+i;
	  for(int dire=0;dire<6;dire++)
	  {
	    fieldout<<Field[index][dire]<<" \t";
	  }
	  fieldout<<endl;
	}
      }
    }
  }
  else if(type==2)
  {
    vector<string> filedex={".bsz",".bsr"};
    for(int dire=0;dire<filedex.size();dire++)
    {
      string filepath=path+filedex[dire];
      fin.open(filepath);
      if(!fin.is_open())
      {
	cout<<"field type2 file error : "<<filepath<<endl;
	int a;
	cin>>a;
	return;
      }
      if(dire==0)
      {
	fin>>nz>>zmax;
	fin>>nr>>rmax;
	Field.resize((nz+1)*(nr+1),vector<double>(2));
	locationMax=location+zmax;
	spacer=(rmax	 )/nr;
	spacez=(zmax	 )/nz;
	getline(fin,str);
	getline(fin,str);
      }
      else
      {
	for(int i = 0;i<3;++i)
	{
	  //to skip several lines
	  getline(fin,str);
	}
      }
      if(stoi(str)!=1)	cout<<"mark error : "<<filepath<<"  "<<str<<endl;


      double temp=0;
      int index=0;
      for(int j=0;j<=nz;++j)
      {
	for(int i=0;i<=nr;++i)
	{
	  fin>>temp;
	  index=j*(nr+1)+i;
	  Field[index][dire]=temp*kb;
	}
      }
      fin.close();
      gridNum=Field.size();
    }
    string fout=filename+"_fieldout"+".dat";
    ofstream fieldout(fout);

    for(int j=0;j<=nz;++j)
    {
      for(int i=0;i<=nr;++i)
      {
	fieldout<<i*spacer<<"  ";
	fieldout<<j*spacez<<"  ";
	double index=j*(nr+1)+i;
	for(int dire=0;dire<2;dire++)
	  fieldout<<Field[index][dire]<<"  ";
	fieldout<<endl;
      }
    }
    /*
       ofstream fieldout("fieldout.dat");

       for(int k=0;k<=nz;++k)
       {
       for(int j=0;j<=ny;++j)
       {
       for(int i=0;i<=nx;++i)
       {
       fieldout<<i*spacex+xmin<<"  ";
       fieldout<<j*spacey+ymin<<"  ";
       fieldout<<k*spacez<<"  ";
       double index=k*(ny+1)*(nx+1)+j*(nx+1)+i;
       for(int dire=0;dire<6;dire++)
       fieldout<<Field[index][dire]<<"  ";
       fieldout<<endl;
       }
       }
       }

*/
  }
  else
  {
    cout<<"Field type error!"<<endl;
  }
}

const vector<double> Field3D::getfield(const vector<double> &xyz,double time)
{
  if(type==0)
  {
    if(	xyz[0]<gridLocation[0][0] || xyz[0]>gridLocation.back()[0]
	||xyz[1]<gridLocation[0][1] || xyz[1]>gridLocation.back()[1]
	||xyz[2]<gridLocation[0][2] || xyz[2]>gridLocation.back()[2]
      )
    {
      cout<<"particle out of range at Field";
      for(int i=0;i<xyz.size();++i)
      {
	cout<<"  "<<xyz[i];
      }
      cout<<endl;
      return {0,0,0,0,0,0};
    }
    //if the particle is out of range , return 0


    vector<int> count(3,1);
    for(int i=0;i<3;++i)
    {
      while(xyz[2*i]>gridLocation[count[i]][i])
      {
	++count[i];
      }
    }
    //get the location where the particle is



    for(int i=0;i<3;++i)
    {
      d1=xyz[2*i]				-gridLocation[count[i]-1][i];
      d2=gridLocation[count[i]][i]	-xyz[2*i];
      if(d1<=0||d2<=0)
      {
	cout<<"particle count error in Field3D"<<endl;
	return zero;
      }
      fieldofParticle[2*i]	=	(eField[count[i]-1][i]*d2	+eField[count[i]][i]*d1)
	/(d1+d2);
      fieldofParticle[2*i+1]	=	(mField[count[i]-1][i]*d2	+mField[count[i]][i]*d1)
	/(d1+d2);
    }
    return fieldofParticle;
  }
  else if(type==1)
  {
    if(	xyz[0]<xmin  || xyz[0]>xmax
	||xyz[2]<ymin|| xyz[2]>ymax
	||xyz[4]<location || xyz[4]>locationEnd)
    {
      if(0)
      {
      cout<<"particle out of range at Field3D type1 "<<filename;
      for(int i=0;i<xyz.size();++i)
      {
	cout<<"  "<<xyz[i];
      }
      cout<<"     "<<location<<" and "<<locationEnd;
      cout<<endl;
      }
      return zero;
    }
    else if(xyz[4]>locationMax)
    {
      return zero;
    }

    double xn,yn,zn;
    int xGrid,yGrid,zGrid,indexGrid;
    vd v(8),fieldtemp(6);
    v2d fieldOnGrid(6,vd(8));
    fieldtemp=zero; 
    xn=(xyz[0]-xmin)	/spacex;
    yn=(xyz[2]-ymin)	/spacey;
    zn=(xyz[4]-location)/spacez;
    for(int i=0;i<8;++i)
    {
      xGrid=((i%2< 1)?ceil(xn):floor(xn));
      yGrid=((i%4< 2)?ceil(yn):floor(yn));
      zGrid=((i  < 4)?ceil(zn):floor(zn));
      v[i]=abs(xn-xGrid) * abs(yn-yGrid) * abs(zn-zGrid);

      xGrid=((i%2==1)?ceil(xn):floor(xn));
      yGrid=((i%4>=2)?ceil(yn):floor(yn));
      zGrid=((i  >=4)?ceil(zn):floor(zn));
      indexGrid=zGrid*(ny+1)*(nx+1)+yGrid*(nx+1)+xGrid;
#pragma unroll
      for(int dim=0;dim<6;++dim)
      {
	fieldOnGrid[dim][i]=v[i]*Field[indexGrid][dim];
	fieldtemp[dim]+=fieldOnGrid[dim][i];
      }
    }
    double pha=2*M_PI*frq*time+synPhase+offsetPhase;
    double phaseCos=cos(2*M_PI*frq*time+synPhase+offsetPhase);
    for(int dim=0;dim<6;++dim)
    {
      fieldtemp[dim]*=phaseCos;
    }
    return fieldtemp;
  }
  else if(type==2)
  {

    double r	=sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2]);
    double theta=atan2(xyz[2],xyz[0]);
    if(	r>rmax
	||xyz[4]<location || xyz[4]>locationEnd)
    {
      cout<<"particle out of range at Field3D type2  "<<filename;
      for(int i=0;i<xyz.size();++i)
      {
	cout<<"  "<<xyz[i];
      }
      cout<<"     "<<location<<" and "<<locationEnd;
      cout<<endl;
      return zero;
    }
    else if(xyz[4]>locationMax)
    {
      return zero;
    }

    double rn,zn;
    int rGrid,zGrid,indexGrid;
    vd v(4),fieldtemp(6,0.0);
    v2d fieldOnGrid(2,vd(4));
    fieldtemp=zero; 
    vd fieldzr(2,0.0); 
    rn=r		/spacer;
    zn=(xyz[4]-location)/spacez;
    for(int i=0;i<4;++i)
    {
      rGrid=((i%2< 1)?ceil(rn):floor(rn));
      zGrid=((i  < 2)?ceil(zn):floor(zn));
      v[i]=abs(rn-rGrid) * abs(zn-zGrid);

      rGrid=((i%2==1)?ceil(rn):floor(rn));
      zGrid=((i  >=2)?ceil(zn):floor(zn));
      indexGrid=zGrid*(nr+1)+rGrid;
#pragma unroll
      for(int dim=0;dim<2;++dim)
      {
	fieldOnGrid[dim][i]=v[i]*Field[indexGrid][dim];
	fieldzr[dim]+=fieldOnGrid[dim][i];
      }
    }

    fieldtemp[3]=fieldzr[1]*cos(theta);
    fieldtemp[4]=fieldzr[1]*sin(theta);
    fieldtemp[5]=fieldzr[0];
    return fieldtemp;
  }

}

void Field3D::synchro(double z_ave,double v_ave,double offset)
{
  double distance	=(location+length/2)-z_ave;
  double time		=distance/v_ave;
  offsetPhase	=-2*M_PI*frq*time;
  //  cout<<location<<"  "<<length<<"  "<<offsetPhase<<endl;
  offsetPhase+=offset/180.0*M_PI;

}
void Field3D::ScanOffset(double offset)
{
  offsetPhase=offset/180.0*M_PI-synPhase;
}
void Field3D::SetOffset(double offset)
{
  offsetPhase=offset/180.0*M_PI;
}

void Field3D::print(){}
void Field3D::Map(){}
