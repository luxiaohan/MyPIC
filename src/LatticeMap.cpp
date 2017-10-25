#include "LatticeMap.h"
#include "MyFunc.h"

#include <fstream>
#include <iostream>
#include <stdio.h> 
const double FloatTolerance=1.0e-8;
LatticeMap::~LatticeMap()
{
}
LatticeMap::LatticeMap()
{
  //cout<<"LatticeMap Created"<<endl;
}
LatticeMap::LatticeMap(const LatticeMap &e)
{
  ELE=e.ELE;
  //  cout<<"copy"<<endl;
  Mat=e.Mat;
  fitness=e.fitness;
  EBeta=e.EBeta;
  EGamma=e.EGamma;
  this->elemtemp=e.elemtemp;
  this->FQ=e.FQ;
  this->BM=e.BM;
  this->DL=e.DL;
  this->SN=e.SN;
  this->RQ=e.RQ;
  this->RM=e.RM;
  this->FIELD3D=e.FIELD3D;
  this->rfqtypepre=e.rfqtypepre;
  this->rfqtypenow=e.rfqtypenow;
  this->rfqtypenext=e.rfqtypenext;
  // this->infile=e.infile;
  this->NDivide=e.NDivide;
  this->Mfail=e.Mfail;
  this->particle_count=e.particle_count;
  this->totalparticalnumber=e.totalparticalnumber;
  this->i=e.i;
  this->j=e.j;
  this->k=e.k;
  this->myid=e.myid;
  this->numprocs=e.numprocs;
  this->positionpre=e.positionpre;
  this->positionnow=e.positionnow;
  this->positionnext=e.positionnext;

}

LatticeMap& LatticeMap::operator =(const LatticeMap& e)
{
  ELE=e.ELE;
  //  cout<<"assign"<<endl;
  Mat=e.Mat;
  fitness=e.fitness;
  EBeta=e.EBeta;
  EGamma=e.EGamma;
  this->elemtemp=e.elemtemp;
  this->FQ=e.FQ;
  this->BM=e.BM;
  this->DL=e.DL;
  this->SN=e.SN;
  this->RQ=e.RQ;
  this->RM=e.RM;
  this->FIELD3D=e.FIELD3D;
  this->rfqtypepre=e.rfqtypepre;
  this->rfqtypenow=e.rfqtypenow;
  this->rfqtypenext=e.rfqtypenext;
  // this->infile=e.infile;
  this->NDivide=e.NDivide;
  this->Mfail=e.Mfail;
  this->particle_count=e.particle_count;
  this->totalparticalnumber=e.totalparticalnumber;
  this->i=e.i;
  this->j=e.j;
  this->k=e.k;
  this->myid=e.myid;
  this->numprocs=e.numprocs;
  this->positionpre=e.positionpre;
  this->positionnow=e.positionnow;
  this->positionnext=e.positionnext;    
  return *this;
}

bool LatticeMap::operator==(const LatticeMap& l)
{
  if(ELE.size()!=l.ELE.size())
    return false;
  for(int i=0;i<ELE.size();i++)
  {
    if(ELE[i]!=l.ELE[i])
      return false;
  }
  return true;
}


void LatticeMap::callelement(Element &v,ifstream &infile)
{
  if(v.name=="drift")
  {
    DriftMap drift(v.p[2],EGamma.back());
    pushback(drift,  DL,  EGamma.back(),driftorder,v);
  }
  else if(v.name=="quadm")
  {
    QuadrupoleMap qf(v.p[2],v.p[4],EGamma.back());
    pushback(qf,  FQ,  EGamma.back(),quadmorder,v);
  }
  else if(v.name=="bendm")
  {   
    BendingMap bendm(v.p[2],v.p[4],EGamma.back());
    pushback(bendm,  BM,  EGamma.back(),bendmorder,v);
  }
  else if(v.name=="solen")
  {
    SolenoidMap solenoid(v.p[2],v.p[4],EGamma.back());
    pushback(solenoid,  SN,  EGamma.back(),solenorder,v);
  }
  else if(v.name=="field")
  {
    Field3D field(v.path,v.position,v.length,v.p[4],v.p[5],v.p[6],int(v.p[1]));
    pushback(field,  FIELD3D,  EGamma.back(),fieldorder,v);
  }
  else if(v.name=="gap")
  {
    //            if(myid==0)
    //              cout<<"BEFORE   "<<EGamma.back()<<endl;
    addEGamma((EGamma.back()*_MASS*C_light*C_light+v.p[4]*1e6*cos(v.p[5]*M_PI/180.0)*abs(_Q))/(_MASS*C_light*C_light));
    //            if(myid==0)
    //              cout<<"AFTER    "<<EGamma.back()<<endl;
    if(v.p[1]==1)
    {
      GapMap gap(v.p[2],v.p[4],v.p[5],v.p[6],EGamma.back());
      gap.Map();
      Mat.push_back(gap.getmap());
    }
    else if(v.p[1]==2)
    {
      char temp[99];
      strcpy(temp,v.path.c_str());
      //cout<<temp<<endl;
      GapMap gap(v.p[2],temp,v.p[5],v.p[6],EGamma.back());
      Mat.push_back(MatrixXd::Identity(6,6));
      gap.getEf();
    }
    else
    {
      cout<<"GAP TYPE error!"<<endl;
      Mat.push_back(MatrixXd::Identity(6,6));
    }
    //cout<<i<<"  drift \n"<<Mat[i]<<endl;
  }
  else if(v.name=="rfq")
  {
    if(v.p[1]==0)
    {
      RFQ rfq(v.path,v.position,v.p[2],v.p[4]);
      pushback(rfq,RQ,EGamma.back(),rfqorder,v);
    }
    else if(v.p[1]==1)
    {
      //(rfqce type length radius voltage phase )
      RFQMap rfq;
      rfq.clear();
      addEGamma(EGamma.back());
      rfq.initial(v.p[2],v.p[4],v.p[5],v.p[6],EGamma.back());
      string rfqtemp;
      infile.seekg(positionpre);
      getline(infile,rfqtemp);
      vector<string> parainput;
      if(rfqtemp.size()>5)
      {
	StringSplit(rfqtemp,parainput);
	if(parainput[0]=="rfq")
	{
	  rfqtypepre=atof(parainput[1].c_str());
	  //cout<<"typepre=\t  "<<para[1].c_str()<<endl;
	}
      }
      infile.seekg(positionnow);
      getline(infile,rfqtemp);
      if(rfqtemp.size()>5)
      {
	StringSplit(rfqtemp,parainput);
	if(parainput[0]=="rfq")
	{
	  rfqtypenow=atof(parainput[1].c_str());
	  //cout<<"typenow=\t  "<<para[1].c_str()<<endl;
	}
      }
      infile.seekg(positionnext);
      getline(infile,rfqtemp);
      if(rfqtemp.size()>5)
      {
	StringSplit(rfqtemp,parainput);

	if(parainput[0]=="rfq")
	{
	  rfqtypenext=atof(parainput[1].c_str());
	  //cout<<"typenext=\t  "<<para[1].c_str()<<endl;
	}
      }
      infile.seekg(positionnext);

      rfq.settype(rfqtypepre,rfqtypenow,rfqtypenext);
      rfq.Map();
      Mat.push_back(rfq.getmap());
      //cout<<i<<"  rfq \n"<<Mat[i]<<endl;
      //rfq.print();
    }
    else
    {
      cout<<"rfq error type"<<endl;
    }
  }
  else
  {
    ELE.pop_back();
    i--;
    Mfail++;
  }

}

void LatticeMap::modify(const int &num,const Element &v)
{
  elemtemp=ELE[num];
  if(num<Mat.size())
  {
    ELE[num]=v;
    if(v.name=="drift")
    {
      DriftMap drift;
      drift.setlength(v.p[2]);
      drift.setgamma(EGamma[num]);
      drift.Map();
      Mat[num]=drift.getmap();
    }
    else if(v.name=="quadm")
    {
      QuadrupoleMap qf;
      qf.setlength(v.p[2]);
      qf.setgradient(v.p[4]);
      qf.setgamma(EGamma[num]);
      qf.Map();
      Mat[num]=qf.getmap();
    }
    else if(v.name=="bendm")
    {
      BendingMap bendm;
      bendm.setlength(v.p[2]);
      bendm.setradius(v.p[4]);
      bendm.Map();
      Mat[num]=bendm.getmap();
    }
    else if(v.name=="solen")
    {
      SolenoidMap solenoid;
      solenoid.setlength(v.p[2]);
      solenoid.setfield(v.p[4]);
      solenoid.Map();
      Mat[num]=solenoid.getmap();
    }
    else if(v.name=="rfq")
    {/*
      //(rfqce type length radius voltage phase )
      RFQMap rfq;
      rfq.clear();
      rfq.initial(v.p[2],v.p[4],v.p[5],v.p[6],EGamma);
      string rfqtemp;
      infile.seekg(positionpre);
      getline(infile,rfqtemp);
      vector<string> parainput;
      if(rfqtemp.size()>5)
      {
      StringSplit(rfqtemp,parainput);
      if(parainput[0]=="rfqce")
      {
      rfqtypepre=atof(parainput[1].c_str());
      //cout<<"typepre=\t  "<<para[1].c_str()<<endl;
      }
      }
      infile.seekg(positionnow);
      getline(infile,rfqtemp);
      if(rfqtemp.size()>5)
      {
      StringSplit(rfqtemp,parainput);
      if(parainput[0]=="rfqce")
      {
      rfqtypenow=atof(parainput[1].c_str());
      //cout<<"typenow=\t  "<<para[1].c_str()<<endl;
      }
      }
      infile.seekg(positionnext);
      getline(infile,rfqtemp);
      if(rfqtemp.size()>5)
      {
      StringSplit(rfqtemp,parainput);

      if(parainput[0]=="rfqce")
      {
      rfqtypenext=atof(parainput[1].c_str());
      //cout<<"typenext=\t  "<<para[1].c_str()<<endl;
      }
      }
      infile.seekg(positionnext);

      rfq.settype(rfqtypepre,rfqtypenow,rfqtypenext);
      rfq.Map();
      Mat[num]=rfq.getmap();
      */
    }
    else
    {
      ELE[num]=elemtemp;
    }
  }
}
double LatticeMap::read(const char *p,double EGam,int numberdivided,int numberrepeat)
{
  if(numberdivided<1) {cout<<"Number of Divided ERROR!"<<endl; return 0;}
  addEGamma(EGam);
  ifstream infile;
  ofstream outlat("latticereadin.dat");
  NDivide=numberdivided;
  string line;
  positionpre=0;
  infile.open(p);
  i=0;
  Mfail=0;
  vector<string> parainput;
  double posi=0;
  for(int repeat=0;repeat<numberrepeat;repeat++)
  {
    infile.seekg(0,infile.beg);
    while(infile.peek()!=EOF){
      positionpre=positionnow;
      positionnow=infile.tellg();
      getline(infile,line);
      positionnext=infile.tellg();
      if(line.size()>1)
      {
	for(int divideCount=0;divideCount<numberdivided;++divideCount)
	{
	  StringSplit(line,parainput);
	  if(parainput.size()==0) continue;
	  elemtemp.name=parainput[0];
	  elemtemp.p.clear();
	  elemtemp.p.push_back(0);
	  elemtemp.length=stod(parainput[2])/numberdivided;
	  elemtemp.position=posi;
	  posi+=elemtemp.length;
	  for(int r=1; r<parainput.size()&&r<12;r++)
	  {
	    if(elemtemp.name=="rfq"
		||elemtemp.name=="gap")
	    {
	      if(r==2)
	      {
		elemtemp.p.push_back(elemtemp.length);
	      }
	      else if(r==5)
	      {
		elemtemp.path=parainput[5];
		elemtemp.p.push_back(0);
	      }
	      else
	      {
		elemtemp.p.push_back(stod(parainput[r]));
	      }
	    }
	    else if(elemtemp.name=="field")
	    {
	      if(r==2)
	      {
		elemtemp.p.push_back(elemtemp.length);
	      }
	      else if(r==7)
	      {
		elemtemp.path=parainput[7];
		elemtemp.p.push_back(0);
	      }
	      else
	      {
		elemtemp.p.push_back(stod(parainput[r]));
	      }
	    }
	    else
	    {
	      if(r==2)
		elemtemp.p.push_back(elemtemp.length);
	      else
		elemtemp.p.push_back(stod(parainput[r]));
	    }
	  }
	  ELE.push_back(elemtemp);
	  //cout<<ELE.back().name<<endl;
	  callelement(ELE.back(),infile);
	  outlat<<elemtemp.name<<"\t ";
	  for(int r=1; r<parainput.size()&&r<12;r++)
	  {
	    if(r==5)
	      outlat<<parainput[5]<<"\t ";
	    else
	      outlat<<elemtemp.p.at(r)<<"\t ";
	  }
	  outlat<<EGamma.back()<<endl;
	}
      }else 
      {
	i--;
	Mfail++;
      }
      i++;
    }//while end

  }
  latticeLength=posi;
    cout<<"Lattice file    : "<<p<<endl;
    cout<<"Lattice length  : "<<posi<<endl;
    cout<<"Element number  : "<<Mat.size()<<endl;
  infile.close();
  return posi;
  //remove(lp);
}

void LatticeMap::prelattice(const char *in,const char *out,int n)
{
  if(myid==0)
  {
    ifstream infilep(in);
    ofstream outfile(out);
    string line;
    vector<string> parainput;
    while(infilep.peek()!=EOF)
    {
      getline(infilep,line);
      string str = line;
      remove(str.begin(), str.end(), ' ');
      remove(str.begin(), str.end(), '\t');
      if(str.size()>1)
      {
	StringSplit(line,parainput);
	double length_prelattice = atof(parainput[2].c_str());
	if((!(parainput[0]=="gap"))&&length_prelattice>0)
	{
	  for(int lnum=0;lnum<n;lnum++)
	  {
	    for(int i1=0;i1<parainput.size();i1++)
	    {
	      if(i1!=2)
		outfile<<parainput[i1]<<"  ";
	      else
		outfile<<length_prelattice/n<<"  ";
	    }
	    outfile<<endl;
	  }
	}
	else
	{
	  for(int i1=0;i1<parainput.size();i1++)
	  {
	    outfile<<parainput[i1]<<"  ";
	  }
	  outfile<<endl;
	}
      }
    }
  }
}

void LatticeMap::print()
{
  for(int ii=0;ii<ELE.size();ii++)
  {
    cout<<ii<<"/"<< ELE.size()<<"  ";
    cout<<ELE[ii].name<<"  ";
    for(int i2=1;i2<ELE[ii].p.size();i2++)
      cout<<ELE[ii].p[i2]<<"  ";
    cout<<endl;
    //cout<<Mat[ii]<<endl;
  }
}
void LatticeMap::addEGamma(const double num)
{
  if(num<1)
  {
    cout<<"EGamma Wrong in lattice"<<endl;
    EGamma.push_back(1);
    EBeta.push_back(0);
    return;
  }
  EGamma.push_back(num);
  EBeta.push_back(sqrt(1-1/num/num));
}

MatrixXd LatticeMap::getMat(int L,double &x,double &dx,double &y,double &dy,double &z,double &dp)
{
  if(ELE[L].name=="gap")
  {
    if(ELE[L].p[1]==1)
    {
      GapMap gap(ELE[L].p[2],ELE[L].p[4],ELE[L].p[5],ELE[L].p[6],EGamma[L]);
      return(gap.getmap(x,dx,y,dy,z,dp));
    }
    else if(ELE[L].p[1]==2)
    {
      char temp[99];
      strcpy(temp,ELE[L].path.c_str());
      GapMap gap(ELE[L].p[2],temp,ELE[L].p[5],ELE[L].p[6],EGamma[L]);
      //cout<<z<<" BEF "<<dp<<endl;
      gap.FieldMethod(x,dx,y,dy,z,dp);
      //cout<<z<<" AFT "<<dp<<endl;
      return MatrixXd::Identity(6,6);
    }
    else
    {
      cout<<"GAP TYPE error!"<<endl;
      return MatrixXd::Identity(6,6);
    }
  }
  else
    return Mat[L];
}
//20151106
const vector<double> LatticeMap::getFieldPre(const vector<double> &xyz,double time)
{ 
  return getFieldBase(xyz,FloatTolerance,time);
}
const vector<double> LatticeMap::getField(const vector<double> &xyz,double time)
{
  return getFieldBase(xyz,-FloatTolerance,time);
}
const vector<double> LatticeMap::getFieldBase(const vector<double> &xyz,double tolerance,double time)
{
  int count=0;

  if(xyz[4]<0||xyz[4]>ELE.back().position+ELE.back().length)
    return vector<double>{0,0,0,0,0,0};
  double particleposition=fmod(xyz.at(4),(ELE.back().position+ELE.back().length));
  if(abs( particleposition-  (ELE.back().position+ELE.back().length) )  < abs(tolerance)) particleposition=0;
  while(particleposition  - ( ELE[count].position+ELE[count].length ) >=tolerance)
  {
    ++count;
  }
  if(ELE[count].name=="drift")
  {
    return DL[ELE[count].ordernumber].getfield(xyz);
  }
  else if(ELE[count].name=="quadm")
  {
    return FQ[ELE[count].ordernumber].getfield(xyz);
  }
  else if(ELE[count].name=="bendm")
  {
    return BM[ELE[count].ordernumber].getfield(xyz);
  }
  else if(ELE[count].name=="solen")
  {
    return SN[ELE[count].ordernumber].getfield(xyz);
  }
  else if(ELE[count].name=="rfq")
  {
    return RQ[ELE[count].ordernumber].getfield(xyz,time);
  }
  else if(ELE[count].name=="field")
  {
    return FIELD3D[ELE[count].ordernumber].getfield(xyz,time);
  }
  else
  {
    return {0,0,0,0,0,0};
  }
}
const int LatticeMap::getEleNum(double z)
{
  int count=0;
  double particleposition=fmod(z,ELE.back().position+ELE.back().length);
  if(abs( particleposition-  (ELE.back().position+ELE.back().length) )  < FloatTolerance) particleposition=0;
  while(particleposition - ( ELE[count].position+ELE[count].length ) >=-FloatTolerance)
  {
    ++count;
  }
  /*
  if(ELE[count].name=="field") 
    cout<<"IN getEleNum   "<<count<<"  "<<z<<"   "<<ELE[count].position<<"  "<<ELE[count].position+ELE[count].length<<"  "<<particleposition<<endl;
    */
  return count;
}
double LatticeMap::getRadius(double z)
{
  int N=getEleNum(z);
  if(ELE[N].name=="rfq")
  {
    return RQ[ELE[N].ordernumber].getradius(z);
  }
  else
  {
    return ELE[N].p[3];
  }
}
  template <typename T1,typename T2>
void LatticeMap::pushback(const T1 &ele,T2 &vecEle,const double ene,int &order,Element &struc)
{
  addEGamma(ene);
  Mat.push_back(ele.getmap());
  vecEle.push_back(ele);
  struc.ordernumber=order++;
}


