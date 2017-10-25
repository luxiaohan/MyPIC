#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
using std::vector;

void CKV_Beam(const vector<double> &param,int particlenumber=10000);


void KVdistributionA(const vector<double> &param,int particlenumber=10000,char *p="particles.dat");


void KVdistributionS(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,int particlenumber=10000,char *p="particles.dat");

void KVdistributionT(const vector<double> &twiss,int particlenumber=10000,char *p="particles.dat");

void KVdistribution6A(double emitx,double betax,double emity,double betay,double emitz,double betaz,double alpx,double alpy,double alpz,int particlenumber=10000,char *p="particles.dat");

void KVdistribution6S(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,double sig_z,double sig_pz,double sig_zpz,int particlenumber=10000,char *p="particles.dat");

void KVdistribution6T(const vector<double> &twiss,int particlenumber,char *p);

void distribution(double a,double b,double c,double d,double rx=0,double ry=0,double rxy=0,double rxtyt=0,double rxyt=0,double rxty=0,int particlenumber=10000,char *p="particles.dat");


void KV2D_T(const vector<double> &twiss,vector<vector<double> > &xy);
void KV2D_S(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,vector<vector<double> > &xy);
void KV2D_A(const vector<double> &param,vector<vector<double> > &xy);
void PB2D_A(const vector<double> &param,vector<vector<double> > &xy);
void WB2D_A(const vector<double> &param,vector<vector<double> > &xy);
void GS2D_A(const vector<double> &param,vector<vector<double> > &xy);

void KV3D_T(const vector<double> &param,vector<vector<double> > &xy);

vector<double> distribution_base1(const vector<double> &param);
void distribution_base2(const vector<double> &par,vector<double> &xy);
#endif
