#ifndef MYFUNC_H
#define MYFUNC_H

#include <vector>
#include <string>
#include <ostream>
#include <Eigen/Dense>//!!! how is this one related to string.remove???

using std::string;
using std::vector;
using Eigen::MatrixXd;
const double min(const double *num,const int n);
const double max(const double *num,const int n);
const double ave(const double *num,const int n);
void StringSplit(const string &s,vector<string> &vec);
int getlinenumber(const string &p);
double gaussrand();
double bessi0(double x);
double bessi1(double x);
double Bessin(int n,double x);
MatrixXd twissmapx(const MatrixXd &R);//beta alpha gamma,in x direction test
MatrixXd twissmapy(const MatrixXd &R);//beta alpha gamma,in y direction test
MatrixXd twissmapz(const MatrixXd &R);//beta alpha gamma,in z direction test
void printvector(std::ostream &par,const vector<vector<double>  > &particle);
#endif
