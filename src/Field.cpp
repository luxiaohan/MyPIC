#include "Field.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ifstream;
double Field::getfield(int z)
{
  return field[z];
}

void Field::read(const string &p)
{
  ifstream infile(p);
  string line;
  while(infile.peek()!=EOF){
    getline(infile,line);
    if(line.size()>0)
    {
      field.push_back(stod(line)*1e0);
    }
  }
}

void Field::addField(double a)
{
  field.push_back(a);
}



