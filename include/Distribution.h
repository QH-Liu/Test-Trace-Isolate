#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#include <utility>
#include <vector>

using namespace std;

const double kShape=2.08;           //shape parameter of gamma distribution of epsilon
const double kRate=0.33;            //rate parameter of gamma distribution of epsilon
const double kHeter=0.3;

const int seed1=10;
const int seed2=2;
const int seed3=3;
const int seed4=4;
const int seed5=5;
const int seed6=6;
const int seed7=7;

int Uniform(int size);               //return random int in (0,size-1)
int Gamma();                        //gamma distribution for incubation period
int Exp(double kMu);                //exponential distribution for remove period
double Prand();                     //return random number in (0.0,1.0)
int Poisson(int m);                 //poisson distribution
int Poisson2(int m);                 //poisson distribution
double Heter();                     //heterogeneity

#endif