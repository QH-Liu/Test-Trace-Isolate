#ifndef TEST_H
#define TEST_H
#include <string>

using namespace std;

const int agentNum=100000;
const int simu=100;

const double r[]={1.3,1.6,1.9,2.2,2.5,2.8};
const string r0[]={"1.3","1.6","1.9","2.2","2.5","2.8"};
const double beta[]={0.316,0.417,0.515,0.640,0.771,0.885};
const double mu[]={0.35,0.35,0.32,0.32,0.32,0.30};
const string kMu[]={"0.35","0.35","0.32","0.32","0.32","0.30"};

// const double beta[]={0.467,0.640,0.810,0.993,1.200,1.430};          //relative infectiousness 50%
// const double mu[]={0.32,0.32,0.32,0.32,0.30,0.30};                  //relative infectiousness 50%
// const string kMu[]={"0.32","0.32","0.32","0.32","0.30","0.30"};     //relative infectiousness 50%

void initial_net();                         //initial network
void calibration(int i);                    //parameter calibration
void degree();                              //network degree statistics
void npi_type(int type);                    //simulation under different npi_type setting
void npi_5_para();                          //npi_5 simulation under different Threshold and Probility setting
void npi_6_para();                          //npi_6 simulation under different Threshold and Probility setting
void test_tsc();                            //simulation under different Tsc setting
void test_tcr();                            //simulation under different Tcr setting
void tcr();                                 //simulation under different Tsc and npi_type setting
void test_ksym();                           //simulation under different Ksym setting
void test_ptest();                          //simulation under different Ptest setting
void test_ptrace();                         //simulation under different Ptrace setting
void test_dtrace();                         //simulation under different Dtrace setting
void test_kseed();                          //simulation under different Kseed setting
void test_RelativeI();                      //simulation under different Infectiousness setting

#endif