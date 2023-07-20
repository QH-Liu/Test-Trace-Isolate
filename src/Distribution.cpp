#include "../include/Distribution.h"
#include <random>
#include <algorithm>
#include <cmath>
#include <utility>

using namespace std;

int Uniform(int size){                        
    static default_random_engine e(seed1);
    static uniform_int_distribution<int> u(1000000);
    int rand=u(e);
    return rand%size;
}

int Gamma(){
    static default_random_engine e(seed2);
    static gamma_distribution<double> u(kShape,1/kRate);
    return round(u(e));
}

int Exp(double kMu){
    static default_random_engine e(seed3);
    static exponential_distribution<double> u(kMu);
    return round(u(e));
}

double Prand(){
    static default_random_engine e(seed4);
    static uniform_real_distribution<double> u(0.0,1);
    return u(e);
}

int Poisson(int m){
    static default_random_engine e(seed5);
    static poisson_distribution<int> u(m);
    return u(e);
}

int Poisson2(int m){
    static default_random_engine e(seed7);
    static poisson_distribution<int> u(m);
    return u(e);
}

double Heter(){
    static default_random_engine e(seed6);
    static gamma_distribution<double> u(kHeter,kHeter);
    return u(e);
}
