#ifndef BASELINE_CALIBRATION_H
#define BASELINE_CALIBRATION_H
#include "Agent.h"
#include "Population.h"
#include "Network.h"

using namespace std;

const int kTestNum=1000;
const int kTestDay=300;

double Fit(vector<int> &data);
void Calibrate(Network &net, double kBeta, string kMu);

#endif