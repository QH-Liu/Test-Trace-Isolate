#ifndef SIMULATION_H
#define SIMULATION_H
#include "../include/Agent.h"
#include "../include/Population.h"
#include "../include/Network.h"
#include "../include/Distribution.h"
#include "../include/Calibration.h"

using namespace std;

class Simulation{
    public:
        int npi_type=2;              //type  of NPIs: 0,1,...,4 
        bool reset=true;

        int endTime=-1;             //time of simulation ends

        int kLimit=50;              //threshold of symptomatic confirmed cases to start gathering limiting
        double kPro=0.5;            //propotion of groups shut down if limiting
        int tLimit=-1;              //time of limiting start

        int kCut=100;               //threshold of symptomatic confirmed cases to start contacts cutting
        double pCut=0.2;            //propotion of contacts to be cut
        int tCut=-1;                //time of cutting indicidual contacts start

        int kDay=365;               //simulation period
        int kSeed=5;                //initial infected seed number
        int kTraceDay=4;            //day to trace
        double kPtrace=0.5;         //probability of sucessfully traced

        double kBeta;               //transmit probability
        double kMu;                 //recovery period
        
        vector<vector<bool>> groupState;    //false=group shut down

        //simulation results
        vector<int> dailyInfection;         //newly infection at each day
        vector<int> dailyCase;              //newly cases at each day
        int maxToTrace=-1;                  //max num to trace
        int maxQua=-1;                      //max num at quarantine

        double symTestRate=0;               //symptom-based true positive test rate
        double traTestRate=0;               //tracing-based true positive test rate
        double allTestRate=0;               
        double avg_trans_test=0;    

        double presymRate=0;                //pre-symptomatic infection count
        double symRate=0;                           
        double asymRate=0;
        
        int totalI=0;                       //total infections num
        int peakN=0;                        //peak num of daily new infection
        int tPeakN=-1;                      //time of peakN
        int peakI=0;                        //peak num of individuals that have been infected and not recovered
        int tPeakI=-1;                      //time of peakI
        int infec_i=0;                      //number of individuals infected through individual contacts
        int infec_g=0;                      //number of individuals indected through group contacts
        int test_qua=0;
        
        //disease burden
        int totalC=0;                       //total comfired infection
        int symptoms=0;                     
        int hospitalisation=0;
        int icu=0; 
        int death=0;
        
        Simulation();
        ~Simulation()=default;

        void transmission(Network &net);        //transmission and npi
        void outCome(Network &net);             //disease burden
        void simuInfo(ofstream &f);             //output simulation result
};

#endif