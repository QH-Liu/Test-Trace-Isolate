#ifndef AGENT_H
#define AGENT_H
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

const double k80=0.5007;                                    //propotion of 80+ in age group 75+
const double kSym=1.0;                                      //relative probability of being symptomatic
const double Ptest=0.8;                                     //probability of being tested if symptomatic
const double kDelta[]={0.58,1.0,1.65};                      //susceptibility to infection by age a
const double Psym[]={0.181,0.224,0.305,0.355,0.646};        //pro of developing symptoms  
const double Phos[]={0.025,0.074,0.053,0.13,0.172,0.275,0.43};
const double Picu[]={0,0.004,0.009,0.026,0.072,0.184};
const double Pdeath[]={0,0.004,0.009,0.056,0.081}; 

const int Lambda=2;                                         //pre-symptomatic infectious stage
const int kTsc=2;                                           //mean delay
const int kTcr=1;                                           //mean delay
const int Tq=14;                                            //duration of quarantine 
const double Ppcr[]={0.979,0.686,0.363,0.300,0};            //pro of testing positive depends on delay between sym-sample
                                    

class Agent
{
public:
    /* NODE INFORMATION */
    int id;                     //0,1,...,n-1
    int age;                    //age_group:0,1,..,15
    int iDegree;                //num of individual contacts

    double heterogeneity;       //heterogeneity
    double susceptibility;      //susceptibility to infection by age a  
    int incubation;             //incubation period
    int recovery;               //recovery period 
    
    int Tsc;                    //Delay between symptom set on and collection of the sample
    int Tcr;                    //Delay between collection of the sample and PCR results
      
    vector<int> iCnt;           //individual contacts
    vector<pair<int,int>> gCnt; //pair=<day_has_gc,index_group_at_day>
    vector<bool> cntFlag;       //flag of contacts, false means this contacts is banned
    
    int secondary=0;            //number of secondary infection
    int trans_presym=0;         //number of transmission pre-symptom
    int trans_aftstm=0;         //number of transmission after symptom onset
    int trans_test=0;           //number of transmission while waiting for test
    int infector=-1;            //infection from infector           
    int transmission=-1;        //transmission=0, individual transmission; transmission=1, group transmission

    /* STATE SYMBOL */
    bool infected=false;        //be infected
    bool symptomatic=false;     //have developed symptom
    bool recoverd=false;        //have recoverd
    
    bool willSymp=false;        //will develop symptom after incubation
    bool willTest=false;        //will be tested if willsymp
    int symTestFlag=-1;         //flag=-1,:have not start sym-test; flag=1: have sampled

    bool isTesting=false;       //is waiting for test results
    bool confirmed=false;       //be confirmed

    int stateIsolation=0;       //0:free; 1:without group contacts; 2:without individual and group contacts
    bool traced=false;          //true if found by contact tracing

    /* TIME POINT */
    int timeInfected=-1;
    int timeSymptomOn=-1;
    int timeRecover=-1;
    int timeConfirmd=-1;
    int timeTraced=-1;

    vector<int> sampleCollected;            //time point when sample collected
    vector<bool> testResult;                //result of test
    vector<int> testOut;                    //time when test result is out
    vector<int> testType;                   //1:symptom driven or 0:by contact tracing

    vector<int> startQ;                     //time of quarantine start
    vector<int> endQ;                       //time of quarantine end
    vector<int> typeQ;                      //type of quarantine                

    /* FUNCTION */
    Agent();
    ~Agent()=default;

    void initial(int i, int a, int d, double h);        //initial parameter
    bool is_infectious(int t);                          //if is infectious
    void check_incubation(int t);                       //if incubation is over
    void update(int t);                                 //update recovery state
    void symTest(int t);                                //add symptom-driven test
    void traceTest(int t);                              //add tracing-driven test
    void checkTest(int t);                              //update test state, return true when result is out
    void reset();                                       //reset state and timepoint
};

#endif