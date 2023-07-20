#include "../include/Agent.h"
#include "../include/Distribution.h"

using namespace std;

Agent::Agent(){
    incubation=Gamma();            
    Tsc=Poisson(kTsc);
    if(Tsc==0) Tsc=1;      //Tsc >=1
    Tcr=Poisson2(kTcr); 
}

void Agent::initial(int i, int a, int d, double h){
    id=i;
    age=a;
    iDegree=d;
    heterogeneity=h;
    
    //initial susceptibility
    if(age<3) susceptibility=kDelta[0];
    else if(age<13) susceptibility=kDelta[1];
    else susceptibility=kDelta[2];
    susceptibility/=kDelta[2];                  

    //initial probability of developing symptoms
    double p;
    if(age<4) p=Psym[0];
    else if(age<8) p=Psym[1];
    else if(age<12) p=Psym[2];
    else if(age<15 || Prand()<k80) p=Psym[3];
    else p=Psym[4];

    p*=kSym;                    
    if(Prand()<p){
        willSymp=true;
        if(Prand()<Ptest)
            willTest=true;
    }
}

bool Agent::is_infectious(int t){
    //infectious
    if(Lambda>=incubation){
        return true;
    }
    else if(t>=(timeInfected+incubation-Lambda)){
        return true;
    }
    else return false;
}

void Agent::check_incubation(int t){
    //symptom onset
    if(t==(timeInfected+incubation)){
        timeSymptomOn=t;
        if(willSymp){
            symptomatic=true;   
        }          
    }
}

void Agent::update(int t){
    //recover state
    if(t==(timeInfected+incubation+recovery)){
        recoverd=true;
        timeRecover=t;
    }
}

void Agent::symTest(int t){
    isTesting=true;
    sampleCollected.push_back(t);
    testType.push_back(1);
    symTestFlag=1;
}

void Agent::traceTest(int t){
    isTesting=true;
    sampleCollected.push_back(t);
    testType.push_back(0);
    traced=true;
    timeTraced=t;
}

void Agent::checkTest(int t){
    //cheak test result after Tcr days of delay since sample collected 
    isTesting=false;
    bool result=false;
    double p=0;

    /*------ null model  ------------*/
    // if(infected && timeSymptomOn==-1 && testType.back()==0 ){
    //     p=0.25;
    //     if(Prand()<p){                  //confirmed
    //         result=true;
    //         confirmed=true;
    //         timeConfirmd=t;
    //         stateIsolation=2;
    //         startQ.push_back(t);
    //         typeQ.push_back(2);
    //     }              
    // }
    
    //only infected and sample collected after incubation individuals will have positive test result
    if(infected && timeSymptomOn>-1 && sampleCollected.back()>=timeSymptomOn){
        int tau=sampleCollected.back()-timeSymptomOn;        //time from symptom onset to sample collected       
        if(tau<=7) p=Ppcr[0];
        else if(tau<=14) p=Ppcr[1];
        else if(tau<=21) p=Ppcr[2];
        else if(tau<=28) p=Ppcr[3];
        else p=Ppcr[4];

        if(Prand()<p){                  //confirmed
            result=true;
            confirmed=true;
            timeConfirmd=t;
            stateIsolation=2;
            startQ.push_back(t);
            typeQ.push_back(2);
        }              
    }

    testResult.push_back(result);
    testOut.push_back(t); 
}

void Agent::reset(){
    secondary=0;
    trans_presym=0;
    trans_aftstm=0;
    trans_test=0;
    infector=-1;
    transmission=-1;
    symTestFlag=-1;  
    
    infected=false;
    symptomatic=false;
    recoverd=false;
 
    isTesting=false;
    confirmed=false;

    stateIsolation=0;
    traced=false;

    timeInfected=-1;
    timeSymptomOn=-1;
    timeRecover=-1;
    timeConfirmd=-1;
    timeTraced=-1;

    sampleCollected.clear();
    testResult.clear();
    testType.clear();

    startQ.clear();
    endQ.clear();
    typeQ.clear();
}

