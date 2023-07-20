#include "../include/Simulation.h"
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>
using namespace std;

Simulation::Simulation(){
    for(int i=0;i<kDay;i++){
        dailyInfection.push_back(0);
        dailyCase.push_back(0);
    }
}

void Simulation::transmission(Network &net){
    vector<int> oldInfection,newInfection;          //record infection
    vector<int> oldSym,newSym;                      //update sym-driven test
    vector<int> oldTesting,newTesting;              //record testing
    vector<int> oldQua,newQua;                      //record quarantine
    vector<int> oldCon,newCon;
    
    int t=0;
    int symcase=0;                                  

    /*-----------------initial---------------------------------*/
    //initial groupState: all groups saved
    for(int i=0;i<kDay;i++){
        groupState.push_back({});
        for(int j=0;j<net.groups[i].size();j++)            
            groupState[i].push_back(true);
    }
    
    //initial seeds
    while(oldInfection.size()<kSeed){
        int seed=Uniform(net.agentNum);
        if(find(oldInfection.begin(),oldInfection.end(),seed)==oldInfection.end()){     //no repeat seed
            net.nodes[seed].infected=true;          //infected
            net.nodes[seed].timeInfected=t;
            oldInfection.push_back(seed);
            if(net.nodes[seed].willTest)            //if sym-test
                oldSym.push_back(seed);
        }
    }
             
    /*----------------transmission-------------------------*/
    while(t<kDay && oldInfection.size()>0){         //transmit till no more infection nodes or last more than kDay
        for(int i=0;i<oldInfection.size();i++){
            int infector=oldInfection[i];

            //if node is not isolated without individual contacts and is infectious
            if(net.nodes[infector].stateIsolation<2 && net.nodes[infector].is_infectious(t) && net.nodes[infector].heterogeneity>0){
                //infect individual contacts that are susceptible and have individual contacts
                for(int j=0;j<net.nodes[infector].iDegree;j++){
                    if(tCut==-1 || net.nodes[infector].cntFlag[j]){                     //if this contact exist or cutting contacts is not triggerd
                        int infectee=net.nodes[infector].iCnt[j];
                        if(!net.nodes[infectee].infected && net.nodes[infectee].stateIsolation<2){
                            if(Prand()<net.nodes[infectee].susceptibility*wi*kBeta*net.nodes[infector].heterogeneity){
                                net.nodes[infectee].infected=true;
                                net.nodes[infectee].timeInfected=t;
                                net.nodes[infectee].infector=infector;
                                net.nodes[infectee].transmission=0;
                                newInfection.push_back(infectee);

                                net.nodes[infector].secondary++;

                                if(net.nodes[infector].timeSymptomOn==-1)
                                    net.nodes[infector].trans_presym++;
                                else net.nodes[infector].trans_aftstm++;

                                if(net.nodes[infector].stateIsolation==1)
                                    net.nodes[infector].trans_test++;
                            }      
                        }
                    }  
                }

                //infect group contacts if not at isolation and no gathering limit
                if(net.nodes[infector].stateIsolation==0 && net.nodes[infector].gCnt.size()>0){
                    int gs=net.nodes[infector].gCnt.size();
                    for(int j=0;j<gs;j++){
                        if(net.nodes[infector].gCnt[j].first==t){               //group at day t
                            int g=net.nodes[infector].gCnt[j].second;           //group index
                            if(groupState[t][g]){                               //this group is not limited
                                int ng=net.groups[t][g].size();
                                for(int k=0;k<ng;k++){ 
                                    int infectee=net.groups[t][g][k];
                                    if(infectee!=infector && !net.nodes[infectee].infected && net.nodes[infectee].stateIsolation==0){
                                        if(Prand()<net.nodes[infectee].susceptibility*wg*kBeta*net.nodes[infector].heterogeneity){
                                            net.nodes[infectee].infected=true;
                                            net.nodes[infectee].timeInfected=t;
                                            net.nodes[infectee].infector=infector;
                                            net.nodes[infectee].transmission=1;
                                            newInfection.push_back(infectee);

                                            net.nodes[infector].secondary++; 

                                            if(net.nodes[infector].timeSymptomOn==-1)
                                                net.nodes[infector].trans_presym++;
                                            else net.nodes[infector].trans_aftstm++;   
                                        }      
                                    }
                                }
                            }
                            break;
                        }
                    }                       
                }       
            }
        }

        //outcome
        dailyInfection[t]=newInfection.size();                  //num of new infection at day t
        if(peakN<newInfection.size()){                          //peak new infection
            peakN=newInfection.size();
            tPeakN=t;
        }

        //update symptom
        for(int i=0;i<oldInfection.size();i++){
            int index=oldInfection[i];
            if(net.nodes[index].timeSymptomOn==-1)
                net.nodes[index].check_incubation(t);
        }
        for(int i=0;i<newInfection.size();i++){                 //if incubation==0
            int index=newInfection[i];
            if(net.nodes[index].timeSymptomOn==-1)
                net.nodes[index].check_incubation(t);
        }  

        //npi
        int tracedNum=0;
        if(npi_type>0){
            //new Sym test
            for(int i=0;i<newInfection.size();i++){
                int index=newInfection[i];
                if(net.nodes[index].willTest){
                    newSym.push_back(index);                //Tsc>=1
                }       
            }

            //symptom-driven test
            for(int i=0;i<oldSym.size();i++){
                int index=oldSym[i];
                if(!net.nodes[index].confirmed){
                    if(net.nodes[index].symptomatic){
                        int time=net.nodes[index].timeSymptomOn+net.nodes[index].Tsc;
                        if(t<time) newSym.push_back(index);     //check later
                        if(t==time && !net.nodes[index].isTesting){
                            net.nodes[index].symTest(t);
                            if(net.nodes[index].Tcr>0){
                                if(npi_type==2 || npi_type==5 || npi_type==6){
                                    net.nodes[index].stateIsolation=1;
                                    net.nodes[index].startQ.push_back(t);
                                    net.nodes[index].typeQ.push_back(1);
                                }
                                if(npi_type==4){
                                    net.nodes[index].stateIsolation=2;
                                    net.nodes[index].startQ.push_back(t);
                                    net.nodes[index].typeQ.push_back(2);
                                }
                                newTesting.push_back(index);
                            }                               
                            else{//check result at this day
                                net.nodes[index].checkTest(t);
                                if(net.nodes[index].confirmed){
                                    oldCon.push_back(index);
                                    newQua.push_back(index);
                                    if(net.nodes[index].symptomatic) symcase++;
                                } 
                            }
                        } 
                    }
                    else newSym.push_back(index);               //check later
                }
            }

            //check test result
            for(int i=0;i<oldTesting.size();i++){
                int index=oldTesting[i];
                int time=net.nodes[index].sampleCollected.back()+net.nodes[index].Tcr;
                if(t<time) newTesting.push_back(index);         //check later
                if(t==time) {
                    net.nodes[index].checkTest(t);
                    if(net.nodes[index].confirmed){
                        oldCon.push_back(index);
                        newQua.push_back(index);
                        if(net.nodes[index].symptomatic) symcase++;
                    }
                    else{
                        if(net.nodes[index].stateIsolation>0){
                            net.nodes[index].stateIsolation=0;
                            net.nodes[index].endQ.push_back(t);
                        }
                    }         
                }
            }
            
            if(npi_type>1){
                for(int i=0;i<oldCon.size();i++){
                    int index=oldCon[i];

                    //trace all individual contacts
                    for(int j=0;j<net.nodes[index].iDegree;j++){
                        if(tCut==-1 || net.nodes[index].cntFlag[j]){                     //if this contact exist or cutting contacts is not triggerd
                            int infectee=net.nodes[index].iCnt[j];
                            if(!net.nodes[infectee].confirmed && !net.nodes[infectee].isTesting){
                                net.nodes[infectee].traceTest(t);
                                tracedNum++;
                                if(net.nodes[infectee].Tcr>0){
                                    newTesting.push_back(infectee);
                                    if(npi_type==2 || npi_type==5 || npi_type==6){
                                        net.nodes[infectee].stateIsolation=1;     //no group contacts when waiting for test results
                                        net.nodes[infectee].startQ.push_back(t);
                                        net.nodes[infectee].typeQ.push_back(1);
                                    }                            
                                    if(npi_type==4){
                                        net.nodes[infectee].stateIsolation=2;     //no contacts when waiting for test results
                                        net.nodes[infectee].startQ.push_back(t);
                                        net.nodes[infectee].typeQ.push_back(2);
                                    }       
                                }                    
                                else{//check result now
                                    net.nodes[infectee].checkTest(t);
                                    if(net.nodes[infectee].confirmed){
                                        newCon.push_back(infectee);
                                        newQua.push_back(infectee);
                                        if(net.nodes[infectee].symptomatic) symcase++;
                                    }
                                }
                            }                             
                        }                                     
                    }

                    //trace group contacts in kTraceday since sample collected with probability of kPtrace
                    if(net.nodes[index].gCnt.size()>0){
                        for(int j=0;j<net.nodes[index].gCnt.size();j++){
                            int day=net.nodes[index].gCnt[j].first;         //day of group contacts
                            int tr=net.nodes[index].sampleCollected.back(); //day of tracing start
                            if(npi_type==3) tr=t;                           //day of being confirmed

                            if(day>tr-kTraceDay && day<=tr){                //group of this day is within tracing time window
                                int g=net.nodes[index].gCnt[j].second;      //group index
                                if(groupState[day][g]){                     //this group is not limited
                                    for(int k=0; k<net.groups[day][g].size();k++){    
                                        if(Prand()<kPtrace){
                                            int infectee=net.groups[day][g][k];
                                            if(infectee!=index && !net.nodes[infectee].confirmed && !net.nodes[infectee].isTesting){
                                                net.nodes[infectee].traceTest(t);
                                                tracedNum++;
                                                if(net.nodes[infectee].Tcr>0){
                                                    if(npi_type==2 || npi_type==5 || npi_type==6){
                                                        net.nodes[infectee].stateIsolation=1;     //no group contacts when waiting for test results
                                                        net.nodes[infectee].startQ.push_back(t);
                                                        net.nodes[infectee].typeQ.push_back(1);
                                                    }
                                                    
                                                    if(npi_type==4){
                                                        net.nodes[infectee].stateIsolation=2;     //no contacts when waiting for test results
                                                        net.nodes[infectee].startQ.push_back(t);
                                                        net.nodes[infectee].typeQ.push_back(2);
                                                    }
                                                    newTesting.push_back(infectee);    
                                                }                                                   
                                                else{//check result now
                                                    net.nodes[infectee].checkTest(t);
                                                    if(net.nodes[infectee].confirmed){
                                                        newCon.push_back(infectee);
                                                        newQua.push_back(infectee);
                                                        if(net.nodes[infectee].symptomatic) symcase++;
                                                    } 
                                                }
                                            }        
                                        }
                                    }
                                }                                
                            }
                        }                               
                    }
                } 
            }

            if(npi_type==5 || npi_type==6){
                //if over threshold, start gathering limit
                if(tLimit==-1 && symcase>kLimit){
                    tLimit=t;
                    for(int j=t;j<kDay;j++){
                        for(int k=0;k<groupState[j].size();k++){
                            if(Prand()<kPro){
                                groupState[j][k]=false;
                            }
                        }
                    }
                }
            }

            if(npi_type==6){
                //if over threshold, start contact cutting
                if(tCut==-1 && symcase>kCut){
                    tCut=t;
                }
            }

            //check quarantine state
            for(int i=0;i<oldQua.size();i++){
                int index=oldQua[i];
                int time=net.nodes[index].timeConfirmd+Tq;
                if(t<time) newQua.push_back(index);
                if(t==time){
                    net.nodes[index].stateIsolation=0;
                    net.nodes[index].endQ.push_back(t);
                } 
            }
            
            //update sym
            oldSym.swap(newSym);
            newSym.clear();
            //update testing
            oldTesting.swap(newTesting);
            newTesting.clear();
            //update quarantim
            int oqs=oldQua.size();
            maxQua=maxQua>oqs?maxQua:oqs;
            oldQua.swap(newQua);
            newQua.clear();
            //update confirm
            int ocs=oldCon.size();
            maxToTrace=maxToTrace>ocs?maxToTrace:ocs;
            oldCon.swap(newCon);
            newCon.clear();

        }
        
        //update oldInfection and newInfection
        for(int i=0;i<oldInfection.size();i++){
            int infector=oldInfection[i];
            net.nodes[infector].update(t);       
            if(!net.nodes[infector].recoverd){
                newInfection.push_back(infector);
            }
        }
        oldInfection.swap(newInfection);
        newInfection.clear();

        if(peakI<oldInfection.size()){      //peak num of individuals that are in infection state
            peakI=oldInfection.size();
            tPeakI=t;
        }
        t++;     
    }
    endTime=t;      
}

void Simulation::outCome(Network &net){
    double symTest=0,traTest=0,avgTest=0,avgWait=0;

    //infection and disease burden
    for(int i=0;i<net.agentNum;i++){
        if(net.nodes[i].infected){
            totalI++;                                       //cumulative infections
            if(net.nodes[i].confirmed)
                totalC++;                                   //cumulative confirmed
            if(net.nodes[i].symptomatic){
                int ts=net.nodes[i].timeSymptomOn;
                dailyCase[ts]++;
            }
            //disease burden
            int a=net.nodes[i].age;
            double ph,pi,pd;
            if(a<3){
                ph=Phos[0];
                pi=Picu[0];
                pd=Pdeath[0];
            }
            else if(a<4){
                ph=Phos[1];
                pi=Picu[0];
                pd=Pdeath[0];
            }
            else if(a<8){
                ph=Phos[2];
                pi=Picu[1];
                pd=Pdeath[0];
            }
            else if(a<10){
                ph=Phos[3];
                pi=Picu[2];
                pd=Pdeath[0];
            }
            else if(a<12){
                ph=Phos[3];
                pi=Picu[2];
                pd=Pdeath[1];
            }
            else if(a<14){
                ph=Phos[4];
                pi=Picu[3];
                pd=Pdeath[2];
            }
            else if(a<15 || Prand()<k80){
                ph=Phos[5];
                pi=Picu[4];
                pd=Pdeath[3];
            }
            else{
                ph=Phos[6];
                pi=Picu[5];
                pd=Pdeath[4];
            }

            if(net.nodes[i].symptomatic) symptoms++;
            if(Prand()<ph) hospitalisation++;
            if(Prand()<pi) icu++;
            if(Prand()<pd) death++; 

            if(net.nodes[i].transmission==0) infec_i++;
            if(net.nodes[i].transmission==1) infec_g++;
        }

        if(npi_type==4 && net.nodes[i].Tcr>0){
           int testNum=net.nodes[i].sampleCollected.size();
           test_qua+=(testNum*net.nodes[i].Tcr);        //quarantine days while waiting for results
        }

        if(net.nodes[i].sampleCollected.size()>0 && net.nodes[i].infected){
            avgTest++;
            avgWait+=net.nodes[i].sampleCollected.size();
            avg_trans_test+=net.nodes[i].trans_test;
            
            int first=0;
            for(int j=0;j<net.nodes[i].sampleCollected.size();j++){
                if(net.nodes[i].sampleCollected[j]>=net.nodes[i].timeInfected){
                    first=j;
                    break;
                }
            }

            if(net.nodes[i].testType[first]==0){
                traTest++;
                if(net.nodes[i].testResult.size()>first && net.nodes[i].testResult[first]){
                    traTestRate++;
                }
            } 
            if(net.nodes[i].testType[first]==1){
                symTest++;
                if(net.nodes[i].testResult.size()>first && net.nodes[i].testResult[first]){
                    symTestRate++;
                }
            }          
        }

        if(net.nodes[i].symptomatic){
            symRate+=net.nodes[i].secondary;
            presymRate+=net.nodes[i].trans_presym;
        }
        else asymRate+=net.nodes[i].secondary;

        if(reset)
            net.nodes[i].reset();
    }

    if(traTest+symTest>0) allTestRate=(traTestRate+symTestRate)/(traTest+symTest);
    if(traTest>0) traTestRate/=traTest;
    if(symTest>0) symTestRate/=symTest;
    if(avgTest>0) avg_trans_test/=avgWait;

    if(symRate>0) presymRate/=symRate;
    if(symRate+asymRate>0) asymRate/=(symRate+asymRate);
}

void Simulation::simuInfo(ofstream &f){
    if(!f.is_open())
        cout<<"cannot open file_simuInfo"<<endl;
    
    f<<noskipws<<totalI<<" "<<death<<endl;
    // f<<noskipws<<totalI<<" "<<death<<" "<<maxToTrace<<" "<<maxQua<<" "<<totalC<<" "<<test_qua<<" "<<symptoms<<" "<<hospitalisation<<" "<<icu<<" "<<peakN<<" "<<infec_i<<" "<<infec_g<<" "<<endTime<<" "<<tPeakN<<" "<<tLimit<<" "<<tCut<<endl;
    // f<<noskipws<<setprecision(4)<<avg_trans_test<<" "<<allTestRate<<" "<<symTestRate<<" "<<traTestRate<<" "<<totalI<<endl;  
    // f<<noskipws<<setprecision(4)<<presymRate<<" "<<asymRate<<endl;
}