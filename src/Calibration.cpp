#include "../include/Calibration.h"
#include "../include/Distribution.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

double Fit(vector<int> &data){              //linear fit
    vector<int> x;
    vector<double> y;
    int s;
    double a=-1;
    double avg_x=0.0;
    double avg_y=0.0;
    double sum_xy=0.0;
    double sum_xx=0.0;

    //period that monotonically increasing in [index1,index2]
    int index1=0;
    int index2=data.size()-1;
    if(index2>=0 && data[index2]>exp(9)){
        while(data[index2]>exp(9))
            index2--;
        while(data[index1]<exp(5))
            index1++;
        for(int i=index1;i<=index2;i++){
            y.push_back(log(data[i]));
        }
    }

    s=y.size();
    if(s>0){
        for(int i=0;i<s;i++){
            x.push_back(i);
            avg_x+=x[i];
            avg_y+=y[i];
        }
        avg_x/=s;
        avg_y/=s;
        for(int i=0;i<s;i++){
            sum_xy+=(x[i]-avg_x)*(y[i]-avg_y);
            sum_xx+=(x[i]-avg_x)*(x[i]-avg_x);
        }    
        a=sum_xy/sum_xx;
    }
    return a;
}

void Calibrate(Network &net, double kBeta, string kMu){
    int seed;
    int t;
    double estimateR=0;                             //estimation of R0
    double estimateT=0;                             //estimation of Tg
    double tmpR;                                    //tmp value
    double tmpT;
    vector<double> testR;                           //test record
    vector<double> testT;
    vector<double> freT(26,0.0);                    //Tg distribution

    string name ="../result/calibration/infectNum.txt"; 
    string name1="../result/calibration/generation_time.txt";
    string name2="../result/calibration/growth_rate.txt";
    ofstream file;

    //set nodes's recovery
    net.setRecovery(kMu);

    /*-----relative infectiousness----------------------*/
    // for(int i=0;i<net.agentNum;i++){
    //     if(!net.nodes[i].willSymp)
    //         net.nodes[i].heterogeneity*=0.5;
    // }

    /* ---------------Tg----------------------*/
    for(int m=0; m<kTestNum*10; m++){                  //simulate for kTestNum times
        vector<int> newInfection;                   //new infection list
        vector<int> infectNum;                      //daily new infection
        t=0;                                         
        seed = Uniform(net.agentNum);               //random seed   
        net.nodes[seed].infected=true;
        net.nodes[seed].timeInfected=t;             //at t=0, set seed node infected     
        
        //transmit till seed is recoverd
        while(!net.nodes[seed].recoverd){           
            infectNum.push_back(0);
            if(net.nodes[seed].is_infectious(t) && net.nodes[seed].heterogeneity>0){
                //infect individual contacts
                for(int i=0; i<net.nodes[seed].iDegree;i++){
                    int infectee=net.nodes[seed].iCnt[i];
                    if(!net.nodes[infectee].infected){
                        double p = Prand();
                        if(p<net.nodes[infectee].susceptibility*wi*kBeta*net.nodes[seed].heterogeneity){
                            net.nodes[infectee].infected=true;
                            net.nodes[infectee].timeInfected=t;
                            infectNum[t]++;
                            newInfection.push_back(infectee);
                        }
                    }
                }
            
                //infect group contacts
                if(net.nodes[seed].gCnt.size()>0){
                    int gs=net.nodes[seed].gCnt.size();
                    for(int i=0;i<gs;i++){
                        if(net.nodes[seed].gCnt[i].first==t){                   //has group at day t
                            int g=net.nodes[seed].gCnt[i].second;
                            int ng=net.groups[t][g].size();
                            for(int j=0;j<ng;j++){
                                int infectee=net.groups[t][g][j];
                                if(infectee!=seed && !net.nodes[infectee].infected){
                                    double p = Prand();
                                    if(p<net.nodes[infectee].susceptibility*wg*kBeta*net.nodes[seed].heterogeneity){
                                        net.nodes[infectee].infected=true;
                                        net.nodes[infectee].timeInfected=t;
                                        infectNum[t]++;
                                        newInfection.push_back(infectee);
                                    }
                                }
                            }
                            break;
                        }
                    }
                }     
            }
            net.nodes[seed].update(t);             
            t++;               
        }

        //cumulative new infection
        tmpR=accumulate(infectNum.begin(),infectNum.end(),0.0);
        tmpT=0;
        if(tmpR>0){
            for(int i=0; i<infectNum.size(); i++){
                tmpT+=(infectNum[i]*i);
            }
            tmpT/=tmpR;
            testT.push_back(tmpT);          // add result tmpT to record list testT
        }
       
        //reset node state
        net.nodes[seed].reset();
        for(int k=0;k<newInfection.size();k++){
            int index=newInfection[k];
            net.nodes[index].reset();
        }
    }

    //save Tg distribution
    for(int i=0;i<testT.size();i++)
            freT[floor(testT[i])]++;       //if T=a.b, is in [a,a+1)
    for(int i=0;i<freT.size();i++)
        freT[i]/=testT.size();

    // file.open(name1,ios::out);
    // for(int i=0;i<freT.size();i++)
    //     file<<freT[i]<<endl;
    // file.close();


    /* ----------------------------R0---------------------------------------*/
    // file.open(name,ios::out);

    for(int m=0; m<kTestNum; m++){
        vector<int> oldInfection,newInfection;               //infection list
        vector<int> record;                                  //record list of nodes to reset
        vector<int> infectNum;                               //daily new infection
        seed = Uniform(net.agentNum);                        //random seed
        t=0;                                                 //at t=0, set seed node infected 
        net.nodes[seed].infected=true;
        net.nodes[seed].timeInfected=t;                                         
        oldInfection.push_back(seed);
        record.push_back(seed);

        //transmit till no more infected nodes or time over kTestDay
        while(!oldInfection.empty() && t<kTestDay){
            infectNum.push_back(0);
            for(int i=0;i<oldInfection.size();i++){
                int infector=oldInfection[i];

                if(net.nodes[infector].is_infectious(t) && net.nodes[infector].heterogeneity>0){
                    //infect individual contacts
                    for(int j=0; j<net.nodes[infector].iDegree;j++){
                        int infectee=net.nodes[infector].iCnt[j];
                        if(!net.nodes[infectee].infected){
                            double p = Prand();
                            if(p<net.nodes[infectee].susceptibility*wi*kBeta*net.nodes[infector].heterogeneity){
                                net.nodes[infectee].infected=true;
                                net.nodes[infectee].timeInfected=t;
                                infectNum[t]++;
                                newInfection.push_back(infectee);
                                record.push_back(infectee);
                            }
                        }
                    }

                    //infect group contacts
                    if(net.nodes[infector].gCnt.size()>0){
                        int gs=net.nodes[infector].gCnt.size();
                        for(int j=0;j<gs;j++){
                            if(net.nodes[infector].gCnt[j].first==t){
                                int g=net.nodes[infector].gCnt[j].second;
                                int ng=net.groups[t][g].size();
                                for(int k=0; k<ng;k++){
                                    int infectee=net.groups[t][g][k];
                                    if(infectee!=infector && !net.nodes[infectee].infected){
                                        double p = Prand();
                                        if(p<net.nodes[infectee].susceptibility*wg*kBeta*net.nodes[infector].heterogeneity){
                                            net.nodes[infectee].infected=true;
                                            net.nodes[infectee].timeInfected=t;
                                            infectNum[t]++;
                                            newInfection.push_back(infectee);
                                            record.push_back(infectee);
                                        }
                                    }
                                }
                                break;
                            }
                        }                   
                    }       
                }
            }

            //update infection state
            for(int i=0;i<oldInfection.size();i++){
                int infector=oldInfection[i];
                net.nodes[infector].update(t);        
                if(!net.nodes[infector].recoverd){
                    newInfection.push_back(infector);
                }
            }

            oldInfection.swap(newInfection);
            newInfection.clear();
            t++;               
        }
        
        //count cumulative infection
        for(int k=1;k<infectNum.size();k++)
            infectNum[k]+=infectNum[k-1];
        // for(int k=0;k<infectNum.size();k++){
        //     file<<noskipws<<infectNum[k]<<" ";
        // }
        // file<<endl;

        tmpR=Fit(infectNum);                    //growth rate for one simulation
        if(tmpR>=0) testR.push_back(tmpR);      //if growth rate >=0, add to test record list

        //reset node state
        for(int k=0;k<record.size();k++){
            int p=record[k];
            net.nodes[p].reset();
        }
    }
    // file.close();

    //save growth rate
    // file.open(name2,ios::out);
    // for(int i=0;i<testR.size();i++){
    //     file<<testR[i]<<endl;
    // }
    // file.close();

    /* -----estimate the value of Tg and R0----------------- */
    estimateT=accumulate(testT.begin(),testT.end(),0.0);
    estimateT/=(testT.size());
    estimateR=accumulate(testR.begin(),testR.end(),0.0);
    estimateR/=testR.size();
    cout<<"growth rate:"<<estimateR<<endl;

    double sum=0.0;  
    for(int i=1;i<freT.size();i++)
        sum+=freT[i]*(exp(-estimateR*(i-1))-exp(-estimateR*i));
    estimateR/=sum;
    
    cout<<"estimateT:"<<estimateT<<endl;
    cout<<"estimateR:"<<estimateR<<endl;
}
