#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>

#include "../include/agent.h"
#include "../include/Calibration.h"
#include "../include/Distribution.h"
#include "../include/Network.h"
#include "../include/Population.h"
#include "../include/Simulation.h"
#include "../include/Test.h"

using namespace std;

void initial_net(){
    /*-----------initial network---------------------*/
    Population p(agentNum);
    p.initial();
    Network net(16, agentNum);
    net.individual_layer(p);
    net.group_layer(p);
    net.netInfo(p);
}

void calibration(int i){
    /*------------calibration----------------------*/
    Network net(16, agentNum);
    net.readNet(); 

    cout<<"---parameterization_of_R0_"<<r0[i]<<endl;
    cout<<"beta:"<<beta[i]<<endl;
    cout<<"paraM:"<<mu[i]<<endl;
    Calibrate(net,beta[i],kMu[i]);
}

void degree(){
    /*--------------degree-----------------------------------------------------*/
    Population p(agentNum);
    p.initial();
    Network net(16, agentNum);
    net.readNet();

    double avg_i=0;
    vector<double> avg_a(16,0);
    vector<double> g;
    double avg_tsc=0;
    double avg_tcr=0;
    double symRate=0;

    for(int i=0;i<agentNum;i++){
        int a=net.nodes[i].age;
        avg_a[a]+=net.nodes[i].iDegree;
        avg_i+=net.nodes[i].iDegree;
        avg_tsc+=net.nodes[i].Tsc;
        avg_tcr+=net.nodes[i].Tcr;

        if(net.nodes[i].willSymp){
            symRate++;
        }
    }

    for(int i=0;i<16;i++){
        avg_a[i]/=p.ageSize[i];
        cout<<setprecision(3)<<i<<":"<<avg_a[i]<<endl;
    }

    avg_i/=agentNum;
    for(int i=0;i<simuDay;i++){
        double tmp=0;
        int gs=net.groups[i].size();
        for(int j=0;j<gs;j++){
            int s=net.groups[i][j].size();
            tmp+=(s*(s-1));
        }
        tmp/=agentNum;
        g.push_back(tmp);
    }
    double avg_g=accumulate(g.begin(),g.end(),0.0)/simuDay;
    cout<<"individidual degree: "<<avg_i<<",group degree: "<<avg_g<<endl; 
    cout<<setprecision(3)<<"avg Tsc: "<<avg_tsc/agentNum<<", avg Tcr: "<<avg_tcr/agentNum<<endl;
    cout<<setprecision(3)<<"symRate: "<<symRate/agentNum<<endl;
}

void npi_type(int type){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);
  
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        string file="../result/baseline_npi_"+to_string(type)+"/r0_"+r0[test_r[i]]+".txt";
        f.open(file,ios::out);
        for(int k=0;k<simu;k++){
            Simulation s;
            s.npi_type=type;
            s.kBeta=beta[test_r[i]];
            s.kMu=mu[test_r[i]];
            s.transmission(net);
            s.outCome(net);
            s.simuInfo(f);
        }
        f.close();  
    }  
}

void npi_5_para(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    double pro[]={0.25,0.5,0.75};
    int thr[]={1,50,100};
    string p[]={"0.25","0.5","0.75"};
    string thresh[]={"thresh_1","thresh_50","thresh_100"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    for(int k=1;k<2;k++){
        net.setRecovery(kMu[test_r[k]]);
        // //threshold=50
        for(int i=0;i<3;i++){
            string file="../result/npi_5_para/"+thresh[1]+"/r0_"+r0[test_r[k]]+"_"+"pro_"+p[i]+".txt";
            f.open(file,ios::out);
            for(int m=0;m<simu;m++){
                Simulation s;
                s.npi_type=5;
                s.kBeta=beta[test_r[k]];
                s.kMu=mu[test_r[k]];
                s.kLimit=thr[1];
                s.kPro=pro[i];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close(); 
        }
        //pro=0.5
        for(int i=0;i<3;i++){
            string file="../result/npi_5_para/"+thresh[i]+"/r0_"+r0[test_r[k]]+"_"+"pro_"+p[1]+".txt";
            f.open(file,ios::out);
            for(int m=0;m<simu;m++){
                Simulation s;
                s.npi_type=5;
                s.kBeta=beta[test_r[k]];
                s.kMu=mu[test_r[k]];
                s.kLimit=thr[i];
                s.kPro=pro[1];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close(); 
        }
    }
}

void npi_6_para(){
    Network net(16, agentNum);
    net.readNet();
    
    int thr[]={1,50,100};
    string thresh[]={"thresh_1","thresh_50","thresh_100"};
    double pro[]={0.1,0.3,0.5};
    string p[]={"0.1","0.3","0.5"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    for(int k=0;k<6;k++){
        net.setRecovery(kMu[test_r[k]]);
        //threshold=50
        for(int i=0;i<3;i++){
            string file="../result/npi_6_para/"+thresh[1]+"/r0_"+r0[test_r[k]]+"_"+"pro_"+p[i]+".txt";
            f.open(file,ios::out);
            for(int m=0;m<simu;m++){
                Simulation s;
                s.npi_type=6;
                s.kBeta=beta[test_r[k]];
                s.kMu=mu[test_r[k]];
                s.kCut=thr[1];
                s.pCut=pro[i];
                net.setCnt(pro[i]);
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close(); 
        }
        //pro=0.3
        net.setCnt(pro[1]);
        for(int i=0;i<3;i++){
            string file="../result/npi_6_para/"+thresh[i]+"/r0_"+r0[test_r[k]]+"_"+"pro_"+p[1]+".txt";
            f.open(file,ios::out);
            for(int m=0;m<simu;m++){
                Simulation s;
                s.npi_type=6;
                s.kBeta=beta[test_r[k]];
                s.kMu=mu[test_r[k]];
                s.kCut=thr[i];
                s.pCut=pro[1];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close(); 
        }
    }
}
void test_tsc(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    string p[]={"1"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;
   
    //Tsc
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<1;j++){
            string file="../result/parameter_npi_2/tsc/r0_"+r0[test_r[i]]+"_tsc_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];           
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void test_tcr(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    string p[]={"4"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    //Tcr
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<1;j++){
            string file="../result/parameter_npi_2/tcr/r0_"+r0[test_r[i]]+"_tcr_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.npi_type=2;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void tcr(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    string p[]={"4"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    //Tcr
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<1;j++){
            string file="../result/tcr_test/npi_3/r0_"+r0[test_r[i]]+"_tcr_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.npi_type=3;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
        for(int j=0;j<1;j++){
            string file="../result/tcr_test/npi_4/r0_"+r0[test_r[i]]+"_tcr_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.npi_type=4;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void test_ksym(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    double ks[]={0.6,1.0,1.2,1.4};
    string k[]={"0.6","0.8","1.0","1.2","1.4"};
    int test_r[]={0,1,2,3,4};
    ofstream f;

    //Ksym
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            string file="../result/parameter_npi_2/ksym/r0_"+r0[test_r[i]]+"ksym_"+k[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                for(int m=0;m<agentNum;m++){
                    int age=net.nodes[m].age;
                    double p;
                    if(age<4) p=Psym[0];
                    else if(age<8) p=Psym[1];
                    else if(age<12) p=Psym[2];
                    else if(age<15 || Prand()<k80) p=Psym[3];
                    else p=Psym[4];
                    p*=ks[j];
                    if(Prand()<p){
                        net.nodes[m].willSymp=true;
                        if(Prand()<Ptest)
                            net.nodes[m].willTest=true;
                        else net.nodes[m].willTest=false;
                    }
                    else net.nodes[m].willSymp=false; 
                }
               
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);       
            }
            f.close();
        }
    }
}

void test_ptest(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    double pt[]={0.6,1.0};
    string p[]={"0.6","1.0"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    //Ptest
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<2;j++){
            string file="../result/parameter_npi_2/ptest/r0_"+r0[test_r[i]]+"_ptest_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                for(int m=0;m<agentNum;m++){
                    if(net.nodes[m].willSymp){
                        if(Prand()<pt[j])
                            net.nodes[m].willTest=true;
                        else net.nodes[m].willTest=false;
                    } 
                }
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void test_ptrace(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    double pt[]={0.3,0.7};
    string p[]={"0.3","0.7"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    //Ptrace
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<2;j++){
            string file="../result/parameter_npi_2/ptrace/r0_"+r0[test_r[i]]+"_ptrace_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                s.kPtrace=pt[j];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void test_dtrace(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    int dt[]={2,6};
    string p[]={"2","6"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    //Ptrace
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<2;j++){
            string file="../result/parameter_npi_2/dtrace/r0_"+r0[test_r[i]]+"_dtrace_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                s.kTraceDay=dt[j];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void test_kseed(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);

    int init[]={1,10};
    string p[]={"1","10"};
    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    //kseed
    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        for(int j=0;j<1;j++){
            string file="../result/parameter_npi_2/kseed/r0_"+r0[test_r[i]]+"_kseed_"+p[j]+".txt";
            f.open(file,ios::out);
            for(int k=0;k<simu;k++){
                Simulation s;
                s.kBeta=beta[test_r[i]];
                s.kMu=mu[test_r[i]];
                s.kSeed=init[j];
                s.transmission(net);
                s.outCome(net);
                s.simuInfo(f);
            }
            f.close();
        }
    }
}

void test_RelativeI(){
    Network net(16, agentNum);
    net.readNet();
    net.setCnt(0);
    /*-----relative infectiousness----------------------*/
    for(int i=0;i<net.agentNum;i++){
        if(!net.nodes[i].willSymp)
            net.nodes[i].heterogeneity*=0.5;
    }

    int test_r[]={0,1,2,3,4,5};
    ofstream f;

    for(int i=0;i<6;i++){
        net.setRecovery(kMu[test_r[i]]);
        string file="../result/parameter_npi_2/relative/r0_"+r0[test_r[i]]+".txt";
        f.open(file,ios::out);
        for(int k=0;k<simu;k++){
            Simulation s;
            s.kBeta=beta[test_r[i]];
            s.kMu=mu[test_r[i]];
            s.transmission(net);
            s.outCome(net);
            s.simuInfo(f);
        }
        f.close();

        string file2="../result/parameter_npi_2/relative/npi_0_r0_"+r0[test_r[i]]+".txt";
        f.open(file2,ios::out);
        for(int k=0;k<simu;k++){
            Simulation s;
            s.npi_type=0;
            s.kBeta=beta[test_r[i]];
            s.kMu=mu[test_r[i]];
            s.transmission(net);
            s.outCome(net);
            s.simuInfo(f);
        }
        f.close();
    }
}
