#include "../include/Agent.h"
#include "../include/Population.h"
#include "../include/Network.h"
#include "../include/Distribution.h"
#include <random>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

Network::Network(int nClass, int nAgent){
    ageClass=nClass;
    agentNum=nAgent;
}

void Network::setRecovery(string kMu){
    string name="../network/kMu/kMu_"+kMu+".txt";
    ifstream file;
    file.open(name,ios::in);
    string line;
    int index=0;
    while(getline(file,line)){
        stringstream ss(line);
        int d;
        while(ss>>d)
            nodes[index].recovery=d;
        index++;
    }
    file.close();
}

void Network::individual_layer(Population &popu){  
    //initial node id , age , idegree, heterogeneity    
    for(int i=0; i<agentNum; i++){                          
        Agent a;
        nodes.push_back(a);
        nodes[i].initial(i,popu.listAge[i],popu.listIdegree[i],popu.listH[i]);
    }

    /* ---------CONFIGURATION NETWORK---------------- */
    vector<vector<int>> stubArray;                                  //store free stubs of each age group
    vector<int> enableAge;                                          //age group that still have stubs
    int total=0;                                                    //total number of free stubs
    int failure=0;                                                  //+1 when connect failed to break while loop

    //initial stubArray:store agent id repeat degree times
    for(int i=0; i<ageClass; i++)                                   
        stubArray.push_back({});
    for(int i=0; i<agentNum; i++){                                  
        int id     = nodes[i].id;
        int age    = nodes[i].age;
        int degree = nodes[i].iDegree;
        for (int j=0;j<degree;j++){
            stubArray[age].push_back(id);
        }        
    }

    //initial enableAge   
    for(int i=0;i<ageClass;i++){
        total+=stubArray[i].size();
        enableAge.push_back(i);
    }
            
    //connect free stubs
    while(total>0 && failure<kFailure){
        int age1   = Uniform(enableAge.size());                         //random pick age group from enableAge
            age1   = enableAge[age1];
        int index1 = Uniform(stubArray[age1].size());                   //random pick one agent from the group
        int id1    = stubArray[age1][index1];

        double p = Prand();                                             //pick age gruop to according to pro
        int age2 = 0;
        while(p>popu.disAgei[age1][age2] && age2<ageClass-1)
            age2++;
        if(stubArray[age2].size()<1)                                    //fail when age gruop2 has not enough stubs 
            failure++;
        else{
            int index2 = Uniform(stubArray[age2].size());
            int id2    = stubArray[age2][index2];

            //agent1 has no iCnt with agent2
            if(find(nodes[id1].iCnt.begin(),nodes[id1].iCnt.end(),id2)==nodes[id1].iCnt.end()){
                if(age1==age2 && id1!=id2){
                    nodes[id1].iCnt.push_back(id2);
                    nodes[id2].iCnt.push_back(id1);
                    stubArray[age1].erase(stubArray[age1].begin()+index1);
                    if(index1<index2)
                        stubArray[age2].erase(stubArray[age2].begin()+index2-1);
                    else stubArray[age2].erase(stubArray[age2].begin()+index2);
                    total-=2;
                }
                else if(age1!=age2){
                    nodes[id1].iCnt.push_back(id2);
                    nodes[id2].iCnt.push_back(id1);
                    stubArray[age1].erase(stubArray[age1].begin()+index1);
                    stubArray[age2].erase(stubArray[age2].begin()+index2);
                    total-=2;
                }
                else failure++;

                //update enableAge: delete age group that has not enough stubs
                if(stubArray[age1].size()<1){
                    for(vector<int>::iterator it=enableAge.begin();it!=enableAge.end();){
                        if(*it==age1){
                            it=enableAge.erase(it);
                        }
                        else ++it;              
                    }
                }
                if(stubArray[age2].size()<1){
                    for(vector<int>::iterator it=enableAge.begin();it!=enableAge.end();){
                        if(*it==age2){
                            it=enableAge.erase(it);
                        }
                        else ++it;              
                    }
                }
            }
            else failure++;     
        }
    }

    for(int i=0;i<agentNum;i++)    
        nodes[i].iDegree=nodes[i].iCnt.size();           //update degree=iCnt size     
}

void Network::group_layer(Population &popu){
    vector<int> gSize;
    for(int i=0;i<popu.disNumg.size();i++)
        for(int k=0;k<popu.disNumg[i]*agentNum;k++)
            gSize.push_back(i+20);
    random_shuffle(gSize.begin(),gSize.end());

    //each day construct a new network
    for(int t=0;t<simuDay;t++){
        vector<vector<int>> tmp;                //temporary list of groups at day t, each group with the first element represent the type of this group and a list of members
        vector<int> hasGroup;                   //list of nodes that have groups
        vector<vector<int>> stubArray;
        vector<int> enableAge;
        int gIndex=-1;                          //id of groups
        int total=0;

        /*initial list of nodes that have groups */
        for(int i=0;i<ageClass;i++){
            stubArray.push_back({});
        }
        for(int i=0;i<agentNum;i++){
            //initial group info of each node
            int age=nodes[i].age;
            double p;
            //if this individual has group contacts
            p=Prand();
            if(p<popu.pro_group[age]){
                hasGroup.push_back(nodes[i].id);
                stubArray[age].push_back(nodes[i].id);
            }
        }
        for(int i=0;i<ageClass;i++){
            if(stubArray[i].size()>0){
                enableAge.push_back(i);
            }
            total+=stubArray[i].size();
        }

        /* construct groups */
        while(total>0){
            int age1=Uniform(enableAge.size());
                age1=enableAge[age1];
            int index1=Uniform(stubArray[age1].size());
            int id1=stubArray[age1][index1];
            
            int gDegree;                  
            vector<int> nodeID;                     //list of members in this group
            pair<int,int> tmp_gcnt;

            int index=Uniform(agentNum);      
            gDegree=gSize[index];

            //first node in the group
            gIndex++; 
            nodeID.push_back(id1);
            tmp_gcnt=make_pair(t,gIndex);
            nodes[id1].gCnt.push_back(tmp_gcnt);
            stubArray[age1].erase(stubArray[age1].begin()+index1);
            total--;
            //update enableAge[age1]
            if(stubArray[age1].size()<1){
                for(vector<int>::iterator it=enableAge.begin();it!=enableAge.end();){
                    if(*it==age1){
                        it=enableAge.erase(it);
                    }
                    else ++it;              
                }
            }

            int failure2=0;
            //choose other nodes in groups
            while(gDegree>0 && failure2<1000){
                //initial the age of the member
                int age2=0;
                double p=Prand();
                while(p>popu.disAgeg[age1][age2] && age2<ageClass-1)
                    age2++;
                
                if(stubArray[age2].size()<1)
                    failure2++;
                else{
                    //choose one node in this age group
                    int index2=Uniform(stubArray[age2].size());
                    int id2=stubArray[age2][index2];

                    nodeID.push_back(id2);
                    nodes[id2].gCnt.push_back(tmp_gcnt);
                    gDegree--;
                    total--;
                    stubArray[age2].erase(stubArray[age2].begin()+index2);
                    //update enableAge[age2]
                    if(stubArray[age2].size()<1){
                        for(vector<int>::iterator it=enableAge.begin();it!=enableAge.end();){
                            if(*it==age2){
                                it=enableAge.erase(it);
                            }
                            else ++it;              
                        }
                    }
                }
            }
            tmp.push_back(nodeID);      
        }
        groups.push_back(tmp);
    }    
}

void Network::setCnt(double proi){
    //initial contact flag
    for(int i=0;i<agentNum;i++){
        for(int j=0;j<nodes[i].iDegree;j++){
            bool flag=true;
            nodes[i].cntFlag[j]=flag;
        }
    }

    //cut a proportion of proi contacts
    if(proi>0)
        for(int i=0;i<agentNum;i++){
            for(int j=0;j<nodes[i].iDegree;j++){
                double p=Prand();
                if(p<proi && nodes[i].cntFlag[j]){      //this contact have bot been cut
                    int index=nodes[i].iCnt[j];         //the other node of this contact
                    nodes[i].cntFlag[j]=false;
                    for(int k=0;k<nodes[index].iDegree;k++){
                        if(nodes[index].iCnt[k]==i){
                            nodes[index].cntFlag[k]=false;
                            break;
                        }
                    }
                }
            }
        }
}

void Network::netInfo(Population &popu){
    string name ="../network/agentInfo.txt";
    string name1="../network/individual_contact_Info.txt";
    string name2="../network/group_contact_Info.txt";
    string name3="../network/groups.txt";

    ofstream file;
    //save inf0
    file.open(name,ios::out);
    
    file<<"id,age,iDegree,incubation,Tsc,Tcr,willSym,willTest"<<endl;
    for(int i=0;i<agentNum;i++){
        file<<nodes[i].id<<noskipws<<" "<<nodes[i].age<<" "<<nodes[i].iDegree<<" "<<nodes[i].incubation<<" "
        <<nodes[i].Tsc<<" "<<nodes[i].Tcr<<" "<<nodes[i].willSymp<<" "<<nodes[i].willTest<<endl;
    }
    file.close();

    file.open(name1,ios::out);
    file<<"individual contacts"<<endl;
    for(int i=0;i<agentNum;i++){
        for(int k=0;k<nodes[i].iDegree;k++){
            file<<noskipws<<nodes[i].iCnt[k]<<" ";
        }
        file<<endl;
    }
    file.close();

    file.open(name2,ios::out);
    file<<"id,group_day,index_group"<<endl;
    for(int i=0;i<agentNum;i++){
        file<<noskipws<<nodes[i].id<<" "<<nodes[i].gCnt.size()<<endl;
        if(nodes[i].gCnt.size()>0){   
            for(int j=0;j<nodes[i].gCnt.size();j++){
                file<<noskipws<<nodes[i].gCnt[j].first<<" "<<nodes[i].gCnt[j].second<<endl;      
            }
        }    
    }
    file.close();

    file.open(name3,ios::out);
    file<<"day,group_size,group_id,nodes"<<endl;
    for(int i=0;i<simuDay;i++){
        file<<noskipws<<i<<" "<<groups[i].size()<<endl;
        for(int j=0;j<groups[i].size();j++){
            file<<noskipws<<j<<" ";
            for(int k=0;k<groups[i][j].size();k++){
                file<<noskipws<<groups[i][j][k]<<" ";
            }
            file<<endl;
        }
    }
    file.close();
}

void Network::readNet(){
    string name ="../network/agentInfo.txt";
    string heter="../network/heterogeneity.txt";
    string name1="../network/individual_contact_Info.txt";
    string name2="../network/group_contact_Info.txt";
    string name3="../network/groups.txt";

    ifstream file;
    string line;
    vector<double> h;

    file.open(heter,ios::in);
    for(int i=0;i<agentNum;i++){
        getline(file,line);
        stringstream ss(line);
        double d;
        ss>>d;
        h.push_back(d);
    }
    file.close();

    file.open(name,ios::in);
    getline(file,line);
    while(getline(file,line)){
        stringstream ss(line);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        int d;
        vector<int> info;
        while(ss>>d)
            info.push_back(d);

        Agent a;
        a.id=info[0];
        a.age=info[1];
        a.iDegree=info[2];
        a.incubation=info[3];
        a.Tsc=info[4];
        a.Tcr=info[5];
        a.willSymp=info[6];
        a.willTest=info[7];

        if(a.age<3) a.susceptibility=kDelta[0];
        else if(a.age<13) a.susceptibility=kDelta[1];
        else a.susceptibility=kDelta[2];
        a.susceptibility/=kDelta[2];                  //normalization
        a.heterogeneity=h[a.id];
        nodes.push_back(a);
        info.clear();
    }
    file.close();

    file.open(name1,ios::in);    
    int index=0;
    getline(file,line);
    while(getline(file,line)){
        stringstream ss(line);
        int d;
        vector<int> cnt;
        while(ss>>d)
            cnt.push_back(d);
        for(int k=0;k<cnt.size();k++){
            nodes[index].iCnt.push_back(cnt[k]);
        }   
        index++;
        cnt.clear();
    }
    file.close();

    file.open(name2,ios::in);    
    getline(file,line);
    for(int i=0;i<agentNum;i++){
        getline(file,line);
        stringstream ss(line);
        int d;
        ss>>d;
        int id=d;
        ss>>d;
        int size=d;
        if(size>0){
            for(int j=0;j<size;j++){
                getline(file,line);
                stringstream ss2(line);
                int d2;
                pair<int,int> list;
                ss2>>d2;
                list.first=d2;
                ss2>>d2;
                list.second=d2;
                nodes[id].gCnt.push_back(list);
            }
        }
    }
    file.close();

    file.open(name3,ios::in);    
    getline(file,line);
    for(int i=0;i<simuDay;i++){
        getline(file,line);
        stringstream ss(line);
        int d;
        ss>>d;//simu Day
        ss>>d; 
        int size=d;

        vector<vector<int>> tmp;
        for(int j=0;j<size;j++){
            getline(file,line);
            stringstream ss(line);
            int d;
            vector<int> info;
            ss>>d;//gIndex
            while(ss>>d)
                info.push_back(d);
            tmp.push_back(info);
        }
        groups.push_back(tmp);
    }
    file.close();

    //initial contact flag
    for(int i=0;i<agentNum;i++){
        for(int j=0;j<nodes[i].iDegree;j++){
            bool flag=true;
            nodes[i].cntFlag.push_back(flag);
        }
    }  
}