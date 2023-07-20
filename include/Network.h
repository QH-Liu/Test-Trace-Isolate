#ifndef BASELINE_NETWORK_H
#define BASELINE_NETWORK_H
#include "Agent.h"
#include "Population.h"
#include <utility>

using namespace std;

const int kFailure=100000;
const int simuDay=400;

class Network{
    public:
        int ageClass;                                                   //number of classes of age
        int agentNum;                                                   //number of agents in regular individual layer   

        vector<Agent> nodes;                                            //agents
        vector<vector<vector<int>>> groups;                             //each day, matrix of groups, list of nodes in each group                                              

        Network(int nClass, int nAgent);
        ~Network()=default;

        void setRecovery(string kMu);                                   //set recoverey of each node
        void individual_layer(Population &popu);                        //set individual layer
        void group_layer(Population &popu);                             //set group layer
        void setCnt(double proi);                                       //set individual contact flag
        void netInfo(Population &popu);                                 //save network info
        void readNet();                                                 //initial network from input file
};

#endif