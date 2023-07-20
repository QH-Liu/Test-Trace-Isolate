#ifndef BASE_LINE_POPULATION_H
#define BASE_LINE_POPULATION_H
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "Agent.h"
#include "Distribution.h"

using namespace std;

const double kAvg=18.9;                                                 //avg num of all the contacts of the survey
const double wi=1;                                                      //individual contacts relative risk=2.84;
const double wg=0.352;                                                  //group contact relative risk=1;
// const double wg=0.176;                                               //relative risk=0.5;
// const double wg=0.704;                                               //relative risk=2;

class Population
{
    public:
        int ageClass;                                                   //number of classes of age
        int agentNum;                                                   //number of agents in regular individual layer   

        vector<vector<double>> disNumi;                                 //distribution of the num of the individual contacts
        vector<vector<double>> disAgei;                                 //distribution of the age of the individual contacts
        vector<double> disNumg;                                         //distribution of the num of the group contacts
        vector<vector<double>> disAgeg;                                 //distribution of the age of the group contacts
        vector<double> pro_group;                                       //proportion of each age group that have group contacts

        vector<int> ageSize;                                            //size of each age group of the population
        vector<vector<int>> ageGroup;                                   //id of nodes in each group

        vector<int> listAge;                                            //list of node age
        vector<int> listIdegree;                                        //list of node idegree
        vector<double> listH;                                           //list of node heterogeneity

        Population(int n);
        ~Population()=default;
        void initial();                                                 //initial population               
};

#endif