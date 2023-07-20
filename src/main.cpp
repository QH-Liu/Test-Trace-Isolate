#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "../include/agent.h"
#include "../include/Population.h"
#include "../include/Network.h"
#include "../include/Distribution.h"
#include "../include/Calibration.h"
#include "../include/Simulation.h"
#include "../include/Test.h"

using namespace std;

int main(){
    clock_t start,end;

    start=clock();          //start time

    /* ------------examples------------------

    initial_net();          //generate network

    degree();               //statics: degrees of the network

    calibration(5);         //calibration Tg and R0

    for(int i=0;i<5;i++)    //simulations under npi_type: 0 to 4
        npi_type(i);        

    npi_5_para();           //simulations of TTI_relax

    npi_6_para();           //simulations of TTI_strict

    test_ptrace();          //simulations with different values of probability to trace temporal contacs   
    
    test_dtrace();          //simulations with different values of days to trace temporal contacs
    
    test_tsc();             //simulations with different values of time from symptom to sample collection
    
    test_tcr();             //simulations with different values of time from sample collection to result

    test_kseed();           //simulations with different values of initial infection
    
    test_ptest();           //simulations with different values of probability to test symptomatic infection
    
    test_RelativeI();       ////simulations with different values of relative infectiousness
 */
    end=clock();            //end time

    cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<endl;
    
    return 0;
}