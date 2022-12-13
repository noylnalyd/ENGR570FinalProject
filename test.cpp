#include "Simulator.cpp"
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <chrono>
int main(void){
    double *args = new double[4];
    args[0] = 20+273.15; // TsrmIndoors, indoor surroundings surface temp (K)
    args[1] = 20+273.15; //TsrmOutdoors, outdoor surroundings surface temp (K)
    args[2] = 300; // trecovery, time until arrival of ambulance (seconds)
    args[3] = 7+273.15; // TECMO, temp of of blood-saline mixture from machine (K)
    double *outs = new double[1];

    BODYMODEL::BodyModel* body = BODYMODEL::defaultBody();
    SIMMODEL::SimModel* sim = new SIMMODEL::KonstasCase();

    SIMULATOR::Simulator* simulator = new SIMULATOR::Simulator(body,sim);
    simulator->initializer();
    simulator->runSim(args,outs);

    cout << outs[0] << endl;
    
    delete simulator;
    delete sim;
    //delete body;
    delete [] args;
    delete [] outs;
}