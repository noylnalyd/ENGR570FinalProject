#include <iostream>
#include <math.h>
#include <cstdlib>
#include "Simulator.hpp"

using namespace SIMULATOR;
using namespace BODYMODEL;
using namespace SIMMODEL;

int main(int argc, char *argv[]){

    double t0 = 27;
    double t1 = 28;
    double t2 = 29;
    double t3 = 30;

    Simulator* mySim = new Simulator();
    BodyModel* myBody = defaultBody();
    SimModel* mySimModel = new SimModel();
    mySim->initializer();
    

}