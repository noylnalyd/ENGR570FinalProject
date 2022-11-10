#ifndef _SIMULATOR_
#define _SIMULATOR_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include "BodyModel.hpp"
#include "SimModel.hpp"

using namespace std;

namespace SIMULATOR
{
    enum SimulatorState { undefined, initialized, bodyLoaded, simLoaded, precomputed, simrunning, simrun, output };


    class Simulator
    {
    private:
        SimulatorState _state = undefined;
    public:
        // Simulation attributes
        const double tInitial = 0.0; // s
        double tFinal = 3600.0; // s
        double dt = 1e-2; // s
        int nSteps; // -
        double* times; // -

        // SimModel attributes
        int ICconfigMode; // How to initialize (0=scratch, 1=T from file, 2=T and PBM from file)
        BODYMODEL::BodyModel body; // Body object
        SIMMODEL::SimModel* sim; // Simulator options object


        // Values to record over sym
        double* Tskm; // K, Mean skin temperature
        double* Thy; // K, Hypothalamic temperature
        double* Tblp; // K, Blood pool temperature
        double* H; // W, Heat load

        // Active controls, also recorded
        double* Sh; // W, Shivering power
        double* Cs; // -, Vasoconstriction ratio
        double* Dl; // W/K, Vasodilation capacitance
        double* Sw; // g/min, Sweat output
    };

}

#endif