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
    enum SimulatorState { undefined, initialized, allocated, bodyLoaded, simLoaded, precomputed, simrunning, simrun, output };


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
        BODYMODEL::BodyModel* body; // Body object
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

        Simulator(int simCase){
            _state = initialized;
            nSteps = (int)ceil((tFinal-tInitial)/dt);

            times = new double[nSteps];
            Tskm = new double[nSteps];
            Thy = new double[nSteps];
            Tblp = new double[nSteps];
            H = new double[nSteps];
            Sh = new double[nSteps];
            Cs = new double[nSteps];
            Dl = new double[nSteps];
            Sw = new double[nSteps];

            body = new BODYMODEL::BodyModel();
            
            switch(simCase){
                case 1:{ // Steady
                    sim = new SIMMODEL::SteadyCase();
                }
                case 2:{ // Transient, static
                    sim = new SIMMODEL::TransientCase();
                }
                case 3:{ // Injury
                    sim = new SIMMODEL::InjuryCase();
                }
                case 4:{ // Konstas
                    sim = new SIMMODEL::KonstasCase();
                }
                case 5:{ // Diao
                    sim = new SIMMODEL::DiaoCase();
                }

            }
        }
    };

}

#endif