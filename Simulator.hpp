#ifndef _SIMULATOR_
#define _SIMULATOR_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include "BodyModel.hpp"
#include "SimModel.hpp"

namespace SIMULATOR
{
    enum SimulatorState { undefined, initialized, bodyLoaded, simLoaded, allocated, precomputed, simrunning, simrun, output };


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


        // Values used in whole simulation.
        double Tskm; // K, Mean skin temperature
        double Thy; // K, Hypothalamic temperature
        double Tblp; // K, Blood pool temperature
        double M; // W, Total metabolism
        double H; // W, Heat load
        double Qresp; // W, Respiratory cooling intake
        // Active controls:
        double Sh; // W, Shivering power
        double Cs; // -, Vasoconstriction ratio
        double Dl; // W/K, Vasodilation capacitance
        double Sw; // g/min, Sweat output

        // Holdover past values
        double TskmPrv; // K, Mean skin temperature
        double ThyPrv; // K, Hypothalamic temperature
        double TblpPrv; // K, Blood pool temperature
        double MPrv; // W, Total metabolism
        double HPrv; // W, Heat load
        double QrespPrv; // W, Respiratory cooling intake
        // Active controls:
        double ShPrv; // W, Shivering power
        double CsPrv; // -, Vasoconstriction ratio
        double DlPrv; // W/K, Vasodilation capacitance
        double SwPrv; // g/min, Sweat output

        // Projected future values
        double TskmNxt; // K, Mean skin temperature
        double ThyNxt; // K, Hypothalamic temperature
        double TblpNxt; // K, Blood pool temperature
        double MNxt; // W, Total metabolism
        double HNxt; // W, Heat load
        double QrespNxt; // W, Respiratory cooling intake
        // Active controls:
        double ShNxt; // W, Shivering power
        double CsNxt; // -, Vasoconstriction ratio
        double DlNxt; // W/K, Vasodilation capacitance
        double SwNxt; // g/min, Sweat output

        // Useful values, past and present
        double* T; // K
        double* T0; // K
        double* Tprv; // K
        double Mprv; // W
        double* QRespprv, // W
        double* Tskprv, // K
        double* qprv, // W/m^3
        double* betaprv, // K
        double* Tppprv, // K
        double* Tprv, // K
        double* Tprv, // K
        double* Tprv, // K
        double* Tprv, // K
        double* Tprv, // K

        // Values to record over sim
        double* HistTskm; // K, Mean skin temperature
        double* HistThy; // K, Hypothalamic temperature
        double* HistTblp; // K, Blood pool temperature
        double* HistH; // W, Heat load
        double* HistSh; // W, Shivering power
        double* HistCs; // -, Vasoconstriction ratio
        double* HistDl; // W/K, Vasodilation capacitance
        double* HistSw; // g/min, Sweat output

        void preallocate(
            double* Tprv, // K
            double* Mprv, // W
            double* QRespprv, // W
            double* Tskprv, // K
            double* qprv, // W/m^3
            double* betaprv, // K
            double* Tppprv, // K
            double* Tprv, // K
            double* Tprv, // K
            double* Tprv, // K
            double* Tprv, // K
            double* Tprv, // K
        );
        void runSim();

        Simulator(int simCase, BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr){
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

            body = bodyPtr;
            sim = simPtr;
            
            _state = allocated;
        }
        void computeThermalLoadParameters(double* M, double* H, double*Q);
    };

}

#endif