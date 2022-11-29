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
    enum SimulatorState { undefined, initialized, bodyLoaded, simLoaded, computed, allocated, initializerRun};


    class Simulator
    {
    private:
        SimulatorState _state = undefined;
        // Initial values
        double* T0,beta0,q0,Tpp0;
        double M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0;


    public:
        Simulator();
        Simulator(BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr);
        ~Simulator();

        // Runs a simple sim case to create initial T vector, TPrv vector, etc.
        void initializer();

        // Run the sim with a given args and write to given outs addresses
        void runSim( double args[], double* outs[] );

        
    };

    class SimulationInstance
    {
        public:
            // Simulation attributes
            const double tInitial = 0.0; // s
            double tFinal = 3600.0; // s
            double dt = 1e-2; // s
            int nSteps = (int)ceil((tFinal-tInitial)/dt); // -
            double* times; // -

            // Sim attributes
            int ICconfigMode; // How to initialize (0=scratch, 1=T from file, 2=T and PBM from file)
            BODYMODEL::BodyModel* body; // Body object
            SIMMODEL::SimModel* sim; // Simulator options object
            PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm; // Linear matrix

            
            // Local temporary pointers
            ELEMENT::Element* element;
            SECTOR::Sector* sector;
            WASHER::Washer* washer;

            // Sim initial values
            double Tskm0=NAN; // K, initial average skin temperature
            double Thy0=NAN; // K, initial hypothalamic temperature
            double Viv0=NAN; // m^3, initial blood volume
            double Vrbc0=NAN; // m^3, initial red blood cell volume

            // Sim values per timestep
            // Sim values
            double time; // s, simulation time
            double* rhs; // K, right hand of system
            // Body values
            double Tskm; // K, Mean skin temperature
            double Thy; // K, Hypothalamic temperature
            double Tblp; // K, Blood pool temperature
            double M; // W, Total metabolism
            double H; // W, Heat load
            double Qresp; // W, Respiratory cooling intake
            // Active error signals
            double TskError;
            double ThyError;
            double TskErrorGradient;
            // Active controls:
            double Sh; // W, Shivering power
            double Cs; // -, Vasoconstriction ratio
            double Dl; // W/K, Vasodilation capacitance
            double Sw; // g/min, Sweat output
            // Misc
            double Viv; // Intravenous blood volume
            double Vrbc; // Red blood cell volume
            double DViv; // Added non-blood volume
            double bvr; // Blood volume / initial blood volume
            double svr; // Saline volume / initial blood volume
            double deltaHCT; // Ratio of hematocrit, ability to carry oxygen
            

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
            double* QRespprv; // W
            double* Tskprv; // K
            double* qprv; // W/m^3
            double* betaprv; // K
            double* Tppprv; // K

            void runSim();
            // Linear system maintenance
            void clearSystem();
            // Thermal body loads
            void computeThermalLoadParameters();
            double computeQresp();
            double deltaQMetabolic();
            // Active system
            double computeMeanSkinTemp();
            double computeHypothalamicTemp();
            void computeErrorSignals();
            void computeActiveControls();
            // BCs
            void BCvalues();
            void BVR();
            void bloodParams();
            void deltaHCT();
            void flowECMOSaline();
            void flowECMOBlood();
            // Linear projection
            double project( double cur, double prv );
            void projectValues();
            
            // Matrix construction
            // Element level
            void elemBloodProps( int eleIdx );
            // Node level
            void heat();
    };

}

#endif