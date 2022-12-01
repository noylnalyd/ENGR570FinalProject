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

        // Body and sim
        BODYMODEL::BodyModel* body; // Body object
        SIMMODEL::SimModel* sim; // Simulator options object
        PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm; // Linear matrix

        // Initial values
        double *T0,*beta0,*q0,*Tpp0;
        double M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0;

    public:
        Simulator();
        Simulator(BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr);
        ~Simulator();
        void setBody(BODYMODEL::BodyModel* bodyPtr);
        void setSim(SIMMODEL::SimModel* simPtr);

        // Runs a simple sim case to create initial values and pbm
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
            int ICconfigMode=0; // How to initialize (0=scratch, 1=T from file, 2=T and PBM from file)
            
            BODYMODEL::BodyModel* body; // Body object
            SIMMODEL::SimModel* sim; // Simulator options object
            PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm; // Linear matrix
            
            // Local temporary pointers
            ELEMENT::Element* element;
            SECTOR::Sector* sector;
            WASHER::Washer* washer;
            WASHER::Washer* coreWasher;
            // Temporary indices
            int idx = NAN; // Used to loop thru temperatures
            int elemIdx = NAN;
            int sectIdx = NAN;
            int washIdx = NAN;
            int coreIdx = NAN; // Used to retain index of core node
            int forwardIdx = NAN;
            int backwardIdx = NAN;

            // Sim initial values
            double Tskm0=NAN; // K, initial average skin temperature
            double Thy0=NAN; // K, initial hypothalamic temperature
            double Viv0=NAN; // m^3, initial blood volume
            double Vrbc0=NAN; // m^3, initial red blood cell volume
            double *T0,*beta0,*w0,*q0,*Tpp0;
            double M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0;

            // Sim values per timestep
            // Sim values
            double time; // s, simulation time
            double* rhs; // K, right hand of system

            // Body values
            // Thermal loads
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
            // Whole body blood content
            double Viv; // Intravenous blood volume
            double Vrbc; // Red blood cell volume
            double DViv; // Added non-blood volume
            double bvr; // Blood volume / initial blood volume
            double svr; // Saline volume / initial blood volume
            double DeltaHCT; // Ratio of hematocrit, ability to carry oxygen

            // Element values
            // Blood values
            double *FlowECMOBlood; // m^3/s, Recirculated blood injection rate
            double *FlowECMOSaline; // m^3/s, Saline injection rate
            double *Cp; // J/kg/K, Blood specific heat capacity
            double *Rho; // kg/m^3, Blood density
            double *TblP; // K, experienced blood pool temperature
            double *TblPNxtRatio; // K/K, future experienced blood pool temperature ratio

            // Node values
            double *T; // Temperature (K), to be solved for
            double *Tpp; // Temperature (K), fictional forward skin node temperature
            double *q; // W/m^3, Heat generation in tissue
            double *w; // -, Blood perfusion rate
            double *beta; // W/m^3/K, Blood perfusion rate factor (rho*c*w)

            // Agglomerated element values
            double *BV; // -, beta*volume over all nodes
            double *BVT; // -, beta*volume*T over all nodes
            double *BPRBPCfactor; // Westin eqn 21 term 1 coefficient
            double *TblAoverlayFactor; // Westin eqn 21 term 2 coefficient
            double *TblAoverlay; // Westin eqn 21 term 2
            double *BPRBPCfactor; // Westin eqn 21 term 1 coefficient
            double *TblA; // K, Arterial blood temperature, Westin eqn 21 resultant

            // Agglomerated body values
            double CplC; // -, Blood pool timestep coupling coefficient (Bottom rightmost matrix entry)
            
            // Heats
            double qDm; // W/m^3, change in metabolism
            double qW; // W/m^3, heat gen due to workload
            double qSh; // W/m^3, heat gen due to shivering
            double qResp; // W/m^3, heat gen due to breathing



            void runSim();
            // Linear system maintenance
            void clearSystem();
            // Whole body parameters
            void computeThermalLoadParameters();
            double computeQresp();
            double deltaQMetabolic(double q, double T, double TNxt);
            // Element parameters
            void elementValues();
            void flowECMOSaline(int eleIdx);
            void flowECMOBlood(int eleIdx);
            void cp(int eleIdx);
            void rho(int eleIdx);

            // Node parameters
            void nodeValues();
            void qwbeta( double *qs, double *ws, double *betas, double *cps, double *rhos, double* Ts);
            void skinT();

            // Agglomerated element parameters
            void agglomeratedElementValues();
            // Agglomerated body parameters
            void agglomeratedBodyValues();
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
            
            // Linear projection
            double project( double cur, double prv );
            void projectBodyValues();
            void projectElementValues();
            void projectNodeValues();
            
            // Matrix construction
            void buildSystem();
            void solveSystem();
    };

}

#endif