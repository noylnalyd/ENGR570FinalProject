#ifndef _SIMULATOR_
#define _SIMULATOR_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include "BodyModel.cpp"
#include "SimModel.hpp"

namespace SIMULATOR
{
    enum SimulatorState { undefined, initialized, bodyLoaded, simLoaded, pbmLoaded, computed, allocated, initializerRun};


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
        // Steady values
        double M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0;

    public:
        Simulator();
        Simulator(BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr);
        ~Simulator();
        void setBody(BODYMODEL::BodyModel* bodyPtr);
        void setSim(SIMMODEL::SimModel* simPtr);

        // Allocates space
        void initializer();

        // Finds initial conditions for the current sim
        void findICs( double args[] );

        // Run the sim with a given args and write to given outs addresses
        void runSim( double args[], double outs[] );


    };

    class SimulationInstance
    {
        private:
            // Local temporary pointers
            BODYMODEL::BodyModel* body; // Body object
            SIMMODEL::SimModel* sim; // Simulator options object
            PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm; // Linear matrix
            ELEMENT::Element* element;
            SECTOR::Sector* sector;
            WASHER::Washer* washer;
            WASHER::Washer* coreWasher;
            // Temporary indices
            int idx = -1; // Used to loop thru temperatures
            int elemIdx = -1;
            int sectIdx = -1;
            int washIdx = -1;
            int coreIdx = -1; // Used to retain index of core node
            int forwardIdx = -1;
            int backwardIdx = -1;

            // Previous timestep values

            // Body values
            // Thermal loads
            double TskmPrv; // K, Mean skin temperature
            double MPrv; // W, Total metabolism
            double HPrv; // W, Heat load
            double QrespPrv; // W, Respiratory cooling intake
            // Active controls:
            double ShPrv; // W, Shivering power
            double CsPrv; // -, Vasoconstriction ratio
            double DlPrv; // W/K, Vasodilation capacitance
            double SwPrv; // g/min, Sweat output
            // Node values
            double *TPrv; // Temperature (K), to be solved for
            double *TppPrv; // Temperature (K) at skin surf

            // Heats
            double qDm; // W/m^3, change in metabolism
            double qW; // W/m^3, heat gen due to workload
            double qSh; // W/m^3, heat gen due to shivering
            double qResp; // W/m^3, heat gen due to breathing
        public:
            // Simulation attributes

            // Sim attributes
            int ICconfigMode=0; // How to initialize (0=scratch, 1=T from file, 2=T and PBM from file)

            // Sim initial values
            double Tskm0=NAN; // K, initial average skin temperature
            double Thy0=NAN; // K, initial hypothalamic temperature
            double Viv0=NAN; // m^3, initial blood volume
            double Vrbc0=NAN; // m^3, initial red blood cell volume
            double *T0,*beta0,*w0,*q0,*Tpp0;
            double M0,QResp0,H0,Sh0,Cs0,Dl0,Sw0;

            // Sim values per timestep
            // Sim values
            double time; // s, simulation time
            double* rhs; // K, right hand of system

            // Body values
            // Thermal loads
            double Tskm; // K, Mean skin temperature
            double Thy; // K, Hypothalamic temperature
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
            double *TblA; // K, Arterial blood temperature, Westin eqn 21 resultant

            // Agglomerated body values
            double CplC; // -, Blood pool timestep coupling coefficient (Bottom rightmost matrix entry)

            // Next timestep values
            double timeNxt; // s, simulation time at next timestep

            // Body values
            // Thermal loads
            double TskmNxt; // K, Mean skin temperature
            double MNxt; // W, Total metabolism
            double HNxt; // W, Heat load
            double QrespNxt; // W, Respiratory cooling intake
            // Active error signals
            double TskErrorNxt;
            double ThyErrorNxt;
            double TskErrorGradientNxt;
            // Active controls:
            double ShNxt; // W, Shivering power
            double CsNxt; // -, Vasoconstriction ratio
            double DlNxt; // W/K, Vasodilation capacitance
            double SwNxt; // g/min, Sweat output
            // Whole body blood content
            double VivNxt; // Intravenous blood volume
            double VrbcNxt; // Red blood cell volume
            double DVivNxt; // Added non-blood volume
            double bvrNxt; // Blood volume / initial blood volume
            double svrNxt; // Saline volume / initial blood volume
            double DeltaHCTNxt; // Ratio of hematocrit, ability to carry oxygen

            // Element values
            // Blood values
            double *FlowECMOBloodNxt; // m^3/s, Recirculated blood injection rate
            double *FlowECMOSalineNxt; // m^3/s, Saline injection rate
            double *CpNxt; // J/kg/K, Blood specific heat capacity
            double *RhoNxt; // kg/m^3, Blood density

            // Node values
            double *TNxt; // Temperature (K), to be solved for
            double *TppNxt; // Temperature (K), fictional forward skin node temperature
            double *qNxt; // W/m^3, Heat generation in tissue
            double *wNxt; // -, Blood perfusion rate
            double *betaNxt; // W/m^3/K, Blood perfusion rate factor (rho*c*w)

            // Agglomerated element values
            double *BVNxt; // -, beta*volume over all nodes
            double *BPRBPCfactorNxt; // Westin eqn 21 term 1 coefficient
            double *TblAoverlayFactorNxt; // Westin eqn 21 term 2 coefficient
            double *TblAoverlayNxt; // Westin eqn 21 term 2            
            

            // Initialize
            SimulationInstance(BODYMODEL::BodyModel* tbody,SIMMODEL::SimModel* tsim, PSEUDOBLOCKMATRIX::PseudoBlockMatrix* tpbm);
            void copyToInitials( double *T0tmp, double *beta0tmp,double *q0tmp, double *Tpp0tmp);
            void fillInitials( double *T0tmp, double *beta0tmp, double *q0tmp, double *Tpp0tmp);
            void copyToSteadys( double *T0tmp, double *beta0tmp, double *q0tmp, double *Tpp0tmp,
                    double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
                    double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp);
            void fillSteadys( double *T0tmp, double *beta0tmp, double *q0tmp, double *Tpp0tmp,
                    double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
                    double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp);
            
            void runSim();
            // Linear system maintenance
            void clearSystem();
            // Whole body parameters
            void computeThermalLoadParameters();
            double computeQresp();
            double deltaQMetabolic(double q, double T, double TNxt);
            // Element parameters
            void elementValues();
            void cp(int eleIdx);
            void rho(int eleIdx);
            void tblp( int eleIdx );

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
            void BVRSVR();
            void flowECMOSaline(int eleIdx);
            void flowECMOBlood(int eleIdx);
            void ECMOtreatment();
            void bloodParams();
            void deltaHCT();
            
            // Linear projection
            double project( double cur, double prv );
            void projectBodyValues();
            void projectNodeValues();
            
            // Matrix construction
            void buildSystem();
            void solveSystem();

            // Timestep permute
            void permuteTimestep();
    };

}

#endif