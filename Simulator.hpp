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
    /// @brief Enum for the state of a Simulator.
    enum SimulatorState { undefined, initialized, bodyLoaded, simLoaded, pbmLoaded, computed, allocated, initializerRun};


    /// @brief Simulator factory class
    class Simulator
    {
    private:
        /// @brief State of this Simulator.
        SimulatorState _state = undefined;

        // Body and sim
        /**
         * @brief BodyModel ptr object.
         * Object used to describe the body and its parameters.
         * Also includes firm biology properties like blood density.
         * Initialized once, never modified.
         */
        BODYMODEL::BodyModel* body;
        /**
         * @brief SimModel ptr object.
         * Object used to describe the BC's and inputs to simulation.
         * Should be just one for many simulation instances, UQs changing.
         */
        SIMMODEL::SimModel* sim;
        /**
         * @brief PseudoBlockMatrix ptr object.
         * Object used as the left-hand-side matrix of simulations.
         * Created and allocated just once to avoid reallocating during runs.
         */
        PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm;

        //@{
        /** Values fed to begin initializing steady state.*/
        // Initial values
        double *T0,*beta0,*w0,*q0,*Tpp0;
        double M0=0,QResp0=0,Tskm0=0,H0=0,Sh0=0,Cs0=0,Dl0=0,Sw0=0;
        //@}

    public:
        /**
         * @brief Construct a new Simulator object
         * Creates a blank Simulator object
         */
        Simulator();
        /**
         * @brief Construct a new Simulator object
         * Creates an allocated Simulator object
         * 
         * @param[in] bodyPtr Pointer to already-created body object.
         * @param[in] simPtr Pointer to already-created simmodel object
         */
        Simulator(BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr);
        /**
         * @brief Destroy the Simulator object
         * 
         */
        ~Simulator();
        /**
         * @brief Set the Body object
         * 
         * @param[in] bodyPtr Ptr to bodymodel to copy.
         */
        void setBody(BODYMODEL::BodyModel* bodyPtr);
        /**
         * @brief Set the Sim object
         * 
         * @param simPtr Ptr to simmodel to copy
         */
        void setSim(SIMMODEL::SimModel* simPtr);

        /**
         * @brief Allocates space, declares pre-initial values
         * 
         */
        void initializer();

        /**
         * @brief Finds initial conditions for the current SimModel
         * 
         * @param[in] args List of current uncertainty quantification values
         */
        void findICs( double args[] );

        /**
         * @brief Run the sim with given args and writes to given outs addresses
         * Finds ICs, begins sim, records until end condition, outputs an array of doubles
         * 
         * @param[in] args Uncertainty Quantification values
         * @param[out] outs Addresses of output parameters to fill
         */
        void runSim( double args[], double outs[] );

    };

    /// @brief Class that performs one single simulation.
    class SimulationInstance
    {
        private:
            ///@{
            /// @brief Local temporary pointers
            BODYMODEL::BodyModel* body; //!< Body object
            SIMMODEL::SimModel* sim; //!< Simulator options object
            PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm; //!< Linear matrix
            ELEMENT::Element* element; //!< Element ptr
            SECTOR::Sector* sector; //!< Sector ptr
            WASHER::Washer* washer; //!< Washer ptr
            WASHER::Washer* coreWasher; //!< Core washer ptr (When core out of scope)
            ///@}
            //@{
            // Temporary indices
            
            /// @brief Temporary indices
            int idx = -1; // Used to loop thru temperatures
            int elemIdx = -1;
            int sectIdx = -1;
            int washIdx = -1;
            int coreIdx = -1; // Used to retain index of core node
            int forwardIdx = -1;
            int backwardIdx = -1;
            //@}

            //@{
            /// @brief Previous timestep values
            // Body values
            // Thermal loads
            double TskmPrv; // K, Mean skin temperature
            double MPrv; // W, Total metabolism
            double HPrv; // W, Heat load
            double QRespPrv; // W, Respiratory cooling intake
            // Active controls:
            double ShPrv; // W, Shivering power
            double CsPrv; // -, Vasoconstriction ratio
            double DlPrv; // W/K, Vasodilation capacitance
            double SwPrv; // g/min, Sweat output
            // Node values
            double *TPrv; // Temperature (K), to be solved for
            double *TppPrv; // Temperature (K) at skin surf
            //@}

            //@{ 
            /// @brief Heats
            double qDm; // W/m^3, change in metabolism
            double qW; // W/m^3, heat gen due to workload
            double qSh; // W/m^3, heat gen due to shivering
            double qResp; // W/m^3, heat gen due to breathing
            //@}
        public:
            // Simulation attributes


            /// @brief Initial condition configuration mode. Only mode 0 is implemented.
            int ICconfigMode=0; // How to initialize (0=scratch, 1=T from file, 2=T and PBM from file)

            /**
             * @defgroup Steadys
             * Initial and steady state values use throughtout the sim's run.
             * @{
             */
            double Tskm0=NAN; //!< K, initial average skin temperature
            double Thy0=NAN; //!< K, initial hypothalamic temperature
            double Viv0=NAN; //!< m^3, initial blood volume
            double Vrbc0=NAN; //!< m^3, initial red blood cell volume
            double *T0; //!< K, initial node temperature values
            double *beta0; //!< W/m^3/K, initial perfusion rate factor
            double *w0; //!< 1/s, initial perfusion rate
            double *q0; //!< W/m^3, initial metabolic heat gen
            double *Tpp0; //!< K, initial projected temperatures (skin nodes only)
            double M0; //!< W, initial body-wide metabolism
            double QResp0; //!< W, initial body-wide respiratory heat loss
            double H0; //!< W, initial body-wide heat load
            double Sh0; //!< W, initial body-wide shivering command
            double Cs0; //!< -, initial vasoconstriction command
            double Dl0; //!< W/K, initial vasodilation command
            double Sw0; //!< g/min, initial sweating command
            /// @}
            
            /**
             * @defgroup current Simulation values at current time
             * @{
             */

            // Sim values
            double time; //!< s, simulation time
            double *rhs; //!< K, right hand of system

            // Body values
            // Thermal loads
            double Tskm; //!< K, Mean skin temperature
            double Thy; //!< K, Hypothalamic temperature
            double M; //!< W, Total metabolism
            double H; //!< W, Heat load
            double QResp; //!< W, Respiratory cooling intake
            // Active error signals
            double TskError; //!< K, Mean skin temperature error
            double ThyError; //!< K, Hypothalamic temperature error
            double TskErrorGradient; //!< K/s, Derivative of mean skin temperature error
            // Active controls:
            double Sh; //!< W, Shivering power
            double Cs; //!< -, Vasoconstriction ratio
            double Dl; //!< W/K, Vasodilation capacitance
            double Sw; //!< g/min, Sweat output
            // Whole body blood content
            double Viv; //!< Intravenous blood volume
            double Vrbc; //!< Red blood cell volume
            double DViv; //!< Added non-blood volume
            double bvr; //!< Blood volume / initial blood volume
            double svr; //!< Saline volume / initial blood volume
            double DeltaHCT; //!< Ratio of hematocrit, ability to carry oxygen

            // Element values
            // Blood values
            double *FlowECMOBlood; //!< m^3/s, Recirculated blood injection rate
            double *FlowECMOSaline; //!< m^3/s, Saline injection rate
            double *Cp; //!< J/kg/K, Blood specific heat capacity
            double *Rho; //!< kg/m^3, Blood density
            double *TblP; //!< K, experienced blood pool temperature
            double *TblPNxtRatio; //!< K/K, future experienced blood pool temperature ratio

            // Node values
            double *T; //!< Temperature (K), to be solved for
            double *Tpp; //!< Temperature (K), fictional forward skin node temperature
            double *q; //!< W/m^3, Heat generation in tissue
            double *w; //!< -, Blood perfusion rate
            double *beta; //!< W/m^3/K, Blood perfusion rate factor (rho*c*w)

            // Agglomerated element values
            double *BV; //!< -, beta*volume over all nodes
            double *BVT; //!< -, beta*volume*T over all nodes
            double *BPRBPCfactor; //!< Westin eqn 21 term 1 coefficient
            double *TblAoverlayFactor; //!< Westin eqn 21 term 2 coefficient
            double *TblAoverlay; //!< Westin eqn 21 term 2
            double *TblA; //!< K, Arterial blood temperature, Westin eqn 21 resultant

            // Agglomerated body values
            double CplC; //!< -, Blood pool timestep coupling coefficient (Bottom rightmost matrix entry)
            /// @}

            /**
             * @defgroup future Future timestep values
             * Future values used throughout the sim's run. Projected and actualized.
             * @{
             */
            double timeNxt; //!< s, simulation time at next timestep

            // Thermal loads
            double TskmNxt; //!< K, Mean skin temperature
            double MNxt; //!< W, Total metabolism
            double HNxt; //!< W, Heat load
            double QRespNxt; //!< W, Respiratory cooling intake
            // Active error signals
            double TskErrorNxt; //!< K, Error in skin mean temperature
            double ThyErrorNxt; //!< K, Error in hypothalamic mean temperature
            double TskErrorGradientNxt; //!< K/s, Derivative of mean skin error
            // Active controls:
            double ShNxt; //!< W, Shivering power
            double CsNxt; //!< -, Vasoconstriction ratio
            double DlNxt; //!< W/K, Vasodilation capacitance
            double SwNxt; //!< g/min, Sweat output
            // Whole body blood content
            double VivNxt; //!< Intravenous blood volume
            double VrbcNxt; //!< Red blood cell volume
            double DVivNxt; //!< Added non-blood volume
            double bvrNxt; //!< Blood volume / initial blood volume
            double svrNxt; //!< Saline volume / initial blood volume
            double DeltaHCTNxt; //!< Ratio of hematocrit, ability to carry oxygen

            // Blood values
            double *FlowECMOBloodNxt; //!< m^3/s, Recirculated blood injection rate
            double *FlowECMOSalineNxt; //!< m^3/s, Saline injection rate
            double *CpNxt; //!< J/kg/K, Blood specific heat capacity
            double *RhoNxt; //!< kg/m^3, Blood density

            // Node values
            double *TNxt; //!< Temperature (K), to be solved for
            double *TppNxt; //!< Temperature (K), fictional forward skin node temperature
            double *qNxt; //!< W/m^3, Heat generation in tissue
            double *wNxt; //!< -, Blood perfusion rate
            double *betaNxt; //!< W/m^3/K, Blood perfusion rate factor (rho*c*w)

            // Agglomerated element values
            double *BVNxt; //!< -, beta*volume over all nodes
            double *BPRBPCfactorNxt; //!< Westin eqn 21 term 1 coefficient
            double *TblAoverlayFactorNxt; //!< Westin eqn 21 term 2 coefficient
            double *TblAoverlayNxt; //!< Westin eqn 21 term 2
            /// @}
            

            /// @brief Create simulation instance
            /// @param tbody Initialized bodymodel ptr
            /// @param tsim Initialized simmodel ptr
            /// @param tpbm Initialized pseudoblockmatrix ptr
            SimulationInstance(BODYMODEL::BodyModel* tbody,SIMMODEL::SimModel* tsim, PSEUDOBLOCKMATRIX::PseudoBlockMatrix* tpbm);
            
            /// @brief Destroy simulation instance
            ~SimulationInstance();

            /// @brief Allocate memory.
            void allocate();
            /// @brief Deallocate memory.
            void deallocate();
            /// @brief Copy all steady values into Prv and current values
            void ghostStart();
            /// @brief Copy initial values from external arrays
            /// @param[in] T0tmp 
            /// @param[in] beta0tmp 
            /// @param[in] w0tmp 
            /// @param[in] q0tmp 
            /// @param[in] Tpp0tmp 
            void copyToInitials( double *T0tmp, double *beta0tmp, double *w0tmp, double *q0tmp, double *Tpp0tmp);
            /// @brief Fill external arrays with initial values
            /// @param[out] T0tmp 
            /// @param[out] beta0tmp 
            /// @param[out] w0tmp 
            /// @param[out] q0tmp 
            /// @param[out] Tpp0tmp 
            /// @param[out] Tskm0 
            void fillInitials( double *T0tmp, double *beta0tmp, double *w0tmp, double *q0tmp, double *Tpp0tmp, double Tskm0);
            /**
             * @brief Copy steady values from external arrays
             * 
             * @param T0tmp 
             * @param beta0tmp 
             * @param w0tmp 
             * @param q0tmp 
             * @param Tpp0tmp 
             * @param M0tmp 
             * @param QResp0tmp 
             * @param Tskm0tmp 
             * @param H0tmp 
             * @param Sh0tmp 
             * @param Cs0tmp 
             * @param Dl0tmp 
             * @param Sw0tmp 
             */
            void copyToSteadys( double *T0tmp, double *beta0tmp, double *w0tmp,double *q0tmp, double *Tpp0tmp,
                    double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
                    double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp);
            /**
             * @brief Fill external arrays with steady values
             * 
             * @param T0tmp 
             * @param beta0tmp 
             * @param w0tmp 
             * @param q0tmp 
             * @param Tpp0tmp 
             * @param M0tmp 
             * @param QResp0tmp 
             * @param Tskm0tmp 
             * @param H0tmp 
             * @param Sh0tmp 
             * @param Cs0tmp 
             * @param Dl0tmp 
             * @param Sw0tmp 
             */
            void fillSteadys( double *T0tmp, double *beta0tmp, double *w0tmp, double *q0tmp, double *Tpp0tmp,
                    double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
                    double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp);
            
            /**
             * @brief Run the simulation.
             * 
             * Performs fixed-timestep iteration until the simmodel's end condition hits or 3600s is achieved.
             * Uses PBM's version of Gauss-Seidel solving.
             * May be modified to produce iterative updates to the terminal.
             */
            void runSim();
            /**
             * @brief Clears the matrix and rhs.
             * 
             */
            void clearSystem();
            // Whole body parameters
            /**
             * @brief Computes whole-body load parameters like M,H,QResp
             * 
             */
            void computeThermalLoadParameters();
            /**
             * @brief Computes QResp
             * 
             * @return double QResp
             */
            double computeQresp();
            /**
             * @brief Computes change in metabolism at a node, following Q10 law.
             * 
             * @param q W/m^3, basal metabolism
             * @param T K, Actual node temperature
             * @param Tneut K, Thermoneutral node temperature
             * @return double W/m^3, q_DM, change in q
             */
            double deltaQMetabolic(double q, double T, double Tneut);
            // Element parameters
            /**
             * @brief Computes element-level values like TblP, Rho, Cp
             * 
             */
            void elementValues();
            /**
             * @brief Computes blood density and specific heat capacity at an element
             * 
             * @param eleIdx Element to compute
             */
            void rhocp(int eleIdx);
            /**
             * @brief Computes EXPERIENCED blood pool temperature.
             * Differs from T[body->N-1] when ECMO introduced.
             * 
             * @param eleIdx Element to compute
             */
            void tblp( int eleIdx );

            // Node parameters
            /**
             * @brief Computes all nodal values
             * 
             */
            void nodeValues();
            /**
             * @brief Computes q,w,beta at all nodes at one timestep.
             * 
             * @param[out] qs W/m^3, volumetric heat gen
             * @param[out] ws 1/s, perfusion rate
             * @param[out] betas W/m^3/K, blood perfusion rate factor
             * @param[in] cps W/kg/K, blood experienced specific heat capacity
             * @param[in] rhos kg/m^3, blood experienced density
             * @param[in] Ts K, temperatures at this timestep
             */
            void qwbeta( double *qs, double *ws, double *betas, double *cps, double *rhos, double* Ts);

            // Agglomerated element parameters
            /**
             * @brief Computes agglomerated element values
             * Some element-level values depend on nodal values, so this comes after nodeValues()
             * Computes TblA coefficients, BV/BVT dot products, BPRBPCfactors
             */
            void agglomeratedElementValues();
            // Agglomerated body parameters
            /**
             * @brief Computes agglomerated body values
             * One value CplC depends on agglomerated element values over the body
             */
            void agglomeratedBodyValues();
            // Active system
            /**
             * @brief Computes mean skin temperature
             * 
             * @return double K, Tskm
             */
            double computeMeanSkinTemp();
            /**
             * @brief Computes "hypothalamic temperature"
             * Note: Brain resolution is small, so this is node 0 in the brain
             * 
             * @return double K, Thy
             */
            double computeHypothalamicTemp();
            /**
             * @brief Computes active control signals
             * Computes TskError,ThyError,TskErrorGradient
             * 
             */
            void computeErrorSignals();
            /**
             * @brief Computes active control outputs
             * Computes Sh,Cs,Dl,Sw
             * 
             */
            void computeActiveControls();

            // BCs
            /**
             * @brief Computes everything to do with simmodel.
             * Updates simmodel's temperatures, gets UQ values.
             * 
             */
            void BCvalues();
            /**
             * @brief Gets blood and saline volume ratios (Injury)
             * 
             */
            void BVRSVR();
            /**
             * @brief Gets saline flow at an element
             * 
             * @param eleIdx Element computed
             */
            void flowECMOSaline(int eleIdx);
            /**
             * @brief Gets blood flow at an element
             * 
             * @param eleIdx Element computed
             */
            void flowECMOBlood(int eleIdx);
            /**
             * @brief Computes all ECMO parameters
             * 
             */
            void ECMOtreatment();
            /**
             * @brief Computes all blood volume parameters
             * 
             */
            void bloodParams();
            /**
             * @brief Computes change in hematocrit
             * Hematocrit is really relative concentration of red blood cells to steady state.
             * Affects blood properties as well as oxygenation. Subject to future work.
             * 
             */
            void deltaHCT();
            
            // Linear projection
            /**
             * @brief Finite difference to get future value.
             * Note: Will equal cur if steady state or initial.
             * 
             * @param cur Current value
             * @param prv Previous value
             * @return double Next value
             */
            double project( double cur, double prv );
            /**
             * @brief Get future body values
             * 
             */
            void projectBodyValues();
            /**
             * @brief ESTIMATE future temperature values.
             * Very important for skin eqns and heat gen, both termperature reliant.
             * 
             */
            void projectTemperatures();
            /**
             * @brief Project fictinal outer-skin temperatures.
             * 
             */
            void projectTpp();
            
            // Matrix construction
            /**
             * @brief Fills in linearized system
             * This function is the workhorse of the simulation.
             * Must have computed every parameter beforehand; This only affects pbm and rhs.
             * Fills in, element-by-element, first the 1D relational FDM equations.
             * Then completes the 2D arterial blood eqns from Westin.
             * Throughout, completes the final "blood pool row" for blood temp.
             */
            void buildSystem();
            /**
             * @brief Uses one of the PBM solvers to fill in TNxt with realistic values.
             * At present, uses Gauss-Seidel since it was indeed the fastest.
             * 
             */
            void solveSystem();

            // Timestep permute
            /**
             * @brief Copies values to previous timesteps in memory.
             * 
             */
            void permuteTimestep();
    };

}

#endif