#ifndef _SIMULATOR_CPP_
#define _SIMULATOR_CPP_

#include "Simulator.hpp"

namespace SIMULATOR
{

    Simulator::Simulator(){
        _state = initialized;
    }
    Simulator::Simulator(BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr) : Simulator() {
        setBody(bodyPtr);
        setSim(simPtr);
    }
    Simulator::~Simulator(){
        pbm->~PseudoBlockMatrix();
        delete [] T0;
        delete [] beta0;
        delete [] q0;
        delete [] Tpp0;
    }
    void Simulator::setBody(BODYMODEL::BodyModel* bodyPtr){
        assert(_state==initialized);
        body = bodyPtr;
        _state = bodyLoaded;
    }
    void Simulator::setSim(SIMMODEL::SimModel* simPtr){
        assert(_state==bodyLoaded);
        sim = simPtr;
        _state = simLoaded;
    }

    // Runs a simple sim case to create initial T vector, TPrv vector, etc.
    void initializer();

    // Run the sim with a given args and write to given outs addresses
    void runSim( double args[], double* outs[] );
    void SimulationInstance::computeThermalLoadParameters()
    {
        double Mbas0 = 0; // W
        double MbasDelta = 0; // W
        double MperWork = 0.2; // W/W
        int Tidx = 0;
        for(int eleIdx = 0; eleIdx < body->nElements; ++eleIdx){
            ELEMENT::Element* elem = body->elements[eleIdx];
            WASHER::Washer* cur = elem->washers[0];
            Mbas0 += cur->volume*elem->sumPhi/(2*M_PI)*cur->q_m;
            MbasDelta += cur->volume*elem->sumPhi/(2*M_PI)*sim->thermovariant*deltaQMetabolic(cur->q_m,T[Tidx],Tn[Tidx],sim->Q10);
            Tidx++;
            for(int secIdx = 0; secIdx<elem->nSectors; ++secIdx){
                SECTOR::Sector* sect = elem->sectors[secIdx];
                for(int washIdx = 1; washIdx<elem->nWashers; ++washIdx){
                    cur = elem->washers[washIdx];
                    Mbas0 += cur->volume*sect->phi/(2*M_PI)*cur->q_m;
                    MbasDelta += cur->volume*sect->phi/(2*M_PI)*sim->thermovariant*deltaQMetabolic(cur->q_m,T[Tidx],Tn[Tidx],sim->Q10);
                    Tidx++;
                }
            }
        }
        double act = 1 + (Sh+MbasDelta)/Mbas0*MperWork; // MET
        double etaW;
        if(act < 1.6)
            etaW = 0.05; // W/W
        else
            etaW = 0.2*tanh(body->b1*act+body->b0); // W/W
        M = act*Mbas0/actBas; // W
        H = M*(1-etaW)-Mbas0; // W
        Qresp = computeQresp(); // W

    }
    double SimulationInstance::computeQresp(){
        return 3.45*M*(0.028+6.5e-5*Tair-4.98e-6*sim->Pair)+
            1.44e-3*M*(32.6-0.934*Tair+1.99e-4*sim->Pair);
    }
    double SimulationInstance::deltaQMetabolic( double q, double T, double Tn)
    {
        return q*(pow(sim->Q10,(T-Tn)/10.0)-1);
    }
    double SimulationInstance::project( double cur, double prv ){
        return cur+(cur-prv)*sim->thermovariant*sim->transient;
    }
    void SimulationInstance::projectBodyValues(){
        MNxt = project(M,MPrv);
        HNxt = project(H,HPrv);
        QrespNxt = project(Qresp,QrespPrv);
        ShNxt = project(Sh,ShPrv);
        CsNxt = project(Cs,CsPrv);
        DlNxt = project(Dl,DlPrv);
        SwNxt = project(Sw,SwPrv);
    }

    void SimulationInstance::elemBloodProps(int eleIdx){
        double Tblptmp = sim->caseBloodProps(eleIdx, TblpElem, time, Viv,DViv,Viv0);
        TblpNxtRatio = Tblptmp/TblpElem;
    }

    void SimulationInstance::bloodHeatParams(int eleIdx){

    }
    void SimulationInstance::nodeValues(){
        qAndBeta(&q,&beta,Cp.Rho,T);
        qAndBeta(&qNxt,&betaNxt,CpNxt.RhoNxt,TNxt);

    }
    void SimulationInstance::qAndBeta( double **qs, double **betas, double Cp, double Rho, double* T){
        // T index
        idx = 0;
        // Loop thru elements
        for(int eleIdx = 0; eleIdx < body->nElements; ++eleIdx){
            element = body->elements[eleIdx];
            washer = element->washers[0];
            qDm = deltaQMetabolic(washer->q_m,T[idx],T0[idx])*sim->thermovariant;
            qW = 0;
            qSh = 0;
            qResp = 0;
            if(element->Vmuscle > 0){
                qW = element->a_sed*H/element->Vmuscle*washer->volume/element->Vmuscle*sim->thermovariant;
                qSh = element->a_sh*Sh/element->Vmuscle*washer->volume/element->Vmuscle*sim->thermovariant;
            }
            if(element->Vresp > 0){
                qResp = Qresp*washer->a_resp/(washer->volume*element->sumPhi/(2*M_PI))*sim->thermovariant;
            }
            *qs[idx] = washer->q_m+qDm+qW+qSh+qResp;
            // Loop thru washers
            for(int washIdx = 0; washIdx < element->nWashers; ++washIdx){
                
            }
        
            // Update body-level blood props


        }
    }
    void SimulationInstance::BVR(){
        sim->BVR(&bvr, &svr, time);
        sim->SVR(&bvrNxt,&svrNxt,timeNxt);
    }
    void SimulationInstance::bloodParams(){
        Viv = body->Viv0*bvr;
        VivNxt = body->Viv0*bvrNxt;

        Vrbc = body->Vrbc0*bvr;
        VrbcNxt = body->Vrbc0*bvrNxt;

        Dviv = body->Viv0*svr;
        DvivNxt = body->Viv0*svrNxt;
    }
    void SimulationInstance::deltaHCT(){
        deltaHCT = (body->p*Vrbc*DViv)/(Viv*(Viv+body->p*DViv));
    }
    void SimulationInstance::flowECMOSaline(){
        sim->ecmoSaline(&flowECMOSaline,time);
        sim->ecmoSaline(&flowECMOSalineNxt,timeNxt);
    }
    void SimulationInstance::flowECMOBlood(){
        sim->ecmoBlood(&flowECMOBlood,time);
        sim->ecmoBlood(&flowECMOBloodNxt,timeNxt);
    }
    void SimulationInstance::BCvalues(){
        // Shock
        BVR();
        // Hemodilution
        bloodParams();
        // Change in hematocrit
        deltaHCT();
        // Saline flows
        flowECMOSaline();
        // Recirculated blood flow
        flowECMOBlood();
        // Environment values
        environmentParams();
    }
    double SimulationInstance::computeMeanSkinTemp(){
        double Tskcur = 0;
        int idx=1;
        for(int eleIdx = 0;eleIdx<body->nElements;eleIdx++){
            element = body->elements[eleIdx];
            for(int secIdx=0;secIdx<element->nSectors;secIdx++){
                sector = element->sectors[secIdx];
                idx+=element->nWashers-2;
                Tskcur += element->a_sk*sector->phi/element->sumPhi*
                        T[idx];
            }
            idx++;
        }
        return Tskcur;
    }
    double SimulationInstance::computeHypothalamicTemp(){
        return T[0];
    }
    void SimulationInstance::computeErrorSignals(){
        // Tsk,m
        Tskm = computeMeanSkinTemp();
        Thy = computeHypothalamicTemp();
        TskError = (Tskm-Tskm0)*sim->transient*sim->thermovariant;
        ThyError = (Thy-Thy0)*sim->transient*sim->thermovariant;
        TskErrorGradient = (Tskm-TskmPrv)/dt/3600.0*sim->transient*sim->thermovariant;
    }

    void SimulationInstance::computeActiveControls(){
        double TskErrorDot = 0;
        if(TskError <=0 && TskErrorGradient <=0){
            TskErrorDot = TskError*TskErrorGradient;
        }
        Sh = 10*(tanh(.48*TskError+3.62)-1)*TskError
            -27.9*ThyError
            +1.7*TskErrorDot
            -28.6; // W
        Cs = 35*(tanh(.34*TskError+1.07)-1)*TskError
            +3.9*TskErrorDot; // -
        Dl = 21*(tanh(.79*TskError-.70)+1)*TskError
            +32*(tanh(3.29*ThyError-1.46)+1)*ThyError; // W/K
        Sw = (.8*tanh(.59*TskError-.19)+1.2)*TskError
            +(5.7*tanh(1.98*ThyError-1.03)+6.3)*ThyError; // g/min
        Sh = min(350.0*(svr+bvr),max(0.0,Sh)); // W      
        Sw = min(30.0,max(Sw,0.0)); // g/min
    }

    void SimulationInstance::clearSystem(){
        pbm->fill(0.0);
        for(int i=0;i<body->N;i++)
            rhs[i] = 0;
    }

    void SimulationInstance::runSim(){
        
        // Check status of body and simmodel. Make sure both completed.
        assert(body->getState() == BODYMODEL::computed);
        assert(sim->getState() == SIMMODEL::allValuesAssigned);

        // Provide initial values

        time = 0;

        // Iterate
        for(int timestep=1;timestep<nSteps;++timestep){
            
            // Compute temporary BC values and properties from sim
            timeNxt = time + dt;
            BCvalues();

            // Compute whole-body thermal load parameters
            computeThermalLoadParameters();

            // Compute error inputs
            computeErrorSignals();

            // Compute active controls
            computeActiveControls();

            // Project to get future body values
            projectBodyValues();

            // Compute node thermal load parameters
            nodeValues();

            // Project to get future node values
            projectNodeValues();

            // Empty PBM and RHS
            clearSystem();

            

            // Loop thru elements
            for(int eleIdx = 0; eleIdx < body->nElements; ++eleIdx){
                element = body->elements[eleIdx];
                // Compute element-level blood pool properties
                elemBloodProps(eleIdx);
                
                // Loop thru sectors
                for(int sectIdx = 0; sectIdx < element->nSectors; ++eleIdx){
                    sector = element->sectors[sectIdx];
                    // Loop thru washers
                }
                // Update body-level blood props


            }

        }

    }
}
#endif