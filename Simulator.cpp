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

    void SimulationInstance::computeThermalLoadParameters()
    {
        double Mbas0 = 0; // W
        double MbasDelta = 0; // W
        double MperWork = 0.2; // W/W
        idx = 0;
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            ELEMENT::Element* elem = body->elements[elemIdx];
            WASHER::Washer* cur = elem->washers[0];
            Mbas0 += cur->volume*elem->sumPhi/(2*M_PI)*cur->q_m;
            MbasDelta += cur->volume*elem->sumPhi/(2*M_PI)*sim->thermovariant*deltaQMetabolic(cur->q_m,T[idx],Tn[idx],sim->Q10);
            idx++;
            for(sectIdx = 0; sectIdx<elem->nSectors; ++sectIdx){
                SECTOR::Sector* sect = elem->sectors[sectIdx];
                for(washIdx = 1; washIdx<elem->nWashers; ++washIdx){
                    cur = elem->washers[washIdx];
                    Mbas0 += cur->volume*sect->phi/(2*M_PI)*cur->q_m;
                    MbasDelta += cur->volume*sect->phi/(2*M_PI)*sim->thermovariant*deltaQMetabolic(cur->q_m,T[idx],Tn[idx],sim->Q10);
                    idx++;
                }
            }
        }
        double act = 1 + (Sh+MbasDelta)/Mbas0*MperWork; // MET
        double etaW;
        if(act < 1.6)
            etaW = 0.05; // W/W
        else
            etaW = 0.2*tanh(body->b1*act+body->b0); // W/W
        M = act*Mbas0/body->actBas; // W
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

    void SimulationInstance::elementValues(){
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];

            flowECMOBlood(elemIdx);
            flowECMOSaline(elemIdx);
            cp(elemIdx);
            rho(elemIdx);
            tblp(elemIdx);
        }
    }

    void SimulationInstance::nodeValues(){
        for(idx=0;idx<body->N;idx++){
            q[idx] = 0; // W/m^3, Heat generation in tissue
            w[idx] = 0; // -, Blood perfusion rate
            beta[idx] = 0; // W/m^3/K, Blood perfusion rate factor (rho*c*w)
        }
        
        // Heat generation, blood perfusion, heats
        qwbeta(q,w,beta,Cp,Rho,T);
        qwbeta(qNxt,wNxt,betaNxt,CpNxt,RhoNxt,TNxt);


    }
    void SimulationInstance::agglomeratedElementValues(){
        idx = 0;
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            BV[elemIdx] = 0; // -, beta*volume over all nodes
            BVT[elemIdx] = 0; // -, beta*volume*T over all nodes
            BVTfactor[elemIdx] = 0; // K, Countercurrent heat exchange component of TblA
            BPRBPCfactor[elemIdx] = 0; // -, Complement of BVT factor, body component of TblA
            TblA[elemIdx] = 0; // K, Arterial blood temperature

            element = body->elements[elemIdx];

            // Core washer
            washer = element->washers[0];
            BV[elemIdx] += beta[idx]*washer->volume*element->sumPhi/(2*M_PI);
            BVNxt[elemIdx] += betaNxt[idx]*washer->volume*element->sumPhi/(2*M_PI);
            BVT[elemIdx] += BV[elemIdx]*T[idx];
            idx++;
            
            // Interior/skin washers
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
                for(washIdx = 1; washIdx<element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];
                    
                    BV[elemIdx] += beta[idx]*washer->volume*sector->phi/(2*M_PI);
                    BVNxt[elemIdx] += betaNxt[idx]*washer->volume*sector->phi/(2*M_PI);
                    BVT[elemIdx] += BV[elemIdx]*T[idx];
                    idx++;
                }
            }
            
            // Derived values
            BVTfactor[elemIdx] = element->hx*(BVT[elemIdx]/BV[elemIdx])/(element->hx+BV[elemIdx]); // K
            BPRBPCfactor[elemIdx] = BVNxt[elemIdx]/(element->hx+BVNxt[elemIdx])*TblPNxtRatio[elemIdx]; // -
            TblA[elemIdx] = TblP[elemIdx]*BV[elemIdx]/(element->hx+BV[elemIdx])+BVTfactor[elemIdx]; // K
        }
    }

    void SimulationInstance::agglomeratedBodyValues(){
        CplC = 0;
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];
            CplC += -(BVNxt[elemIdx]*BVNxt[elemIdx])/(element->hx+BVnxt[elemIdx])*TblPNxtRatio;
        }
    }

    void SimulationInstance::qwbeta( double *qs, double *ws, double *betas, double *cps, double *rhos, double* Ts){
        // T index
        idx = 0;
        // Loop thru elements
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];

            // Core washer
            washer = element->washers[0];
            qDm = deltaQMetabolic(washer->q_m,Ts[idx],T0[idx])*sim->thermovariant;
            qW = 0;
            qSh = 0;
            qResp = 0;
            if(element->Vmuscle > 0){
                qW = element->a_sed*H/element->Vmuscle*sim->thermovariant;
                qSh = element->a_sh*Sh/element->Vmuscle*sim->thermovariant;
            }
            if(element->Vresp > 0){
                qResp = -Qresp*washer->a_resp/(washer->volume*element->sumPhi/(2*M_PI))*sim->thermovariant;
            }
            qs[idx] = washer->q_m+qDm+qW+qSh+qResp;
            ws[idx] = max(w0[idx]*pow(sim->KonstasAlpha,sim->KonstasBeta*(Ts[idx]-T0[idx])*sim->thermovariant)*(1-sim->KonstasGamma*DeltaHCT),0.0);
            betas[idx] = (beta0[idx]+0.932*(qDm+qSh+qW))
                    *ws[idx]/w0[idx]
                    *rhos[idx]/body->rhoBlood
                    *cps[idx]/body->cpBlood;
            idx++;
            
            // Loop thru washers
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
                for(washIdx = 1; washIdx < element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];
                    qDm = deltaQMetabolic(washer->q_m,Ts[idx],T0[idx])*sim->thermovariant;
                    qW = 0;
                    qSh = 0;
                    qResp = 0;
                    if(element->Vmuscle > 0){
                        qW = element->a_sed*H/element->Vmuscle*sim->thermovariant;
                        qSh = element->a_sh*Sh/element->Vmuscle*sim->thermovariant;
                    }
                    if(element->Vresp > 0){
                        qResp = -Qresp*washer->a_resp/(washer->volume*sector->phi/(2*M_PI))*sim->thermovariant;
                    }
                    qs[idx] = washer->q_m+qDm+qW+qSh+qResp;
                    ws[idx] = max(w0[idx]*pow(sim->KonstasAlpha,sim->KonstasBeta*(Ts[idx]-T0[idx])*sim->thermovariant)*(1-sim->KonstasGamma*DeltaHCT),0.0);
                    betas[idx] = (beta0[idx]+0.932*(qDm+qSh+qW))
                            *ws[idx]/w0[idx]
                            *rhos[idx]/body->rhoBlood
                            *cps[idx]/body->cpBlood;
                    // Add vasodilation, vasoconstriction for outermost skin
                    if(washIdx == element->nWashers-1){
                        // Westin eqn 84
                        betas[idx] = (betas[idx]+element->a_dl*Dl)/(1+element->a_cs*Cs*exp(-Dl/80));
                        // Westin eqn 85
                        betas[idx] = min((386.9-.32*.932*H)*(Viv+DViv)/Viv0,betas[idx]);
                    }
                    idx++;
                }
            }
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
    void SimulationInstance::flowECMOSaline(elemIdx){
        sim->ecmoSaline(&(flowECMOSaline[elemIdx]), elemIdx, time);
        sim->ecmoSaline(&(flowECMOSalineNxt[elemIdx]), elemIdx, timeNxt);
    }
    void SimulationInstance::flowECMOBlood(elemIdx){
        sim->ecmoBlood(&(flowECMOBlood[elemIdx]), elemIdx, time);
        sim->ecmoBlood(&(flowECMOBloodNxt[elemIdx]), elemIdx, timeNxt);
    }
    void SimulationInstance::BCvalues(){
        // Shock
        BVR();
        // Hemodilution
        bloodParams();
        // Change in hematocrit
        deltaHCT();
        // Environment values
        environmentParams();
    }
    double SimulationInstance::computeMeanSkinTemp( double* T ){
        double Tskcur = 0;
        int idx=1;
        for(elemIdx = 0;elemIdx<body->nElements;elemIdx++){
            element = body->elements[elemIdx];
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
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

    void SimulationInstance::buildSystem(){
        idx = 0;
        // Loop thru elements
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];
            
            // Handle core
            coreIdx = idx;
            coreWasher = element->washers[0];
            idx++;

            // Loop thru sectors
            for(sectIdx = 0; sectIdx < element->nSectors; ++elemIdx){
                sector = element->sectors[sectIdx];
                // Loop thru washers
                for(washIdx = 1; washIdx < element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];

                    // Core adjacent nodes have the same previous node, the core
                    if(washIdx == 2){
                        // Westin eqn 57 term 2
                        pbm->add(coreIdx,idx,-element->theta*coreWasher->AForwardNxt*element->sumPhi);
                    }
                    // Other nodes have an interior node as a previous node
                    else{
                        // Westin eqn 55 term 1
                        rhs[idx] += (1-washer->gamma)*washer->ABackwardPrv*TPrv[idx-1];
                    }
                    // Forwards


            }
            // Update body-level blood props


        }
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

            // Compute body values
            // Compute whole-body thermal load parameters
            computeThermalLoadParameters();
            // Compute error inputs
            computeErrorSignals();
            // Compute active controls
            computeActiveControls();
            projectBodyValues();

            // Compute element values
            elementValues();
            projectElementValues();

            // Compute node thermal load parameters
            nodeValues();
            projectNodeValues();

            // Compute agglomerated element values
            agglomeratedElementValues();
            projectAgglomeratedElementValues();

            // Compute agglomerated body values
            agglomeratedBodyValues();
            projectAgglomeratedBodyValues();

            // Empty PBM and RHS
            clearSystem();

            // Build system
            buildSystem();

            // Solve system

            // Overwrite values with "Prv"


           

        }

    }
}
#endif