#ifndef _SIMULATOR_CPP_
#define _SIMULATOR_CPP_

#include "Simulator.hpp"

namespace SIMULATOR
{

    Simulator::Simulator(){
        _state = initialized;
    }
    Simulator::Simulator( BODYMODEL::BodyModel* bodyPtr, SIMMODEL::SimModel* simPtr ) : Simulator() {
        setBody(bodyPtr);
        setSim(simPtr);
        int *blockN = new int[body->nElements],
            *blockM = new int[body->nElements];
        for(int i=0;i<body->nElements;i++){
            blockN[i] = body->elements[i]->N;
            blockM[i] = body->elements[i]->N;
        }
        PSEUDOBLOCKMATRIX::PseudoBlockMatrix* pbm = new PSEUDOBLOCKMATRIX::PseudoBlockMatrix(body->nElements,blockN,blockM);
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

    void Simulator::initializer(){
        // Allocate initial matrices
        T0 = new double[body->N];
        beta0 = new double[body->N];
        q0 = new double[body->N];
        Tpp0 = new double[body->N];
    }
    void Simulator::findICs( double args[]){
        for(int i=0;i<body->N;i++){
            T0[i] = 33+273.15; // K
            Tpp0[i] = 27+273.15; // K
        }
        // Populate beta, q, and prvs for steady case
        SIMMODEL::InitialCase* simInit = new SIMMODEL::InitialCase();
        simInit->setUQs(args);
        SimulationInstance* siInit = new SimulationInstance(body,simInit,pbm);
        siInit->runSim();
        siInit->copyToInitials(T0,beta0,q0,Tpp0);

        // Now find the steady case
        SIMMODEL::SteadyCase* simSteady = new SIMMODEL::SteadyCase();
        simSteady->setUQs(args);
        SimulationInstance* siSteady = new SimulationInstance(body,simSteady,pbm);
        siSteady->fillInitials(T0,beta0,q0,Tpp0);
        siSteady->runSim();
        siSteady->copyToSteadys(T0,beta0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);

        // Now normalize by using the transient, no injury case with active controls
        SIMMODEL::TransientCase* simTransient = new SIMMODEL::TransientCase();
        simTransient->setUQs(args);
        SimulationInstance* siTransient = new SimulationInstance(body,simTransient,pbm);
        siTransient->fillSteadys(T0,beta0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        siTransient->runSim();
        siTransient->copyToSteadys(T0,beta0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);

        // Free vars
        delete simInit;
        delete simSteady;
        delete simTransient;
        delete siInit;
        delete siSteady;
        delete siTransient;
    }

    void Simulator::runSim( double args[], double outs[] )
    {
        // Modify simmodel with args
        // UQ attributes worth exploring for all cases
        sim->setUQs(args);
        // Create instance
        SimulationInstance* si = new SimulationInstance(body,sim,pbm);
        // Determine ICs
        findICs(args);
        // Fill initial values
        si->fillSteadys(T0,beta0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        // Run si
        si->runSim();
        // Fill outputs
        outs[0] = si->time;
    }


    SimulationInstance::SimulationInstance(BODYMODEL::BodyModel* tbody,SIMMODEL::SimModel* tsim, PSEUDOBLOCKMATRIX::PseudoBlockMatrix* tpbm)
    {
        // Copy over ptrs
        body = tbody;
        sim = tsim;
        pbm = tpbm;

        
    }
    void SimulationInstance::copyToInitials( double *T0tmp, double *beta0tmp,double *q0tmp, double *Tpp0tmp){
        for(idx=0; idx<body->N; idx++){
            T0tmp[idx] = T0[idx];
            beta0tmp[idx] = beta0[idx];
            q0tmp[idx] = q0[idx];
            Tpp0tmp[idx] = Tpp0[idx];
        }
    }
    void SimulationInstance::fillInitials( double *T0tmp, double *beta0tmp, double *q0tmp, double *Tpp0tmp){
        for(idx=0; idx<body->N; idx++){
            T0[idx] = T0tmp[idx];
            beta0[idx] = beta0tmp[idx];
            q0[idx] = q0tmp[idx];
            Tpp0[idx] = Tpp0tmp[idx];
        }
    }
    void SimulationInstance::copyToSteadys( double *T0tmp, double *beta0tmp, double *q0tmp, double *Tpp0tmp,
            double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
            double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp){
        for(idx=0; idx<body->N; idx++){
            T0tmp[idx] = T0[idx];
            beta0tmp[idx] = beta0[idx];
            q0tmp[idx] = q0[idx];
            Tpp0tmp[idx] = Tpp0[idx];
        }
        M0tmp = M0;
        QResp0tmp = QResp0;
        Tskm0tmp = Tskm0;
        H0tmp = H0;
        Sh0tmp = Sh0;
        Cs0tmp = Cs0;
        Dl0tmp = Dl0;
        Sw0tmp = Sw0;
    }
    void SimulationInstance::fillSteadys( double *T0tmp, double *beta0tmp, double *q0tmp, double *Tpp0tmp,
            double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
            double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp){
        for(idx=0; idx<body->N; idx++){
            T0[idx] = T0tmp[idx];
            beta0[idx] = beta0tmp[idx];
            q0[idx] = q0tmp[idx];
            Tpp0[idx] = Tpp0tmp[idx];
        }
        M0 = M0tmp;
        QResp0 = QResp0tmp;
        Tskm0 = Tskm0tmp;
        H0 = H0tmp;
        Sh0 = Sh0tmp;
        Cs0 = Cs0tmp;
        Dl0 = Dl0tmp;
        Sw0 = Sw0tmp;
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
            Mbas0 += body->V[idx]*cur->q_m;
            MbasDelta += body->V[idx]*sim->thermovariant*deltaQMetabolic(cur->q_m,T[idx],T0[idx]);
            idx++;
            for(sectIdx = 0; sectIdx<elem->nSectors; ++sectIdx){
                SECTOR::Sector* sect = elem->sectors[sectIdx];
                for(washIdx = 1; washIdx<elem->nWashers; ++washIdx){
                    cur = elem->washers[washIdx];
                    Mbas0 += body->V[idx]*cur->q_m;
                    MbasDelta += body->V[idx]*sim->thermovariant*deltaQMetabolic(cur->q_m,T[idx],T0[idx]);
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
        return 3.45*M*(0.028+6.5e-5*sim->Tair-4.98e-6*sim->Pair)+
            1.44e-3*M*(32.6-0.934*sim->Tair+1.99e-4*sim->Pair);
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
        TskmNxt = project(Tskm,TskmPrv);
        ShNxt = project(Sh,ShPrv);
        CsNxt = project(Cs,CsPrv);
        DlNxt = project(Dl,DlPrv);
        SwNxt = project(Sw,SwPrv);
        
    }
    void SimulationInstance::projectNodeValues(){
        for(idx=0;idx<body->N;idx++){
            // Project nodal temperature only to compute q,w,beta
            TNxt[idx] = project(T[idx],TPrv[idx]); // K, Nodal temperature.
        }
        
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
    void SimulationInstance::cp( int eleIdx ){
        Cp[eleIdx] = (sim->salineRho*sim->salineCp*(DViv/(Viv+DViv)*FlowECMOBlood[eleIdx]+FlowECMOSaline[eleIdx])
            +body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)*FlowECMOBlood[eleIdx]))/
            (sim->salineRho*(DViv/(Viv+DViv)*FlowECMOBlood[eleIdx]+FlowECMOSaline[eleIdx])+
            body->rhoBlood*(Viv/(Viv+DViv)*FlowECMOBlood[eleIdx]));
        CpNxt[eleIdx] = (sim->salineRho*sim->salineCp*(DVivNxt/(VivNxt+DVivNxt)*FlowECMOBloodNxt[eleIdx]+FlowECMOSalineNxt[eleIdx])
            +body->rhoBlood*body->cpBlood*(VivNxt/(VivNxt+DVivNxt)*FlowECMOBloodNxt[eleIdx]))/
            (sim->salineRho*(DVivNxt/(VivNxt+DVivNxt)*FlowECMOBloodNxt[eleIdx]+FlowECMOSalineNxt[eleIdx])+
            body->rhoBlood*(VivNxt/(VivNxt+DVivNxt)*FlowECMOBloodNxt[eleIdx]));
    }
    void SimulationInstance::rho( int eleIdx ){
        Rho[eleIdx] = (sim->salineRho*(DViv/(Viv+DViv)*FlowECMOBlood[eleIdx]+FlowECMOSaline[eleIdx])
            +body->rhoBlood*(Viv/(Viv+DViv)*FlowECMOBlood[eleIdx]))/
            (FlowECMOBlood[eleIdx]+FlowECMOSaline[eleIdx]);
        RhoNxt[eleIdx] = (sim->salineRho*(DVivNxt/(VivNxt+DVivNxt)*FlowECMOBloodNxt[eleIdx]+FlowECMOSalineNxt[eleIdx])
            +body->rhoBlood*(VivNxt/(VivNxt+DVivNxt)*FlowECMOBloodNxt[eleIdx]))/
            (FlowECMOBloodNxt[eleIdx]+FlowECMOSalineNxt[eleIdx]);
    }

    void SimulationInstance::tblp( int eleIdx ){
        TblP[eleIdx] = T[body->N-1];
        /*
        TblP[eleIdx] = (sim->salineRho*sim->salineCp*(DViv/(Viv+DViv)*FlowECMOBlood[eleIdx]*T[body->N-1]+FlowECMOSaline[eleIdx]*sim->TECMO)+
            body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)*FlowECMOBlood[eleIdx]*T[body->N-1]))/
            (sim->salineRho*sim->salineCp*(DViv/(Viv+DViv)*FlowECMOBlood[eleIdx]+FlowECMOSaline[eleIdx])+
            body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)*FlowECMOBlood[eleIdx]));
        */
        TblPNxtRatio[eleIdx] = TblP[eleIdx]/T[body->N-1];
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
            BVNxt[elemIdx] = 0; // -, beta*volume over all nodes
            element = body->elements[elemIdx];

            // Core washer
            washer = element->washers[0];
            BV[elemIdx] += beta[idx]*body->V[idx];
            BVNxt[elemIdx] += betaNxt[idx]*body->V[idx];
            BVT[elemIdx] += BV[elemIdx]*T[idx];
            idx++;
            
            // Interior/skin washers
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
                for(washIdx = 1; washIdx<element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];
                    
                    BV[elemIdx] += beta[idx]*body->V[idx];
                    BVNxt[elemIdx] += betaNxt[idx]*body->V[idx];
                    BVT[elemIdx] += BV[elemIdx]*T[idx];
                    idx++;
                }
            }
            
            // Derived values
            //  See Westin eqn 21
            BPRBPCfactor[elemIdx] = BV[elemIdx]/(element->hx+BV[elemIdx]); // -
            BPRBPCfactorNxt[elemIdx] = BVNxt[elemIdx]/(element->hx+BVNxt[elemIdx])*TblPNxtRatio[elemIdx]; // -

            TblAoverlayFactor[elemIdx] = element->hx/BV[elemIdx] / (element->hx + BV[elemIdx]); // -
            TblAoverlayFactorNxt[elemIdx] = element->hx/BVNxt[elemIdx] / (element->hx + BVNxt[elemIdx]); // -
            TblAoverlay[elemIdx] = BVT[elemIdx] * TblAoverlayFactor[elemIdx]; // K
            
            TblA[elemIdx] = TblP[elemIdx]*BPRBPCfactor[elemIdx] + TblAoverlay[elemIdx]; // K
        }
    }

    void SimulationInstance::agglomeratedBodyValues(){
        CplC = 0;
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];
            CplC += -(BVNxt[elemIdx]*BVNxt[elemIdx])/(element->hx+BVNxt[elemIdx])*TblPNxtRatio[elemIdx];
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
                qResp = -Qresp*washer->a_resp/(body->V[idx])*sim->thermovariant;
            }
            qs[idx] = washer->q_m+qDm+qW+qSh+qResp;
            ws[idx] = max(w0[idx]*pow(sim->KonstasAlpha,sim->KonstasBeta*(Ts[idx]-T0[idx])*sim->thermovariant)*(1-sim->KonstasGamma*DeltaHCT),0.0);
            betas[idx] = (beta0[idx]+0.932*(qDm+qSh+qW))
                    *ws[idx]/w0[idx]
                    *rhos[elemIdx]/body->rhoBlood
                    *cps[elemIdx]/body->cpBlood;
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
                        qResp = -Qresp*washer->a_resp/(body->V[idx])*sim->thermovariant;
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
    void SimulationInstance::BVRSVR(){
        sim->BVRSVR(&bvr, &svr, time);
        sim->BVRSVR(&bvrNxt,&svrNxt,timeNxt);
    }
    void SimulationInstance::bloodParams(){
        Viv = body->Viv0*bvr;
        VivNxt = body->Viv0*bvrNxt;

        Vrbc = body->Vrbc0*bvr;
        VrbcNxt = body->Vrbc0*bvrNxt;

        DViv = body->Viv0*svr;
        DVivNxt = body->Viv0*svrNxt;
    }
    void SimulationInstance::deltaHCT(){
        DeltaHCT = (body->p*Vrbc*DViv)/(Viv*(Viv+body->p*DViv));
    }
    void SimulationInstance::flowECMOSaline(int elemIdx){
        sim->ecmoSaline(&(FlowECMOSaline[elemIdx]), elemIdx, time);
        sim->ecmoSaline(&(FlowECMOSalineNxt[elemIdx]), elemIdx, timeNxt);
    }
    void SimulationInstance::flowECMOBlood(int elemIdx){
        sim->ecmoBlood(&(FlowECMOBlood[elemIdx]), elemIdx, time);
        sim->ecmoBlood(&(FlowECMOBloodNxt[elemIdx]), elemIdx, timeNxt);
    }
    void SimulationInstance::ECMOtreatment(){
        for(elemIdx=0;elemIdx<body->nElements;elemIdx++){
            flowECMOBlood(elemIdx);
            flowECMOSaline(elemIdx);
        }
    }
    void SimulationInstance::BCvalues(){
        // Shock
        BVRSVR();
        // ECMO treatment
        ECMOtreatment();
        // Hemodilution
        bloodParams();
        // Change in hematocrit
        deltaHCT();
    }
    double SimulationInstance::computeMeanSkinTemp(){
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
        TskErrorGradient = (Tskm-TskmPrv)/sim->dt/3600.0*sim->transient*sim->thermovariant;
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
            
            // Core node eqn (Westin eqn 57)
            coreIdx = idx;
            // Current (Core) node
            washer = element->washers[0];
            coreWasher = washer;
            // Westin eqn 57 rhs term 1
            rhs[idx] += (
                washer->zeta/sim->dt
                -washer->del*beta[idx]
                +element->theta*(washer->AForwardCur-1)*element->sumPhi
            ) * T[idx];
            // Westin eqn 57 lhs term 1
            pbm->add(idx,idx,
                washer->zeta/sim->dt
                +washer->del*betaNxt[idx]
                +element->theta*(1-washer->AForwardCur)*element->sumPhi
            );
            // Core node heat gen
            // Westin eqn 57 rhs term 3
            rhs[idx] += washer->del*(q[idx] + qNxt[idx]);
            
            // Core node blood perfusion warming
            // Westin eqn 57 rhs term 4
            rhs[idx] += washer->del*beta[idx]*TblA[idx];
            // Westin eqn 57 lhs term 3
            //      Westin eqns 21 term 2 (Also called TblA overlay)
            for(int tIdx=coreIdx; tIdx<coreIdx+element->N; tIdx++){
                pbm->add(idx,tIdx,
                    -washer->del
                    *betaNxt[idx]
                    *betaNxt[tIdx]
                    *body->V[tIdx]
                    *TblAoverlayFactorNxt[elemIdx]
                );
            }
            // Westin eqn 57 lhs term 3
            //      Westin eqns 69
            pbm->add(idx,body->N-1,
                -washer->del*betaNxt[idx]*BPRBPCfactorNxt[elemIdx]
            );
            //      Westin eqns 70
            pbm->add(body->N-1,idx,
                betaNxt[idx]*body->V[idx]*BPRBPCfactorNxt[elemIdx]
            );

            idx++;
            // Loop thru sectors
            for(sectIdx = 0; sectIdx < element->nSectors; ++elemIdx){
                sector = element->sectors[sectIdx];
                // Loop thru washers
                for(washIdx = 1; washIdx < element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];


                    // Backwards nodes
                    backwardIdx = -1;
                    // Core adjacent nodes have the same previous node, the core
                    if(washIdx == 1){
                        backwardIdx = coreIdx;

                        // Add this node to core node eqn
                        // Westin eqn 57 lhs term 2
                        pbm->add(coreIdx,idx,-element->theta*coreWasher->AForwardNxt*sector->phi);
                        // Westin eqn 57 rhs term 2
                        rhs[coreIdx] += element->theta*coreWasher->AForwardNxt*sector->phi*T[idx];
                    }
                    // Other nodes have an interior node as a previous node
                    else{
                        backwardIdx = idx-1;
                    }
                    // Westin eqn 55 rhs term 1
                    rhs[idx] += (1-washer->gamma)*washer->ABackwardPrv*T[backwardIdx];
                    // Westin eqn 55 lhs term 1
                    pbm->add(idx,backwardIdx,(washer->gamma-1)*washer->ABackwardPrv);

                    // Forwards nodes
                    forwardIdx = -1;
                    // Internal nodes have forward nodes
                    if(washIdx < element->nWashers-1){
                        forwardIdx = idx+1;
                        // Westin eqn 55 rhs term 3
                        rhs[idx] += (1+washer->gamma)*washer->AForwardNxt*T[forwardIdx];
                        // Westin eqn 55 lhs term 3
                        pbm->add(idx,forwardIdx,-(washer->gamma+1)*washer->AForwardNxt);
                    }
                    // Skin nodes use the virtual temperature Tpp instead
                    else{
                        // Westin eqn 68 rhs term 3
                        rhs[idx] += (1+washer->gamma)*washer->AForwardNxt*(Tpp[idx]+TppNxt[idx]);
                    }

                    // Current node coupling
                    // Westin eqn 55 rhs term 2 / Westin eqn 68 rhs term 2 (The same thing)
                    rhs[idx] += (
                        (1-washer->gamma)*washer->ABackwardCur
                        + washer->zeta/sim->dt
                        -2
                        -washer->del*beta[idx]
                        +(1+washer->gamma)*washer->AForwardCur
                    ) * T[idx];
                    // Westin eqn 55 lhs term 2 / Westin eqn 68 lhs term 2 (The same thing)
                    pbm->add(idx,idx,
                        (washer->gamma-1)*washer->ABackwardCur
                        +washer->zeta/sim->dt
                        +2
                        +washer->del*betaNxt[idx]
                        -(1+washer->gamma)*washer->AForwardCur
                    );

                    // Current node heat gen
                    // Westin eqn 55 rhs term 4 / Westin eqn 68 rhs term 4 (The same thing)
                    rhs[idx] += washer->del*(q[idx] + qNxt[idx]);

                    // Current node blood perfusion warming
                    // Westin eqn 55 rhs term 5 / Westin eqn 68 rhs term 5 (The same thing)
                    rhs[idx] += washer->del*beta[idx]*TblA[idx];
                    // Westin eqn 55 lhs term 4 / Westin eqn 68 lhs term 3 (The same thing)
                    //      Westin eqns 21 term 2 (Also called TblA overlay)
                    for(int tIdx=coreIdx; tIdx<coreIdx+element->N; tIdx++){
                        pbm->add(idx,tIdx,
                            -washer->del
                            *betaNxt[idx]
                            *betaNxt[tIdx]
                            *body->V[tIdx]
                            *TblAoverlayFactorNxt[elemIdx]
                        );
                    }
                    //      Westin eqns 69
                    pbm->add(idx,body->N-1,
                        -washer->del*betaNxt[idx]*BPRBPCfactorNxt[elemIdx]
                    );
                    //      Westin eqns 70
                    pbm->add(body->N-1,idx,
                        betaNxt[idx]*body->V[idx]*BPRBPCfactorNxt[elemIdx]
                    );
                }
            }
        }

        // Westin eqn 71, CplC
        pbm->add(body->N-1,body->N-1,
            CplC
        );
    }

    void SimulationInstance::solveSystem(){
        pbm->GaussSeidel(rhs,T,1e-5,1e-5,TNxt);
    }
    void SimulationInstance::permuteTimestep(){
        // Initial values
        double *T0,*beta0,*q0,*Tpp0;
        // Steady values
        double M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0;

        // Overwrite Prv
        for(idx=0;idx<body->N;idx++){
            TPrv[idx] = T[idx];
        }
        MPrv = M;
        QrespPrv = Qresp;
        TskmPrv = Tskm;
        HPrv = H;
        ShPrv = Sh;
        CsPrv = Cs;
        DlPrv = Dl;
        SwPrv = Sw;

        // Overwrite Cur
        for(idx=0;idx<body->N;idx++){
            T[idx] = TNxt[idx];
        }
    }

    void SimulationInstance::runSim(){
        
        // Check status of body and simmodel. Make sure both completed.
        assert(body->getState() == BODYMODEL::computed);
        assert(sim->getState() == SIMMODEL::allValuesAssigned);

        // Provide initial values
        time = 0;

        // Iterate
        for(int timestep=1;timestep<sim->nSteps;++timestep){
            
            // Compute temporary BC values and properties from sim
            timeNxt = time + sim->dt;
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

            // Compute node thermal load parameters
            projectNodeValues();
            nodeValues();

            // Compute agglomerated element values
            agglomeratedElementValues();

            // Compute agglomerated body values
            agglomeratedBodyValues();

            // Empty PBM and RHS
            clearSystem();

            // Build system
            buildSystem();

            // Solve system
            solveSystem();

            // Overwrite values with "Prv"
            permuteTimestep();

            // End condition
            if(sim->endCondition(Thy))
                break;

        }

    }
}
#endif