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
        initializer();
        pbm = new PSEUDOBLOCKMATRIX::PseudoBlockMatrix(body->nElements,blockN,blockM);
    }
    Simulator::~Simulator(){
        delete pbm;
        delete [] T0;
        delete [] beta0;
        delete [] q0;
        delete [] Tpp0;
        delete [] w0;
        cout << "Destroyed?" << endl;
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
        T0 = new double[body->N]();
        beta0 = new double[body->N]();
        w0 = new double[body->N]();
        q0 = new double[body->N]();
        Tpp0 = new double[body->N]();
    }
    void Simulator::findICs( double args[]){
        for(int i=0;i<body->N;i++){
            T0[i] = 34.8+273.15; // K
            Tpp0[i] = NAN;// 30+273.15; // K
            q0[i] = NAN; // W/m^3
            w0[i] = NAN; // m^3/m^3
            beta0[i] = NAN; // W/kgK
        }
        Sh0 = NAN; // W
        H0 = NAN; // W
        M0 = NAN;
        QResp0 = NAN;
        Tskm0 = NAN;
        Cs0 = NAN;
        Dl0 = NAN;
        Sw0 = NAN;

        Tskm0 = 23+273.15; // K
        // Populate beta, q, and prvs for steady case
        SIMMODEL::InitialCase* simInit = new SIMMODEL::InitialCase();
        simInit->setUQs(args);
        SimulationInstance* siInit = new SimulationInstance(body,simInit,pbm);
        siInit->fillSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        siInit->runSim();
        siInit->copyToSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        
        for(int i=0;i<body->N-1;i++)
            cout << T0[i] <<" " << siInit->Tpp[i] << " " << i << " "<< pbm->blockPicker[i] << endl;
        
        // Now find the steady case
        SIMMODEL::SteadyCase* simSteady = new SIMMODEL::SteadyCase();
        simSteady->setUQs(args);
        SimulationInstance* siSteady = new SimulationInstance(body,simSteady,pbm);
        siSteady->fillSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        siSteady->runSim();
        siSteady->copyToSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);

        
        
        
        // // Now normalize by using the transient, no injury case with active controls
        // SIMMODEL::TransientCase* simTransient = new SIMMODEL::TransientCase();
        // simTransient->setUQs(args);
        // SimulationInstance* siTransient = new SimulationInstance(body,simTransient,pbm);
        // siTransient->fillSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        // siTransient->runSim();
        // //siTransient->copyToSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        
        

        // // Free vars
        delete simInit;
        delete simSteady;
        // delete simTransient;
        
        delete siInit;
        delete siSteady;
        // delete siTransient;
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
        si->fillSteadys(T0,beta0,w0,q0,Tpp0,M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0);
        // Run si
        si->runSim();
        // Fill outputs
        
        outs[0] = si->time;
        si->~SimulationInstance();
    }


    SimulationInstance::SimulationInstance(BODYMODEL::BodyModel* tbody,SIMMODEL::SimModel* tsim, PSEUDOBLOCKMATRIX::PseudoBlockMatrix* tpbm)
    {
        // Copy over ptrs
        body = tbody;
        sim = tsim;
        pbm = tpbm;
        
        allocate();
    }
    SimulationInstance::~SimulationInstance(){
        deallocate();
    }

    void SimulationInstance::allocate(){

        // Allocate system rhs
        rhs = new double[body->N]();

        // Allocate element arrays
        // Current
        FlowECMOBlood = new double[body->nElements];
        FlowECMOSaline = new double[body->nElements];
        Cp = new double[body->nElements];
        Rho = new double[body->nElements];
        TblP = new double[body->nElements];
        TblPNxtRatio = new double[body->nElements];
        // Next
        FlowECMOBloodNxt = new double[body->nElements];
        FlowECMOSalineNxt = new double[body->nElements];
        CpNxt = new double[body->nElements];
        RhoNxt = new double[body->nElements];

        // Allocate node arrays
        // Initial
        T0 = new double[body->N];
        beta0 = new double[body->N];
        w0 = new double[body->N];
        q0 = new double[body->N];
        Tpp0 = new double[body->N];
        // Previous
        TPrv = new double[body->N];
        TppPrv = new double[body->N];
        // Current
        T = new double[body->N];
        beta = new double[body->N];
        w = new double[body->N];
        q = new double[body->N];
        Tpp = new double[body->N];
        // Next
        TNxt = new double[body->N];
        betaNxt = new double[body->N];
        wNxt = new double[body->N];
        qNxt = new double[body->N];
        TppNxt = new double[body->N];

        // Allocate agglomerated element arrays
        // Current
        BV = new double[body->nElements];
        BVT = new double[body->nElements];
        BPRBPCfactor = new double[body->nElements];
        TblAoverlayFactor = new double[body->nElements];
        TblAoverlay = new double[body->nElements];
        TblA = new double[body->nElements];
        // Next
        BVNxt = new double[body->nElements];
        BPRBPCfactorNxt = new double[body->nElements];
        TblAoverlayFactorNxt = new double[body->nElements];
        TblAoverlayNxt = new double[body->nElements];
    }
    void SimulationInstance::deallocate(){
        // Allocate system rhs
        delete [] rhs;

        // Allocate element arrays
        // Current
        delete [] FlowECMOBlood;
        delete [] FlowECMOSaline;
        delete [] Cp;
        delete [] Rho;
        delete [] TblP;
        delete [] TblPNxtRatio;
        // Next
        delete [] FlowECMOBloodNxt;
        delete [] FlowECMOSalineNxt;
        delete [] CpNxt;
        delete [] RhoNxt;

        // Allocate node arrays
        // Initial
        delete [] T0;
        delete [] beta0;
        delete [] w0;
        delete [] q0;
        delete [] Tpp0;
        // Previous
        delete [] TPrv;
        delete [] TppPrv;
        // Current
        delete [] T;
        delete [] beta;
        delete [] w;
        delete [] q;
        delete [] Tpp;
        // Next
        delete [] TNxt;
        delete [] betaNxt;
        delete [] wNxt;
        delete [] qNxt;
        delete [] TppNxt;

        // Allocate agglomerated element arrays
        // Current
        delete [] BV;
        delete [] BVT;
        delete [] BPRBPCfactor;
        delete [] TblAoverlayFactor;
        delete [] TblAoverlay;
        delete [] TblA;
        // Next
        delete [] BVNxt;
        delete [] BPRBPCfactorNxt;
        delete [] TblAoverlayFactorNxt;
        delete [] TblAoverlayNxt;
    }
    void SimulationInstance::ghostStart(){
        for(idx=0;idx<body->N;idx++){
            // Previous
            TPrv[idx] = T0[idx];
            TppPrv[idx] = Tpp0[idx];
            // Current
            T[idx] = T0[idx];
            beta[idx] = beta0[idx];
            w[idx] = w0[idx];
            q[idx] = q0[idx];
            Tpp[idx] = Tpp0[idx];
            // Nxt
            TppNxt[idx] = Tpp[idx];
        }
        // Thermal loads
        TskmPrv = Tskm =Tskm0; // K, Mean skin temperature
        Thy0 = Thy = T[0];
        MPrv = M = M0; // W, Total metabolism
        HPrv=H=H0; // W, Heat load
        QRespPrv=QResp=QResp0; // W, Respiratory cooling intake
        // Active controls:
        ShPrv=Sh=Sh0; // W, Shivering power
        CsPrv=Cs=Cs0; // -, Vasoconstriction ratio
        DlPrv=Dl=Dl0; // W/K, Vasodilation capacitance
        SwPrv=Sw=Sw0; // g/min, Sweat output
    }
    void SimulationInstance::copyToInitials( double *T0tmp, double *beta0tmp,double *w0tmp,double *q0tmp, double *Tpp0tmp){
        for(idx=0; idx<body->N; idx++){
            T0tmp[idx] = T[idx];
            beta0tmp[idx] = beta[idx];
            w0tmp[idx] = w[idx];
            q0tmp[idx] = q[idx];
            Tpp0tmp[idx] = Tpp[idx];
        }
    }
    void SimulationInstance::fillInitials( double *T0tmp, double *beta0tmp, double *w0tmp, double *q0tmp, double *Tpp0tmp, double Tskm0tmp){
        for(idx=0; idx<body->N; idx++){
            T0[idx] = T0tmp[idx];
            beta0[idx] = beta0tmp[idx];
            w0[idx] = w0tmp[idx];
            q0[idx] = q0tmp[idx];
            Tpp0[idx] = Tpp0tmp[idx];
        }
        Tskm0 = Tskm0tmp;
    }
    void SimulationInstance::copyToSteadys( double *T0tmp, double *beta0tmp, double *w0tmp, double *q0tmp, double *Tpp0tmp,
            double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
            double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp){
        for(idx=0; idx<body->N; idx++){
            T0tmp[idx] = T[idx];
            beta0tmp[idx] = beta[idx];
            w0tmp[idx] = w[idx];
            q0tmp[idx] = q[idx];
            Tpp0tmp[idx] = Tpp[idx];
        }
        M0tmp = M;
        QResp0tmp = QResp;
        Tskm0tmp = Tskm;
        H0tmp = H;
        Sh0tmp = Sh;
        Cs0tmp = Cs;
        Dl0tmp = Dl;
        Sw0tmp = Sw;
    }
    void SimulationInstance::fillSteadys( double *T0tmp, double *beta0tmp,double *w0tmp, double *q0tmp, double *Tpp0tmp,
            double M0tmp, double QResp0tmp, double Tskm0tmp, double H0tmp,
            double Sh0tmp, double Cs0tmp, double Dl0tmp, double Sw0tmp){
        for(idx=0; idx<body->N; idx++){
            T0[idx] = T0tmp[idx];
            beta0[idx] = beta0tmp[idx];
            w0[idx] = w0tmp[idx];
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
            MbasDelta += body->V[idx]*sim->thermovariant*sim->transient*deltaQMetabolic(cur->q_m,T[idx],T0[idx]);
            idx++;
            for(sectIdx = 0; sectIdx<elem->nSectors; ++sectIdx){
                SECTOR::Sector* sect = elem->sectors[sectIdx];
                for(washIdx = 1; washIdx<elem->nWashers; ++washIdx){
                    cur = elem->washers[washIdx];
                    Mbas0 += body->V[idx]*cur->q_m;
                    MbasDelta += body->V[idx]*sim->thermovariant*sim->transient*deltaQMetabolic(cur->q_m,T[idx],T0[idx]);
                    idx++;
                }
            }
        }
        assert(!isnan(abs(Sh)));
        assert(!isnan(abs(MbasDelta)));
        Mbas0 = 58.2;
        double act = 1 + (Sh+MbasDelta)/Mbas0*MperWork; // MET
        double etaW;
        if(act < 1.6)
            etaW = 0.05; // W/W
        else
            etaW = 0.2*tanh(body->b1*act+body->b0); // W/W
        assert(!isnan(abs(act)));
        M = act*Mbas0/body->actBas; // W
        assert(!isnan(abs(M)));
        H = M*(1-etaW)-Mbas0; // W
        assert(!isnan(abs(H)));
        QResp = computeQresp(); // W
        assert(!isnan(abs(QResp)));
    }
    double SimulationInstance::computeQresp(){
        assert(!isnan(abs(sim->Pair)));
        assert(!isnan(abs(sim->Tair)));
        return 3.45*M*(0.028+6.5e-5*sim->Tair-4.98e-6*sim->Pair)+
            1.44e-3*M*(32.6-0.934*sim->Tair+1.99e-4*sim->Pair);
    }
    double SimulationInstance::deltaQMetabolic( double q, double T, double Tn)
    {
        return q*(pow(sim->Q10,(T-Tn)*sim->thermovariant*sim->transient/10.0)-1);
    }
    double SimulationInstance::project( double cur, double prv ){
        assert(!isnan(abs(cur)));
        if(sim->transient ==0 || sim->thermovariant==0){
            return cur;
        }
        else{
            assert(!isnan(abs(prv)));
            return cur+(cur-prv);
        }
        
    }
    void SimulationInstance::projectBodyValues(){
        // MNxt = project(M,MPrv);
        // HNxt = project(H,HPrv);
        // QRespNxt = project(QResp,QRespPrv);
        // TskmNxt = project(Tskm,TskmPrv);
        // ShNxt = project(Sh,ShPrv);
        // CsNxt = project(Cs,CsPrv);
        // DlNxt = project(Dl,DlPrv);
        // SwNxt = project(Sw,SwPrv);
    }
    void SimulationInstance::projectTemperatures(){
        for(idx=0;idx<body->N;idx++){
            // Project nodal temperature only to compute q,w,beta
            TNxt[idx] = project(T[idx],TPrv[idx]); // K, Nodal temperature.
        }
        
    }
    void SimulationInstance::projectTpp(){
        for(idx=0;idx<body->N;idx++){
            // Project nodal temperature for qsk
            TppNxt[idx] = project(Tpp[idx],TppPrv[idx]);
        }
        
    }

    void SimulationInstance::elementValues(){
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];
            rhocp(elemIdx);
            sim->tblp(TblP,TblPNxtRatio,body,elemIdx,T[body->N-1],Viv,DViv,FlowECMOBlood[elemIdx],FlowECMOSaline[elemIdx],time);
        }
    }
    void SimulationInstance::rhocp( int eleIdx ){
        assert(!isnan(abs(DViv)));
        assert(!isnan(abs(Viv)));
        assert(!isnan(abs(FlowECMOBlood[eleIdx])));
        assert(!isnan(abs(FlowECMOSaline[eleIdx])));
        sim->rhocp(Rho,Cp,body,eleIdx,Viv,DViv,FlowECMOBlood[eleIdx],FlowECMOSaline[idx],time);
        sim->rhocp(RhoNxt,CpNxt,body,eleIdx,VivNxt,DVivNxt,FlowECMOBloodNxt[eleIdx],FlowECMOSalineNxt[idx],time);
        assert(!isnan(abs(Rho[eleIdx])));
        assert(!isnan(abs(RhoNxt[eleIdx])));
    }

    void SimulationInstance::nodeValues(){
        for(idx=0;idx<body->N;idx++){
            q[idx] = 0; // W/m^3, Heat generation in tissue
            w[idx] = 0; // -, Blood perfusion rate
            beta[idx] = 0; // W/m^3/K, Blood perfusion rate factor (rho*c*w)
        }
        
        // Skin values
        sim->skinBC(Tpp,body,T,T0,Sw,time);
        // Heat generation, blood perfusion, heats
        qwbeta(q,w,beta,Cp,Rho,T);

        // Future terms
        if(sim->thermovariant > 0 && sim->transient > 0 && 1){
            qwbeta(qNxt,wNxt,betaNxt,CpNxt,RhoNxt,TNxt);
            sim->skinBC(TppNxt,body,TNxt,T0,Sw,timeNxt);
        }
        else{
            for(idx=0;idx<body->N;idx++){
                qNxt[idx] = q[idx]; // W/m^3, Heat generation in tissue
                wNxt[idx] = w[idx]; // -, Blood perfusion rate
                betaNxt[idx] = beta[idx]; // W/m^3/K, Blood perfusion rate factor (rho*c*w)
                TppNxt[idx] = Tpp[idx];
            }
        }

        
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
            assert(!isnan(abs(body->V[idx])));
            assert(!isnan(abs(beta[idx])));
            assert(!isnan(abs(betaNxt[idx])));
            idx++;
            
            // Interior/skin washers
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
                for(washIdx = 1; washIdx<element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];
                    
                    BV[elemIdx] += beta[idx]*body->V[idx];
                    BVNxt[elemIdx] += betaNxt[idx]*body->V[idx];
                    BVT[elemIdx] += beta[idx]*body->V[idx]*T[idx];
                    assert(!isnan(abs(body->V[idx])));
                    assert(!isnan(abs(beta[idx])));
                    assert(!isnan(abs(betaNxt[idx])));
                    idx++;
                }
            }
            // Derived values
            //  See Westin eqn 21
            BPRBPCfactor[elemIdx] = BV[elemIdx]/(element->hx+BV[elemIdx]); // -
            BPRBPCfactorNxt[elemIdx] = BVNxt[elemIdx]/(element->hx+BVNxt[elemIdx]); // -

            TblAoverlayFactor[elemIdx] = element->hx/BV[elemIdx] / (element->hx + BV[elemIdx]); // -
            TblAoverlayFactorNxt[elemIdx] = element->hx/BVNxt[elemIdx] / (element->hx + BVNxt[elemIdx]); // -
            TblAoverlay[elemIdx] = BVT[elemIdx] * TblAoverlayFactor[elemIdx]; // K
            
            TblA[elemIdx] = TblP[elemIdx]*BPRBPCfactor[elemIdx] + TblAoverlay[elemIdx]; // K
            assert(!isnan(abs(TblAoverlayFactor[elemIdx])));
            assert(!isnan(abs(TblAoverlayFactorNxt[elemIdx])));
            assert(!isnan(abs(TblAoverlay[elemIdx])));
        }
    }

    void SimulationInstance::agglomeratedBodyValues(){
        CplC = 0;
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];
            CplC += -(BVNxt[elemIdx]*BVNxt[elemIdx])/(element->hx+BVNxt[elemIdx])/TblPNxtRatio[elemIdx];
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
            qDm = deltaQMetabolic(washer->q_m,Ts[idx],T0[idx])*sim->thermovariant*sim->transient;
            qW = 0;
            qSh = 0;
            qResp = 0;
            if(element->Vmuscle > 0 && washer->muscle > 0){
                qW = element->a_sed*H/element->Vmuscle;
                qSh = element->a_sh*Sh/element->Vmuscle*sim->thermovariant*sim->transient;
            }
            if(element->Vresp > 0){
                qResp = QResp*washer->a_resp/(body->V[idx]);
            }
            qs[idx] = washer->q_m+qDm+qW+qSh+qResp;
            ws[idx] = max(pow(sim->KonstasAlpha,sim->KonstasBeta*(Ts[idx]-T0[idx])*sim->thermovariant*sim->transient)*(1-sim->KonstasGamma*DeltaHCT*sim->thermovariant*sim->transient),0.0);
            betas[idx] = (washer->w_bl*body->rhoBlood*body->cpBlood+0.932*(qDm+qSh+qW))
                    *ws[idx]
                    *rhos[elemIdx]/body->rhoBlood
                    *cps[elemIdx]/body->cpBlood;
            ws[idx] *= washer->w_bl;
            
            assert(!isnan(abs(betas[idx])));
            idx++;
            
            // Loop thru washers
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
                for(washIdx = 1; washIdx < element->nWashers; ++washIdx){
                    washer = element->washers[washIdx];

                    qDm = deltaQMetabolic(washer->q_m,Ts[idx],T0[idx])*sim->thermovariant*sim->transient;
                    
                    qW = 0;
                    qSh = 0;
                    qResp = 0;
                    if(element->Vmuscle > 0 && washer->muscle > 0){
                        qW = element->a_sed*H/element->Vmuscle;
                        qSh = element->a_sh*Sh/element->Vmuscle*sim->thermovariant*sim->transient;
                    }
                    if(element->Vresp > 0){
                        qResp = QResp*washer->a_resp/washer->volume;
                    }

                    assert(!isnan(abs(qDm)));
                    assert(!isnan(abs(qW)));
                    assert(!isnan(abs(qSh)));
                    assert(!isnan(abs(qResp)));
                    qs[idx] = washer->q_m+qDm+qW+qSh+qResp;
                    ws[idx] = max(pow(sim->KonstasAlpha,sim->KonstasBeta*(Ts[idx]-T0[idx])*sim->thermovariant*sim->transient)*(1-sim->KonstasGamma*DeltaHCT),0.0);
                    assert(!isnan(abs(ws[idx])));    
                    assert(!isnan(abs(rhos[elemIdx])));
                    assert(!isnan(abs(cps[elemIdx])));
                    betas[idx] = (washer->w_bl*body->rhoBlood*body->cpBlood+0.932*(qDm+qSh+qW))
                            *ws[idx]
                            *rhos[elemIdx]/body->rhoBlood
                            *cps[elemIdx]/body->cpBlood;
                    ws[idx] *= washer->w_bl;
                    assert(!isnan(abs(betas[idx])));
                    // Add vasodilation, vasoconstriction for outermost skin
                    if(washIdx == element->nWashers-1){
                        // Westin eqn 84
                        betas[idx] = (betas[idx]+element->a_dl*Dl)/(1+element->a_cs*Cs*exp(-Dl/80));
                        // Westin eqn 85
                        assert(!isnan(abs(H)));
                        assert(!isnan(abs(Viv)));
                        assert(!isnan(abs(DViv)));
                        betas[idx] = min((386.9-.32*.932*H)*(Viv+DViv)/body->Viv0,betas[idx]);
                    }
                    // cout<<idx<<" "<<(washer->w_bl*body->rhoBlood*body->cpBlood+0.932*(qDm+qSh+qW))<<" "<<ws[idx]<<" "<<qs[idx]<<endl;
                    assert(!isnan(abs(betas[idx])));
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
        sim->ecmoSaline(&(FlowECMOSaline[elemIdx]), body, elemIdx, Viv, bvr, time);
        sim->ecmoSaline(&(FlowECMOSalineNxt[elemIdx]), body, elemIdx, Viv, bvr, timeNxt);
    }
    void SimulationInstance::flowECMOBlood(int elemIdx){
        sim->ecmoBlood(&(FlowECMOBlood[elemIdx]), body, elemIdx, Viv, time);
        sim->ecmoBlood(&(FlowECMOBloodNxt[elemIdx]), body, elemIdx, VivNxt, timeNxt);
    }
    void SimulationInstance::ECMOtreatment(){
        for(elemIdx=0;elemIdx<body->nElements;elemIdx++){
            flowECMOBlood(elemIdx);
            flowECMOSaline(elemIdx);
        }
    }
    void SimulationInstance::BCvalues(){
        // Surroundings BCs
        sim->setTairTsrm(time);
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
        int idx=0;
        for(elemIdx = 0;elemIdx<body->nElements;elemIdx++){
            element = body->elements[elemIdx];
            for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                sector = element->sectors[sectIdx];
                idx+=element->nWashers-1;
                Tskcur += element->a_sk*sector->phi/element->sumPhi*
                        T[idx];
            }
            idx++;
        }
        assert(!isnan(abs(Tskcur)));
        return Tskcur;
    }
    double SimulationInstance::computeHypothalamicTemp(){
        return T[0];
    }
    void SimulationInstance::computeErrorSignals(){
        // Tsk,m
        Tskm = computeMeanSkinTemp();
        Thy = computeHypothalamicTemp();
        assert(!isnan(abs(Tskm0)));
        assert(!isnan(abs(Thy0)));
        TskError = (Tskm-Tskm0)*sim->thermovariant;
        ThyError = (Thy-Thy0)*sim->thermovariant;
        TskErrorGradient = (Tskm-TskmPrv)/sim->dt/3600.0*sim->thermovariant;
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
        assert(!isnan(abs(Dl)));
        Sw = (.8*tanh(.59*TskError-.19)+1.2)*TskError
            +(5.7*tanh(1.98*ThyError-1.03)+6.3)*ThyError; // g/min
        Sh = min(350.0*(svr+bvr),max(0.0,Sh)); // W
        Sw = min(30.0,max(Sw,0.0)); // g/min
        sim->shiverDrugs(&Sh,time);
        assert(!isnan(abs(Dl)));
        assert(!isnan(abs(Cs)));
        assert(!isnan(abs(Sh)));
        assert(!isnan(abs(Sw)));
    }

    void SimulationInstance::clearSystem(){
        pbm->fill(0.0);
        for(int i=0;i<body->N;i++)
            rhs[i] = 0;
    }

    void safeAdd(double* arr, int idx, double val){
        assert(!isnan(abs(val)));
        arr[idx] += val;
    }

    void SimulationInstance::buildSystem(){
        idx = 0;
        double tval = NAN;
        // Loop thru elements
        for(elemIdx = 0; elemIdx < body->nElements; ++elemIdx){
            element = body->elements[elemIdx];
            
            // Core node eqn (Westin eqn 57)
            coreIdx = idx;
            // Current (Core) node
            washer = element->washers[0];
            coreWasher = washer;
            // Westin eqn 57 rhs term 1
            
            safeAdd(rhs,idx,
                (
                    washer->zeta/sim->dt*sim->thermovariant*sim->transient
                    -washer->del*beta[idx]
                    +element->theta*(washer->AForwardCur-1)*element->sumPhi
                ) * T[idx]*sim->transient
            );
            // Westin eqn 57 lhs term 1
            pbm->add(idx,idx,
                washer->zeta/sim->dt*sim->thermovariant*sim->transient
                +washer->del*betaNxt[idx]
                +element->theta*(1-washer->AForwardCur)*element->sumPhi
            );
            // Core node heat gen
            // Westin eqn 57 rhs term 3
            safeAdd(rhs,idx,
                washer->del*(q[idx] + qNxt[idx]*sim->thermovariant*sim->transient)
            );
            
            // Core node blood perfusion warming
            // Westin eqn 57 rhs term 4
            safeAdd(rhs,idx,
                washer->del*beta[idx]*TblA[elemIdx]*sim->thermovariant*sim->transient
            );
            // Westin eqn 57 lhs term 3
            //      Westin eqns 21 term 2 (Also called TblA overlay)
            for(int tIdx=coreIdx; tIdx<coreIdx+element->N; tIdx++){
                pbm->add(idx,tIdx,
                    -washer->del
                    *betaNxt[idx]
                    *TblAoverlayFactorNxt[elemIdx]
                    *betaNxt[tIdx]
                    *body->V[tIdx]
                );
            }
            // Westin eqn 57 lhs term 3
            //      Westin eqns 69
            pbm->add(idx,body->N-1,
                -washer->del*betaNxt[idx]*BPRBPCfactorNxt[elemIdx]*TblPNxtRatio[elemIdx]
            );
            //      Westin eqns 70
            pbm->add(body->N-1,idx,
                betaNxt[idx]*body->V[idx]*BPRBPCfactorNxt[elemIdx]*TblPNxtRatio[elemIdx]
            );

            idx++;
            // Loop thru sectors
            for(sectIdx = 0; sectIdx < element->nSectors; ++sectIdx){
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
                        safeAdd(rhs,coreIdx,
                            element->theta*coreWasher->AForwardNxt*sector->phi*T[idx]*sim->thermovariant*sim->transient
                        );
                    }
                    // Other nodes have an interior node as a previous node
                    else{
                        backwardIdx = idx-1;
                    }
                    // Westin eqn 55 rhs term 1 / Westin eqn 68 rhs term 1
                    safeAdd(rhs,idx,(1-washer->gamma)*washer->ABackwardPrv*T[backwardIdx]*sim->thermovariant*sim->transient);
                    // Westin eqn 55 lhs term 1 / Westin eqn 68 lhs term 1
                    pbm->add(idx,backwardIdx,(washer->gamma-1)*washer->ABackwardPrv);

                    // Forwards nodes
                    forwardIdx = -1;
                    // Internal nodes have forward nodes
                    if(washIdx < element->nWashers-1){
                        forwardIdx = idx+1;
                        // Westin eqn 55 rhs term 3
                        safeAdd(rhs,idx,(1+washer->gamma)*washer->AForwardNxt*T[forwardIdx]*sim->thermovariant*sim->transient);
                        
                        // Westin eqn 55 lhs term 3
                        pbm->add(idx,forwardIdx,-(washer->gamma+1)*washer->AForwardNxt);
                    }
                    // Skin nodes use the virtual temperature Tpp instead
                    else{
                        // Westin eqn 68 rhs term 3
                        // cout << Tpp[idx] << " " << idx << endl;
                        assert(!isnan(abs(Tpp[idx])));
                        safeAdd(rhs,idx,(1+washer->gamma)*washer->AForwardNxt*(Tpp[idx]+TppNxt[idx]*sim->thermovariant*sim->transient));
                    }

                    // Current node coupling
                    // Westin eqn 55 rhs term 2 / Westin eqn 68 rhs term 2 (The same thing)
                    safeAdd(rhs,idx,(
                        (1-washer->gamma)*washer->ABackwardCur
                        + washer->zeta/sim->dt*sim->thermovariant*sim->transient
                        -2
                        -washer->del*beta[idx]
                        +(1+washer->gamma)*washer->AForwardCur
                    ) * T[idx]*sim->thermovariant*sim->transient);
                    // Westin eqn 55 lhs term 2 / Westin eqn 68 lhs term 2 (The same thing)
                    assert(!isnan(abs(washer->AForwardCur)));
                    pbm->add(idx,idx,
                        (washer->gamma-1)*washer->ABackwardCur
                        + washer->zeta/sim->dt*sim->thermovariant*sim->transient
                        +2
                        +washer->del*betaNxt[idx]
                        -(1+washer->gamma)*washer->AForwardCur
                    );

                    // Current node heat gen
                    // Westin eqn 55 rhs term 4 / Westin eqn 68 rhs term 4 (The same thing)
                    safeAdd(rhs,idx,washer->del*(q[idx] + qNxt[idx]*sim->thermovariant*sim->transient));
                    

                    // Current node blood perfusion warming
                    // Westin eqn 55 rhs term 5 / Westin eqn 68 rhs term 5 (The same thing)
                    safeAdd(rhs,idx,washer->del*beta[idx]*TblA[elemIdx]*sim->thermovariant*sim->transient);
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
                        -washer->del*betaNxt[idx]*BPRBPCfactorNxt[elemIdx]*TblPNxtRatio[elemIdx]
                    );
                    //      Westin eqns 70
                    pbm->add(body->N-1,idx,
                        betaNxt[idx]*body->V[idx]*BPRBPCfactorNxt[elemIdx]*TblPNxtRatio[elemIdx]
                    );

                    idx++;
                }
            }
        }

        // Westin eqn 71, CplC
        pbm->add(body->N-1,body->N-1,
            CplC
        );
        rhs[body->N-1] = 0.0;
    }

    void SimulationInstance::solveSystem(){
        // pbm->directsolve(rhs,TNxt);
        pbm->GaussSeidel(rhs,T,1e-5,1e-5,TNxt);
    }
    void SimulationInstance::permuteTimestep(){
        // Initial values
        // *T0,*beta0,*w0,*q0,*Tpp0;
        // Steady values
        // M0,QResp0,Tskm0,H0,Sh0,Cs0,Dl0,Sw0;

        // Overwrite Prv
        for(idx=0;idx<body->N;idx++){
            TPrv[idx] = T[idx];
            TppPrv[idx] = Tpp[idx];
        }
        MPrv = M;
        QRespPrv = QResp;
        TskmPrv = Tskm;
        HPrv = H;
        ShPrv = Sh;
        CsPrv = Cs;
        DlPrv = Dl;
        SwPrv = Sw;

        // Overwrite Nxt
        for(idx=0;idx<body->N;idx++){
            T[idx] = TNxt[idx];
            //Tpp[idx] = TppNxt[idx];
        }

    }

    void SimulationInstance::runSim(){
        
        // Check status of body and simmodel. Make sure both completed.
        assert(body->getState() == BODYMODEL::computed);
        assert(sim->getState() == SIMMODEL::allValuesAssigned);
        
        // Provide initial values
        time = 0;
        ghostStart();

        // Iterate
        for(int timestep=1;timestep<=sim->nSteps;++timestep){
            // Compute temporary BC values and properties from sim
            timeNxt = time + sim->dt;
            
            // No dependencies
            // Linearly project TNxt for thermal parameters
            projectTemperatures();
            // Compute error inputs
            computeErrorSignals();
            // Find properties of sim and blood
            BCvalues();
            
            
            // Compute active controls
            //      Must run computeErrorSignals
            computeActiveControls();
            
            // Compute whole-body thermal load parameters
            //      Must run computeActiveControls
            computeThermalLoadParameters();
            
            // Compute body values
            projectBodyValues();

            // Compute element values
            //      Must run BCvalues
            elementValues();

            // Compute node thermal load parameters
            nodeValues();

            // Compute agglomerated element values
            agglomeratedElementValues();

            // Compute agglomerated body values
            agglomeratedBodyValues();

            // Empty PBM and RHS
            clearSystem();

            // Build system
            buildSystem();
            // pbm->Print();
            //cout << "rhs " << rhs[210] << endl;
            
            // Solve system
            //pbm->Print();
            // for(int i=25;i<35;++i){
            //     cout << rhs[i]<<endl;
            // }
            solveSystem();
            
            // Overwrite values with "Prv"
            permuteTimestep();
            time += sim->dt;

            // Occasional updates
            if(false){
                cout << "ITER " << time << endl;
                cout << "Thy "<< timestep <<" " << Thy << endl;
                cout << "Tskcur "<< timestep <<" " << Tskm << endl;
                cout << "M "<< timestep <<" " << M << endl;
                cout << "Tskeg" << TskErrorGradient << endl;
                cout << "Tair" << sim->Tair << endl;

                int maxb = 0, maxwidx = 0;
                double maxt = 0;
                idx = 0;
                for(elemIdx = 0;elemIdx<body->nElements;elemIdx++){
                    element = body->elements[elemIdx];
                    if(T[idx]>maxt){
                        maxt=T[idx];
                        maxb=elemIdx;
                        maxwidx=0;
                    }
                    idx++;
                    for(sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                        sector = element->sectors[sectIdx];
                        for(washIdx=1;washIdx<element->nWashers;washIdx++){
                            washer = element->washers[washIdx];
                            if(T[idx]>maxt){
                                maxt=T[idx];
                                maxb=elemIdx;
                                maxwidx=washIdx;
                            }
                            // cout << idx << "\tgamma"<< washer->AForwardNxt << "\tbeta" << beta[idx]<<"\tq" << q[idx]<< "\tqNxt" << qNxt[idx]  << "\tSrDelta" << pbm->SumRow(idx)*T[idx]/rhs[idx] <<"\tT0"<<T0[idx] << "Tpp" << Tpp[idx] <<endl;
                            idx++;
                        }
                        
                    }
                }
                cout << "Max T: " << maxt <<" at " << maxb << " " << maxwidx << endl;
                
            }

            // End condition
            if(sim->endCondition(Thy))
                break;
        }
        // cout << "TblA " << TblA[0] << endl;
        // cout << "TblP " << TblP[0] << endl;
        // for(int i=0;i<body->N;++i){
        //     cout << i << "\tT"<<T[i] << "\trhs" << rhs[i]<<"\tq" << q[i]<< "\tw" << w[i]  << "\tbeta" << beta[i] <<"\tT0"<<T0[i] << "Tpp" << Tpp[i] <<endl;
        // }
        // pbm->Print();
    }
}
#endif