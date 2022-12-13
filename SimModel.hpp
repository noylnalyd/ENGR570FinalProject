#ifndef _SIMMODEL_
#define _SIMMODEL_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>

using namespace std;

namespace SIMMODEL
{
    enum SimModelState { undefined, initialized, allValuesAssigned};


    class SimModel
    {
    private:
        SimModelState _state = undefined;

        // Local pointers
        ELEMENT::Element* element;
        SECTOR::Sector* sector;
        WASHER::Washer* washer;

        
    public:
        // Time attributes
        const double tInitial = 0.0; // s
        double tFinal; // s
        double dt; // s
        int nSteps;

        // Model attributes
        int modelCase; // Number of case (if automating)
        double thermovariant; // 0 for thermoneutral, 1 for variant
        double transient; // 0 for steady, 1 for transient
        double Q10 = 3; // Q10 constant. Often 2, sometimes 3
        const double KonstasBeta = 0.08401; // K^-1
        const double KonstasGamma =  2.245; // -
        const double KonstasAlpha = 2.961; // -
        
        // Saline constants
        const double salineRho = 1004.6; // kg/m^3, saline density
        const double salineCp = 4180; // J/kg/K, saline specific heat capacity

        // Environmental constants
        const double Pair = 2300; // Pa, partial vapor air pressure
        const double vAirEff = 0.05; // m/s, effective air velocity
        const double echelonSrOutdoors = 0.46; // -
        const double echelonSrIndoors = 0.93; // -
        const double alphaSrOutdoors = 0.3; // -
        const double alphaSrIndoors = 0.8; // -
        const double sOutdoors = 10752/126.7; // W/m^2
        const double sIndoors = 500/50; // W/m^2
        const double Tmelt = 273.15; // K, melting point of water

        // UQ attributes worth exploring for all cases
        double TsrmIndoors; // K, mean surroundings temperature inside
        double TsrmOutdoors; // K, mena surroundings temperature outside
        double trecovery; // s, time at which treatment begins
        double TECMO; // K, temperature of reincorporated ECMO blood/saline mix
        // UQ-dependencies
        double TairIndoors; // K, inside / treatment air temperature
        double TairOutdoors; // K, outside / pretreatment air temperature
        // BC values to update at each iteration
        double Tsrm; // K, surrounding temperature mean
        double Tair; // K, air temperature
        double echelonSrm; // -, Surroundings emissivity
        double alphaSrm; // --, Surroundings absorptivity
        double s; // W/m^2, Light intensity

        SimModel(){
            _state = initialized;
        };
        SimModelState getState();
        double skinHeat(
                BODYMODEL::BodyModel* body,
                ELEMENT::Element* elem,
                SECTOR::Sector* sect,
                int elemIdx,
                int sectIdx,
                double Sw, // g/min, sweat command
                double TsfAssumed, // K, assumed skin surface temperature
                double Tsf0, // K, skin thermoneutral temperature
                int njActual, // -, number of equivalent clothings on this sector
                double psiActual // rad^2, exposure angle
        );
        virtual void rhocp(
            double *Rho,
            double *Cp,
            BODYMODEL::BodyModel *body,
            int eleIdx,
            double Viv,
            double DViv,
            double fEB,
            double fES,
            double time
        );
        virtual void skinBC(
                double *Tpp,
                BODYMODEL::BodyModel* body,
                double *T,
                double *T0,
                double Sw,
                double time
        );
        void setUQs( double args[] );
        virtual void setTairTsrm( double time );
        virtual void ecmoBlood(
                double *fEB,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double time
            );
        virtual void ecmoSaline(
                double *fES,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double bvr,
                double time
            );
        virtual void tblp(
                double *tblp,
                double *tblpNxtRatio,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double tblpCur,
                double Viv,
                double DViv,
                double fEB,
                double fES,
                double time
            );
        virtual void shiverDrugs(
                double *Sh,
                double time
            );
        virtual void BVRSVR( double *bvr, double *svr, double time );
        virtual bool endCondition( double Thy );
    };
    
    void SimModel::setUQs(double args[])
    {
        TsrmIndoors = args[0];
        TsrmOutdoors = args[1];
        trecovery = args[2];
        TECMO = args[3];

        // Dependencies!
        TairIndoors = TsrmIndoors;
        TairOutdoors = TsrmOutdoors;

        // Other dependencies
        nSteps = ceil((tFinal-tInitial)/dt);

        _state = allValuesAssigned;
    }
    void SimModel::setTairTsrm(double time)
    {
        Tair = TairOutdoors;
        Tsrm = TsrmOutdoors;
        echelonSrm = echelonSrOutdoors;
        alphaSrm = alphaSrOutdoors;
        s = sOutdoors;
    }
    void SimModel::rhocp(
        double *Rho,
        double *Cp,
        BODYMODEL::BodyModel *body,
        int eleIdx,
        double Viv,
        double DViv,
        double fEB,
        double fES,
        double time
    ){
        Rho[eleIdx] = body->rhoBlood;
        Cp[eleIdx] = body->cpBlood;
    }
    void SimModel::tblp(
            double *tblp,
            double *tblpNxtRatio,
            BODYMODEL::BodyModel *body,
            int eleIdx,
            double tblpCur,
            double Viv,
            double DViv,
            double fEB,
            double fES,
            double time
    ){
        tblp[eleIdx] = tblpCur;
        tblpNxtRatio[eleIdx] = 1.0;
    }
    void SimModel::shiverDrugs(
            double *Sh,
            double time
    ){
        
    };
    void SimModel::skinBC(
            double *Tpp,
            BODYMODEL::BodyModel* body,
            double *T,
            double *T0,
            double Sw,
            double time)
    {
        int idx=0;
        for(int elemIdx=0;elemIdx<body->nElements;elemIdx++){
            element = body->elements[elemIdx];
            washer = element->washers[element->nWashers-1];
            for(int sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                idx+=element->nWashers-1;
                sector = element->sectors[sectIdx];
                
                double qsk = skinHeat(body,element,sector,elemIdx,sectIdx,Sw,(T[idx]),T0[idx],element->n_j,sector->psiStnd);
                assert(!isnan(abs(qsk)));
                // Fiala eqn A6
                if(element->isCylinder){
                    Tpp[idx] = T[idx]-
                            qsk*(washer->r+washer->deltaR/2.0)/washer->k
                            *log((washer->r+washer->deltaR)/washer->r)
                            ;
                    assert(!isnan(abs(Tpp[idx])));
                }
                // Fiala eqn A9
                else{
                    Tpp[idx] = T[idx]-
                            qsk*(
                                (washer->r+washer->deltaR/2.0)*(washer->r+washer->deltaR/2.0)*
                                (1/washer->r-1/(washer->r+washer->deltaR)))/
                                washer->k;
                            ;
                    assert(!isnan(abs(Tpp[idx])));
                }
            }
            idx++;
        }
    }
    void SimModel::BVRSVR( double *bvr, double *svr, double time )
    {
        *bvr = 1.0;
        *svr = 0.0;
    }
    void SimModel::ecmoBlood(
                double *fEB,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double time
            )
    {
        *fEB = 0;
    }
    void SimModel::ecmoSaline(
                double *fES,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double bvr,
                double time
            )
    {
        *fES = 0;
    }
    bool SimModel::endCondition(double Thy)
    {
        return false;
    }
    // Default heat flux of the skin at a select elment and sector at a particular moment in time
    double SimModel::skinHeat(
        BODYMODEL::BodyModel* body,
        ELEMENT::Element* elem,
        SECTOR::Sector* sect,
        int elemIdx,
        int sectIdx,
        double Sw, // g/min, sweat command
        double TsfAssumed, // K, assumed skin surface temperature
        double Tsf0, // K, skin thermoneutral temperature
        int njActual, // -, number of equivalent clothings on this sector
        double psiActual // rad^2, exposure angle
        ){
        
        double qsk = 0, qC = 0, qR = 0, qsR = 0, qe = 0;
        
        double f = 1;
        double echelonSf = body->elements[elemIdx]->sectors[sectIdx]->echelon;
        double echelonSr = echelonSrm;
        if(njActual > 0){
            f = 0.9;
            echelonSr = 1.0;
        }

        double TsfAbs = TsfAssumed; // K
        double TsrmAbs = Tsrm; // K
        double TsfC = TsfAssumed-273.15; // C
        double TsrmC = TsrmAbs - 273.15; // C
        double TairAbs = Tair; // K
        double TairC = Tair-273.15; // C
        
        double hcmix = sqrt(abs(
                elem->a_nat*sqrt(abs(TsfC-TairC))
                +elem->a_frc*vAirEff
                +elem->a_mix
        )); // W/m^2/K
        
        double hR = body->StefanBoltzmann*echelonSf*echelonSrm
                *(TsfAbs*TsfAbs + TsrmAbs*TsrmAbs)
                *(TsfAbs + TsrmAbs); // W/m^2/K
        
        double U = 1/(njActual*body->Icl+1/(f*(hcmix+hR))); // W/m^2/K
        qsR = alphaSrm*psiActual*s;
        

        double Tweight = (hcmix*TairAbs + hR*(TsrmAbs) + qsR)/(hcmix + hR);
        double Posksat = 100.0*exp(18.956-4030.0/(TsfC+235));
        
        double Ue = body->LewisConstant/(
            njActual*body->Icl/body->icl
            +1/(f*hcmix)
        );
        double dsweat = elem->a_sw*sect->areaShare*Sw*pow(2,(TsfAssumed-Tsf0)*transient/Q10) * 1e-3/60;
        double Psk = (body->LambdaH20*dsweat + Posksat*body->Rinverse + Ue*Pair)/
                (Ue + body->Rinverse);
        double ans = U*(TsfAbs-Tweight) + Ue*(Psk-Pair);
        assert(!isnan(abs(ans)));
        return ans;
    }

    SimModelState SimModel::getState()
    {
        return _state;
    }
    class InitialCase : public SimModel
    {
        public:
            InitialCase() : SimModel() {
                thermovariant = 0;
                transient = 0;
                tFinal = 1;
                dt = 1;
            };

    };
    class SteadyCase : public SimModel
    {
        public:
            SteadyCase() : SimModel() {
                thermovariant = 1;
                transient = 0;
                tFinal = 10;
                dt = 1;
            };
    };
    class TransientCase : public SimModel
    {
        public:
            TransientCase() : SimModel() {
                thermovariant = 1;
                transient = 1;
                tFinal = 20.0;
                dt = .1;
            };
    };
    class InjuryCase : public SimModel
    {
        public:
            double tinjury = 300; // s
            double bvrRate = -.4/tinjury; // m^3/m^3 / s
            InjuryCase() : SimModel() {
                thermovariant = 1;
                transient = 1;
                tFinal = 3600.0;
                dt = 0.1;
            };
            void BVRSVR( double *bvr, double *svr, double time ) override
            {
                if(time <= tinjury){
                    *bvr = 1.0 + bvrRate*(time);
                    *svr = 0;
                }else{
                    *bvr = 1.0 + bvrRate*(tinjury);
                    *svr = 0;
                }
            };
            bool endCondition(double Thy)
            {
                return Thy<=10+273.15;
            };
    };
    class KonstasCase : public SimModel
    {
        public:
            double tinjury = 300; // s
            double bvrRate = -.4/tinjury; // m^3/m^3 / s
            double trestore = 120; // s, time to restore blood volume
            KonstasCase() : SimModel() {
                thermovariant = 1;
                transient = 1;
                tFinal = 3600.0;
                dt = 1.0;
            };
            void rhocp(
                double *Rho,
                double *Cp,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double DViv,
                double fEB,
                double fES,
                double time
            ) override {
                if(eleIdx<3 && time>=trecovery){
                    Rho[eleIdx] = (salineRho*(DViv/(Viv+DViv)*fEB+fES)+
                            body->rhoBlood*(Viv/(Viv+DViv)*fEB))/
                            (fES+fEB);
                    Cp[eleIdx] = (salineRho*salineCp*(DViv/(Viv+DViv)*fEB+fES)+
                            body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)*fEB))/
                            (salineRho*(DViv/(Viv+DViv)*fEB+fES)+
                            body->rhoBlood*(Viv/(Viv+DViv)*fEB));
                } else {
                    Rho[eleIdx] = salineRho*(DViv/(Viv+DViv))+body->rhoBlood*(Viv/(Viv+DViv));
                    Cp[eleIdx] = (salineRho*salineCp*(DViv/(Viv+DViv))+body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)))/
                            salineRho*(DViv/(Viv+DViv))+body->rhoBlood*(Viv/(Viv+DViv));
                }
                // cout << "ElemIdx" << eleIdx << " Rho " << Rho[eleIdx] << " Cp " << Cp[eleIdx] << endl;
            }
            void setTairTsrm(double time) override
            {
                if(time < trecovery){
                    Tair = TairOutdoors;
                    Tsrm = TsrmOutdoors;
                    echelonSrm = echelonSrOutdoors;
                    alphaSrm = alphaSrOutdoors;
                    s = sOutdoors;
                }else{
                    Tair = TairIndoors;
                    Tsrm = TsrmIndoors;
                    echelonSrm = echelonSrIndoors;
                    alphaSrm = alphaSrIndoors;
                    s = sIndoors;
                }
            }
            void BVRSVR( double *bvr, double *svr, double time ) override
            {
                if(time < trecovery){
                    *bvr = 1+bvrRate*min(tinjury,time);
                    *svr = 0.0;
                }
                else{
                    *bvr = 1+bvrRate*min(tinjury,trecovery);
                    *svr = (1-*bvr)*min(1.0,(time-trecovery)/trecovery);
                }
            };

            void ecmoBlood(
                double *fEB,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double time
            ) override {
                if(eleIdx<3 && time>=trecovery){
                    *fEB = Viv/body->Viv0 * 250 / 60 * 1e-6;
                } else {
                    *fEB = 0;
                }
            }
            void ecmoSaline(
                double *fES,
                BODYMODEL::BodyModel *body,
                int eleIdx,
                double Viv,
                double bvr,
                double time
            ) override {
                if(eleIdx<3 && time>=trecovery && time < trestore+trecovery){
                    *fES = body->Viv0*(1.0-bvr)/trestore;
                } else if(eleIdx<3 && time>=trecovery && time >= trestore+trecovery){
                    *fES = 0;
                } else{
                    *fES = 0;
                }
            }

            void tblp(
                    double *tblp,
                    double *tblpNxtRatio,
                    BODYMODEL::BodyModel *body,
                    int eleIdx,
                    double tblpCur,
                    double Viv,
                    double DViv,
                    double fEB,
                    double fES,
                    double time
            ) override {
                if(eleIdx<3 && time>trecovery)
                    tblp[eleIdx] = (salineRho*salineCp*(DViv/(Viv+DViv)*fEB*TECMO+fES*TECMO)+
                    body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)*fEB*TECMO))/
                    (salineRho*salineCp*(DViv/(Viv+DViv)*fEB+fES)+
                    body->rhoBlood*body->cpBlood*(Viv/(Viv+DViv)*fEB));
                else
                    tblp[eleIdx] = tblpCur;
                tblpNxtRatio[eleIdx] = tblp[eleIdx]/tblpCur;
                assert(!isnan(abs(tblp[eleIdx])));
            }
            void shiverDrugs(
                    double *Sh,
                    double time
            ) override {
                if(time > trecovery+30)
                    *Sh = 0;
            }
            bool endCondition(double Thy) override
            {
                return Thy<=10+273.15;
            };

    };
}

#endif