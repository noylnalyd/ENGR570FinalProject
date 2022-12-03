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
        // Private variables (set by BC)
        double s; // W/m^2, Light intensity
        double echelonSf; // -, Skin surface emissivity
        double echelonSrm; // -, Surroundings emissivity
        double alphaSf; // -, Skin surface absorptivity
        double psi; // -, View factor

        // UQ-dependencies
        double TairIndoors; // K, inside / treatment air temperature
        double TairOutdoors; // K, outside / pretreatment air temperature
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
        double Q10; // Q10 constant. Often 2, sometimes 3
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

        // BC values to update at each iteration
        double Tsrm; // K, surrounding temperature mean
        double Tair; // K, air temperature

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
        virtual void skinBC(
                double *Tpp,
                BODYMODEL::BodyModel* body,
                double *T,
                double *T0,
                double Sw,
                double time
        );
        void setUQs( double args[] );
        virtual void ecmoBlood( double *fEB, int elemIdx, double time);
        virtual void ecmoSaline( double *fES, int elemIdx, double time);
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
            idx++;
            washer = element->washers[element->nWashers-1];
            for(int sectIdx=0;sectIdx<element->nSectors;sectIdx++){
                idx+=element->nWashers-2;
                sector = element->sectors[sectIdx];
                double qsk = skinHeat(body,element,sector,elemIdx,sectIdx,Sw,Tpp[idx],T0[idx],element->n_j,sector->psiStnd);
                // Westin eqn 64
                Tpp[idx] = qsk*(washer->r+washer->deltaR/2.0)/washer->k
                        *log((washer->r+washer->deltaR)/washer->r)
                        +T[idx];
            }
        }
    }
    void SimModel::BVRSVR( double *bvr, double *svr, double time )
    {
        *bvr = 1.0;
        *svr = 0.0;
    }
    void SimModel::ecmoBlood( double *fEB, int elemIdx, double time)
    {
        *fEB = 0;
    }
    void SimModel::ecmoSaline( double *fES, int elemIdx, double time )
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
        double TairC = TairC; // C

        double hcmix = sqrt(abs(
                elem->a_nat*sqrt(abs(TsfC-TairC))
                +elem->a_frc*vAirEff
                +elem->a_mix
        )); // W/m^2/K

        double hR = body->StefanBoltzmann*echelonSf*echelonSrm
                *(TsfAbs*TsfAbs + TsrmAbs*TsrmAbs)
                *(TsfAbs + TsrmAbs); // W/m^2/K
        
        double U = 1/(njActual*body->Icl+1/(f*(hcmix+hR))); // W/m^2/K

        qsR = alphaSf*psiActual*s;

        double Tweight = (hcmix*TairC + hR*(TsrmC) + qsR)/(hcmix + hR);
        double Posksat = 100.0*exp(18.956-4030.0/(TsfC+235));
        
        double Ue = body->LewisConstant/(
            njActual*body->Icl/body->icl
            +1/(f*hcmix)
        );
        double dsweat = elem->a_sw*sect->areaShare*Sw*pow(2,(TsfAssumed-Tsf0)/Q10) * 1e-3/60;
        double Psk = (body->LambdaH20*dsweat + Posksat*body->Rinverse + Ue*Pair)/
                (Ue + body->Rinverse);
        
        return U*(TsfC-Tweight) + Ue*(Psk-Pair);
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
                tFinal = 1.0;
                dt = 1.0;
            };

    };
    class SteadyCase : public SimModel
    {
        public:
            SteadyCase() : SimModel() {
                thermovariant = 1;
                transient = 0;
                tFinal = 1.0;
                dt = 1.0;
            };
    };
    class TransientCase : public SimModel
    {
        public:
            TransientCase() : SimModel() {
                thermovariant = 1;
                transient = 1;
                tFinal = 200.0;
                dt = 1.0;
            };
    };
    class InjuryCase : public SimModel
    {
        public:
            double tinjury = 0;
            double trecovery = 300;
            double bvrFinal = .6;
            InjuryCase() : SimModel() {
                thermovariant = 1;
                transient = 1;
                tFinal = 3600.0;
                dt = 1.0;
            };
            void BVRSVR( double *bvr, double *svr, double time ) override
            {
                if(time > tinjury && time < trecovery)
                    *bvr = bvrFinal + (1-bvrFinal)*(time-tinjury)/(trecovery-tinjury);
                    *svr = 0;
            };
    };
}

#endif