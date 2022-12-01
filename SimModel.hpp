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
    public:
        // Model attributes
        int modelCase; // Number of case (if automating)
        double thermovariant; // 0 for thermoneutral, 1 for variant
        double transient; // 0 for steady, 1 for transient
        double Q10; // Q10 constant. Often 2, sometimes 3
        const double KonstasBeta = 0.08401; // K^-1
        const double KonstasGamma =  2.245; // -
        const double KonstasAlpha = 2.961; // -
        
        // Saline constants
        const double salineRho; // kg/m^3, saline density
        const double salineCp; // J/kg/K, saline specific heat capacity

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
        double TairIndoors; // K, inside / treatment air temperature
        double TairOutdoors; // K, outside / pretreatment air temperature
        double TsrmIndoors; // K, mean surroundings temperature inside
        double TsrmOutdoors; // K, mean surroundings temperature outside

        // BC values
        SimModelState getState();
        void skinHeat(BODYMODEL::BodyModel* body,
            int elemIdx,
            int sectIdx,
            double Twash, // K, outermost node temperature
            double Sw, // g/min, sweat command
            double TsfAssumed, // K, assumed skin surface temperature
            int njActual // -, number of equivalent clothings on this sector
            );
        virtual void skinBC( int element, int sector, double time, double* T, double* Amatrixii, double* Trhs);
        virtual void BVRSVR( double *bvr, double *svr, double time );
    };

    void SimModel::BVRSVR( double *bvr, double *svr, double time )
    {
        *bvr = 1.0;
        *svr = 1.0;
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

        qsR = alphaSr*psi*s;

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
    class InitialCase : SimModel
    {
        thermovariant = 0;
        transient = 0;

    }
    class SteadyCase : SimModel
    {
        thermovariant = 0;
        transient = 1;
    }
    class TransientCase : SimModel
    {

    }
    class KonstasCase : SimModel
    {

    };


}

#endif