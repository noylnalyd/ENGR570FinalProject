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

        SimModelState getState();
        virtual void skinBC( int element, int sector, double time, double* T, double* Amatrixii, double* Trhs)=0;
    };

    SimModelState SimModel::getState()
    {
        return _state;
    }
    class KonstasCase : SimModel
    {

    };


}

#endif