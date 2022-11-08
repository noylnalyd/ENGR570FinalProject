#ifndef _SIMMODEL_
#define _SIMMODEL_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>

using namespace std;

namespace SIMMODEL
{
    enum SimState { undefined, initialized, elementsAdded, elementsComputed, staticPBMAssembled};


    class SimModel
    {
    private:
        SimState _state = undefined;
    public:
        // Model attributes
        int modelCase; // Number of case (if automating)

        // Saline constants
        const double salineRho; // kg/m^3, saline density
        const double salineCp; // J/kg/K, saline specific heat capacity

        // Environmental constants
        const double Pair; // Pa, air pressure
        const double vAirEff; // m/s, effective air velocity

        // UQ attributes worth exploring
        double TairIndoors; // K, inside / treatment air temperature
        double TairOutdoors; // K, outside / pretreatment air temperature
        double TsrmIndoors; // K, mean surroundings temperature inside
        double TsrmOutdoors; // K, mean surroundings temperature outside

    };


}

#endif