#ifndef _SIMULATOR_CPP_
#define _SIMULATOR_CPP_

#include "Simulator.hpp"

namespace SIMULATOR{
    void Simulator::computeThermalLoadParameters(double* M, double* H, double*Q)
    {
        double Mbas0 = 0; // W
        double MbasDelta = 0; // W
        double MperWork = 0.2; // W/W
        int Tidx = 0;
        for(int eleIdx = 0; eleIdx < body->nElements; ++eleIdx){
            ELEMENT::Element* elem = body->elements[eleIdx];
            WASHER::Washer* cur = elem->washers[0];
            Mbas0 += cur
        }
    }
}

#endif