#ifndef _SIMULATOR_CPP_
#define _SIMULATOR_CPP_

#include "Simulator.hpp"

namespace SIMULATOR
{

    void Simulator::preallocate(){
        assert(_state == SIMULATOR::simLoaded);
        
    }
    void Simulator::runSim(){
        
        // Check status of body and simmodel. Make sure both completed.
        assert(body->getState() == BODYMODEL::staticPBMAssembled);
        assert(sim->getState() == SIMMODEL::allValuesAssigned);

        

    }
}

int main(){
    BODYMODEL::BodyModel* body = BODYMODEL::defaultBody();
    return 0;
}
#endif