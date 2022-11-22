#ifndef _SIMULATOR_CPP_
#define _SIMULATOR_CPP_

#include "Simulator.hpp"

namespace SIMULATOR
{

    void Simulator::preallocate(){
        assert(_state == SIMULATOR::simLoaded);

    }
    void Simulator::computeThermalLoadParameters()
    {
        double Mbas0 = 0; // W
        double MbasDelta = 0; // W
        double MperWork = 0.2; // W/W
        int Tidx = 0;
        for(int eleIdx = 0; eleIdx < body->nElements; ++eleIdx){
            ELEMENT::Element* elem = body->elements[eleIdx];
            WASHER::Washer* cur = elem->washers[0];
            Mbas0 += cur->volume*elem->sumPhi/(2*M_PI)*cur->q_m;
            MbasDelta += cur->volume*elem->sumPhi/(2*M_PI)*sim->thermoneutral*deltaQMetabolic(cur->q_m,T[Tidx],Tn[Tidx],sim->Q10);
            Tidx++;
            for(int secIdx = 0; secIdx<elem->nSectors; ++secIdx){
                SECTOR::Sector* sect = elem->sectors[secIdx];
                for(int washIdx = 1; washIdx<elem->nWashers; ++washIdx){
                    cur = elem->washers[washIdx];
                    Mbas0 += cur->volume*sect->phi/(2*M_PI)*cur->q_m;
                    MbasDelta += cur->volume*sect->phi/(2*M_PI)*sim->thermoneutral*deltaQMetabolic(cur->q_m,T[Tidx],Tn[Tidx],sim->Q10);
                    Tidx++;
                }
            }
        }
        double act = 1 + (Sh+MbasDelta)/Mbas0*MperWork; // MET
        double etaW;
        if(act < 1.6)
            etaW = 0.05; // W/W
        else
            etaW = 0.2*tanh(body->b1*act+body->b0); // W/W
        M = act*Mbas0/actBas; // W
        H = M*(1-etaW)-Mbas0; // W
        Qresp = computeQresp(); // W

    }
    double computeQresp(){
        return 3.45*M*(0.028+6.5e-5*sim->Tair-4.98e-6*sim->Pair)+
            1.44e-3*M*(32.6-0.934*sim->Tair+1.99e-4*sim->Pair);
    }
    double deltaQMetabolic( double q, double T, double Tn, double Q10 )
    {
        return q*(pow(Q10,(T-Tn)/10.0)-1);
    }

    void Simulator::runSim(){
        
        // Check status of body and simmodel. Make sure both completed.
        assert(body->getState() == BODYMODEL::staticPBMAssembled);
        assert(sim->getState() == SIMMODEL::allValuesAssigned);

        // Declare local temporary values
        ELEMENT::Element* element;
        SECTOR::Sector* sector;
        WASHER::Washer* washer;

        // Provide initial values

        // Iterate
        for(int timestep=1;timestep<nSteps;++timestep){
            
            // Compute temporary BC values and properties from sim

            // Skin handling, Brain handling

            // Compute error inputs

            // Compute active controls

            // Compute whole-body thermal load parameters

            // Project to get future values

            // Empty PBM

            // Loop thru elements
            for(int eleIdx = 0; eleIdx < body->nElements; ++eleIdx){
                element = body->elements[eleIdx];
                // Compute element-level blood pool properties
                
                // Initialize core node values, including cusp-less projection

                // Loop thru sectors
                for(int sectIdx = 0; sectIdx < elem->nSectors; ++eleIdx){
                    sector = element->sectors[sectIdx];
                    // Loop thru washers
                }

            }

        }

    }
}

int main(){
    BODYMODEL::BodyModel* body = BODYMODEL::defaultBody();
    return 0;
}
#endif