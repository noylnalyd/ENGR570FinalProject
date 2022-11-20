#ifndef _ELEMENT_CPP_
#define _ELEMENT_CPP_

#include "element.hpp";

namespace ELEMENT{
    void Element::subCompute( double* totalSkinArea ){
        // Assert all w/s addded
        assert(_state==wsAdded);
        sumPhi = 0;
        double rskin = washers[nWashers-1]->r+washers[nWashers-1]->deltaR/2.0;
        double coeff = 4*M_PI*rskin*rskin/(2*M_PI);
        for(int i=0;i<nSectors;i++){
            sumPhi += sectors[i]->phi;
            sectors[i]->areaSkin = coeff*sectors[i]->phi;
            *totalSkinArea += sectors[i]->areaSkin;
        }
    };
    void Element::compute( double totalSkinArea ){
        // Assert all w/s addded
        assert(_state==wsAdded);
        // Compute skin share
        for(int i=0;i<nSectors;i++)
            sectors[i]->areaShare = sectors[i]->areaSkin/totalSkinArea;

        // Compute basic attributes
        for(int i=0;i<nWashers;i++){
            washers[i]->compute();
        }

        // Compute relational attributes
        for(int i=0;i<nWashers-1;i++){
            computeAForward(washers[i], washers[i+1]);
        }
        for(int i=1;i<nWashers;i++){
            computeABackward(washers[i],washers[i-1]);
        }

        computeN();
        _state = computed;
    };

    void Element::computeN()
    {
        N = 1 + (nWashers-1)*nSectors;
    }
    ElementState Element::getState()
    {
        return _state;
    }
    void Element::addWashers(
        int nAdd,   // Number of washers to add
        double r0,  // m, inner radius
        double rf,  // m, outer radius
        double k,   // W/m^2, conductivity
        double rho, // kg/m^3, density
        double c,   // J/kg/K, specific heat capacity
        double w_bl,// 1/s, tissue permeability
        double q_m) // W/m^3, specific basal metabolism
    {
        // Must have been allocated!
        assert(_state==allocated);
        // Must have length configured!
        assert(length>-2);
        // Must have space remaining!
        assert(nAdd+washerIdx<=nWashers);
        // Must have defined omega!
        assert(omega>0);

        double drr = (rf-r0)/nAdd;
        for(int i=0;i<nAdd;i++){
            double r = r0+(i+.5)*drr;
            if(isCylinder){
                WASHER::CylinderWasher* tmp = new WASHER::CylinderWasher();
                tmp->r = r;
                tmp->k = k;
                tmp->rho = rho;
                tmp->c = c;
                tmp->w_bl = w_bl;
                tmp->q_m = q_m;
                tmp->deltaR = drr;
                washers[washerIdx++] = tmp;
            }
            else{
                WASHER::SphereWasher* tmp = new WASHER::SphereWasher();
                tmp->r = r;
                tmp->k = k;
                tmp->rho = rho;
                tmp->c = c;
                tmp->w_bl = w_bl;
                tmp->q_m = q_m;
                tmp->deltaR = drr;
                washers[washerIdx++] = tmp;
            }
        }
        if(sectorIdx==nSectors && washerIdx==nWashers)
            _state = wsAdded;
    }
    void Element::addSector(SECTOR::Sector* sec){
        // Must have been allocated!
        assert(_state==allocated);
        // Must have space remaining!
        assert(sectorIdx<nSectors);

        sectors[sectorIdx++] = sec;
        if(sectorIdx==nSectors && washerIdx==nWashers)
            _state = wsAdded;
    }
}

#endif