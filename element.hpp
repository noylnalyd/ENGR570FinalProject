#ifndef _ELEMENT_
#define _ELEMENT_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include "washer.hpp"
#include "sector.hpp"

using namespace std;

namespace ELEMENT
{
    enum ElementState { undefined, initialized, allocated, wsAdded, computed};

    class Element
    {
    private:
        ElementState _state = undefined;
        void computeN();
        void addWashers(
            int nAdd,   // Number of washers to add
            double r0,  // m, inner radius
            double rf,  // m, outer radius
            double k,   // W/m^2, conductivity
            double rho, // kg/m^3, density
            double c,   // J/kg/K, specific heat capacity
            double w_bl,// 1/s, tissue permeability
            double q_m); // W/m^3, specific basal metabolism

    public:
        bool isCylinder;
        int N; // Number of distinct T nodes
        int nWashers; // Number of washers
        int washerIdx=0; // index of free washer
        int nSectors; // Number of sectors
        double length = -2; // m, length of element iff cylinder
        double hx; // W/K, countercurrent heat transfer coefficient
        int n_j; // -, number of equivalent clothing layers
        double omega=-1; // -, 1 for cylinder / 2 for sphere
        double sumPhi; // rad, sum of sector angles

        void compute();
        

        SECTOR::Sector **sectors; // Array of sectors
        WASHER::Washer **washers; // Array of washers

        Element()
        {
            assert(_state==undefined);

        }
        Element(int nWash, int nSect)
        {
            nWashers = nWash;
            nSectors = nSect;
            washers = new WASHER::Washer*[nWashers];
            sectors = new SECTOR::Sector*[nSectors];
        }
        ~Element()
        {
            delete [] washers;
            delete [] sectors;
        }
        
        


    };

    void Element::compute(){
        // Set core sumphi
        washers[0].

        // Compute basic attributes
        for(int i=0;i<nWashers;i++){
            washers[i].compute();
        }

        // Compute relational attributes
        for(int i=0;i<nWashers-1;i++){
            washers[i].computeAForward(&(washers[i+1]));
        }
        for(int i=1;i<nWashers;i++){
            washers[i].computeABackward(&(washers[i-1]));
        }

        computeN();
        _state = computed;
    }

    void Element::computeN()
    {
        N = 1 + (nWashers-1)*nSectors;
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
        // Must have core already loaded!
        assert(washerIdx>0);
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
                WASHER::CylinderWasher tmp;
                tmp.r = r;
                tmp.k = k;
                tmp.rho = rho;
                tmp.c = c;
                tmp.w_bl = w_bl;
                tmp.q_m = q_m;
                tmp.deltaR = drr;
                tmp.computeGamma(omega);
                tmp.computeVolume(length);
                washers[washerIdx++] = tmp;
            }
            else{
                WASHER::SphereWasher tmp;
                tmp.r = r;
                tmp.k = k;
                tmp.rho = rho;
                tmp.c = c;
                tmp.w_bl = w_bl;
                tmp.q_m = q_m;
                tmp.deltaR = drr;
                tmp.computeGamma(omega);
                tmp.computeVolume();
                washers[washerIdx++] = tmp;
            }
        }
    }

}

#endif