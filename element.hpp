#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include <vector>
#include <memory>
#include "washer.hpp"
#include "sector.hpp"

using namespace std;

namespace ELEMENT
{
    enum ElementState { undefined, initialized, allocated, wsAdded, subcomputed, computed};

    class Element
    {
    private:
        ElementState _state = undefined;
        void computeN();
        

    public:
        // Class attributes
        bool isCylinder;
        int N; // Number of distinct T nodes
        int nWashers; // Number of washers
        int washerIdx=0; // index of free washer
        int nSectors; // Number of sectors
        int sectorIdx=0; // index of free sector

        // Body level attributes
        double length = -2; // m, length of element iff cylinder
        double hx; // W/K, countercurrent heat transfer coefficient
        int n_j; // -, number of equivalent clothing layers
        double omega = -1; // - , shape parameter (-1 for undefined)
        double sumPhi = 0; // rad, sum of sector angles
        double theta = 0; // 1/rad, geometry parameter
        double Vmuscle = 0; // m^3, volume of muscle

        // Fractional coefficients
        double a_nat; // Share of natural convection
        double a_frc; // Share of forced convection
        double a_mix; // Share of mixed convection
        double a_sk; // Share of skin area
        double a_sed; // Share of sedentary load
        double a_stnd; // Share of standard load
        double a_sw; // Share of sweat
        double a_dl; // Share of vasodilation
        double a_cs; // Share of vasoconstriction
        double a_sh; // Share of shivering

        // Precomputation methods
        void subCompute( double* totalSkinArea );
        void compute( double totalSkinArea );
        virtual void computeVolume( WASHER::Washer* cur ) = 0;
        virtual void computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt) = 0;
        virtual void computeABackward( WASHER::Washer* cur, WASHER::Washer* prv) = 0;
        void computeGamma( WASHER::Washer* cur );
        virtual void computeTheta() = 0;

        void allocate(int nWash, int nSect){
            nWashers = nWash;
            nSectors = nSect;
            washers = new WASHER::Washer*[nWashers];
            sectors = new SECTOR::Sector*[nSectors];
            _state = allocated;
        }
        ElementState getState();
        void addWashers(
            int nAdd,   // Number of washers to add
            double r0,  // m, inner radius
            double rf,  // m, outer radius
            double k,   // W/m^2, conductivity
            double rho, // kg/m^3, density
            double c,   // J/kg/K, specific heat capacity
            double w_bl,// 1/s, tissue permeability
            double q_m, // W/m^3, specific basal metabolism
            double muscle);
        void addSector(SECTOR::Sector* sec);

        SECTOR::Sector **sectors; // Array of sectors
        WASHER::Washer **washers; // Array of washers

        Element()
        {
            assert(_state==undefined);
            _state = initialized;
        }
        
        Element(int nWash, int nSect)
        {
            assert(_state==undefined);
            _state = initialized;
            allocate(nWash, nSect);
        }
        ~Element()
        {
            delete [] washers;
            delete [] sectors;
        }
    };
    class CylinderElement : public Element{
        public:
            CylinderElement() : Element()
            {
                isCylinder = true;
                omega = 1;
            }
            CylinderElement( int nWash, int nSect ) : CylinderElement()
            {
                allocate(nWash,nSect);
            }
            void computeVolume( WASHER::Washer* cur ) override;
            void computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt) override;
            void computeABackward( WASHER::Washer* cur, WASHER::Washer* prv) override;
            void computeTheta() override;
    };
    class SphereElement : public Element{
        public:
            SphereElement() : Element()
            {
                isCylinder = false;
                omega = 2;
            }
            SphereElement( int nWash, int nSect ) : SphereElement()
            {
                allocate(nWash,nSect);
            }
            void computeVolume( WASHER::Washer* cur ) override;
            void computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt) override;
            void computeABackward( WASHER::Washer* cur, WASHER::Washer* prv) override;
            void computeTheta() override;
    };
}

#endif