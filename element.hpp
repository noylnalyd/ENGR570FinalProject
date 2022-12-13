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
    /// @brief Enum for element's state
    enum ElementState { undefined, initialized, allocated, wsAdded, computed};

    /// @brief Abstract class for a Fiala element.
    class Element
    {
    private:
        /// @brief This element's state.
        ElementState _state = undefined;
        /// @brief Compute number of thermal nodes.
        void computeN();
        

    public:
        // Class attributes
        bool isCylinder; //!< Is a cylinder (Used sparingly!)
        int N; //!< Number of distinct T nodes
        int nWashers; //!< Number of washers
        int washerIdx=0; //!< index of free washer
        int nSectors; //!< Number of sectors
        int sectorIdx=0; //!< index of free sector

        // Body level attributes
        double length = -2; //!< m, length of element iff cylinder
        double hx; //!< W/K, countercurrent heat transfer coefficient
        int n_j; //!< -, number of equivalent clothing layers
        double omega = -1; //!< - , shape parameter (-1 for undefined)
        double sumPhi = 0; //!< rad, sum of sector angles
        double theta = 0; //!< 1/rad, geometry parameter
        double Vmuscle = 0; //!< m^3, volume of muscle
        double Vresp = 0; //!< m^3, volume of respiratory exposure

        // Fractional coefficients
        double a_nat; //!< Share of natural convection
        double a_frc; //!< Share of forced convection
        double a_mix; //!< Share of mixed convection
        double a_sk; //!< Share of skin area
        double a_sed; //!< Share of sedentary load
        double a_stnd; //!< Share of standard load
        double a_sw; //!< Share of sweat
        double a_dl; //!< Share of vasodilation
        double a_cs; //!< Share of vasoconstriction
        double a_sh; //!< Share of shivering

        // Precomputation methods

        /// @brief Compute skin area, add to totalSkinArea.
        /// @param[out] totalSkinArea Ptr to body skin area.
        virtual void subCompute( double* totalSkinArea ) = 0;

        /// @brief Compute element subparameters.
        /// @param[in] totalSkinArea Body total skin area.
        void compute( double totalSkinArea );

        /// @brief Compute this washer's volume.
        /// @param cur Washer to evaluate
        virtual void computeVolume( WASHER::Washer* cur ) = 0;

        /// @brief Compute this washer's AForward coefficients
        /// @param cur Washer to evaluate
        /// @param nxt Washer to compare to
        virtual void computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt) = 0;

        /// @brief Compute this washer's ABackward coefficients
        /// @param cur Washer to evaluate
        /// @param prv Washer to compare to
        virtual void computeABackward( WASHER::Washer* cur, WASHER::Washer* prv) = 0;

        /// @brief Compute this washer's AForward coefficients (skin only!)
        /// @param cur Washer to evaluate
        virtual void computeAForwardSkin( WASHER::Washer* cur ) = 0;

        /// @brief Find resp share of washer vs. element share.
        /// @param tresp Element resp share
        /// @param ricur inner radius of washer
        /// @param rocur outer radius of washer
        /// @param ri inner radius of washer set
        /// @param ro outer radius of washer set
        /// @return a_resp for one washer
        virtual double handlearesp(double tresp, double ricur, double rocur, double ri, double ro) = 0;

        /// @brief Computes gamma coefficient, coefficient of shearflow
        /// @param cur Washer to evaluate
        void computeGamma( WASHER::Washer* cur );

        /// @brief Computes theta, core washer coefficient
        virtual void computeTheta() = 0;

        /// @brief Allocates element arrays.
        /// @param nWash Number of washers
        /// @param nSect Number of sectors
        void allocate(int nWash, int nSect){
            nWashers = nWash;
            nSectors = nSect;
            washers = new WASHER::Washer*[nWashers];
            sectors = new SECTOR::Sector*[nSectors];
            _state = allocated;
        }

        /// @brief Allows access to protected __state
        /// @return __state
        ElementState getState();

        /// @brief Adds a block of washers to avoid redundancy.
        void addWashers(
            int nAdd,   //!< Number of washers to add
            double r0,  //!< m, inner radius of block
            double rf,  //!< m, outer radius of block
            double k,   //!< W/m^2, conductivity
            double rho, //!< kg/m^3, density
            double c,   //!< J/kg/K, specific heat capacity
            double w_bl,//!< 1/s, tissue permeability
            double q_m, //!< W/m^3, specific basal metabolism
            double muscle, //!< -, ratio of muscle (usually 0 or 1)
            double resp); //!< -, ratio of respiratory exposure

        
        /// @brief Adds a Sector to this element
        /// @param sec Sector to add.
        void addSector(SECTOR::Sector* sec);

        /// @brief Array of Sector ptrs.
        SECTOR::Sector **sectors; // Array of sectors

        /// @brief Array of Washer ptrs.
        WASHER::Washer **washers; // Array of washers

        /// @brief Blank constructor
        Element()
        {
            assert(_state==undefined);
            _state = initialized;
        }
        
        /// @brief Constructor with allocated arrays.
        /// @param nWash Number of washers
        /// @param nSect Number of sectors
        Element(int nWash, int nSect)
        {
            assert(_state==undefined);
            _state = initialized;
            allocate(nWash, nSect);
        }
        /// @brief Destructor
        ~Element()
        {
            delete [] washers;
            delete [] sectors;
        }
    };
    /// @brief Cylindrical element concrete subclass
    class CylinderElement : public Element{
        public:
            /*! @copydoc Element() */
            CylinderElement() : Element()
            {
                isCylinder = true;
                omega = 1;
            }
            /*! @copydoc Element(int,int) */
            CylinderElement( int nWash, int nSect ) : CylinderElement()
            {
                allocate(nWash,nSect);
            }
            void subCompute( double* totalSkinArea );
            double handlearesp(double tresp, double ricur, double rocur, double ri, double ro) override;
            void computeVolume( WASHER::Washer* cur ) override;
            void computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt) override;
            void computeAForwardSkin( WASHER::Washer* cur ) override;
            void computeABackward( WASHER::Washer* cur, WASHER::Washer* prv) override;
            void computeTheta() override;
    };

    /// @brief Spherical element concrete subclass
    class SphereElement : public Element{
        public:
            /*! @copydoc Element() */
            SphereElement() : Element()
            {
                isCylinder = false;
                omega = 2;
            }
            /*! @copydoc Element(int,int) */
            SphereElement( int nWash, int nSect ) : SphereElement()
            {
                allocate(nWash,nSect);
            }
            void subCompute( double* totalSkinArea );
            double handlearesp(double tresp, double ricur, double rocur, double ri, double ro) override;
            void computeVolume( WASHER::Washer* cur ) override;
            void computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt) override;
            void computeAForwardSkin( WASHER::Washer* cur ) override;
            void computeABackward( WASHER::Washer* cur, WASHER::Washer* prv) override;
            void computeTheta() override;
    };
}

#endif