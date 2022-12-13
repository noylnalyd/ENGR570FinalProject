#ifndef _WASHER_
#define _WASHER_

#include <iostream>
#include <cstddef>
#include <math.h>


namespace WASHER
{
    /// @brief Enum for state of washer.
    enum WasherState { undefined, initialized, computed };

    /// @brief Class for all washers.
    class Washer
    {
    private:
        /// @brief State of this washer.
        WasherState _state = undefined;
        /// @brief Compute this washer's zeta.
        void computeZeta();
        /// @brief Compute this washer's del.
        void computeDel();
        
    public:
        double r; //!< m, Mean radius
        double k; //!< W/m/K, specific conductivity
        double rho; //!< kg/m^3, density
        double c; //!< J/kg/K, specific constant pressure heat capacity
        double w_bl; //!< 1/s, blood perfusivity
        double q_m; //!< W/m^3, basal metabolic heat gen
        double a_resp=0; //!< -, respiratory power share
        double zeta; //!< s, Time-dependence coefficient
        double gamma; //!< -, volume shear coefficient
        double del; //!< m^4/W, flux to perfusion ratio
        double deltaR; //!< m, outer minus inner radius
        double volume; //!< m^3, volume for whole torus
        double muscle; //!< -, muscle ratio (1 or 0)
        double resp; //!< -, resp ratio per washer
        double AForwardCur=0; //!< -, Forward current coefficient
        double AForwardNxt=0; //!< -, Forward relative coefficient
        double ABackwardCur=0; //!< -, Backward current coefficient
        double ABackwardPrv=0; //!< -, Backward relative coefficient

        /// @brief Constructor
        Washer(){
            _state = initialized;
        };
        /// @brief Destructor
        ~Washer();

        /// @brief Compute all subparameters.
        void compute();

    };

    void Washer::computeZeta()
    {
        zeta = 2*deltaR*deltaR*rho*c/k;
    }
    void Washer::computeDel()
    {
        del = deltaR*deltaR/k;
    }
    void Washer::compute()
    {
        assert(_state == initialized);
        computeZeta();
        computeDel();
        _state = computed;
    }


}

#endif