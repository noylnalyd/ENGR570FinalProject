#ifndef _WASHER_
#define _WASHER_

#include <iostream>
#include <cstddef>
#include <math.h>


namespace WASHER
{
    enum WasherState { undefined, initialized, computed };

    class Washer
    {
    private:
        WasherState _state = undefined;
        void computeZeta();
        void computeDel();
        
    public:
        double r; // m
        double k; // W/m^2
        double rho; // kg/m^3
        double c; // J/kg/K
        double w_bl; // 1/s
        double q_m; // W/m^3
        double a_resp=0; // -
        double zeta; // s
        double gamma; // -
        double del; // m^4/W
        double deltaR; // m
        double volume; // m^3
        double muscle; // -
        double resp; // -
        double AForwardCur=0; // -
        double AForwardNxt=0; // -
        double ABackwardCur=0; // -
        double ABackwardPrv=0; // -

        Washer(){
            _state = initialized;
        };
        ~Washer();

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