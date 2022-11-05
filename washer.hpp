#ifndef _WASHER_
#define _WASHER_

#include <iostream>
#include <cstddef>
#include <math.h>

using namespace std;

namespace WASHER
{
    class Washer
    {
    private:

    public:
        double r; // m
        double k; // W/m^2
        double rho; // kg/m^3
        double c; // J/kg/K
        double w_bl; // 1/s
        double q_m; // W/m^3
        double a_resp; // -
        double zeta; // s
        double gamma; // -
        double del; // m^4/W
        double deltaR; // m
        double volume; // m^3
        double AForwardCur=0; // -
        double AForwardNxt=0; // -
        double ABackwardCur=0; // -
        double ABackwardPrv=0; // -

        Washer();
        virtual ~Washer();

        void computeZeta()
        {
        zeta = 2*deltaR*deltaR*rho*c/k;
        }
        void computeGamma(double omega)
        {
        gamma = deltaR/r/omega;
        }
        void computeDel()
        {
        del = deltaR*deltaR/k;
        }
    };

    class CylinderWasher : public Washer
    {
        public:
            void computeAForward( CylinderWasher nxt )
            {
                double rIFC = r+deltaR/2;
                double L = k/log((r+deltaR)/r);
                double Lnxt = nxt.k/log(nxt.r/(nxt.r-nxt.deltaR));
                double D = log(rIFC/r)/log((r+deltaR)/r);
                double Ds = log(rIFC/(r+deltaR))/log((r+deltaR)/r);
                double Dnxt = log(rIFC/nxt.r)/log(nxt.r/(nxt.r+nxt.deltaR));

                AForwardCur = (Ds*Lnxt-Dnxt*L)/(D*Lnxt-Dnxt*L);
                AForwardNxt = Lnxt/(D*Lnxt-Dnxt*L);
            }
            void computeABackward( CylinderWasher prv )
            {
                double rIFC = r-deltaR/2;
                double L = k/log((r)/(r-deltaR));
                double Lprv = prv.k/log((prv.r+prv.deltaR)/prv.r);
                double D = log(rIFC/r)/log(r/(r-deltaR));
                double Ds = log(rIFC/(r-deltaR))/log(r/(r-deltaR));
                double Dprv = log(rIFC/prv.r)/log((prv.r+prv.deltaR)/prv.r);

                ABackwardCur = (Dprv*L-Ds*Lprv)/(Dprv*L-D*Lprv);
                ABackwardPrv = Lprv/(Dprv*L-D*Lprv);
            }
            void computeVolume( double length /*m*/ ){
                volume = length*M_PI*((r+deltaR/2)*(r+deltaR/2)-(r-deltaR/2)*(r-deltaR/2));
            }
    };

    class CylinderCoreWasher : public CylinderWasher
    {
        public:
            double theta; // 1/rad
            void computeTheta( double sumPhi /*rad*/ )
            {
                theta = 8.0*(1.0-1.0/sqrt(2.0))*(1.0-1.0/sqrt(2.0))/(log(2.0*sqrt(2.0)-1.0)*sumPhi);
            }
    };

    class SphereWasher : public Washer
    {
        public:
            void computeAForward( CylinderWasher nxt )
            {
                double L = k/(1/r-1/(r+deltaR));
                double Lnxt = nxt.k/(1/(nxt.r-nxt.deltaR)-1/nxt.r);

                AForwardCur = (L-Lnxt*r/(nxt.r-nxt.deltaR))/(L+Lnxt*((r+deltaR)/(nxt.r-nxt.deltaR)));
                AForwardNxt = (Lnxt*(1+nxt.r/(nxt.r-nxt.deltaR)))/(L+Lnxt*((r+deltaR)/(nxt.r-nxt.deltaR)));
            }
            void computeABackward( CylinderWasher prv )
            {
                double L = k/(1/(r-deltaR)-1/r);
                double Lprv = prv.k/(1/prv.r-1/(prv.r+prv.deltaR));

                ABackwardCur = (L-Lprv*(r/(prv.r+prv.deltaR)))/(L+Lprv*((r-deltaR)/prv.r+prv.deltaR));
                ABackwardPrv = (Lprv*(1+prv.r/(prv.r+prv.deltaR)))/(L+Lprv*((r-deltaR)/(prv.r+prv.deltaR)));
            }
            void computeVolume(){
                volume = 4.0*M_PI/3.0*((r+deltaR/2)*(r+deltaR/2)*(r+deltaR/2)-(r-deltaR/2)*(r-deltaR/2)*(r-deltaR/2));
            }
    };

    class SphereCoreWasher : public SphereWasher
    {
        public:
            double theta; // 1/rad
            void computeTheta( double sumPhi /*rad*/ )
            {
                theta = 12*(cbrt(16)-1)*(1-1/cbrt(2))*(1-1/cbrt(2))/((cbrt(32)-cbrt(16))*sumPhi);
            }
    };


}

#endif