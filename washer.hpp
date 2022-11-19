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
        void computeGamma();
        void computeDel();
        virtual void computeVolume() = 0;
        
    public:
        double r; // m
        double k; // W/m^2
        double rho; // kg/m^3
        double c; // J/kg/K
        double w_bl; // 1/s
        double q_m; // W/m^3
        double a_resp; // -
        double omega=0; // -, False value
        double zeta; // s
        double gamma; // -
        double del; // m^4/W
        double deltaR; // m
        double volume; // m^3
        double AForwardCur=0; // -
        double AForwardNxt=0; // -
        double ABackwardCur=0; // -
        double ABackwardPrv=0; // -
        
        virtual void computeAForward( Washer* nxt ) = 0;
        virtual void computeABackward( Washer* prv ) = 0;

        Washer(){
            _state = initialized;
        };
        ~Washer();
        void sharedCompute();
        void compute();

    };

    void Washer::computeZeta()
    {
        zeta = 2*deltaR*deltaR*rho*c/k;
    }
    void Washer::computeGamma()
    {
        gamma = deltaR/r/omega;
    }
    void Washer::computeDel()
    {
        del = deltaR*deltaR/k;
    }
    void Washer::sharedCompute()
    {
        computeDel();
        computeGamma();
        computeZeta();
        computeVolume();
        _state = computed;
    }
    void Washer::compute()
    {
        sharedCompute();
    }

    class CylinderWasher : public Washer
    {
        private:
            void computeVolume(){
                volume = length*M_PI*((r+deltaR/2)*(r+deltaR/2)-(r-deltaR/2)*(r-deltaR/2));
            }
        public:
            double length;
            CylinderWasher() : Washer() {
                omega = 1.0;
            }

            void computeAForward( Washer* nxt )
            {
                double rIFC = r+deltaR/2;
                double L = k/log((r+deltaR)/r);
                double Lnxt = nxt->k/log(nxt->r/(nxt->r-nxt->deltaR));
                double D = log(rIFC/r)/log((r+deltaR)/r);
                double Ds = log(rIFC/(r+deltaR))/log((r+deltaR)/r);
                double Dnxt = log(rIFC/nxt->r)/log(nxt->r/(nxt->r+nxt->deltaR));

                AForwardCur = (Ds*Lnxt-Dnxt*L)/(D*Lnxt-Dnxt*L);
                AForwardNxt = Lnxt/(D*Lnxt-Dnxt*L);
            }
            void computeABackward( Washer* prv )
            {
                double rIFC = r-deltaR/2;
                double L = k/log((r)/(r-deltaR));
                double Lprv = prv->k/log((prv->r+prv->deltaR)/prv->r);
                double D = log(rIFC/r)/log(r/(r-deltaR));
                double Ds = log(rIFC/(r-deltaR))/log(r/(r-deltaR));
                double Dprv = log(rIFC/prv->r)/log((prv->r+prv->deltaR)/prv->r);

                ABackwardCur = (Dprv*L-Ds*Lprv)/(Dprv*L-D*Lprv);
                ABackwardPrv = Lprv/(Dprv*L-D*Lprv);
            }
            
    };

    class CylinderCoreWasher : public CylinderWasher
    {
        private:
            void computeTheta()
            {
                theta = 8.0*(1.0-1.0/sqrt(2.0))*(1.0-1.0/sqrt(2.0))/(log(2.0*sqrt(2.0)-1.0)*sumPhi);
            }
        public:
            double theta; // 1/rad
            double sumPhi; // rad
            
            void compute(){
                computeTheta();
                sharedCompute();
            }
    };

    class SphereWasher : public Washer
    {
        public:
            SphereWasher() : Washer() {
                omega = 2.0;
            }
            void computeAForward( Washer* nxt )
            {
                double L = k/(1/r-1/(r+deltaR));
                double Lnxt = nxt->k/(1/(nxt->r-nxt->deltaR)-1/nxt->r);

                AForwardCur = (L-Lnxt*r/(nxt->r-nxt->deltaR))/(L+Lnxt*((r+deltaR)/(nxt->r-nxt->deltaR)));
                AForwardNxt = (Lnxt*(1+nxt->r/(nxt->r-nxt->deltaR)))/(L+Lnxt*((r+deltaR)/(nxt->r-nxt->deltaR)));
            }
            void computeABackward( Washer* prv )
            {
                double L = k/(1/(r-deltaR)-1/r);
                double Lprv = prv->k/(1/prv->r-1/(prv->r+prv->deltaR));

                ABackwardCur = (L-Lprv*(r/(prv->r+prv->deltaR)))/(L+Lprv*((r-deltaR)/prv->r+prv->deltaR));
                ABackwardPrv = (Lprv*(1+prv->r/(prv->r+prv->deltaR)))/(L+Lprv*((r-deltaR)/(prv->r+prv->deltaR)));
            }
            void computeVolume(){
                volume = 4.0*M_PI/3.0*((r+deltaR/2)*(r+deltaR/2)*(r+deltaR/2)-(r-deltaR/2)*(r-deltaR/2)*(r-deltaR/2));
            }
    };

    class SphereCoreWasher : public SphereWasher
    {
        private:
            void computeTheta()
            {
                theta = 12*(cbrt(16)-1)*(1-1/cbrt(2))*(1-1/cbrt(2))/((cbrt(32)-cbrt(16))*sumPhi);
            }
        public:
            double theta; // 1/rad
            double sumPhi; // rad
            
            void compute(){
                computeTheta();
                sharedCompute();
            }
    };


}

#endif