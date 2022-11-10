#ifndef _BODYMODEL_
#define _BODYMODEL_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include "element.hpp"
#include "PseudoBlockMatrix.hpp"

using namespace std;

namespace BODYMODEL
{
    enum BodyState { undefined, initialized, elementsAdded, elementsComputed, staticPBMAssembled};


    class BodyModel
    {
    private:
        BodyState _state = undefined;
    public:
        // Model attributes
        int nElements; // Number of body elements
        ELEMENT::Element* elements; // Body elements
        int N; // Number of nodes in whole-body model
        PSEUDOBLOCKMATRIX::PseudoBlockMatrix staticPBM; // Constant entries in PBM

        // Basal constants and attributes
        double Mbas0; // W, Whole body basal metabolism
        double actBas; // MET, basal activity fraction of effort
        const double b0; // 1/MET, linear slope of actBas
        const double b1; // MET/MET, y intercept of actBas
        double etaWork; // J/J, work efficiency

        // Skin constants and attributes
        const double StefanBoltzmann = 5.670374419e-8; // W/m^2/K^4
        double Icl; // W/m^2/K, clothing thermal resistance
        double icl; // -, effective clothing resistance
        const double Rinverse; // W*Pa/m^2, vapor resistivity
        const double LambdaH20; // J/kg, latent heat of water
        const double LewisConstant; // K/Pa
        
        // Blood constants and attributes
        const double p; // -, RBC fraction of blood
        double Viv0; // m^3, initial blood plasma volume
        double Vrbc0; // m^3, initial RBC volume
        const double rhoBlood; // kg/m^3, nominal blood density
        const double cpBlood; // J/kg/K, blood constant pressure specific heat capacity

        BodyModel( int tnElements )
        {
            nElements = tnElements;
            _state = initialized;
        }
        ~BodyModel()
        {
            delete [] elements;
        }
        void computeN()
        {
            N = 1;
            for(int i=0;i<nElements;i++){
                elements[i].computeN();
                N += elements[i].N;
            }
        }
        void buildStaticPBM(){

        }
        void addWashers(
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


    };


}

#endif