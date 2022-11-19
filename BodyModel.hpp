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
    enum BodyState { undefined, initialized, allocated, elementsAdded, computed, staticPBMAssembled};

    class BodyModel
    {
    private:
        BodyState _state = undefined;

        void computeN();

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

        void compute();
        

        BodyModel()
        {
            _state = initialized;
        }
        BodyModel( int tnElements ) : BodyModel()
        {
            nElements = tnElements;
        }
        ~BodyModel()
        {
            delete [] elements;
        }

        

    };

    void BodyModel::computeN()
    {
        assert(_state==computed);

        N = 1;
        for(int i=0;i<nElements;i++){
            N += elements[i].N;
        }
    }

    void BodyModel::compute()
    {
        for(int i=0;i<nElements;i++)
            elements[i].compute();
        
        computeN();
    }



}

#endif