#ifndef _BODYMODEL_CPP_
#define _BODYMODEL_CPP_

#include "BodyModel.hpp"
using namespace std;

namespace BODYMODEL
{
    BodyModel* defaultBody(){
        BodyModel* body = new BodyModel(12);
    }
    ELEMENT::Element* head(){
        ELEMENT::Element* head = new ELEMENT::Element();
        // Class attributes
        head->isCylinder = false;
        head->omega = 1;
        head->nWashers = 13;
        head->nSectors = 2;

        // Other attributes
        head->length = -1;
        head->hx = 0;
        head->a_nat = 3.0;
        head->a_frc = 113.0;
        head->a_mix = -5.70;
        head->a_sk = 0.0835;
        head->a_sed = 0.00;
        head->a_stnd = 0;
        head->a_sw = 0.095;
        head->a_dl = 0.055;
        head->a_cs = 0.030;
        head->a_sh = 0.0;
        head->n_j = 0;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        head->addWashers(
            5,   // Number of washers to add
            0,  // m, inner radius
            0.0860,  // m, outer radius
            0.49,   // W/m^2, conductivity
            1080, // kg/m^3, density
            3850,   // J/kg/K, specific heat capacity
            10.1320e-3,// 1/s, tissue permeability
            13400); // W/m^3, specific basal metabolism
        head->addWashers(
            2,   // Number of washers to add
            0.0860,  // m, inner radius
            0.1005,  // m, outer radius
            1.16,   // W/m^2, conductivity
            1500, // kg/m^3, density
            1591,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        head->addWashers(
            2,   // Number of washers to add
            0.1005,  // m, inner radius
            0.1020,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            .0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        head->addWashers(
            4,   // Number of washers to add
            0.1020,  // m, inner radius
            0.1040,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            5.4800e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
    }

}

#endif