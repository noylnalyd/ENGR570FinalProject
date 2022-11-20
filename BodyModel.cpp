#ifndef _BODYMODEL_CPP_
#define _BODYMODEL_CPP_

#include "BodyModel.hpp"
using namespace std;

namespace BODYMODEL
{
    BodyModel* defaultBody(){
        BodyModel* body = new BodyModel(10);
        
        // Check allocation status
        assert(body->getState() == BODYMODEL::allocated);

        // Add elements
        body->addElement(head());
        body->addElement(face());
        body->addElement(neck());
        body->addElement(shoulders());
        body->addElement(thorax());
        body->addElement(abdomen());
        body->addElement(arms());
        body->addElement(hands());
        body->addElement(legs());
        body->addElement(feet());

        // Check status now
        assert(body->getState() == BODYMODEL::elementsAdded);

        // Compute body parameters
        body->compute();

        // Check status again
        assert(body->getState() == BODYMODEL::computed);

        
    }
    
    void controlSignals(
        double* Sh,
        double* Cs,
        double* Dl,
        double* Sw,
        double TskError,
        double ThyError,
        double TskErrorGradient,
        double bvr,
        double svr)
        {
            double TskErrorDot = 0;
            if(TskError <=0 && TskErrorGradient <=0){
                TskErrorDot = TskError*TskErrorGradient;
            }
            *Sh = 10*(tanh(.48*TskError+3.62)-1)*TskError
                -27.9*ThyError
                +1.7*TskErrorDot
                -28.6; // W
            *Cs = 35*(tanh(.34*TskError+1.07)-1)*TskError
                +3.9*TskErrorDot; // -
            *Dl = 21*(tanh(.79*TskError-.70)+1)*TskError
                +32*(tanh(3.29*ThyError-1.46)+1)*ThyError; // W/K
            *Sw = (.8*tanh(.59*TskError-.19)+1.2)*TskError
                +(5.7*tanh(1.98*ThyError-1.03)+6.3)*ThyError; // g/min
            *Sh = min(350.0*(svr+bvr),max(0.0,*Sh)); // W      
            *Sw = min(30.0,max(*Sw,0.0)); // g/min
        }

    ELEMENT::Element* head(){
        ELEMENT::Element* head = new ELEMENT::SphereElement(13,2);

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
        
        SECTOR::Sector *forehead = new SECTOR::Sector();
        forehead->phi = 10*M_PI/180.0; // rad, interior angle of sector
        forehead->psiSed = 1.0; // -, sedentary view factor
        forehead->psiStnd = 1.0; // -, standard view factor
        forehead->echelon = 0.99; // -, emission coefficient
        head->addSector(forehead);
        SECTOR::Sector *headBack = new SECTOR::Sector();
        headBack->phi = 170*M_PI/180.0; // rad, interior angle of sector
        headBack->psiSed = 0.90; // -, sedentary view factor
        headBack->psiStnd = 0.90; // -, standard view factor
        headBack->echelon = 0.80; // -, emission coefficient
        head->addSector(headBack);

        // Make sure totally filled!
        assert(head->getState()==ELEMENT::wsAdded);
        return head;
    }
    ELEMENT::Element* face(){
        ELEMENT::Element* face = new ELEMENT::CylinderElement(7,1);

        // Other attributes
        face->length = 0.0984;
        face->hx = 0;
        face->a_nat = 3.0;
        face->a_frc = 113.0;
        face->a_mix = -5.70;
        face->a_sk = 0.0418;
        face->a_sed = 0.00;
        face->a_stnd = 0;
        face->a_sw = 0.054;
        face->a_dl = 0.046;
        face->a_cs = 0.033;
        face->a_sh = 0.002;
        face->n_j = 0;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        face->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0268,  // m, outer radius
            0.420,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/K, specific heat capacity
            0.538e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        face->addWashers(
            1,   // Number of washers to add
            0.0268,  // m, inner radius
            0.0542,  // m, outer radius
            1.16,   // W/m^2, conductivity
            1500, // kg/m^3, density
            1591,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        face->addWashers(
            1,   // Number of washers to add
            0.0542,  // m, inner radius
            0.0680,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/K, specific heat capacity
            0.538e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        face->addWashers(
            2,   // Number of washers to add
            0.0680,  // m, inner radius
            0.0760,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            588); // W/m^3, specific basal metabolism
        face->addWashers(
            2,   // Number of washers to add
            0.0760,  // m, inner radius
            0.0780,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            11.17e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *faceFront = new SECTOR::Sector();
        faceFront->phi = 210*M_PI/180.0; // rad, interior angle of sector
        faceFront->psiSed = 0.9; // -, sedentary view factor
        faceFront->psiStnd = 0.9; // -, standard view factor
        faceFront->echelon = 0.99; // -, emission coefficient
        face->addSector(faceFront);

        // Make sure totally filled!
        assert(face->getState()==ELEMENT::wsAdded);
        return face;
    }
    ELEMENT::Element* neck(){
        ELEMENT::Element* neck = new ELEMENT::CylinderElement(11,2);

        // Other attributes
        neck->length = 0.0842;
        neck->hx = 0;
        neck->a_nat = 1.6;
        neck->a_frc = 130.0;
        neck->a_mix = -6.50;
        neck->a_sk = 0.0417;
        neck->a_sed = 0.03;
        neck->a_stnd = 0.01;
        neck->a_sw = 0.042;
        neck->a_dl = 0.031;
        neck->a_cs = 0.025;
        neck->a_sh = 0.002;
        neck->n_j = 0;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        neck->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0190,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        neck->addWashers(
            4,   // Number of washers to add
            0.0190,  // m, inner radius
            0.0546,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/K, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        neck->addWashers(
            2,   // Number of washers to add
            0.0546,  // m, inner radius
            0.0556,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        neck->addWashers(
            4,   // Number of washers to add
            0.0556,  // m, inner radius
            0.0567,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            6.8e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *neckAnterior = new SECTOR::Sector();
        neckAnterior->phi = 180*M_PI/180.0; // rad, interior angle of sector
        neckAnterior->psiSed = 0.7; // -, sedentary view factor
        neckAnterior->psiStnd = 0.7; // -, standard view factor
        neckAnterior->echelon = 0.99; // -, emission coefficient
        neck->addSector(neckAnterior);
        SECTOR::Sector *neckPosterior = new SECTOR::Sector();
        neckPosterior->phi = 180*M_PI/180.0; // rad, interior angle of sector
        neckPosterior->psiSed = 0.75; // -, sedentary view factor
        neckPosterior->psiStnd = 0.75; // -, standard view factor
        neckPosterior->echelon = 0.99; // -, emission coefficient
        neck->addSector(neckPosterior);

        // Make sure totally filled!
        assert(neck->getState()==ELEMENT::wsAdded);
        return neck;
    }
    ELEMENT::Element* shoulders(){
        ELEMENT::Element* shoulders = new ELEMENT::CylinderElement(7,1);

        // Other attributes
        shoulders->length = 0.3200;
        shoulders->hx = 0.80;
        shoulders->a_nat = 5.9;
        shoulders->a_frc = 216.0;
        shoulders->a_mix = -10.80;
        shoulders->a_sk = 0.03;
        shoulders->a_sed = 0.05;
        shoulders->a_stnd = 0.02;
        shoulders->a_sw = 0.037;
        shoulders->a_dl = 0.020;
        shoulders->a_cs = 0.010;
        shoulders->a_sh = 0.0002;
        shoulders->n_j = 1;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        shoulders->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0370,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        shoulders->addWashers(
            2,   // Number of washers to add
            0.0370,  // m, inner radius
            0.0390,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/K, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        shoulders->addWashers(
            2,   // Number of washers to add
            0.0390,  // m, inner radius
            0.0440,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        shoulders->addWashers(
            2,   // Number of washers to add
            0.0440,  // m, inner radius
            0.0460,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            1.01e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *shouldersTop = new SECTOR::Sector();
        shouldersTop->phi = 130*M_PI/180.0; // rad, interior angle of sector
        shouldersTop->psiSed = 0.90; // -, sedentary view factor
        shouldersTop->psiStnd = 0.90; // -, standard view factor
        shouldersTop->echelon = 0.99; // -, emission coefficient
        shoulders->addSector(shouldersTop);

        // Make sure totally filled!
        assert(shoulders->getState()==ELEMENT::wsAdded);
        return shoulders;
    }
    ELEMENT::Element* thorax(){
        ELEMENT::Element* thorax = new ELEMENT::CylinderElement(19,3);

        // Other attributes
        thorax->length = 0.3060;
        thorax->hx = 0;
        thorax->a_nat = 0.5;
        thorax->a_frc = 180.0;
        thorax->a_mix = -7.40;
        thorax->a_sk = 0.1290;
        thorax->a_sed = 0.12;
        thorax->a_stnd = 0.07;
        thorax->a_sw = 0.101;
        thorax->a_dl = 0.141;
        thorax->a_cs = 0.0005;
        thorax->a_sh = 0.6305;
        thorax->n_j = 1;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        thorax->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0773,  // m, outer radius
            0.28,   // W/m^2, conductivity
            550, // kg/m^3, density
            3718,   // J/kg/K, specific heat capacity
            (4.9e-3)/60/5.744209713,// 1/s, tissue permeability. Note: Respiratory!
            600); // W/m^3, specific basal metabolism
        thorax->addWashers(
            3,   // Number of washers to add
            0.0773,  // m, inner radius
            0.0891,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        thorax->addWashers(
            3,   // Number of washers to add
            0.0891,  // m, inner radius
            0.1234,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/K, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        thorax->addWashers(
            6,   // Number of washers to add
            0.1234,  // m, inner radius
            0.1268,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        thorax->addWashers(
            6,   // Number of washers to add
            0.1268,  // m, inner radius
            0.1290,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            1.58e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *thoraxAnterior = new SECTOR::Sector();
        thoraxAnterior->phi = 150*M_PI/180.0; // rad, interior angle of sector
        thoraxAnterior->psiSed = 0.80; // -, sedentary view factor
        thoraxAnterior->psiStnd = 0.90; // -, standard view factor
        thoraxAnterior->echelon = 0.99; // -, emission coefficient
        thorax->addSector(thoraxAnterior);
        SECTOR::Sector *thoraxPosterior = new SECTOR::Sector();
        thoraxPosterior->phi = 150*M_PI/180.0; // rad, interior angle of sector
        thoraxPosterior->psiSed = 0.95; // -, sedentary view factor
        thoraxPosterior->psiStnd = 0.95; // -, standard view factor
        thoraxPosterior->echelon = 0.99; // -, emission coefficient
        thorax->addSector(thoraxPosterior);
        SECTOR::Sector *thoraxInferior = new SECTOR::Sector();
        thoraxInferior->phi = 60*M_PI/180.0; // rad, interior angle of sector
        thoraxInferior->psiSed = 0.05; // -, sedentary view factor
        thoraxInferior->psiStnd = 0.10; // -, standard view factor
        thoraxInferior->echelon = 0.99; // -, emission coefficient
        thorax->addSector(thoraxInferior);
        

        // Make sure totally filled!
        assert(thorax->getState()==ELEMENT::wsAdded);
        return thorax;
    }
    ELEMENT::Element* abdomen(){
        ELEMENT::Element* abdomen = new ELEMENT::CylinderElement(19,3);

        // Other attributes
        abdomen->length = 0.5520;
        abdomen->hx = 0;
        abdomen->a_nat = 0.5;
        abdomen->a_frc = 180.0;
        abdomen->a_mix = -9.0;
        abdomen->a_sk = 0.1210;
        abdomen->a_sed = 0.46;
        abdomen->a_stnd = 0.20;
        abdomen->a_sw = 0.181;
        abdomen->a_dl = 0.161;
        abdomen->a_cs = 0.0205;
        abdomen->a_sh = 0.24;
        abdomen->n_j = 1;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        abdomen->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0785,  // m, outer radius
            0.53,   // W/m^2, conductivity
            1000, // kg/m^3, density
            3697,   // J/kg/K, specific heat capacity
            4.31e-3,// 1/s, tissue permeability. Note: Respiratory!
            4100); // W/m^3, specific basal metabolism
        abdomen->addWashers(
            3,   // Number of washers to add
            0.0785,  // m, inner radius
            0.0834,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        abdomen->addWashers(
            3,   // Number of washers to add
            0.0834,  // m, inner radius
            0.1090,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/PosteriorK, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        abdomen->addWashers(
            6,   // Number of washers to add
            0.1090,  // m, inner radius
            0.1244,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        abdomen->addWashers(
            6,   // Number of washers to add
            0.1244,  // m, inner radius
            0.1260,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            1.44e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *abdomenAnterior = new SECTOR::Sector();
        abdomenAnterior->phi = 150*M_PI/180.0; // rad, interior angle of sector
        abdomenAnterior->psiSed = 0.80; // -, sedentary view factor
        abdomenAnterior->psiStnd = 0.90; // -, standard view factor
        abdomenAnterior->echelon = 0.99; // -, emission coefficient
        abdomen->addSector(abdomenAnterior);
        SECTOR::Sector *abdomenPosterior = new SECTOR::Sector();
        abdomenPosterior->phi = 150*M_PI/180.0; // rad, interior angle of sector
        abdomenPosterior->psiSed = 0.95; // -, sedentary view factor
        abdomenPosterior->psiStnd = 0.95; // -, standard view factor
        abdomenPosterior->echelon = 0.99; // -, emission coefficient
        abdomen->addSector(abdomenPosterior);
        SECTOR::Sector *abdomenInferior = new SECTOR::Sector();
        abdomenInferior->phi = 60*M_PI/180.0; // rad, interior angle of sector
        abdomenInferior->psiSed = 0.20; // -, sedentary view factor
        abdomenInferior->psiStnd = 0.30; // -, standard view factor
        abdomenInferior->echelon = 0.99; // -, emission coefficient
        abdomen->addSector(abdomenInferior);
        

        // Make sure totally filled!
        assert(abdomen->getState()==ELEMENT::wsAdded);
        return abdomen;
    }
    ELEMENT::Element* arms(){
        ELEMENT::Element* arms = new ELEMENT::CylinderElement(16,3);

        // Other attributes
        arms->length = 1.274;
        arms->hx = 4.13;
        arms->a_nat = 8.3;
        arms->a_frc = 216.0;
        arms->a_mix = -10.8;
        arms->a_sk = 0.1800;
        arms->a_sed = 0.19;
        arms->a_stnd = 0.08;
        arms->a_sw = 0.133;
        arms->a_dl = 0.095;
        arms->a_cs = 0.1945;
        arms->a_sh = 0.04;
        arms->n_j = 0;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        arms->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0153,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        arms->addWashers(
            3,   // Number of washers to add
            0.0153,  // m, inner radius
            0.0343,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/PosteriorK, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        arms->addWashers(
            6,   // Number of washers to add
            0.0343,  // m, inner radius
            0.0401,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        arms->addWashers(
            6,   // Number of washers to add
            0.0401,  // m, inner radius
            0.0418,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            1.1e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *armsAnterior = new SECTOR::Sector();
        armsAnterior->phi = 135*M_PI/180.0; // rad, interior angle of sector
        armsAnterior->psiSed = 0.75; // -, sedentary view factor
        armsAnterior->psiStnd = 0.85; // -, standard view factor
        armsAnterior->echelon = 0.99; // -, emission coefficient
        arms->addSector(armsAnterior);
        SECTOR::Sector *armsPosterior = new SECTOR::Sector();
        armsPosterior->phi = 135*M_PI/180.0; // rad, interior angle of sector
        armsPosterior->psiSed = 0.80; // -, sedentary view factor
        armsPosterior->psiStnd = 0.90; // -, standard view factor
        armsPosterior->echelon = 0.99; // -, emission coefficient
        arms->addSector(armsPosterior);
        SECTOR::Sector *armsInferior = new SECTOR::Sector();
        armsInferior->phi = 90*M_PI/180.0; // rad, interior angle of sector
        armsInferior->psiSed = 0.10; // -, sedentary view factor
        armsInferior->psiStnd = 0.20; // -, standard view factor
        armsInferior->echelon = 0.99; // -, emission coefficient
        arms->addSector(armsInferior);
        

        // Make sure totally filled!
        assert(arms->getState()==ELEMENT::wsAdded);
        return arms;
    }
    ELEMENT::Element* hands(){
        ELEMENT::Element* hands = new ELEMENT::CylinderElement(9,2);

        // Other attributes
        hands->length = .62;
        hands->hx = .57;
        hands->a_nat = 8.3;
        hands->a_frc = 216.0;
        hands->a_mix = -10.8;
        hands->a_sk = 0.0900;
        hands->a_sed = 0.02;
        hands->a_stnd = 0.01;
        hands->a_sw = 0.049;
        hands->a_dl = 0.121;
        hands->a_cs = 0.1100;
        hands->a_sh = 0.002;
        hands->n_j = 0;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        hands->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0070,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        hands->addWashers(
            2,   // Number of washers to add
            0.0070,  // m, inner radius
            0.0174,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/PosteriorK, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        hands->addWashers(
            2,   // Number of washers to add
            0.0174,  // m, inner radius
            0.0204,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        hands->addWashers(
            4,   // Number of washers to add
            0.0204,  // m, inner radius
            0.0226,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            4.54e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *handsHandback = new SECTOR::Sector();
        handsHandback->phi = 180*M_PI/180.0; // rad, interior angle of sector
        handsHandback->psiSed = 0.80; // -, sedentary view factor
        handsHandback->psiStnd = 0.80; // -, standard view factor
        handsHandback->echelon = 0.99; // -, emission coefficient
        hands->addSector(handsHandback);
        SECTOR::Sector *handsPalm = new SECTOR::Sector();
        handsPalm->phi = 180*M_PI/180.0; // rad, interior angle of sector
        handsPalm->psiSed = 0.10; // -, sedentary view factor
        handsPalm->psiStnd = 0.20; // -, standard view factor
        handsPalm->echelon = 0.99; // -, emission coefficient
        hands->addSector(handsPalm);

        // Make sure totally filled!
        assert(hands->getState()==ELEMENT::wsAdded);
        return hands;
    }
    ELEMENT::Element* legs(){
        ELEMENT::Element* legs = new ELEMENT::CylinderElement(19,3);

        // Other attributes
        legs->length = 1.39;
        legs->hx = 6.90;
        legs->a_nat = 5.3;
        legs->a_frc = 220.0;
        legs->a_mix = -11;
        legs->a_sk = 0.2080;
        legs->a_sed = 0.11;
        legs->a_stnd = 0.60;
        legs->a_sw = 0.261;
        legs->a_dl = 0.230;
        legs->a_cs = 0.2;
        legs->a_sh = 0.0813;
        legs->n_j = 1;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        legs->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0220,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        legs->addWashers(
            6,   // Number of washers to add
            0.0220,  // m, inner radius
            0.0480,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/PosteriorK, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        legs->addWashers(
            6,   // Number of washers to add
            0.0480,  // m, inner radius
            0.0533,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        legs->addWashers(
            6,   // Number of washers to add
            0.0533,  // m, inner radius
            0.0553,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            1.05e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *legsAnterior = new SECTOR::Sector();
        legsAnterior->phi = 150*M_PI/180.0; // rad, interior angle of sector
        legsAnterior->psiSed = 0.85; // -, sedentary view factor
        legsAnterior->psiStnd = 0.90; // -, standard view factor
        legsAnterior->echelon = 0.99; // -, emission coefficient
        legs->addSector(legsAnterior);
        SECTOR::Sector *legsPosterior = new SECTOR::Sector();
        legsPosterior->phi = 150*M_PI/180.0; // rad, interior angle of sector
        legsPosterior->psiSed = 0.95; // -, sedentary view factor
        legsPosterior->psiStnd = 0.90; // -, standard view factor
        legsPosterior->echelon = 0.99; // -, emission coefficient
        legs->addSector(legsPosterior);
        SECTOR::Sector *legsInferior = new SECTOR::Sector();
        legsInferior->phi = 60*M_PI/180.0; // rad, interior angle of sector
        legsInferior->psiSed = 0.10; // -, sedentary view factor
        legsInferior->psiStnd = 0.65; // -, standard view factor
        legsInferior->echelon = 0.99; // -, emission coefficient
        legs->addSector(legsInferior);

        // Make sure totally filled!
        assert(legs->getState()==ELEMENT::wsAdded);
        return legs;
    }
    ELEMENT::Element* feet(){
        ELEMENT::Element* feet = new ELEMENT::CylinderElement(11,2);

        // Other attributes
        feet->length = 1.39;
        feet->hx = 6.90;
        feet->a_nat = 5.3;
        feet->a_frc = 220.0;
        feet->a_mix = -11.0;
        feet->a_sk = 0.0750;
        feet->a_sed = 0.1100;
        feet->a_stnd = 0.60;
        feet->a_sw = 0.047;
        feet->a_dl = 0.100;
        feet->a_cs = 0.3765;
        feet->a_sh = 0.002;
        feet->n_j = 3;
         // addWashers(obj,N,r0,rf,k,rho,c,w_bl,q_m,index)
        feet->addWashers(
            1,   // Number of washers to add
            0,  // m, inner radius
            0.0200,  // m, outer radius
            0.75,   // W/m^2, conductivity
            1357, // kg/m^3, density
            1700,   // J/kg/K, specific heat capacity
            0e-3,// 1/s, tissue permeability
            0); // W/m^3, specific basal metabolism
        feet->addWashers(
            2,   // Number of washers to add
            0.0200,  // m, inner radius
            0.0250,  // m, outer radius
            0.42,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3768,   // J/kg/PosteriorK, specific heat capacity
            0.5380e-3,// 1/s, tissue permeability
            684); // W/m^3, specific basal metabolism
        feet->addWashers(
            4,   // Number of washers to add
            0.0250,  // m, inner radius
            0.0326,  // m, outer radius
            0.16,   // W/m^2, conductivity
            850, // kg/m^3, density
            2300,   // J/kg/K, specific heat capacity
            0.0036e-3,// 1/s, tissue permeability
            58); // W/m^3, specific basal metabolism
        feet->addWashers(
            4,   // Number of washers to add
            0.0326,  // m, inner radius
            0.0350,  // m, outer radius
            0.47,   // W/m^2, conductivity
            1085, // kg/m^3, density
            3680,   // J/kg/K, specific heat capacity
            1.5e-3,// 1/s, tissue permeability
            368); // W/m^3, specific basal metabolism
        
        SECTOR::Sector *feetInstep = new SECTOR::Sector();
        feetInstep->phi = 180*M_PI/180.0; // rad, interior angle of sector
        feetInstep->psiSed = 0.90; // -, sedentary view factor
        feetInstep->psiStnd = 0.90; // -, standard view factor
        feetInstep->echelon = 0.99; // -, emission coefficient
        feet->addSector(feetInstep);
        SECTOR::Sector *feetSole = new SECTOR::Sector();
        feetSole->phi = 180*M_PI/180.0; // rad, interior angle of sector
        feetSole->psiSed = 1.0; // -, sedentary view factor
        feetSole->psiStnd = 1.0; // -, standard view factor
        feetSole->echelon = 0.99; // -, emission coefficient
        feet->addSector(feetSole);

        // Make sure totally filled!
        assert(feet->getState()==ELEMENT::wsAdded);
        return feet;
    }
}

#endif