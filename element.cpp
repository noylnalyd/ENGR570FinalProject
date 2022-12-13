#ifndef _ELEMENT_CPP_
#define _ELEMENT_CPP_

#include "element.hpp"

namespace ELEMENT{
    
    void Element::compute( double totalSkinArea ){
        // Assert left DFS'd
        // Compute skin share
        for(int i=0;i<nSectors;i++)
            sectors[i]->areaShare = sectors[i]->areaSkin/totalSkinArea;

        // Compute sub attributes (Right DFS)
        for(int i=0;i<nWashers;i++){
            computeGamma( washers[i] );
            computeVolume( washers[i] );
            washers[i]->compute();
            Vmuscle += washers[i]->volume*washers[i]->muscle*sumPhi/(2*M_PI);
            Vresp += washers[i]->volume*washers[i]->resp*sumPhi/(2*M_PI);
        }
        computeAForwardSkin(washers[nWashers-1]);

        // Compute relational attributes
        for(int i=0;i<nWashers-1;i++){
            computeAForward(washers[i], washers[i+1]);
        }
        for(int i=1;i<nWashers;i++){
            computeABackward(washers[i],washers[i-1]);
        }

        // Compute node attributes
        computeTheta();
        computeN();
        _state = computed;
    };

    void Element::computeN()
    {
        N = 1 + (nWashers-1)*nSectors;
    };
    ElementState Element::getState()
    {
        return _state;
    };
    void Element::addWashers(
        int nAdd,   // Number of washers to add
        double r0,  // m, inner radius
        double rf,  // m, outer radius
        double k,   // W/m^2, conductivity
        double rho, // kg/m^3, density
        double c,   // J/kg/K, specific heat capacity
        double w_bl,// 1/s, tissue permeability
        double q_m, // W/m^3, specific basal metabolism
        double muscle, // -, ratio of muscle (usually 0 or 1)
        double resp) // -, ratio of respiratory exposure
    {
        // Must have been allocated!
        assert(_state==allocated);
        // Must have length configured!
        assert(length>-2);
        // Must have space remaining!
        assert(nAdd+washerIdx<=nWashers);
        // Must have defined omega!
        assert(omega>0);

        double drr = (rf-r0)/nAdd;
        for(int i=0;i<nAdd;i++){
            double r = r0+(i+.5)*drr;
            WASHER::Washer* tmp = new WASHER::Washer();
            tmp->r = r;
            tmp->k = k;
            tmp->rho = rho;
            tmp->c = c;
            tmp->w_bl = w_bl;
            tmp->q_m = q_m;
            tmp->deltaR = drr;
            tmp->muscle = muscle;
            tmp->resp = resp;
            tmp->a_resp = handlearesp(resp,r-.5*drr,r+.5*drr,r0,rf);
            washers[washerIdx++] = tmp;
        }
        if(sectorIdx==nSectors && washerIdx==nWashers)
            _state = wsAdded;
    };
    void Element::addSector(SECTOR::Sector* sec){
        // Must have been allocated!
        assert(_state==allocated);
        // Must have space remaining!
        assert(sectorIdx<nSectors);

        sectors[sectorIdx++] = sec;
        if(sectorIdx==nSectors && washerIdx==nWashers)
            _state = wsAdded;
    };

    void Element::computeGamma( WASHER::Washer* cur )
    {
        cur->gamma = cur->deltaR/cur->r/2*omega;
    }
    void CylinderElement::subCompute( double* totalSkinArea ){
        // Assert all w/s addded
        sumPhi = 0;
        double rskin = washers[nWashers-1]->r+washers[nWashers-1]->deltaR/2.0;
        double coeff = 2*M_PI*rskin*length/(2*M_PI);
        for(int i=0;i<nSectors;i++){
            sumPhi += sectors[i]->phi;
            sectors[i]->areaSkin = coeff*sectors[i]->phi;
            *totalSkinArea += sectors[i]->areaSkin;
        }
    };
    void SphereElement::subCompute( double* totalSkinArea ){
        // Assert all w/s addded
        sumPhi = 0;
        double rskin = washers[nWashers-1]->r+washers[nWashers-1]->deltaR/2.0;
        double coeff = 4*M_PI*rskin*rskin/(2*M_PI);
        for(int i=0;i<nSectors;i++){
            sumPhi += sectors[i]->phi;
            sectors[i]->areaSkin = coeff*sectors[i]->phi;
            *totalSkinArea += sectors[i]->areaSkin;
        }
    };
    void CylinderElement::computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt )
    {
        double rIFC = (cur->r+nxt->r)/2.0;
        double L = cur->k/log((cur->r+cur->deltaR)/cur->r);
        double Lnxt = nxt->k/log(nxt->r/(nxt->r-nxt->deltaR));
        double D = log(rIFC/cur->r)/log((cur->r+cur->deltaR)/cur->r);
        double Ds = log(rIFC/(cur->r+cur->deltaR))/log((cur->r+cur->deltaR)/cur->r);
        double Dnxt = log(rIFC/nxt->r)/log(nxt->r/(nxt->r-nxt->deltaR));

        cur->AForwardCur = (Ds*Lnxt-Dnxt*L)/(D*Lnxt-Dnxt*L);
        cur->AForwardNxt = Lnxt/(D*Lnxt-Dnxt*L);
        assert(!isnan(cur->AForwardCur));
        assert(!isnan(cur->AForwardNxt));
    }
    void CylinderElement::computeAForwardSkin(WASHER::Washer* cur){
        double rIFC = cur->r+cur->deltaR/2;
        double L = cur->k/log((cur->r+cur->deltaR)/cur->r);
        double Lnxt = cur->k/log((cur->r+cur->deltaR)/((cur->r+cur->deltaR)-cur->deltaR));
        double D = log(rIFC/cur->r)/log((cur->r+cur->deltaR)/cur->r);
        double Ds = log(rIFC/(cur->r+cur->deltaR))/log((cur->r+cur->deltaR)/cur->r);
        double Dnxt = log(rIFC/(cur->r+cur->deltaR))/log((cur->r+cur->deltaR)/((cur->r+cur->deltaR)-cur->deltaR));
        cur->AForwardCur = (Ds*Lnxt-Dnxt*L)/(D*Lnxt-Dnxt*L);
        cur->AForwardNxt = Lnxt/(D*Lnxt-Dnxt*L);
        assert(!isnan(cur->AForwardCur));
        assert(!isnan(cur->AForwardNxt));
    }
    
    void CylinderElement::computeABackward( WASHER::Washer* cur, WASHER::Washer* prv )
    {
        double rIFC = (cur->r+prv->r)/2.0;
        double L = cur->k/log((cur->r)/(cur->r-cur->deltaR));
        double Lprv = prv->k/log((prv->r+prv->deltaR)/prv->r);
        double D = log(rIFC/cur->r)/log(cur->r/(cur->r-cur->deltaR));
        double Ds = log(rIFC/(cur->r-cur->deltaR))/log(cur->r/(cur->r-cur->deltaR));
        double Dprv = log(rIFC/prv->r)/log((prv->r+prv->deltaR)/prv->r);

        cur->ABackwardCur = (Dprv*L-Ds*Lprv)/(Dprv*L-D*Lprv);
        cur->ABackwardPrv = Lprv/(Dprv*L-D*Lprv);
        assert(!isnan(cur->ABackwardCur));
        assert(!isnan(cur->ABackwardPrv));
    }
    double CylinderElement::handlearesp(double resp, double ricur, double rocur, double ri, double ro)
    {
        return resp*(rocur*rocur-ricur*ricur)/(ro*ro-ri*ri);
    }
    void CylinderElement::computeVolume( WASHER::Washer* cur ){
        double r_i = cur->r-cur->deltaR/2.0;
        double r_o = cur->r+cur->deltaR/2.0;
        cur->volume = length*M_PI*(r_o*r_o-r_i*r_i);
    }
    void CylinderElement::computeTheta(){
        theta = 8.0*(1.0-1.0/sqrt(2.0))*(1.0-1.0/sqrt(2.0))/(log(2.0*sqrt(2.0)-1.0)*sumPhi);
    }
    void SphereElement::computeAForward( WASHER::Washer* cur, WASHER::Washer* nxt )
    {
        double L = cur->k/(1/cur->r-1/(cur->r+cur->deltaR));
        double Lnxt = nxt->k/(1/(nxt->r-nxt->deltaR)-1/nxt->r);

        cur->AForwardCur = (L-Lnxt*cur->r/(nxt->r-nxt->deltaR))/
                (L+Lnxt*((cur->r+cur->deltaR)/(nxt->r-nxt->deltaR)));
        cur->AForwardNxt = (Lnxt*(1+nxt->r/(nxt->r-nxt->deltaR)))/
                (L+Lnxt*((cur->r+cur->deltaR)/(nxt->r-nxt->deltaR)));
        assert(!isnan(cur->AForwardCur));
        assert(!isnan(cur->AForwardNxt));
    }
    void SphereElement::computeAForwardSkin( WASHER::Washer* cur)
    {
        double L = cur->k/(1/cur->r-1/(cur->r+cur->deltaR));
        double Lnxt = cur->k/(1/((cur->r+cur->deltaR)-cur->deltaR)-1/(cur->r+cur->deltaR));

        cur->AForwardCur = (L-Lnxt*cur->r/((cur->r+cur->deltaR)-cur->deltaR))/
                (L+Lnxt*((cur->r+cur->deltaR)/((cur->r+cur->deltaR)-cur->deltaR)));
        cur->AForwardNxt = (Lnxt*(1+(cur->r+cur->deltaR)/((cur->r+cur->deltaR)-cur->deltaR)))/
                (L+Lnxt*((cur->r+cur->deltaR)/((cur->r+cur->deltaR)-cur->deltaR)));
        assert(!isnan(cur->AForwardCur));
        assert(!isnan(cur->AForwardNxt));
    }
    void SphereElement::computeABackward( WASHER::Washer* cur, WASHER::Washer* prv )
    {
        double L = cur->k/(1/(cur->r-cur->deltaR)-1/cur->r);
        double Lprv = prv->k/(1/prv->r-1/(prv->r+prv->deltaR));

        cur->ABackwardCur = (L-Lprv*(cur->r/(prv->r+prv->deltaR)))/(L+Lprv*((cur->r-cur->deltaR)/(prv->r+prv->deltaR)));
        cur->ABackwardPrv = (Lprv*(1+prv->r/(prv->r+prv->deltaR)))/(L+Lprv*((cur->r-cur->deltaR)/(prv->r+prv->deltaR)));
        assert(!isnan(cur->ABackwardCur));
        assert(!isnan(cur->ABackwardPrv));
    }
    double SphereElement::handlearesp(double resp, double ricur, double rocur, double ri, double ro){
        return resp*(rocur*rocur*rocur-ricur*ricur*ricur)/(ro*ro*ro-ri*ri*ri);
    }
    void SphereElement::computeVolume( WASHER::Washer* cur ){
        double r_i = cur->r-cur->deltaR/2.0;
        double r_o = cur->r+cur->deltaR/2.0;
        cur->volume = 4.0*M_PI/3.0*(r_o*r_o*r_o-r_i*r_i*r_i);
    }
    void SphereElement::computeTheta(){
        theta = 12*(cbrt(16)-1)*(1-1/cbrt(2))*(1-1/cbrt(2))/((cbrt(32)-cbrt(16))*sumPhi);
    }
}

#endif