#ifndef _BODYMODEL_HPP_
#define _BODYMODEL_HPP_

#include <iostream>
#include <cstddef>
#include <math.h>
#include <assert.h>
#include "element.cpp"
#include "PseudoBlockMatrix.hpp"

using namespace std;

namespace BODYMODEL
{
    /// @brief Enum for state of BodyModel
    enum BodyState { undefined, initialized, allocated, elementsAdded, computed, staticPBMAssembled};

    /// @brief Class for a body.
    class BodyModel
    {
    private:
        /// @brief The state of this bodymodel.
        BodyState _state = undefined;

        /// @brief Computes number of temperature nodes.
        void computeN();

    public:
        /// @brief Model attributes
        int nElements; //!< Number of body elements
        ELEMENT::Element** elements; //!< Body elements ptr array
        int N; //!< Number of nodes in whole-body model
        double *V; //!< m^3, volume of each node (N entries)
        int elemIdx=0; //!< Index of first free element
        PSEUDOBLOCKMATRIX::PseudoBlockMatrix staticPBM; //!< Constant entries in PBM

        // Basal constants and attributes
        double Mbas0; //!< W, Whole body basal metabolism. 58.2 by default.
        const double actBas = 0.8; //!< MET, basal activity fraction of effort
        const double b0 = -0.60; //!< 1/MET, linear slope of actBas
        const double b1 = 0.39; //!< MET/MET, y intercept of actBas
        double etaWork = 0.05; //!< J/J, work efficiency

        // Skin constants and attributes
        double skinSurfaceArea = 0; //!< m^2, total skin surface area
        const double StefanBoltzmann = 5.670374419e-8; //!< W/m^2/K^4
        const double Icl = 0.093; //!< W/m^2/K, clothing thermal resistance
        const double icl = 0.34; //!< -, effective clothing resistance
        const double Rinverse = 0.003; //!< W*Pa/m^2, vapor resistivity
        const double LambdaH20 = 2256e3; //!< J/kg, latent heat of water
        const double LewisConstant = 0.0165; //!< K/Pa
        
        // Blood constants and attributes
        const double p = 0.33; //!< -, RBC fraction of blood
        const double Viv0 = 5.0e-3; //!< m^3, initial blood plasma volume
        const double Vrbc0 = 2.1e-3; //!< m^3, initial RBC volume
        const double rhoBlood = 1069.0; //!< kg/m^3, nominal blood density
        const double cpBlood = 3650.0; //!< J/kg/K, nominal blood constant pressure specific heat capacity

        /// @brief Access protected state object.
        /// @return The body's _state
        BodyState getState();


        /// @brief Add element to the bodymodel.
        /// @param element Element to add
        void addElement( ELEMENT::Element* element );
        
        /// @brief Computes all subparameters of BodyModel.
        void Compute();
        
        /// @brief Blank constructor
        BodyModel()
        {
            _state = initialized;
        }
        /// @brief Constructor for fixed number of elements.
        /// @param tnElements Total number of elements.
        BodyModel( int tnElements ) : BodyModel()
        {
            nElements = tnElements;
            elements = new ELEMENT::Element*[nElements];
            _state = allocated;
        }
        /// @brief Destructor.
        ~BodyModel()
        {
            for(elemIdx=0;elemIdx<nElements;elemIdx++)
                delete elements[elemIdx];
            delete [] elements;
            delete [] V;
        }
    };

    //! Creates pointer to a default body.
    BodyModel* defaultBody();

    //! Default body helper methods
    /*! Add elements to the body. These are the
        default elements described by Fiala et. al.
    */
    //! Creates default element head
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* head();
    //! Creates default element face
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* face();
    //! Creates default element neck
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* neck();
    //! Creates default element shoulders
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* shoulders();
    //! Creates default element thorax
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* thorax();
    //! Creates default element abdomen
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* abdomen();
    //! Creates default element arms
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* arms();
    //! Creates default element hands
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* hands();
    //! Creates default element legs
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* legs();
    //! Creates default element feet
    /*!
    @return Element ptr.
    */
    ELEMENT::Element* feet();


    /// @brief Access protected state object.
    /// @return The body's _state
    BodyState BodyModel::getState()
    {
        return _state;
    }

}

#endif