#ifndef _ELEMENT_
#define _ELEMENT_

#include <iostream>
#include <cstddef>
#include <math.h>
#include "washer.hpp"
#include "sector.hpp"

using namespace std;

namespace ELEMENT
{
    class Element
    {
    private:

    public:
        bool isCylinder;
        int N; // Number of distinct T nodes
        double length = -1; // m, length of element iff cylinder
        double hx; // W/K, countercurrent heat transfer coefficient
        int n_j; // -, number of equivalent clothing layers
        double omega; // -, 1 for cylinder / 2 for sphere
        double sumPhi; // rad, sum of sector angles
        Sector *sectors; // Array of sectors
        Washer *washers; // Array of washers


        // Derived proportional coefficients in (0,1)
        double a_nat; // -
        double a_frc; // -
        double a_mix; // -
        double a_sk; // -
        double a_sed; // -
        double a_stnd; // -
        double a_sw; // -
        double a_dl; // -
        double a_cs; // -
        double a_sh; // -

    };

}

#endif