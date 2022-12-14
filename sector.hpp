#ifndef _SECTOR_
#define _SECTOR_

#include <iostream>
#include <cstddef>
#include <math.h>

namespace SECTOR
{
    /// @brief Angular sector class
    class Sector
    {
    public:
        double phi; //!< rad, interior angle of sector
        double psiSed; //!< -, sedentary view factor
        double psiStnd; //!< -, standard view factor
        double echelon; //!< -, emission coefficient
        double areaSkin; //!< m^2, skin surface area
        double areaShare; //!< -, skin share proportion
    };
}

#endif