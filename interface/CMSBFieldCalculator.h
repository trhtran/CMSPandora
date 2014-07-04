#ifndef CMS_BFIELD_CALCULATOR_H
#define CMS_BFIELD_CALCULATOR_H 1

#include "Utilities/BFieldCalculator.h"

class CMSBFieldCalculator : public pandora::BFieldCalculator
{
public:
    static float        m_bField;                ///< The bfield, units Tesla

private:
    void InitializeGeometry();
    float GetBField(const pandora::CartesianVector &positionVector) const;
};

#endif 
