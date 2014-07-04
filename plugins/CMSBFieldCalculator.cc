#include "PFCal/runPandora/interface/CMSBFieldCalculator.h"

float CMSBFieldCalculator::m_bField = 2.f;

//------------------------------------------------------------------------------------------------------------------------------------------

void CMSBFieldCalculator::InitializeGeometry()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CMSBFieldCalculator::GetBField(const pandora::CartesianVector &/*positionVector*/) const
{
    return m_bField;
}

