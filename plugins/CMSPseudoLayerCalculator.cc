/**
 * 
 *  $Log: $
 */

#include "PFCal/runPandora/interface/CMSPseudoLayerCalculator.h"

#include <algorithm>

CMSPseudoLayerCalculator::CMSPseudoLayerCalculator() :
    PseudoLayerCalculator(),
    m_barrelInnerEdgeR(0.f),
    m_barrelOuterEdgeR(0.f),
    m_endCapInnerEdgeZ(0.f),
    m_endCapOuterEdgeZ(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CMSPseudoLayerCalculator::InitializeGeometry()
{

  // std::cout<<"begin"<<std::endl;
  try
    {
      this->StoreLayerPositions();
      // std::cout << "Trying to store outer edge information" << std::endl ; 
      this->StoreDetectorOuterEdge();  //this function requires more information to work
      // std::cout << "Success!!!" << std::endl ; 
    }
  catch (pandora::StatusCodeException &statusCodeException)
    {
      std::cout << "CMSPseudoLayerCalculator: Failed to obtain geometry information." << std::endl;
      throw statusCodeException;
    }

  if(pandora::GeometryHelper::IsInitialized())
    {
      std::cout<<"geometry initialized"<<std::endl;
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::PseudoLayer CMSPseudoLayerCalculator::GetPseudoLayer(const pandora::CartesianVector &positionVector) const
{
    const float zCoordinate(std::fabs(positionVector.GetZ()));
    const float rCoordinate(std::sqrt(positionVector.GetX() * positionVector.GetX() + positionVector.GetY() * positionVector.GetY()));

    if ((zCoordinate > m_endCapOuterEdgeZ) || (rCoordinate > std::max(m_barrelOuterEdgeR,m_endcapOuterEdgeR))) {
      std::cout << "You die here because z = " << zCoordinate << " (" << m_endCapOuterEdgeZ << ") and r = " << rCoordinate 
		<< " (" << m_barrelOuterEdgeR << "," << m_endcapOuterEdgeR << ")" << std::endl; 
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }

    pandora::PseudoLayer pseudoLayer(0);

    if (zCoordinate < m_endCapInnerEdgeZ)
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->FindMatchingLayer(rCoordinate, m_barrelLayerPositions, pseudoLayer));
    }
    else if (rCoordinate < m_barrelInnerEdgeR)
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->FindMatchingLayer(zCoordinate, m_endCapLayerPositions, pseudoLayer));
    }
    else
    {
        pandora::PseudoLayer bestBarrelLayer(0);
        const pandora::StatusCode barrelStatusCode(this->FindMatchingLayer(rCoordinate, m_barrelLayerPositions, bestBarrelLayer));

        pandora::PseudoLayer bestEndCapLayer(0);
        const pandora::StatusCode endCapStatusCode(this->FindMatchingLayer(zCoordinate, m_endCapLayerPositions, bestEndCapLayer));

        if ((pandora::STATUS_CODE_SUCCESS != barrelStatusCode) && (pandora::STATUS_CODE_SUCCESS != endCapStatusCode))
            throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

        pseudoLayer = std::max(bestBarrelLayer, bestEndCapLayer);
    }

    // Reserve a pseudo layer for track projections, etc.
    return (1 + pseudoLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CMSPseudoLayerCalculator::StoreLayerPositions()
{

  pandora::GeometryHelper::SubDetectorParameters ecalBarrelParameters = pandora::GeometryHelper::GetECalBarrelParameters() ; 
  pandora::GeometryHelper::SubDetectorParameters hcalBarrelParameters = pandora::GeometryHelper::GetHCalBarrelParameters() ; 
  pandora::GeometryHelper::SubDetectorParameters ecalEndCapParameters = pandora::GeometryHelper::GetECalEndCapParameters() ; 
  pandora::GeometryHelper::SubDetectorParameters hcalEndCapParameters = pandora::GeometryHelper::GetHCalEndCapParameters() ; 
  
  // if ( !ecalBarrelParameters.IsInitialized() || !ecalEndCapParameters.IsInitialized() || 
  //      !hcalBarrelParameters.IsInitialized() || !hcalEndCapParameters.IsInitialized() ) {
  //   std::cout << "CMSPseudoLayerCalculator: Failure to initialize (EB:" 
  // 	      << ecalBarrelParameters.IsInitialized() << ", EE:"
  // 	      << ecalEndCapParameters.IsInitialized() << ", HB:"
  // 	      << hcalBarrelParameters.IsInitialized() << ", HE:"
  // 	      << hcalEndCapParameters.IsInitialized() << ")" << std::endl ; 
  //   throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
  // }

  // if ( !ecalBarrelParameters.IsMirroredInZ() || !ecalEndCapParameters.IsMirroredInZ() || 
  //      !hcalBarrelParameters.IsMirroredInZ() || !hcalEndCapParameters.IsMirroredInZ() ) {
  //   std::cout << "CMSPseudoLayerCalculator: Detector not mirrored in Z (EB:" 
  // 	      << ecalBarrelParameters.IsMirroredInZ() << ", EE:"
  // 	      << ecalEndCapParameters.IsMirroredInZ() << ", HB:"
  // 	      << hcalBarrelParameters.IsMirroredInZ() << ", HE:"
  // 	      << hcalEndCapParameters.IsMirroredInZ() << ")" << std::endl ; 
  //   throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  // }

  this->StoreLayerPositions(ecalBarrelParameters, m_barrelLayerPositions); 
  this->StoreLayerPositions(hcalBarrelParameters, m_barrelLayerPositions); 
  this->StoreLayerPositions(ecalEndCapParameters, m_endCapLayerPositions); 
  this->StoreLayerPositions(hcalEndCapParameters, m_endCapLayerPositions); 

  std::cout << "Number of layer positions: " << m_barrelLayerPositions.size() << " (barrel) " << m_endCapLayerPositions.size() << " (endcap)" << std::endl ; 

  if (m_barrelLayerPositions.empty() || m_endCapLayerPositions.empty())  // no barrel part 
    // if (m_barrelLayerPositions.empty())
    {
      std::cout << "CMSPseudoLayerCalculator: No layer positions specified." << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
    }

  std::sort(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end());
  std::sort(m_endCapLayerPositions.begin(), m_endCapLayerPositions.end());

  // for (LayerPositionList::const_iterator iter=m_endCapLayerPositions.begin(); iter!=m_endCapLayerPositions.end(); iter++) { 
  //   std::cout << "Layer value: " << (*iter) << std::endl ; 
  // }

  LayerPositionList::const_iterator barrelIter = std::unique(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end());
  LayerPositionList::const_iterator endcapIter = std::unique(m_endCapLayerPositions.begin(), m_endCapLayerPositions.end());

  if ((m_barrelLayerPositions.end() != barrelIter) || (m_endCapLayerPositions.end() != endcapIter))
    // if ((m_barrelLayerPositions.end() != barrelIter))
    {
      std::cout << "CMSPseudoLayerCalculator: Duplicate layer position detected." << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CMSPseudoLayerCalculator::StoreLayerPositions(const pandora::GeometryHelper::SubDetectorParameters &subDetectorParameters, 
						   LayerPositionList &LayerPositionList)
{
  if (!subDetectorParameters.IsInitialized()) {
    std::cout << "CMSPseudoLayerCalculator: Error, required geometry information not initialized." << std::endl;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
  }

  if (!subDetectorParameters.IsMirroredInZ()) {
    std::cout << "CMSPseudoLayerCalculator: Error, detector must be symmetrical about z=0 plane." << std::endl;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }

  const pandora::GeometryHelper::LayerParametersList &layerParametersList(subDetectorParameters.GetLayerParametersList());

  for (pandora::GeometryHelper::LayerParametersList::const_iterator iter = layerParametersList.begin(), iterEnd = layerParametersList.end();
       iter != iterEnd; ++iter) {
    LayerPositionList.push_back(iter->m_closestDistanceToIp);
  }
}

void CMSPseudoLayerCalculator::StoreLayerPositions(const std::string &subDetectorName, LayerPositionList &LayerPositionList)
{
    static const pandora::GeometryHelper::SubDetectorParametersMap &parametersMap(pandora::GeometryHelper::GetAdditionalSubDetectors());
    const pandora::GeometryHelper::SubDetectorParametersMap::const_iterator iter(parametersMap.find(subDetectorName));
    
    // std::cout << "Size of parameters map: " << parametersMap.size() << std::endl ; 
    // for (pandora::GeometryHelper::SubDetectorParametersMap::const_iterator i=parametersMap.begin(); i!=parametersMap.end(); i++) { 
    //   std::cout << "Map instance: " << iter->first << std::endl ; 
    // }

    std::cout << "Looking to store positions for " << subDetectorName << std::endl ; 

    if ((parametersMap.end() == iter)  || (!iter->second.IsInitialized()))
    {
      std::cout << "CMSPseudoLayerCalculator: Error, required geometry information not available." << std::endl;
      if ( parametersMap.end() == iter ) std::cout << "CMSPseudoLayerCalculator: Could not find " << subDetectorName << std::endl ; 
      else std::cout << "CMSPseudoLayerCalculator: " << subDetectorName << " not initialized" << std::endl ; 
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }

    if (!iter->second.IsMirroredInZ())
    {
        std::cout << "CMSPseudoLayerCalculator: Error, detector must be symmetrical about z=0 plane." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    const pandora::GeometryHelper::LayerParametersList &layerParametersList(iter->second.GetLayerParametersList());

    for (pandora::GeometryHelper::LayerParametersList::const_iterator iter = layerParametersList.begin(), iterEnd = layerParametersList.end();
        iter != iterEnd; ++iter)
    {
        LayerPositionList.push_back(iter->m_closestDistanceToIp);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
void CMSPseudoLayerCalculator::StoreDetectorOuterEdge()
{
  const pandora::GeometryHelper::SubDetectorParameters ecalBarrelParameters = pandora::GeometryHelper::GetECalBarrelParameters() ; 
  const pandora::GeometryHelper::SubDetectorParameters hcalBarrelParameters = pandora::GeometryHelper::GetHCalBarrelParameters() ; 
  const pandora::GeometryHelper::SubDetectorParameters ecalEndCapParameters = pandora::GeometryHelper::GetECalEndCapParameters() ; 
  const pandora::GeometryHelper::SubDetectorParameters hcalEndCapParameters = pandora::GeometryHelper::GetHCalEndCapParameters() ; 
    
  // if ( (!hcalBarrelParameters.IsInitialized()) || (!hcalEndCapParameters.IsInitialized()) ) { // debugging
  if ( (!ecalBarrelParameters.IsInitialized()) || (!ecalEndCapParameters.IsInitialized()) || 
       (!hcalBarrelParameters.IsInitialized()) || (!hcalEndCapParameters.IsInitialized()) ) { // debugging
  // if ( (!ecalBarrelParameters.IsInitialized()) || (!hcalBarrelParameters.IsInitialized()) || 
  //      (!ecalEndCapParameters.IsInitialized()) ) { 
    std::cout << "CMSPseudoLayerCalculator: Error, required geometry information not available." << std::endl;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
  }

  m_barrelInnerEdgeR = ecalBarrelParameters.GetInnerRCoordinate(); 
  // m_barrelInnerEdgeR = hcalBarrelParameters.GetInnerRCoordinate(); // debugging
  m_barrelOuterEdgeR = hcalBarrelParameters.GetOuterRCoordinate(); 
  // m_barrelOuterEdgeR = ecalBarrelParameters.GetOuterRCoordinate(); // debugging
  m_endcapOuterEdgeR = std::max(ecalEndCapParameters.GetOuterRCoordinate(),hcalEndCapParameters.GetOuterRCoordinate()); 
  // m_endcapOuterEdgeR = ecalEndCapParameters.GetOuterRCoordinate() ; // debugging
  m_endCapInnerEdgeZ = ecalEndCapParameters.GetInnerZCoordinate(); 
  // m_endCapInnerEdgeZ = hcalEndCapParameters.GetInnerZCoordinate(); // debugging
  m_endCapOuterEdgeZ = hcalEndCapParameters.GetOuterZCoordinate(); // debugging
  
  if ((m_barrelLayerPositions.end() != std::upper_bound(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end(), m_barrelOuterEdgeR)) ||
      (m_endCapLayerPositions.end() != std::upper_bound(m_endCapLayerPositions.begin(), m_endCapLayerPositions.end(), m_endCapOuterEdgeZ))) 
  // if ((m_barrelLayerPositions.end() != std::upper_bound(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end(), m_barrelOuterEdgeR)))
    {
      std::cout << "FineGranularityPseudoLayerCalculator: Layers specified outside detector edge." << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }

  m_barrelLayerPositions.push_back(m_barrelOuterEdgeR);
  m_endCapLayerPositions.push_back(m_endCapOuterEdgeZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CMSPseudoLayerCalculator::FindMatchingLayer(const float position, const LayerPositionList &layerPositionList,
    pandora::PseudoLayer &layer) const
{
    LayerPositionList::const_iterator upperIter = std::upper_bound(layerPositionList.begin(), layerPositionList.end(), position);

    if (layerPositionList.end() == upperIter)
    {
        return pandora::STATUS_CODE_NOT_FOUND;
    }

    if (layerPositionList.begin() == upperIter)
    {
        layer = 0;
        return pandora::STATUS_CODE_SUCCESS;
    }

    LayerPositionList::const_iterator lowerIter = upperIter - 1;

    if (std::fabs(position - *lowerIter) < std::fabs(position - *upperIter))
    {
        layer = std::distance(layerPositionList.begin(), lowerIter);
    }
    else
    {
        layer = std::distance(layerPositionList.begin(), upperIter);
    }

    return pandora::STATUS_CODE_SUCCESS;
}
