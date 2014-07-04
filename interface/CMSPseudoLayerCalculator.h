/**
 *  @file   CMSPandora/src/CMSPseudoLayerCalculator.h
 * 
 *  @brief  Header file for the pseudo layer calculator class.
 * 
 *  $Log: $
 */
#ifndef CMS_PSEUDO_LAYER_CALCULATOR_H
#define CMS_PSEUDO_LAYER_CALCULATOR_H 1

#include "Utilities/PseudoLayerCalculator.h"

class CMSPseudoLayerCalculator : public pandora::PseudoLayerCalculator
{
public:
    /**
     *  @brief  Default constructor
     */
    CMSPseudoLayerCalculator();

private:
    void InitializeGeometry();
    pandora::PseudoLayer GetPseudoLayer(const pandora::CartesianVector &positionVector) const;
    pandora::PseudoLayer GetPseudoLayerAtIp() const;

    /**
     *  @brief  Store all revelevant barrel and endcap layer positions upon initialization
     */
    void StoreLayerPositions();

    typedef std::vector<float> LayerPositionList;

    /**
     *  @brief  Store layer positions for a specific subdetector upon initialization
     * 
     *  @param  subDetectorParameters the name of the subdetector
     *  @param  layerParametersList the layer parameters list
     */
    void StoreLayerPositions(const std::string &subDetectorName, LayerPositionList &LayerPositionList);

    /**
     *  @brief  Store layer positions for a specific subdetector upon initialization
     * 
     *  @param  subDetectorParameters the name of the subdetector
     *  @param  layerParametersList the layer parameters list
     */
    void StoreLayerPositions(const pandora::GeometryHelper::SubDetectorParameters &subDetectorParameters, 
			     LayerPositionList &LayerPositionList);

    /**
     *  @brief  Store positions of barrel and endcap outer edges upon initialization
     */
    void StoreDetectorOuterEdge();

    /**
     *  @brief  Find the layer number corresponding to a specified position, via reference to a specified layer position list
     * 
     *  @param  position the specified position
     *  @param  layerPositionList the specified layer position list
     *  @param  layer to receive the layer number
     */
    pandora::StatusCode FindMatchingLayer(const float position, const LayerPositionList &layerPositionList, pandora::PseudoLayer &layer) const;

    LayerPositionList       m_barrelLayerPositions;     ///< List of barrel layer positions
    LayerPositionList       m_endCapLayerPositions;     ///< List of endcap layer positions

    float                   m_barrelInnerEdgeR;         ///< Barrel inner edge r coordinate
    float                   m_barrelOuterEdgeR;         ///< Barrel outer edge r coordinate
    float                   m_endCapInnerEdgeZ;         ///< EndCap inner edge z coordinate
    float                   m_endCapOuterEdgeZ;         ///< Endcap outer edge z coordinate
    float                   m_endcapOuterEdgeR;         ///< Endcap outer edge r coordinate

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::PseudoLayer CMSPseudoLayerCalculator::GetPseudoLayerAtIp() const
{
    static const pandora::PseudoLayer pseudoLayerAtIp(this->GetPseudoLayer(pandora::CartesianVector(0.f, 0.f, 0.f)));
    return pseudoLayerAtIp;
}

#endif // #ifndef CMS_PSEUDO_LAYER_CALCULATOR_H
