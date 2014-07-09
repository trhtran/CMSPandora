#include "PFCal/runPandora/interface/runPandora.h"
#include "PFCal/runPandora/interface/CMSBFieldCalculator.h"
#include "PFCal/runPandora/interface/CMSPseudoLayerCalculator.h"
#include "FineGranularityContent.h"
#include "PFCal/runPandora/interface/CMSTemplateAlgorithm.h"
#include "PandoraMonitoringApi.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

// Addition for HGC geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/FlatTrd.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//We need the speed of light
#include "CLHEP/Units/PhysicalConstants.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/RootAutoLibraryLoader/interface/RootAutoLibraryLoader.h"

#include "TGClient.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TRint.h"
 
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace edm;
using namespace reco;
  

pandora::Pandora * runPandora::m_pPandora = NULL;
//
// constructors and destructor
//
runPandora::runPandora(const edm::ParameterSet& iConfig) 
{
  //now do what ever initialization is needed
  m_pPandora = new pandora::Pandora();
 
  inputTagEcalRecHitsEB_ = iConfig.getParameter<InputTag>("ecalRecHitsEB");
  inputTagHcalRecHitsHBHE_ = iConfig.getParameter<InputTag>("hcalRecHitsHBHE");
  inputTagHGCEErechit_ = iConfig.getParameter<InputTag>("HGCEErechitCollection");
  inputTagHGCHEFrechit_ = iConfig.getParameter<InputTag>("HGCHEFrechitCollection");
  inputTagHGCHEBrechit_ = iConfig.getParameter<InputTag>("HGCHEBrechitCollection");
  inputTagGeneralTracks_ = iConfig.getParameter< std::vector < InputTag > >("generaltracks");
  inputTagtPRecoTrackAsssociation_ = iConfig.getParameter<InputTag>("tPRecoTrackAsssociation");
  inputTagGenParticles_ = iConfig.getParameter<InputTag>("genParticles");
  m_pandoraSettingsXmlFile = iConfig.getParameter<std::string>("inputconfigfile");

  m_calibrationParameterFile = iConfig.getParameter<std::string>("calibrParFile");
  _outputFileName = iConfig.getParameter<std::string>("outputFile");

  
// NS // SHOWER PROFILE CALCULATOR
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetShowerProfileCalculator(*m_pPandora,new FineGranularityShowerProfileCalculator()));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, FineGranularityContent::RegisterAlgorithms(*m_pPandora));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, FineGranularityContent::RegisterHelperFunctions(*m_pPandora));
  
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldCalculator(*m_pPandora, new CMSBFieldCalculator()));	
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerCalculator(*m_pPandora, new CMSPseudoLayerCalculator()));	
  
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora, "Template", new CMSTemplateAlgorithm::Factory));

  // prepareGeometry(iSetup);
}

runPandora::~runPandora()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void runPandora::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "Analyzing events" << std::endl ; 
  if ( firstEvent_ ) { 
    firstEvent_ = false ; 
    std::cout << "At the first event...preparing geometry" << std::endl ; 
    prepareGeometry(iSetup) ; 
    std::cout << "Done with Geometry setup...moving along" << std::endl ; 
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_pandoraSettingsXmlFile));

  }

  // std::cout << "Analyzing events 1 " << std::endl ;

  // Get the primary vertex
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
  reco::Vertex pv = pvHandle->at(0) ; 
  
  // std::cout << "Analyzing events 2 " << std::endl ;
  
  // get the Calorimeter rechit collections
  edm::Handle<EcalRecHitCollection> ecalRecHitHandleEB;
  edm::Handle<HBHERecHitCollection> hcalRecHitHandleHBHE;
  edm::Handle<HGCRecHitCollection> HGCeeRecHitHandle;
  edm::Handle<HGCRecHitCollection> HGChefRecHitHandle;
  edm::Handle<HGCRecHitCollection> HGChebRecHitHandle;
  
  // std::cout << iEvent.getByLabel(inputTagHGCEErechit_, HGCeeRecHitHandle) << " " 
  // 	    << iEvent.getByLabel(inputTagHGCHEFrechit_, HGChefRecHitHandle) << " " 
  // 	    << iEvent.getByLabel(inputTagHGCHEBrechit_, HGChebRecHitHandle) << std::endl ; 

  bool found = iEvent.getByLabel(inputTagEcalRecHitsEB_, ecalRecHitHandleEB) && 
    iEvent.getByLabel(inputTagHcalRecHitsHBHE_, hcalRecHitHandleHBHE) && 
    iEvent.getByLabel(inputTagHGCEErechit_, HGCeeRecHitHandle) && 
    iEvent.getByLabel(inputTagHGCHEFrechit_, HGChefRecHitHandle) && 
    iEvent.getByLabel(inputTagHGCHEBrechit_, HGChebRecHitHandle);

  edm::Handle<reco::RecoToSimCollection > rectosimCollection;
  iEvent.getByLabel(inputTagtPRecoTrackAsssociation_, rectosimCollection);
  const reco::RecoToSimCollection pRecoToSim = *(rectosimCollection.product());
    
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  //Why PF uses global point (0,0,0) for all events?
  math::XYZVector B_(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));

  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart);
  
  if(!found ) {
    std::ostringstream err;
    err<<"cannot find rechits: "<< HGCeeRecHitHandle.isValid() << "," << HGChefRecHitHandle.isValid() << "," << HGChebRecHitHandle.isValid() ;
    LogError("runPandora")<<err.str()<<std::endl;
    throw cms::Exception( "MissingProduct", err.str());
  } 

  prepareTrack(B_,pRecoToSim,iEvent,iSetup);
  prepareHits( ecalRecHitHandleEB,hcalRecHitHandleHBHE,HGCeeRecHitHandle,HGChefRecHitHandle,HGChebRecHitHandle,pv,iEvent,iSetup );
  preparemcParticle(genpart);
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,PandoraApi::ProcessEvent(*m_pPandora));
  preparePFO(iEvent,iSetup);
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,PandoraApi::Reset(*m_pPandora));

}

void runPandora::initPandoraCalibrParameters()
{

   m_Calibr_ADC2GeV_EE  = 0.000012 ; //w/o absorber thickness correction
   m_Calibr_ADC2GeV_HEF = 0.0000176; //w/o absorber thickness correction
   m_Calibr_ADC2GeV_HEB = 0.0003108; //w/o absorber thickness correction

   m_addCalibrEE = 10.;
   m_addCalibrHEF= 10.;
   m_addCalibrHEB= 10.;

   m_eCalThresBarrel    = 0.;
   m_eCalThresEndCap    = 27.55e-6; //EE
   m_hCalThresEndCapHEF = 42.50e-6;
   m_hCalThresEndCapHEB = 742.2e-6;
   m_hCalThresBarrel    = 0.;


   m_eCalMipThresBarrel    = 0.5;
   m_eCalMipThresEndCap    = 0.5;
   m_hCalMipThresEndCapHEF = 0.5;
   m_hCalMipThresEndCapHEB = 0.5;
   m_hCalMipThresBarrel    = 0.5;

   m_eCalToMipBarrel       = 3.3333333;
   m_eCalToMipEndCap       = 18149.;
   m_hCalToMipEndCapHEF    = 11765.;
   m_hCalToMipEndCapHEB    = 667.4;
   m_hCalToMipBarrel       = 3.3333333;

   m_eCalToEMGeVEndCap     = 1.;
   m_eCalToEMGeVBarrel     = 1.;
   m_hCalToEMGeVEndCapHEF  = 1.;
   m_hCalToEMGeVEndCapHEB  = 1.;
   m_hCalToEMGeVBarrel     = 1.;

   m_eCalToHadGeVEndCap    = 1.;
   m_eCalToHadGeVBarrel    = 1.;
   m_hCalToHadGeVEndCapHEF = 1.;
   m_hCalToHadGeVEndCapHEB = 1.;
   m_hCalToHadGeVBarrel    = 1.;
   m_muonToMip             = 1.;

   return;
}

void runPandora::readCalibrParameterFile()
{

   std::ifstream calibrParFile(m_calibrationParameterFile , std::ifstream::in );

   if (!calibrParFile.is_open()) {
      std::cout << "runPandora::readCalibrParameterFile: calibrParFile does not exist ("
         << m_calibrationParameterFile << ")" << std::endl;
      return;
   }

   while ( !calibrParFile.eof() ) {
      std::string linebuf;
      getline( calibrParFile, linebuf );
      if (linebuf.substr(0,1) == "#") continue;
      if (linebuf.substr(0,2) == "//") continue;

      if (linebuf.empty()) continue;

      std::string paraName;
      float paraValue;
      std::stringstream ss(linebuf);
      ss >> paraName >> paraValue;
      std::cout << "reading calibr parameter " << paraName << " ";

      if (paraName=="Calibr_ADC2GeV_EE"     ) {m_Calibr_ADC2GeV_EE       = paraValue; std::cout <<  m_Calibr_ADC2GeV_EE     << std::endl;}
      if (paraName=="Calibr_ADC2GeV_HEF"    ) {m_Calibr_ADC2GeV_HEF      = paraValue; std::cout <<  m_Calibr_ADC2GeV_HEF    << std::endl;}
      if (paraName=="Calibr_ADC2GeV_HEB"    ) {m_Calibr_ADC2GeV_HEB      = paraValue; std::cout <<  m_Calibr_ADC2GeV_HEB    << std::endl;}

      if (paraName=="addCalibrEE"    ) {m_addCalibrEE = paraValue; std::cout <<  m_addCalibrEE << std::endl;}
      if (paraName=="addCalibrHEF"   ) {m_addCalibrHEF= paraValue; std::cout <<  m_addCalibrHEF<< std::endl;}
      if (paraName=="addCalibrHEB"   ) {m_addCalibrHEB= paraValue; std::cout <<  m_addCalibrHEB<< std::endl;}

      if (paraName=="ECalThresBarrel"       ) {m_eCalThresBarrel         = paraValue; std::cout <<  m_eCalThresBarrel       << std::endl;}
      if (paraName=="ECalThresEndCap"       ) {m_eCalThresEndCap         = paraValue; std::cout <<  m_eCalThresEndCap       << std::endl;}
      if (paraName=="HCalThresEndCapHEF"    ) {m_hCalThresEndCapHEF      = paraValue; std::cout <<  m_hCalThresEndCapHEF    << std::endl;}
      if (paraName=="HCalThresEndCapHEB"    ) {m_hCalThresEndCapHEB      = paraValue; std::cout <<  m_hCalThresEndCapHEB    << std::endl;}
      if (paraName=="HCalThresBarrel"       ) {m_hCalThresBarrel         = paraValue; std::cout <<  m_hCalThresBarrel       << std::endl;}


      if (paraName=="ECalMipThresEndCap"    ) {m_eCalMipThresEndCap      = paraValue; std::cout <<  m_eCalMipThresEndCap    << std::endl;}
      if (paraName=="ECalMipThresBarrel"    ) {m_eCalMipThresBarrel      = paraValue; std::cout <<  m_eCalMipThresBarrel    << std::endl;}
      if (paraName=="HCalMipThresEndCapHEF" ) {m_hCalMipThresEndCapHEF   = paraValue; std::cout <<  m_hCalMipThresEndCapHEF << std::endl;}
      if (paraName=="HCalMipThresEndCapHEB" ) {m_hCalMipThresEndCapHEB   = paraValue; std::cout <<  m_hCalMipThresEndCapHEB << std::endl;}
      if (paraName=="HCalMipThresBarrel"    ) {m_hCalMipThresBarrel      = paraValue; std::cout <<  m_hCalMipThresBarrel    << std::endl;}

      if (paraName=="ECalToMipEndCap"       ) {m_eCalToMipEndCap         = paraValue; std::cout <<  m_eCalToMipEndCap       << std::endl;}
      if (paraName=="ECalToMipBarrel"       ) {m_eCalToMipBarrel         = paraValue; std::cout <<  m_eCalToMipBarrel       << std::endl;}
      if (paraName=="HCalToMipEndCapHEF"    ) {m_hCalToMipEndCapHEF      = paraValue; std::cout <<  m_hCalToMipEndCapHEF    << std::endl;}
      if (paraName=="HCalToMipEndCapHEB"    ) {m_hCalToMipEndCapHEB      = paraValue; std::cout <<  m_hCalToMipEndCapHEB    << std::endl;}
      if (paraName=="HCalToMipBarrel"       ) {m_hCalToMipBarrel         = paraValue; std::cout <<  m_hCalToMipBarrel       << std::endl;}

      if (paraName=="ECalToEMGeVEndCap"     ) {m_eCalToEMGeVEndCap       = paraValue; std::cout <<  m_eCalToEMGeVEndCap     << std::endl;}
      if (paraName=="ECalToEMGeVBarrel"     ) {m_eCalToEMGeVBarrel       = paraValue; std::cout <<  m_eCalToEMGeVBarrel     << std::endl;}
      if (paraName=="HCalToEMGeVEndCapHEF"  ) {m_hCalToEMGeVEndCapHEF    = paraValue; std::cout <<  m_hCalToEMGeVEndCapHEF  << std::endl;}
      if (paraName=="HCalToEMGeVEndCapHEB"  ) {m_hCalToEMGeVEndCapHEB    = paraValue; std::cout <<  m_hCalToEMGeVEndCapHEB  << std::endl;}
      if (paraName=="HCalToEMGeVBarrel"     ) {m_hCalToEMGeVBarrel       = paraValue; std::cout <<  m_hCalToEMGeVBarrel     << std::endl;}

      if (paraName=="ECalToHadGeVEndCap"    ) {m_eCalToHadGeVEndCap      = paraValue; std::cout <<  m_eCalToHadGeVEndCap    << std::endl;}
      if (paraName=="ECalToHadGeVBarrel"    ) {m_eCalToHadGeVBarrel      = paraValue; std::cout <<  m_eCalToHadGeVBarrel    << std::endl;}
      if (paraName=="HCalToHadGeVEndCapHEF" ) {m_hCalToHadGeVEndCapHEF   = paraValue; std::cout <<  m_hCalToHadGeVEndCapHEF << std::endl;}
      if (paraName=="HCalToHadGeVEndCapHEB" ) {m_hCalToHadGeVEndCapHEB   = paraValue; std::cout <<  m_hCalToHadGeVEndCapHEB << std::endl;}
      if (paraName=="HCalToHadGeVBarrel"    ) {m_hCalToHadGeVBarrel      = paraValue; std::cout <<  m_hCalToHadGeVBarrel    << std::endl;}

      if (paraName=="MuonToMip"             ) {m_muonToMip               = paraValue; std::cout <<  m_muonToMip             << std::endl;}
   }

   calibrParFile.close();

   return;
}



void runPandora::prepareGeometry(const edm::EventSetup& iSetup){ // function to setup a geometry for pandora

  // std::cout << "I am preparing my geometry!!!" << std::endl ; 

  PandoraApi::Geometry::Parameters geometryParameters;

  // 
  // Add ECAL/HCAL parameters to geometry
  //
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
  // Get the ecal/hcal barrel, endcap geometry
  const CaloSubdetectorGeometry *ebtmp = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  // const CaloSubdetectorGeometry *eetmp = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
  const CaloSubdetectorGeometry *hbtmp = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  // const CaloSubdetectorGeometry *hetmp = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);

  // Additions
  edm::ESHandle<HGCalGeometry> hgceeGeoHandle ; 
  edm::ESHandle<HGCalGeometry> hgchefGeoHandle ; 
  edm::ESHandle<HGCalGeometry> hgchebGeoHandle ; 
  
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle) ; 
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle) ; 
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle) ; 

  const HGCalGeometry &hgceetmp = *hgceeGeoHandle ; 
  const HGCalGeometry &hgcheftmp = *hgchefGeoHandle ; 
  const HGCalGeometry &hgchebtmp = *hgchebGeoHandle ; 

  // std::cout << "Basic check: " << hgceetmp.producerTag() << " " << hgcheftmp.producerTag() << " " << hgchebtmp.producerTag() << std::endl ; 
  
  std::vector<DetId> ecalBarrelCells = geoHandle->getValidDetIds(DetId::Ecal, EcalBarrel);
  std::vector<DetId> ecalEndcapCells = hgceeGeoHandle->getValidDetIds(DetId::Forward, HGCEE);
  // std::vector<DetId> ecalEndcapCells = geoHandle->getValidDetIds(DetId::Ecal, EcalEndcap);
  std::vector<DetId> hcalBarrelCells = geoHandle->getValidDetIds(DetId::Hcal, HcalBarrel);
  std::vector<DetId> hcalEndcapCellsFront = hgchefGeoHandle->getValidDetIds(DetId::Forward, HGCHEF);
  std::vector<DetId> hcalEndcapCellsBack  = hgchebGeoHandle->getValidDetIds(DetId::Forward, HGCHEB);
  
  const EcalBarrelGeometry* ecalBarrelGeometry = dynamic_cast< const EcalBarrelGeometry* > (ebtmp);
  // const EcalEndcapGeometry* ecalEndcapGeometry = dynamic_cast< const EcalEndcapGeometry* > (eetmp);
  const HcalGeometry* hcalBarrelGeometry = dynamic_cast< const HcalGeometry* > (hbtmp);
  // const HcalGeometry* hcalEndcapGeometry = dynamic_cast< const HcalGeometry* > (hetmp);
  
  assert( ecalBarrelGeometry );
  // assert( ecalEndcapGeometry );
  assert( hcalBarrelGeometry );
  // assert( hcalEndcapGeometry );
  
  const HGCalGeometry &HGCEEGeometry  = dynamic_cast< const HGCalGeometry& > (hgceetmp);
  const HGCalGeometry &HGCHEFGeometry = dynamic_cast< const HGCalGeometry& > (hgcheftmp);
  const HGCalGeometry &HGCHEBGeometry = dynamic_cast< const HGCalGeometry& > (hgchebtmp);
  assert( &HGCEEGeometry );
  assert( &HGCHEFGeometry );
  assert( &HGCHEBGeometry );

  PandoraApi::GeometryParameters::SubDetectorParameters *ebParameters;
  PandoraApi::GeometryParameters::SubDetectorParameters *eeParameters;
  PandoraApi::GeometryParameters::SubDetectorParameters *hbParameters;
  PandoraApi::GeometryParameters::SubDetectorParameters *heParameters;

  // Additions
  // PandoraApi::GeometryParameters::SubDetectorParameters *hgceeParameters;
  // PandoraApi::GeometryParameters::SubDetectorParameters *hgchefParameters;
  // PandoraApi::GeometryParameters::SubDetectorParameters *hgchebParameters;

  ebParameters = new PandoraApi::GeometryParameters::SubDetectorParameters();
  eeParameters = new PandoraApi::GeometryParameters::SubDetectorParameters();
  hbParameters = new PandoraApi::GeometryParameters::SubDetectorParameters();
  heParameters = new PandoraApi::GeometryParameters::SubDetectorParameters();

  // hgceeParameters  = new PandoraApi::GeometryParameters::SubDetectorParameters();
  // hgchefParameters = new PandoraApi::GeometryParameters::SubDetectorParameters();
  // hgchebParameters = new PandoraApi::GeometryParameters::SubDetectorParameters();

  PandoraApi::Geometry::Parameters::LayerParameters *ebLayerParameters;
  // PandoraApi::Geometry::Parameters::LayerParameters *eeLayerParameters;
  PandoraApi::Geometry::Parameters::LayerParameters *hbLayerParameters;
  // PandoraApi::Geometry::Parameters::LayerParameters *heLayerParameters;
  ebLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
  // eeLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
  hbLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
  // heLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();

  std::vector<PandoraApi::Geometry::Parameters::LayerParameters*> hgcEELayerParameters;
  std::vector<PandoraApi::Geometry::Parameters::LayerParameters*> hgcHEFLayerParameters;
  std::vector<PandoraApi::Geometry::Parameters::LayerParameters*> hgcHEBLayerParameters;
  // std::vector<PandoraApi::Geometry::Parameters::LayerParameters*> hgcHELayerParameters;

  unsigned int nHGCeeLayers = 32, nHGChefLayers = 32, nHGChebLayers = 21 ; 
  std::vector<double> min_innerR_depth_ee, min_innerZ_depth_ee ; 
  std::vector<double> min_innerR_depth_hef, min_innerZ_depth_hef ; 
  std::vector<double> min_innerR_depth_heb, min_innerZ_depth_heb ; 
  min_innerR_depth_ee.clear() ; min_innerZ_depth_ee.clear() ; 
  for (unsigned int i=0; i<nHGCeeLayers; i++) { 
    PandoraApi::Geometry::Parameters::LayerParameters *eeLayerParameters;
    eeLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
    hgcEELayerParameters.push_back( eeLayerParameters ) ; 
    min_innerR_depth_ee.push_back( 99999.0 ) ; 
    min_innerZ_depth_ee.push_back( 99999.0 ) ; 

    PandoraApi::Geometry::Parameters::LayerParameters *hefLayerParameters;
    hefLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
    hgcHEFLayerParameters.push_back( hefLayerParameters ) ; 
    min_innerR_depth_hef.push_back( 99999.0 ) ; 
    min_innerZ_depth_hef.push_back( 99999.0 ) ; 

    if ( i < nHGChebLayers ) { 
      PandoraApi::Geometry::Parameters::LayerParameters *hebLayerParameters;
      hebLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
      hgcHEBLayerParameters.push_back( hebLayerParameters ) ; 
      min_innerR_depth_heb.push_back( 99999.0 ) ; 
      min_innerZ_depth_heb.push_back( 99999.0 ) ; 
    }
  }

  // To be enabled
  // PandoraApi::Geometry::Parameters::LayerParameters *hgceeLayerParameters;
  // PandoraApi::Geometry::Parameters::LayerParameters *hgchefLayerParameters;
  // PandoraApi::Geometry::Parameters::LayerParameters *hgchebLayerParameters;

  // hgceeLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
  // hgchefLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
  // hgchebLayerParameters = new PandoraApi::Geometry::Parameters::LayerParameters();
  
  // Phi Coordinate is when start drawing the detector, wrt x-axis.  
  // Assuming this is 0 since CMS ranges from -pi to pi
  ebParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  ebParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  eeParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  eeParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  hbParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  hbParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  heParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  heParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 

  // hgceeParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  // hgceeParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  // hgchefParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  // hgchefParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  // hgchebParameters->m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  // hgchebParameters->m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 

  // Symmetry order is how you draw the "polygon" detector.  
  // Circle approximation for now (0), but can be configured to match N(cells)
  ebParameters->m_innerSymmetryOrder = 0 ; 
  ebParameters->m_outerSymmetryOrder = 0 ; 
  eeParameters->m_innerSymmetryOrder = 0 ; 
  eeParameters->m_outerSymmetryOrder = 0 ; 
  hbParameters->m_innerSymmetryOrder = 0 ; 
  hbParameters->m_outerSymmetryOrder = 0 ; 
  heParameters->m_innerSymmetryOrder = 0 ; 
  heParameters->m_outerSymmetryOrder = 0 ; 

  // hgceeParameters->m_innerSymmetryOrder = 0 ; 
  // hgceeParameters->m_outerSymmetryOrder = 0 ; 
  // hgchefParameters->m_innerSymmetryOrder = 0 ; 
  // hgchefParameters->m_outerSymmetryOrder = 0 ; 
  // hgchebParameters->m_innerSymmetryOrder = 0 ; 
  // hgchebParameters->m_outerSymmetryOrder = 0 ; 

  // Determine: inner/outer detector radius
  double min_innerRadius = 99999.0 ; double max_outerRadius = 0.0 ;
  double min_innerZ = 99999.0 ; double max_outerZ = 0.0 ; 
  for (std::vector<DetId>::const_iterator ib=ecalBarrelCells.begin(); ib!=ecalBarrelCells.end(); ib++) {
    const CaloSubdetectorGeometry* geom = ecalBarrelGeometry;
    const CaloCellGeometry *thisCell = geom->getGeometry(*ib);
  
    // Inner radius taken as average magnitude of (x,y) for corners 0-3
    // Outer radius taken as average magnitude of (x,y) for corners 4-7
  
    // Inner Z taken as maximum z coordinate for corners 0-3.  Is this right?
    // Outer Z taken as maximum z coordinate for corners 4-7.  Is this right?
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
  
    double avgX_inner = 0.25 * (corners[0].x() + corners[1].x() + corners[2].x() + corners[3].x()) ;
    double avgY_inner = 0.25 * (corners[0].y() + corners[1].y() + corners[2].y() + corners[3].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
    double avgX_outer = 0.25 * (corners[4].x() + corners[5].x() + corners[6].x() + corners[7].x()) ;
    double avgY_outer = 0.25 * (corners[4].y() + corners[5].y() + corners[6].y() + corners[7].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
    for( unsigned int isubcell = 0; isubcell<8; isubcell++){
        if( fabs(corners[isubcell].z()) < min_innerZ ) min_innerZ = fabs(corners[isubcell].z()) ;
        if ( corners[isubcell].z() > max_outerZ ) max_outerZ = corners[isubcell].z();
    }
  }

  ebParameters->m_innerRCoordinate = min_innerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  ebParameters->m_outerRCoordinate = max_outerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  ebParameters->m_innerZCoordinate = min_innerZ * 10.0 ; // CMS units cm, Pandora expects mm
  ebParameters->m_outerZCoordinate = max_outerZ * 10.0 ; // CMS units cm, Pandora expects mm
  ebParameters->m_isMirroredInZ = true ; // Duplicate detector +/- z
  //ebParameters->m_nLayers = 1 ; // One ECAL layer
  ebParameters->m_nLayers = 1 ; // One ECAL layer

  ebLayerParameters->m_closestDistanceToIp = min_innerRadius * 10.0 ; 
  ebLayerParameters->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  ebLayerParameters->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  ebParameters->m_layerParametersList.push_back(*ebLayerParameters) ; 

  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;

  for (std::vector<DetId>::const_iterator ie=ecalEndcapCells.begin(); ie!=ecalEndcapCells.end(); ie++) {

    const CaloCellGeometry *thisCell = HGCEEGeometry.getGeometry(*ie);
    
    // HGCEEDetId detid = HGCEEDetId(*ie) ; 
    unsigned int layer = (unsigned int) ((HGCEEDetId)(*ie)).layer() ; 

    // Inner radius taken as average magnitude of (x,y) for corners 0,3,4,7 
    // Outer radius taken as average magnitude of (x,y) for corners 1,2,5,6
  
    // Inner Z taken as minimum z coordinate for corners 0-3.  Is this right?
    // Outer Z taken as maximum z coordinate for corners 4-7.  Is this right?
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
  
    double avgX_inner = 0.25 * (corners[0].x() + corners[3].x() + corners[4].x() + corners[7].x()) ;
    double avgY_inner = 0.25 * (corners[0].y() + corners[3].y() + corners[4].y() + corners[7].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
    if ( innerRadius < min_innerR_depth_ee.at(layer) ) min_innerR_depth_ee.at(layer) = innerRadius ; 
    
    // if ( layer == 0 ) std::cout << "R: " << innerRadius << " Z: " ; 
    // else std::cout << layer << std::endl ; 

    for( unsigned int isubcell = 0; isubcell<8; isubcell++){
      // if ( layer == 0 ) std::cout << fabs(corners[isubcell].z()) << " " ; 
        if ( fabs(corners[isubcell].z()) < min_innerZ ) min_innerZ = fabs(corners[isubcell].z()) ;
        if ( fabs(corners[isubcell].z()) < min_innerZ_depth_ee.at(layer) ) 
	  min_innerZ_depth_ee.at(layer) = fabs(corners[isubcell].z()) ;
        if ( corners[isubcell].z() > max_outerZ ) max_outerZ = corners[isubcell].z();
    }
    // if ( layer == 0 ) std::cout << std::endl ;

    double avgX_outer = 0.25 * (corners[1].x() + corners[2].x() + corners[5].x() + corners[6].x()) ;
    double avgY_outer = 0.25 * (corners[1].y() + corners[2].y() + corners[5].y() + corners[6].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
  }
    
  eeParameters->m_innerRCoordinate = min_innerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  eeParameters->m_outerRCoordinate = max_outerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  eeParameters->m_innerZCoordinate = min_innerZ * 10.0 ; // CMS units cm, Pandora expects mm
  eeParameters->m_outerZCoordinate = max_outerZ * 10.0 ; // CMS units cm, Pandora expects mm
  eeParameters->m_isMirroredInZ = true ; // Duplicate detector +/- z
  // eeParameters->m_nLayers = nHGCeeLayers ; // HACK, see below
  
//  std::cout << "Parameters for HGC EE geometry (in mm): " << min_innerRadius * 10.0 << " to " << max_outerRadius * 10.0
//	    << " in R and " << min_innerZ * 10.0 << " to " << max_outerZ * 10.0 << " in Z" << std::endl ; 

  int nLayers = 0 ; 
  for (unsigned int i=0; i<nHGCeeLayers; i++) { 
    double distToIP = 10.0 * sqrt(min_innerR_depth_ee.at(i)*min_innerR_depth_ee.at(i) + min_innerZ_depth_ee.at(i)*min_innerZ_depth_ee.at(i)) ; 
    hgcEELayerParameters.at(i)->m_closestDistanceToIp = distToIP ; 
    hgcEELayerParameters.at(i)->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    hgcEELayerParameters.at(i)->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    if ( distToIP < 10000.0 ) { 
      nLayers++ ; 
      // std::cout << "Storing this value: " << distToIP << " for HGC EE layer " << i << std::endl ;
      eeParameters->m_layerParametersList.push_back(*(hgcEELayerParameters.at(i))) ; 
    }
  }
  std::cout << "EE layers: " << nLayers << std::endl ; 
  eeParameters->m_nLayers = nLayers ; // HACK(?) to account for the invalid layer 0(???)

  // eeLayerParameters->m_closestDistanceToIp = min_innerRadius * 10.0 ; 
  // eeLayerParameters->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  // eeLayerParameters->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  // eeParameters->m_layerParametersList.push_back(*eeLayerParameters) ; 

  // 
  // HCAL
  // 
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  for (std::vector<DetId>::const_iterator ib=hcalBarrelCells.begin(); ib!=hcalBarrelCells.end(); ib++) {
    const CaloSubdetectorGeometry* geom = hcalBarrelGeometry;
    const CaloCellGeometry *thisCell = geom->getGeometry(*ib);
    
    // Inner radius taken as average magnitude of (x,y) for corners 0-3 
    // Outer radius taken as average magnitude of (x,y) for corners 4-7
   
    // Inner Z taken as maximum z coordinate for corners 0-3.  Is this right?
    // Outer Z taken as maximum z coordinate for corners 4-7.  Is this right?
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
  
    double avgX_inner = 0.25 * (corners[0].x() + corners[1].x() + corners[2].x() + corners[3].x()) ;
    double avgY_inner = 0.25 * (corners[0].y() + corners[1].y() + corners[2].y() + corners[3].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
  
    for( unsigned int isubcell = 0; isubcell<8; isubcell++){
        if( fabs(corners[isubcell].z()) < min_innerZ ) min_innerZ = fabs(corners[isubcell].z()) ;
        if ( corners[isubcell].z() > max_outerZ ) max_outerZ = corners[isubcell].z();
    }
    
    double avgX_outer = 0.25 * (corners[4].x() + corners[5].x() + corners[6].x() + corners[7].x()) ;
    double avgY_outer = 0.25 * (corners[4].y() + corners[5].y() + corners[6].y() + corners[7].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
    
  }
    
  hbParameters->m_innerRCoordinate = min_innerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  hbParameters->m_outerRCoordinate = max_outerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  hbParameters->m_innerZCoordinate = min_innerZ * 10.0 ; // CMS units cm, Pandora expects mm
  hbParameters->m_outerZCoordinate = max_outerZ * 10.0 ; // CMS units cm, Pandora expects mm
  hbParameters->m_isMirroredInZ = true ; // Duplicate detector +/- z
  hbParameters->m_nLayers = 1 ; // One ECAL layer --> To be fixed

  hbLayerParameters->m_closestDistanceToIp = min_innerRadius * 10.0 ; 
  hbLayerParameters->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  hbLayerParameters->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  hbParameters->m_layerParametersList.push_back(*hbLayerParameters) ; 

  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  nLayers = 0 ; 
  for (std::vector<DetId>::const_iterator ie=hcalEndcapCellsFront.begin(); ie!=hcalEndcapCellsFront.end(); ie++) {

    const CaloCellGeometry *thisCell = HGCHEFGeometry.getGeometry(*ie);
    unsigned int layer = (unsigned int) ((HGCHEDetId)(*ie)).layer() ; 
    
    // Inner radius taken as average magnitude of (x,y) for corners 0,3,4,7 
    // Outer radius taken as average magnitude of (x,y) for corners 1,2,5,6
  
    // Inner Z taken as minimum z coordinate for corners 0-3.  Is this right?
    // Outer Z taken as maximum z coordinate for corners 4-7.  Is this right?
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
  
    double avgX_inner = 0.25 * (corners[0].x() + corners[3].x() + corners[4].x() + corners[7].x()) ;
    double avgY_inner = 0.25 * (corners[0].y() + corners[3].y() + corners[4].y() + corners[7].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
    if ( innerRadius < min_innerR_depth_hef.at(layer) ) min_innerR_depth_hef.at(layer) = innerRadius ;   
    for( unsigned int isubcell = 0; isubcell<8; isubcell++) {
      if( fabs(corners[isubcell].z()) < min_innerZ ) min_innerZ = fabs(corners[isubcell].z()) ;
      if ( fabs(corners[isubcell].z()) < min_innerZ_depth_hef.at(layer) ) 
	min_innerZ_depth_hef.at(layer) = fabs(corners[isubcell].z()) ;
    }

    double avgX_outer = 0.25 * (corners[1].x() + corners[2].x() + corners[5].x() + corners[6].x()) ;
    double avgY_outer = 0.25 * (corners[1].y() + corners[2].y() + corners[5].y() + corners[6].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
  }

  for (std::vector<DetId>::const_iterator ie=hcalEndcapCellsBack.begin(); ie!=hcalEndcapCellsBack.end(); ie++) {

    const CaloCellGeometry *thisCell = HGCHEBGeometry.getGeometry(*ie);
    unsigned int layer = (unsigned int) ((HGCHEDetId)(*ie)).layer() ; 
    
    // Inner radius taken as average magnitude of (x,y) for corners 0,3,4,7 
    // Outer radius taken as average magnitude of (x,y) for corners 1,2,5,6  
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
  
    double avgX_inner = 0.25 * (corners[0].x() + corners[3].x() + corners[4].x() + corners[7].x()) ;
    double avgY_inner = 0.25 * (corners[0].y() + corners[3].y() + corners[4].y() + corners[7].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
    if ( innerRadius < min_innerR_depth_heb.at(layer) ) min_innerR_depth_heb.at(layer) = innerRadius ;   
    for( unsigned int isubcell = 0; isubcell<8; isubcell++) {
      if ( corners[isubcell].z() > max_outerZ ) max_outerZ = corners[isubcell].z();
      if ( fabs(corners[isubcell].z()) < min_innerZ_depth_heb.at(layer) ) 
	min_innerZ_depth_heb.at(layer) = fabs(corners[isubcell].z()) ;
    }

    double avgX_outer = 0.25 * (corners[1].x() + corners[2].x() + corners[5].x() + corners[6].x()) ;
    double avgY_outer = 0.25 * (corners[1].y() + corners[2].y() + corners[5].y() + corners[6].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
  }
  
  std::cout << "Parameters for HGC HCAL F/B geometry (in mm): " 
	    << min_innerRadius * 10.0 << " to " << max_outerRadius * 10.0
	    << " in R and " << min_innerZ * 10.0 << " to " << max_outerZ * 10.0 << " in Z" << std::endl ; 

  heParameters->m_innerRCoordinate = min_innerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  heParameters->m_outerRCoordinate = max_outerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  heParameters->m_innerZCoordinate = min_innerZ * 10.0 ; // CMS units cm, Pandora expects mm
  heParameters->m_outerZCoordinate = max_outerZ * 10.0 ; // CMS units cm, Pandora expects mm
  heParameters->m_isMirroredInZ = true ; // Duplicate detector +/- z
  // heParameters->m_nLayers = nHGChefLayers + nHGChebLayers - 1 ; // HACK...see below

  nLayers = 0 ; 
  for (unsigned int i=0; i<nHGChefLayers; i++) { 
    double distToIP = 10.0 * sqrt(min_innerR_depth_hef.at(i)*min_innerR_depth_hef.at(i) + min_innerZ_depth_hef.at(i)*min_innerZ_depth_hef.at(i)) ; 
    hgcHEFLayerParameters.at(i)->m_closestDistanceToIp = distToIP ; 
    hgcHEFLayerParameters.at(i)->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    hgcHEFLayerParameters.at(i)->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    if ( distToIP < 10000.0 ) { 
      nLayers++ ; 
      // std::cout << "Storing this value: " << distToIP << " for HGC HEF layer " << i << std::endl ;
      heParameters->m_layerParametersList.push_back(*(hgcHEFLayerParameters.at(i))) ; 
    }
  }

  for (unsigned int i=0; i<nHGChebLayers; i++) { 
    double distToIP = 10.0 * sqrt(min_innerR_depth_heb.at(i)*min_innerR_depth_heb.at(i) + min_innerZ_depth_heb.at(i)*min_innerZ_depth_heb.at(i)) ; 
    hgcHEBLayerParameters.at(i)->m_closestDistanceToIp = distToIP ; 
    hgcHEBLayerParameters.at(i)->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    hgcHEBLayerParameters.at(i)->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    if ( distToIP < 10000.0 ) { 
      nLayers++ ; 
      // std::cout << "Storing this value: " << distToIP << " for HGC HEB layer " << i << std::endl ;
      heParameters->m_layerParametersList.push_back(*(hgcHEBLayerParameters.at(i))) ; 
    }
  }

  std::cout << "HEF+B layers: " << nLayers << std::endl ; 
  heParameters->m_nLayers = nLayers ; // HACK

  // heLayerParameters->m_closestDistanceToIp = min_innerRadius * 10.0 ; 
  // heLayerParameters->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  // heLayerParameters->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  // heParameters->m_layerParametersList.push_back(*heLayerParameters) ; 
  
  geometryParameters.m_eCalBarrelParameters = *ebParameters ;
  geometryParameters.m_eCalEndCapParameters = *eeParameters ;
  geometryParameters.m_hCalBarrelParameters = *hbParameters ;
  geometryParameters.m_hCalEndCapParameters = *heParameters ;

  // std::cout << "before set GEO" << std::endl;
  // std::cout << "Idle check: " << geometryParameters.m_InnerRCoordinate << std::endl ; 
  PandoraApi::Geometry::Create(*m_pPandora, geometryParameters);
  // std::cout << "after set GEO" << std::endl;

}

void runPandora::prepareTrack( math::XYZVector B_, const reco::RecoToSimCollection pRecoToSim, const edm::Event& iEvent, const edm::EventSetup& iSetup){ // function to setup tracks in an event for pandora
  PandoraApi::Track::Parameters trackParameters;
  //We need the speed of light
  double speedoflight = (CLHEP::c_light*CLHEP::mm)/CLHEP::ns;
  std::cout<< speedoflight << " mm/ns" << std::endl;

 std::cout << "prepareTrack 1 " << inputTagGeneralTracks_.size() << std::endl ;

  for (unsigned int istr=0; istr<inputTagGeneralTracks_.size();istr++){
    
 std::cout << "prepareTrack 2 " << std::endl;

    //Track collection
    // edm::Handle<reco::TrackCollection> tkRefCollection;
    edm::Handle<edm::View<reco::Track> > tkRefCollection;
    bool found1 = iEvent.getByLabel(inputTagGeneralTracks_[istr], tkRefCollection);
    if(!found1 ) {
      std::ostringstream err;
      err<<"cannot find generalTracks: "<< inputTagGeneralTracks_[istr];
      LogError("runPandora")<<err.str()<<std::endl;
      throw cms::Exception( "MissingProduct", err.str());
    } 
        

 std::cout << "prepareTrack 3 " << tkRefCollection->size() << std::endl;

    for(reco::TrackCollection::size_type i=0; i<tkRefCollection->size(); i++) {
      
      std::cout << "prepareTrack 5 " << std::endl;

      const reco::Track * track = &(*tkRefCollection)[i];
	
      //For the d0 = -dxy
      trackParameters.m_d0 = track->d0() * 10. ; //in mm
      //For the z0
      trackParameters.m_z0 = track->dz() * 10. ; //in mm
      //For the Track momentum at the 2D distance of closest approach
      //For tracks reconstructed in the CMS Tracker, the reference position is the point of closest approach to the centre of CMS. (math::XYZPoint posClosest = track->referencePoint();)
      // According to TrackBase.h the momentum() method returns momentum vector at innermost (reference) point on track
      const pandora::CartesianVector momentumAtClosestApproach(track->momentum().x(),track->momentum().y(),track->momentum().z()); //in GeV
      trackParameters.m_momentumAtDca = momentumAtClosestApproach;
 
      //For the track of the state at the start in mm and GeV
      const pandora::CartesianVector positionAtStart(track->innerPosition().x()* 10.,track->innerPosition().y()* 10., track->innerPosition().z() * 10. );
      const pandora::CartesianVector momentumAtStart(track->innerMomentum().x(),track->innerMomentum().y(), track->innerMomentum().z() );
      trackParameters.m_trackStateAtStart = pandora::TrackState(positionAtStart,momentumAtStart);
      //For the track of the state at the end in mm and GeV
      const pandora::CartesianVector positionAtEnd(track->outerPosition().x() * 10.,track->outerPosition().y() * 10., track->outerPosition().z() * 10.);
      const pandora::CartesianVector momentumAtEnd(track->outerMomentum().x(),track->outerMomentum().y(), track->outerMomentum().z() );
      trackParameters.m_trackStateAtEnd = pandora::TrackState(positionAtEnd,momentumAtEnd);
      //For the charge
      double charge = track->charge();
      // std::cout << "charge " << charge << std::endl;
      trackParameters.m_charge = charge;
      //Associate the reconstructed Track (in the Tracker) with the corresponding MC true simulated particle

      edm::RefToBase<reco::Track> tr(tkRefCollection, i);
      std::vector<std::pair<TrackingParticleRef, double> > tp;
      TrackingParticleRef tpr; 

      if(pRecoToSim.find(tr) != pRecoToSim.end()){
	tp = pRecoToSim[tr];
	std::cout << "Reco Track pT: "  << track->pt() <<  " matched to " << tp.size() << " MC Tracks" << " associated with quality: " << tp.begin()->second << std::endl;
	tpr = tp.begin()->first;
      }
      

      trackParameters.m_particleId = 211; // INITIALIZATION // NS

      //So the pdg code of the track is 
      if(pRecoToSim.find(tr) != pRecoToSim.end()) {

       std::cout << " EDW  KAI PALI tpr->pdgId() " << tpr->pdgId() << std::endl;

       trackParameters.m_particleId = (tpr->pdgId());
       std::cout << "the pdg id of this track is " << (tpr->pdgId()) << std::endl;
       //The parent vertex (from which this track was produced) has daughter particles.
       //These are the desire siblings of this track which we need. 
       TrackingParticleRefVector simSiblings = getTpSiblings(tpr);

       const TrackingParticle * sib; 
       int numofsibs = 0;
       std::vector<int> pdgidofsibs; pdgidofsibs.clear();

	    
      if (simSiblings.isNonnull()) {
	for(TrackingParticleRefVector::iterator si = simSiblings.begin(); si != simSiblings.end(); si++){
	  //Check if the current sibling is the track under study
	  if ( (*si) ==  tpr  ) {continue;}
	  sib = &(**si);
	  pdgidofsibs.push_back(sib->pdgId());
	  ++numofsibs;
	  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora, track, sib)); 
	}
	std::cout<< "This track has " << numofsibs << " sibling tracks with pdgids:" << std::endl; 
	for (std::vector<int>::iterator sib_pdg_it = pdgidofsibs.begin(); sib_pdg_it != pdgidofsibs.end(); sib_pdg_it++){
	  std::cout << (*sib_pdg_it) << std::endl;
	}
      } else {
	std::cout << "Particle pdgId = "<< (tpr->pdgId()) << " produced at rho = " << (tpr->vertex().Rho()) << ", z = " << (tpr->vertex().Z()) << ", has NO siblings!"<< std::endl;
      }
     
      //Now the track under study has daughter particles. To find them we study the decay vertices of the track
      TrackingParticleRefVector simDaughters = getTpDaughters(tpr);
      const TrackingParticle * dau; 
      int numofdaus = 0;
      std::vector<int> pdgidofdaus; pdgidofdaus.clear();

      if (simDaughters.isNonnull()) {
	for(TrackingParticleRefVector::iterator di = simDaughters.begin(); di != simDaughters.end(); di++){
	  //We have already checked that simDaughters don't contain the track under study
	  dau = &(**di);
	  pdgidofdaus.push_back(dau->pdgId());
	  ++numofdaus;
	  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora, track, dau)); 
	}
	std::cout<< "This track has " << numofdaus << " daughter tracks with pdgids:" << std::endl; 
	for (std::vector<int>::iterator dau_pdg_it = pdgidofdaus.begin(); dau_pdg_it != pdgidofdaus.end(); dau_pdg_it++){
	  std::cout << (*dau_pdg_it) << std::endl;
	}
      } else {
	std::cout << "Particle pdgId = "<< (tpr->pdgId()) << " produced at rho = " << (tpr->vertex().Rho()) << ", z = " << (tpr->vertex().Z()) << ", has NO daughters!"<< std::endl;
      }

} // KLEINW TO IF pRecoToSim.find(tr) != pRecoToSim.end() NS FIX


      //The mass 
      trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackParameters.m_particleId.Get());
 

      //For the ECAL entrance
      // Starting from outermost hit position of the track and propagating to ECAL Entrance
      float pfoutenergy = track->outerMomentum().Mag2();

      //the input BaseParticlePropagator needs to be cm and GeV
      BaseParticlePropagator theOutParticle = BaseParticlePropagator( RawParticle(XYZTLorentzVector(track->outerMomentum().x(),
												    track->outerMomentum().y(),
												    track->outerMomentum().z(),
												    pfoutenergy),
										  XYZTLorentzVector(track->outerPosition().x(),
												    track->outerPosition().y(),
												    track->outerPosition().z(),
												    0.)),
								      0.,0.,B_.z());

      theOutParticle.setCharge(track->charge());
      math::XYZPoint theOutParticle_position = math::XYZPoint(theOutParticle.vertex());
      math::XYZTLorentzVector theOutParticle_momentum = theOutParticle.momentum();
      std::cout << "magnetic field z " << B_.z() << std::endl;
      std::cout << "theOutParticle x position before propagation in cm "<< theOutParticle_position.x()<< std::endl;
      std::cout << "theOutParticle x momentum before propagation in cm "<< theOutParticle_momentum.x()<< std::endl;
      theOutParticle.propagateToEcalEntrance(false);
      bool reachesCalorimeter = false;
      bool isonendcap = false;

      //Set position and momentum after propagation to ECAL
      theOutParticle_position = math::XYZPoint(theOutParticle.vertex());
      theOutParticle_momentum = theOutParticle.momentum();
 
      std::cout << "theOutParticle x position after propagation to ECAL "<< theOutParticle_position.x()<< std::endl;
      std::cout << "theOutParticle x momentum after propagation to ECAL "<< theOutParticle_momentum.x()<< std::endl;
 	
      if(theOutParticle.getSuccess()!=0){
	// std::cout<< "!!!Reached ECAL!!! "<< std::endl;
	reachesCalorimeter = true;
      }
      if( abs(theOutParticle.getSuccess()) == 2){
	// std::cout<< "It is on the endcaps "<< std::endl;
	isonendcap = true;
      }

      trackParameters.m_reachesCalorimeter = reachesCalorimeter;
      if (reachesCalorimeter){
	const pandora::CartesianVector positionAtCalorimeter(theOutParticle_position.x() * 10.,theOutParticle_position.y() * 10.,theOutParticle_position.z() * 10.);//in mm
	const pandora::CartesianVector momentumAtCalorimeter(theOutParticle_momentum.x(),theOutParticle_momentum.y(),theOutParticle_momentum.z());
	trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(positionAtCalorimeter, momentumAtCalorimeter);
	// For the time at calorimeter we need the speed of light
	//This is in BaseParticlePropagator c_light() method in mm/ns but is protected (299.792458 mm/ns)
	//So we take it from CLHEP
	trackParameters.m_timeAtCalorimeter = positionAtCalorimeter.GetMagnitude() / speedoflight; // in ns
      
      } else { 
	trackParameters.m_trackStateAtCalorimeter = trackParameters.m_trackStateAtEnd.Get();
	trackParameters.m_timeAtCalorimeter = std::numeric_limits<float>::max();
      }

      trackParameters.m_isProjectedToEndCap = isonendcap; 

      bool canFormPfo = false; 
      bool canFormClusterlessPfo = false;
 
      //Add more criteria here
      if (trackParameters.m_reachesCalorimeter.Get()){
	canFormPfo = true;
	canFormClusterlessPfo = true;
	std::cout<< "Yes, this track can form pfo" << std::endl;
      }
      trackParameters.m_canFormPfo = canFormPfo;
      trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;


     //The parent address
      trackParameters.m_pParentAddress =  (void *) track;

      //Some cout
      // std::cout <<  track->innerDetId() << std::endl;
      // std::cout <<  track->outerPx() << std::endl;
      // std::cout <<  track->d0() << std::endl;
      // std::cout <<  track->dz() << std::endl;
      // std::cout <<  pfrt->pdgCode() << std::endl;
 
       PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
    }

    // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraMonitoringApi::VisualizeTracks(  &(*tkRefCollection)  , "currentTrackList", AUTO, false, true  ) );
    // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraMonitoringApi::VisualizeTracks(  &(*tkRefCollection)  , "currentTrackList",  true  ) );
    // PANDORA_MONITORING_API(VisualizeTracks(  &(*tkRefCollection)  , "currentTrackList", AUTO, false, true  ) );
    // PANDORA_MONITORING_API(ViewEvent() );


  }

}

void runPandora::prepareHits( edm::Handle<EcalRecHitCollection> ecalRecHitHandleEB,
			      edm::Handle<HBHERecHitCollection> hcalRecHitHandleHBHE, 
			      edm::Handle<HGCRecHitCollection> HGCeeRecHitHandle,
			      edm::Handle<HGCRecHitCollection> HGChefRecHitHandle,
			      edm::Handle<HGCRecHitCollection> HGChebRecHitHandle,
			      reco::Vertex& pv, 
			      const edm::Event& iEvent, const edm::EventSetup& iSetup){

  PandoraApi::RectangularCaloHitParameters caloHitParameters;

  double speedoflight = (CLHEP::c_light/CLHEP::cm)/CLHEP::ns;
  std::cout<< speedoflight << " cm/ns" << std::endl;

  double sumCaloEnergy = 0.;
  double sumCaloEnergyEM = 0.;
  double sumCaloEnergyHAD = 0.;
  double sumCaloECALEnergyEM = 0.;
  double sumCaloHCALEnergyEM  = 0.;
  double sumCaloECALEnergyHAD= 0.;
  double sumCaloHCALEnergyHAD = 0.;



  // 
  // Process ECAL barrel rechits 
  // 
  for(unsigned i=0; i<ecalRecHitHandleEB->size(); i++) {
	  
    const EcalRecHit * erh = &(*ecalRecHitHandleEB)[i];
    const DetId& detid = erh->detid();
    double energy = erh->energy();

    if (energy < m_eCalThresBarrel) continue;

    double time = erh->time();
    // std::cout << "energy " << energy <<  " time " << time <<std::endl;
    EcalSubdetector esd=(EcalSubdetector)detid.subdetId();
    if (esd != 1) continue;
	  
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
    // get the ecalBarrel geometry
    const CaloSubdetectorGeometry *ebtmp = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  
    const EcalBarrelGeometry* ecalBarrelGeometry = dynamic_cast< const EcalBarrelGeometry* > (ebtmp);
    assert( ecalBarrelGeometry );

    // get the ecalBarrel topology
    EcalBarrelTopology ecalBarrelTopology(geoHandle);

    const CaloSubdetectorGeometry* geom = ecalBarrelGeometry;
    
    const CaloCellGeometry *thisCell = geom->getGeometry(detid);
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraECAL") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      continue;
    }
  	  
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
 //   for (unsigned int i=0; i<8; i++) { 
 //     std::cout << "Corners " << i << ": x " << corners[i].x() << " y " << corners[i].y() << " z " << corners[i].z() << std::endl ; 
 //   }
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    // dist = cm, c = cm/nsec, rechit t in psec
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; 

    caloHitParameters.m_hitType = pandora::ECAL;
    caloHitParameters.m_detectorRegion = pandora::BARREL;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_eCalToEMGeVBarrel * energy; 
    caloHitParameters.m_mipEquivalentEnergy = m_eCalToMipBarrel * energy; // = energy; // HTAN 0 NS AP


    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_eCalMipThresBarrel) {
       //std::cout << "EcalBarrel MIP threshold rejected" << std::endl;
       continue;
    }

    sumCaloEnergy += energy;
    sumCaloEnergyEM += energy  ;//FIXME * absorberCorrectionEM;
    sumCaloEnergyHAD += energy ;//FIXME * absorberCorrectionHAD;
    sumCaloECALEnergyEM  += energy  ;//* absorberCorrectionEM;
    sumCaloECALEnergyHAD += energy  ;//* absorberCorrectionHAD;

    caloHitParameters.m_hadronicEnergy = m_eCalToHadGeVBarrel * energy; // = energy; 
    caloHitParameters.m_layer = 1.;//PFLayer::ECAL_BARREL;
    caloHitParameters.m_nCellRadiationLengths = 0.0; // 6.;
    caloHitParameters.m_nCellInteractionLengths = 0.0; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) erh;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 
    
  }	

  //
  // process HCAL Barrel Hits
  //
  for(unsigned i=0; i<hcalRecHitHandleHBHE->size(); i++) {
	  
    const HBHERecHit * hrh = &(*hcalRecHitHandleHBHE)[i];
    const DetId& detid = hrh->detid();
    double energy = hrh->energy();
    if (energy < m_hCalThresBarrel) continue;
    double time = hrh->time();

    HcalSubdetector hsd=(HcalSubdetector)detid.subdetId();
    if (hsd != 1) continue; // HCAL BARREL ONLY

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
    // get the hcal geometry
    const CaloSubdetectorGeometry *hbtmp = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
    const HcalGeometry* hcalBarrelGeometry = dynamic_cast< const HcalGeometry* > (hbtmp);
    assert( hcalBarrelGeometry );

    const CaloSubdetectorGeometry* hbgeom = hcalBarrelGeometry;
    const CaloCellGeometry *thisCell = hbgeom->getGeometry(detid) ;
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHCAL") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      continue;
    }
  	  
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    // dist = cm, c = cm/nsec, rechit t in psec
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; 

    caloHitParameters.m_hitType = pandora::HCAL;
    caloHitParameters.m_detectorRegion = pandora::BARREL ;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVBarrel * energy; 
    caloHitParameters.m_mipEquivalentEnergy =  m_hCalToMipBarrel * energy; // = energy; // HTAN 0 NS AP  

    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_hCalMipThresBarrel) {
       //std::cout << "HcalBarrel MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy += energy;
    sumCaloEnergyEM  += energy;// * absorberCorrectionEM ;
    sumCaloEnergyHAD += energy;// * absorberCorrectionHAD;
    sumCaloHCALEnergyEM  += energy  ;//* absorberCorrectionEM;
    sumCaloHCALEnergyHAD += energy  ;//* absorberCorrectionHAD;


    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVBarrel * energy; // = energy; 

    caloHitParameters.m_layer = 1.;
    caloHitParameters.m_nCellRadiationLengths = 0.0; // 6.;
    caloHitParameters.m_nCellInteractionLengths = 0.0; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) hrh;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
  }

  //
  // Process HGC EE rec hits 
  // 
  int nNotFoundEE = 0, nCaloHitsEE = 0 ; 
  for(unsigned i=0; i<HGCeeRecHitHandle->size(); i++) {
	  
    const HGCRecHit *eerh = &(*HGCeeRecHitHandle)[i];
    const HGCEEDetId& detid = eerh->id();
    //double energy = eerh->energy();
    double energy = eerh->energy() * m_Calibr_ADC2GeV_EE;
    if (energy<m_eCalThresEndCap) continue;
    energy *= m_addCalibrEE; //this must be after threshold rejection
    double time = eerh->time();

    float nCellInteractionLengths = 0.;
    float nCellRadiationLengths = 0.;
    float absorberCorrectionEM = 1.;
    float absorberCorrectionHAD = 1.;
    getLayerPropertiesEE(eerh, detid.layer(),
          nCellInteractionLengths, nCellRadiationLengths,
          absorberCorrectionEM, absorberCorrectionHAD);


     //std::cout << "HGC EE rechit cell " << detid.cell() << ", sector " << detid.sector() 
     //	      << ", subsector " << detid.subsector() << ", layer " << detid.layer() << ", z " << detid.zside() 
     //	      << " energy " << energy <<  " and time " << time << std::endl;

    ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
    if (thesubdet != 3) continue; // HGC EE
	  
    math::XYZVector axis;

    edm::ESHandle<HGCalGeometry> hgceeGeoHandle ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle) ; 
    const HGCalGeometry &hgceetmp = *hgceeGeoHandle ; 

    const HGCalGeometry &hgceeGeometry = dynamic_cast< const HGCalGeometry& > (hgceetmp);
    assert( &hgceeGeometry );

    // Get the HGC topology
    // edm::ESHandle<HGCalTopology> hgceeTopoHandle ; 
    // iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeTopoHandle) ; 
    // const HGCalTopology &hgceeTopology = *hgceeTopoHandle ; 
    
    // const HGCalGeometry geom = hgceeGeometry;
    const FlatTrd *thisCell = static_cast<const FlatTrd*>(hgceeGeometry.getGeometry(detid));	
    // const CaloCellGeometry *thisCell = hgceeGeometry.getGeometry(detid);
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHGCEE") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      nNotFoundEE++ ; 
      continue;
    }
  	  
    const HGCalGeometry::CornersVec corners = ( std::move( hgceeGeometry.getCorners( detid ) ) );
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    //Hieu
    //caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    //caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; // dist = cm, c = cm/nsec, rechit t in psec
    
    caloHitParameters.m_hitType = pandora::ECAL;
    caloHitParameters.m_detectorRegion = pandora::ENDCAP;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_eCalToEMGeVEndCap * energy * absorberCorrectionEM; 
    caloHitParameters.m_mipEquivalentEnergy =  m_eCalToMipEndCap * energy / m_addCalibrEE; // = energy; // HTAN 0 NS AP  
    
    float angleCorrectionMIP(1.);
    float hitR = distToFrontFace; 
    float hitZ = zf;
    angleCorrectionMIP = hitR/hitZ;
    h_MIP_EE -> Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
    h_MIP_Corr_EE -> Fill(caloHitParameters.m_mipEquivalentEnergy.Get()*angleCorrectionMIP);

    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_eCalMipThresEndCap) {
       //std::cout << "EE MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy    += energy;
    sumCaloEnergyEM  += energy * absorberCorrectionEM ;
    sumCaloEnergyHAD += energy * absorberCorrectionHAD;
    sumCaloECALEnergyEM  += energy  * absorberCorrectionEM;
    sumCaloECALEnergyHAD += energy  * absorberCorrectionHAD;



    caloHitParameters.m_hadronicEnergy = m_eCalToHadGeVEndCap * energy * absorberCorrectionHAD;
    caloHitParameters.m_layer = detid.layer() ;//PFLayer::ECAL_BARREL;
    caloHitParameters.m_nCellRadiationLengths = nCellRadiationLengths; // 6.;
    caloHitParameters.m_nCellInteractionLengths = nCellInteractionLengths; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) eerh;

    //std::cout << "TEST input caloHitParameters: "
    //   << ", hadronicEnergy = " << caloHitParameters.m_hadronicEnergy.Get()
    //   << ", hitType = " << caloHitParameters.m_hitType.Get()
    //   << ", m_detectorRegion = " << caloHitParameters.m_detectorRegion.Get()
    //   << std::endl;

    // std::cout << "Parameters for input: " << std::endl ; 
    // std::cout << "position vector X,Y,Z: " << caloHitParameters.m_positionVector.GetX() << " " 
    // 	      << caloHitParameters.m_positionVector.GetY() << " " << caloHitParameters.m_positionVector.GetZ() << std::endl ; 
    // std::cout << "expected direction X,Y,Z: " << caloHitParameters.m_expectedDirection.GetX() << " " 
    // 	      << caloHitParameters.m_expectedDirection.GetY() << " " << caloHitParameters.m_expectedDirection.GetZ() << std::endl ; 
    

    // std::cout << "Processing HGC EE rec hit with position (in cm) " << x << " " << y << " " << z << std::endl ; 
    
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 
    nCaloHitsEE++ ; 
     //std::cout << "BMD: This calo hit has been created" << std::endl ;   
  }


  //
  // Process HGC HEF rec hits 
  // 
  int nNotFoundHEF = 0, nCaloHitsHEF = 0 ; 
  for(unsigned i=0; i<HGChefRecHitHandle->size(); i++) {
	  
    const HGCRecHit *hefrh = &(*HGChefRecHitHandle)[i];
    const HGCHEDetId& detid = hefrh->id();
    double energy = hefrh->energy() * m_Calibr_ADC2GeV_HEF;
    if (energy<m_hCalThresEndCapHEF)
       continue;
    energy *= m_addCalibrHEF;
    double time = hefrh->time();

    float nCellInteractionLengths = 0.;
    float nCellRadiationLengths = 0.;
    float absorberCorrectionEM = 1.;
    float absorberCorrectionHAD = 1.;
    getLayerPropertiesHEF(hefrh, detid.layer(),
      nCellInteractionLengths, nCellRadiationLengths,
      absorberCorrectionEM, absorberCorrectionHAD);


    // std::cout << "HGC HEF rechit cell " << detid.cell() << ", sector " << detid.sector() 
    // 	      << ", subsector " << detid.subsector() << ", layer " << detid.layer() << ", z " << detid.zside() 
    // 	      << " energy " << energy <<  " and time " << time << std::endl;

    ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
    if (thesubdet != 4) continue; // HGC HEF
	  
    math::XYZVector axis;

    edm::ESHandle<HGCalGeometry> hgchefGeoHandle ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle) ; 
    const HGCalGeometry &hgcheftmp = *hgchefGeoHandle ; 

    const HGCalGeometry &hgchefGeometry = dynamic_cast< const HGCalGeometry& > (hgcheftmp);
    assert( &hgchefGeometry );

    // Get the HGC topology
    // edm::ESHandle<HGCalTopology> hgceeTopoHandle ; 
    // iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeTopoHandle) ; 
    // const HGCalTopology &hgceeTopology = *hgceeTopoHandle ; 
    
    // const HGCalGeometry geom = hgceeGeometry;
    const FlatTrd *thisCell = static_cast<const FlatTrd*>(hgchefGeometry.getGeometry(detid));	
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHGCHEF") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      nNotFoundHEF++ ; 
      continue;
    }
  	  
    const HGCalGeometry::CornersVec corners = ( std::move( hgchefGeometry.getCorners( detid ) ) );
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; // dist = cm, c = cm/nsec, rechit t in psec
    
    caloHitParameters.m_hitType = pandora::HCAL;
    caloHitParameters.m_detectorRegion = pandora::ENDCAP;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVEndCapHEF * energy; 
    caloHitParameters.m_mipEquivalentEnergy =  m_hCalToMipEndCapHEF * energy / m_addCalibrHEF; // = energy; // HTAN 0 NS AP  

    h_MIP_HEF->Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
    float angleCorrectionMIP(1.);
    float hitR = distToFrontFace; 
    float hitZ = zf;
    angleCorrectionMIP = hitR/hitZ;
    h_MIP_Corr_HEF->Fill(caloHitParameters.m_mipEquivalentEnergy.Get()*angleCorrectionMIP);

    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_hCalMipThresEndCapHEF) {
       //std::cout << "HEF MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy += energy;
    sumCaloEnergyEM  += energy * absorberCorrectionEM ;
    sumCaloEnergyHAD += energy * absorberCorrectionHAD;
    sumCaloHCALEnergyEM  += energy  * absorberCorrectionEM;
    sumCaloHCALEnergyHAD += energy  * absorberCorrectionHAD;



    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVEndCapHEF * energy ; // = energy; 
    caloHitParameters.m_layer = detid.layer() ;
    caloHitParameters.m_nCellRadiationLengths = 0.0; // 6.;
    caloHitParameters.m_nCellInteractionLengths = 0.0; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) hefrh;

    // std::cout << "Parameters for input: " << std::endl ; 
    // std::cout << "position vector X,Y,Z: " << caloHitParameters.m_positionVector.GetX() << " " 
    // 	      << caloHitParameters.m_positionVector.GetY() << " " << caloHitParameters.m_positionVector.GetZ() << std::endl ; 
    // std::cout << "expected direction X,Y,Z: " << caloHitParameters.m_expectedDirection.GetX() << " " 
    // 	      << caloHitParameters.m_expectedDirection.GetY() << " " << caloHitParameters.m_expectedDirection.GetZ() << std::endl ; 
    

    // std::cout << "Processing HGC HEF rec hit with position (in cm) " << x << " " << y << " " << z << std::endl ; 
    
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 

    // std::cout << "BMD: This calo hit has been created" << std::endl ;   
    nCaloHitsHEF++ ; 
  }


  //
  // Process HGC HEB rec hits 
  // 
  int nNotFoundHEB = 0, nCaloHitsHEB = 0 ; 
  for(unsigned i=0; i<HGChebRecHitHandle->size(); i++) {
	  
    const HGCRecHit *hebrh = &(*HGChebRecHitHandle)[i];
    const HGCHEDetId& detid = hebrh->id();
    double energy = hebrh->energy() * m_Calibr_ADC2GeV_HEB;
    if (energy<m_hCalThresEndCapHEB) continue;
    energy *= m_addCalibrHEB;
    double time = hebrh->time();
    float nCellInteractionLengths = 0.;
    float nCellRadiationLengths = 0.;
    float absorberCorrectionEM = 1.;
    float absorberCorrectionHAD = 1.;
    getLayerPropertiesHEB(hebrh, detid.layer(),
      nCellInteractionLengths, nCellRadiationLengths,
      absorberCorrectionEM, absorberCorrectionHAD);


     //std::cout << "HGC HEB rechit cell " << detid.cell() << ", sector " << detid.sector() 
     //	      << ", subsector " << detid.subsector() << ", layer " << detid.layer() << ", z " << detid.zside() 
     //	      << " energy " << energy <<  " and time " << time << std::endl;

    ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
    if (thesubdet != 5) continue; // HGC HEB
	  
    math::XYZVector axis;

    edm::ESHandle<HGCalGeometry> hgchebGeoHandle ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle) ; 
    const HGCalGeometry &hgchebtmp = *hgchebGeoHandle ; 

    const HGCalGeometry &hgchebGeometry = dynamic_cast< const HGCalGeometry& > (hgchebtmp);
    assert( &hgchebGeometry );

    // Get the HGC topology
    // edm::ESHandle<HGCalTopology> hgceeTopoHandle ; 
    // iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeTopoHandle) ; 
    // const HGCalTopology &hgceeTopology = *hgceeTopoHandle ; 
    
    // const HGCalGeometry geom = hgceeGeometry;
    const FlatTrd *thisCell = static_cast<const FlatTrd*>(hgchebGeometry.getGeometry(detid));	
    // const CaloCellGeometry *thisCell = hgchebGeometry.getGeometry(detid);
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHGCHEB") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      nNotFoundHEB++ ; 
      continue;
    }
  	  
    const HGCalGeometry::CornersVec corners = ( std::move( hgchebGeometry.getCorners( detid ) ) );
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; // dist = cm, c = cm/nsec, rechit t in psec
    
    caloHitParameters.m_hitType = pandora::HCAL;
    caloHitParameters.m_detectorRegion = pandora::ENDCAP;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVEndCapHEB * energy; 
    caloHitParameters.m_mipEquivalentEnergy = m_hCalToMipEndCapHEB * energy / m_addCalibrHEB; // = energy; // HTAN 0 NS AP  
    h_MIP_HEB->Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_hCalMipThresEndCapHEB) {
       //std::cout << "HEB MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy += energy;
    sumCaloEnergyEM  += energy * absorberCorrectionEM ;
    sumCaloEnergyHAD += energy * absorberCorrectionHAD;
    sumCaloHCALEnergyEM  += energy  * absorberCorrectionEM;
    sumCaloHCALEnergyHAD += energy  * absorberCorrectionHAD;


    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVEndCapHEB * energy; // = energy; 
    caloHitParameters.m_layer = 31 + detid.layer();//PFLayer::ECAL_BARREL;
    caloHitParameters.m_nCellRadiationLengths = 0.0; // 6.;
    caloHitParameters.m_nCellInteractionLengths = 0.0; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) hebrh;

    // std::cout << "Parameters for input: " << std::endl ; 
    // std::cout << "position vector X,Y,Z: " << caloHitParameters.m_positionVector.GetX() << " " 
    // 	      << caloHitParameters.m_positionVector.GetY() << " " << caloHitParameters.m_positionVector.GetZ() << std::endl ; 
    // std::cout << "expected direction X,Y,Z: " << caloHitParameters.m_expectedDirection.GetX() << " " 
    // 	      << caloHitParameters.m_expectedDirection.GetY() << " " << caloHitParameters.m_expectedDirection.GetZ() << std::endl ; 
    

    // std::cout << "Processing HGC HEB rec hit with position (in cm) " << x << " " << y << " " << z << std::endl ; 
    
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 

    nCaloHitsHEB++ ; 
    // std::cout << "BMD: This calo hit has been created" << std::endl ;   
  }

  h_sumCaloE->Fill(sumCaloEnergy);
  h_sumCaloEM->Fill(sumCaloEnergyEM);
  h_sumCaloHad->Fill(sumCaloEnergyHAD);

  h2_Calo_EM_hcalEecalE->Fill(sumCaloECALEnergyEM, sumCaloHCALEnergyEM);
  h2_Calo_Had_hcalEecalE->Fill(sumCaloECALEnergyHAD, sumCaloHCALEnergyHAD);

  std::cout << "sumCaloEnergy = " << sumCaloEnergy << std::endl;
  std::cout << "sumCaloEnergyEM  = " << sumCaloEnergyEM  << std::endl;
  std::cout << "sumCaloEnergyHAD = " << sumCaloEnergyHAD << std::endl;

  std::cout << "prepareHits HGC summary: " << std::endl ; 
  std::cout << "HGC Calo Hits               : " << nCaloHitsEE << " (HGC EE) " 
	    << nCaloHitsHEF << " (HGC HEF) " << nCaloHitsHEB << " (HGC HEB) " << std::endl ;
  std::cout << "DetIDs not found in geometry: " << nNotFoundEE << " (HGC EE) " 
	    << nNotFoundHEF << " (HGC HEF) " << nNotFoundHEB << " (HGC HEB) " << std::endl ;

}

void runPandora::preparemcParticle(edm::Handle<std::vector<reco::GenParticle> > genpart){ // function to setup a mcParticle for pandora
  // PandoraPFANew/v00-09/include/Api/PandoraApi.h
  //class MCParticleParameters
  //{
  //public:
  //    pandora::InputFloat             m_energy;                   ///< The energy of the MC particle, units GeV
  //    pandora::InputCartesianVector   m_momentum;                 ///< The momentum of the MC particle, units GeV
  //    pandora::InputCartesianVector   m_vertex;                   ///< The production vertex of the MC particle, units mm
  //    pandora::InputCartesianVector   m_endpoint;                 ///< The endpoint of the MC particle, units mm
  //    pandora::InputInt               m_particleId;               ///< The MC particle's ID (PDG code)
  //    pandora::InputAddress           m_pParentAddress;           ///< Address of the parent MC particle in the user framework
  //};
  // for(std::vector<reco::GenParticle>::const_iterator cP = genpart->begin();  cP != genpart->end(); cP++ ) {
  
  for(size_t i = 0; i < genpart->size(); ++ i) {
    const GenParticle * pa = &(*genpart)[i];
    PandoraApi::MCParticle::Parameters parameters;
    parameters.m_energy = pa->energy();
    parameters.m_momentum = pandora::CartesianVector(pa->px() , pa->py(),  pa->pz() );
    parameters.m_vertex = pandora::CartesianVector(pa->vx() * 10. , pa->vy() * 10., pa->vz() * 10. ); //in mm
    // parameters.m_endpoint = pandora::CartesianVector(position.x(), position.y(), position.z());
    // Definition of the enpoint depends on the application that created the particle, e.g. the start point of the shower in a calorimeter. 
    // If the particle was not created as a result of a continuous process where the parent particle continues, i.e.
    // hard ionization, Bremsstrahlung, elastic interactions, etc. then the vertex of the daughter particle is the endpoint.
    parameters.m_endpoint = pandora::CartesianVector(pa->vx() * 10. , pa->vy() * 10., pa->vz() * 10. ); //IS THIS CORRECT?!
    parameters.m_particleId = pa->pdgId();
    parameters.m_pParentAddress = (void*) pa;
    if( abs(pa->pdgId()) == 211 ){Egenpart->Fill( pa->energy() ); }; //AP NS // 11 electrons, 13 muons, 211 pions 
    std::cout << "The mc particle pdg id " << pa->pdgId() << " with energy " << pa->energy() << std::endl;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, parameters));
    
    // Create parent-daughter relationships
    // HepMC::GenParticle * theParticle = new HepMC::GenParticle(HepMC::FourVector(0.,0.,0.,0.),12,1);
    size_t n = pa->numberOfDaughters();
    //std::cout << "The mc particle pdg id " << pa->pdgId() << " with energy " << pa->energy() << " and " << n << " daughters " <<  std::endl;
    for(size_t j = 0; j < n; ++ j) {
      const Candidate * d = pa->daughter( j );
      //if we want to keep it also in GenParticle uncomment here
      const GenParticle * da = NULL;
      //We need to check if this daughter has an integer charge
      bool integercharge = ( ( (int) d->charge() ) - (d->charge()) ) == 0 ? true : false;
      std::cout << "checking integer charge: Real Charge " << d->charge() << " int part " << ( (int) d->charge() ) << " bool " << std::endl;
      da = new GenParticle( d->charge(), d->p4() , d->vertex() , d->pdgId() , d->status() , integercharge); 
      
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora, pa , da));
      
    }


  }

}

void runPandora::preparePFO(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	// PandoraPFANew/v00-09/include/Pandora/PandoraInternal.h
	// typedef std::set<ParticleFlowObject *> PfoList;  
	// 	PandoraPFANew/v00-09/include/Api/PandoraContentApi.h
    	//	class ParticleFlowObject
    	//	{
    	//	public:
    	//	    /**
    	//	     *  @brief  Parameters class
    	//	     */
    	//	    class Parameters
    	//	    {
    	//	    public:
    	//	        pandora::InputInt               m_particleId;       ///< The particle flow object id (PDG code)
    	//	        pandora::InputInt               m_charge;           ///< The particle flow object charge
    	//	        pandora::InputFloat             m_mass;             ///< The particle flow object mass
    	//	        pandora::InputFloat             m_energy;           ///< The particle flow object energy
    	//	        pandora::InputCartesianVector   m_momentum;         ///< The particle flow object momentum
    	//	        pandora::ClusterList            m_clusterList;      ///< The clusters in the particle flow object
    	//	        pandora::TrackList              m_trackList;        ///< The tracks in the particle flow object
    	//	    };
    	//	    /**
    	//	     *  @brief  Create a particle flow object
    	//	     * 
    	//	     *  @param  algorithm the algorithm creating the particle flow object
    	//	     *  @param  particleFlowObjectParameters the particle flow object parameters
    	//	     */
    	//	    static pandora::StatusCode Create(const pandora::Algorithm &algorithm, const Parameters &parameters);
    	//	};
    	//	typedef ParticleFlowObject::Parameters ParticleFlowObjectParameters;

	// const pandora::CartesianVector momentum(1., 2., 3.);
	// for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){
	//   (*itPFO)->SetParticleId();
	//   (*itPFO)->SetCharge();
	//   (*itPFO)->SetMass();
	//   (*itPFO)->SetEnergy();
	//   (*itPFO)->SetMomentum();
	// }
  
  const pandora::PfoList *pPfoList = NULL;
  // PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));

  const pandora::StatusCode statusCode(PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));  
  if (pandora::STATUS_CODE_SUCCESS != statusCode){throw pandora::StatusCodeException(statusCode);}

  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart);
//NS ADD

  for(size_t i = 0; i < genpart->size(); ++ i) {

    const GenParticle * pa = &(*genpart)[i];
    double ene_true    = pa->energy();
    double charge_true = 0;
    if(pa->pdgId()>0) charge_true = 1;
    if(pa->pdgId()<0) charge_true = -1;
    std::cout << " iN GENPART LOOP " << std::endl;
     if(abs(pa->pdgId())== 211){
      double diff_min     = 1e9;
      double energy_match = -1;

      for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){ // 4
        double charge = (*itPFO)->GetCharge() ;
        double energy = (*itPFO)->GetEnergy() ;
        double diff = 1e10;
        std::cout << " charge true " << charge_true << " ene_true " << ene_true << std::endl;
        std::cout << " charge " << charge << " energy " << energy << std::endl;
        if(charge_true == charge){
          diff = abs(energy-ene_true);
          if(diff<diff_min){ // 3
           diff_min = diff;
           energy_match = energy;
          }    // 3
        }
      } // 4

      if(energy_match>0){ // 2
       Energy_res->Fill((energy_match-ene_true)/ene_true);
       std::cout << " FOUND MATCH " << energy_match << " to " << ene_true << std::endl;
      } // 2
   } // pdg
 } // 1

  std::cout << "Starting reading PFOs" << std::endl;
  float _sumPFOEnergy(0.);
  double sumClustEMEcalE(0.); //PFO cluster energy in Ecal
  double sumClustEMHcalE(0.); //PFO cluster energy in Hcal
  double sumClustHADEcalE(0.); //PFO cluster energy in Ecal
  double sumClustHADHcalE(0.); //PFO cluster energy in Hcal

  int nbPFOs(0);
  for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){
    nbPFOs++;
    std::cout << "Particle Id: " << (*itPFO)->GetParticleId() << std::endl;
    //std::cout << "Charge: " << (*itPFO)->GetCharge() << std::endl;
    //std::cout << "Mass: " << (*itPFO)->GetMass() << std::endl;
    std::cout << "Energy: " << (*itPFO)->GetEnergy() << std::endl;
    if ( (*itPFO)->GetParticleId() > 0  ){Epfos->Fill( (*itPFO)->GetEnergy() );} // AP NS//

    //For the cluster we will deal with it after the finishing of the calo hit
    const pandora::ClusterAddressList clusterAddressList((*itPFO)->GetClusterAddressList());
    const pandora::TrackAddressList trackAddressList((*itPFO)->GetTrackAddressList());
    const pandora::TrackList trackList((*itPFO)->GetTrackList());
   
    //TGClient *gclient  = NULL;
    //PandoraMonitoringApi::VisualizeTracks(  &trackList  , "currentTrackList", AUTO);
    //PANDORA_MONITORING_API(ViewEvent() );

    PANDORA_MONITORING_API(VisualizeTracks(  trackAddressList  , "currentTrackList", AUTO, false, true  ) );    
    PANDORA_MONITORING_API(ViewEvent() );

   std::cout << " ARE YOU GEETING IN HERE 1  ns, PFO(" << nbPFOs << "): " << clusterAddressList.size() << std::endl;

    for (pandora::ClusterAddressList::const_iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end(); itCluster != itClusterEnd; ++itCluster){
      const unsigned int nHitsInCluster((*itCluster).size());

      std::cout << "# of hits in cluster " << nHitsInCluster << std::endl; 
      
      //const pandora::CaloHitAddressList &caloHitAddressList(*itCluster);

      //EcalRecHit * hgcrh = NULL;
      //HBHERecHit * hrh = NULL;

      
      //int nbNonEHcalHit = 0;
      for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
      {

         const HGCRecHit *hgcHit = (HGCRecHit*)((*itCluster)[iHit]);
         //const HGCEEDetId& detidEE = hgcHit->id();
         float nCellInteractionLengths = 0.;
         float nCellRadiationLengths = 0.;
         float absorberCorrectionEM = 1.;
         float absorberCorrectionHAD = 1.;

         const DetId& detid = hgcHit->id();
         if (!detid)
            continue;

         //const HGCHEDetId& detidHE = hgcHit->id();

         ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
         //ForwardSubdetector thesubdetHE = (ForwardSubdetector)detidHE.subdetId();
         if (thesubdet == 3) {
            int layer = (int) ((HGCEEDetId)(detid)).layer() ;
            getLayerPropertiesEE(hgcHit, layer,
                  nCellInteractionLengths, nCellRadiationLengths,
                  absorberCorrectionEM, absorberCorrectionHAD);
            sumClustEMEcalE += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * absorberCorrectionEM * m_addCalibrEE;
            sumClustHADEcalE += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * absorberCorrectionHAD * m_addCalibrEE;
         }
         else if (thesubdet == 4) {
            sumClustEMHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToEMGeVEndCapHEF * absorberCorrectionEM * m_addCalibrHEF;
            sumClustHADHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToHadGeVEndCapHEF * absorberCorrectionEM * m_addCalibrHEF;
         }
         else if (thesubdet == 5) {
            sumClustEMHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * absorberCorrectionEM * m_addCalibrHEB;
            sumClustHADHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * absorberCorrectionEM * m_addCalibrHEB;
         }
         else {
         }

         //std::cout << "TEST ecal energy: " << sumClustEcalE << std::endl;
         //std::cout << "TEST hcal energy: " << sumClustHcalE << std::endl;
      }

//      for (pandora::CaloHitAddressList::const_iterator hIter = caloHitAddressList.begin(), hIterEnd = caloHitAddressList.end(); hIter != hIterEnd; ++hIter){
	
//	pandora::CaloHit * ch = (pandora::CaloHit *) (*hIter);
//	EcalRecHit * hgcrh = NULL;
//	HBHERecHit * hrh = NULL;
//
//	if (ch->GetHitType() ==  pandora::ECAL) { 
//	  hgcrh = (EcalRecHit *) (*hIter);
//	  std::cout << "EcalRecHit energy " << hgcrh->energy() <<  std::endl;
//     sumClustEcalE += hgcrh->energy();
//	} else if (ch->GetHitType() ==  pandora::HCAL) {  
//	  hrh = (HBHERecHit *) (*hIter); 
//	  std::cout << "HcalRecHit energy " << hrh->energy() <<  std::endl;		  
//     sumClustHcalE += hrh->energy();
//	}
//	else {
//	  std::cout << " No ECAL or HCAL??? What is this? " << ch->GetHitType() << std::endl;
//      nbNonEHcalHit++;
//	}
	
//      }
//      std::cout << "nbNonEHcalHit: " << nbNonEHcalHit << std::endl;
 

      // for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit){
      // 	EVENT::CalorimeterHit *pCalorimeterHit = (CalorimeterHit*)((*itCluster)[iHit]);
      // }
 
    }


    for (pandora::TrackAddressList::const_iterator itTrack = trackAddressList.begin(), itTrackEnd = trackAddressList.end();itTrack != itTrackEnd; ++itTrack){
      reco::Track * track =  (reco::Track *) (*itTrack);
      std::cout<< "Track from pfo charge " << track->charge() << std::endl;
      std::cout<< "Track from pfo transverse momentum " << track->pt() << std::endl;

    }

  }

  h2_EM_hcalEecalE->Fill(sumClustEMEcalE,sumClustEMHcalE);
  h2_Had_hcalEecalE->Fill(sumClustHADEcalE,sumClustHADHcalE);
  h_sumPfoE->Fill(_sumPFOEnergy);
  //h_sumPfoEEM->Fill(_sumPFOEnergyEM);
  //h_sumPfoEHad->Fill(_sumPFOEnergyHad);
  h_nbPFOs->Fill(nbPFOs);

	
}

//Get the track siblings
TrackingParticleRefVector runPandora::getTpSiblings(TrackingParticleRef tp){

  if (tp.isNonnull() && tp->parentVertex().isNonnull() && !tp->parentVertex()->daughterTracks().empty()) {
    return tp->parentVertex()->daughterTracks();
  } else {
    return TrackingParticleRefVector();
  }

}
//Get the track daughters
TrackingParticleRefVector runPandora::getTpDaughters(TrackingParticleRef tp){

  TrackingVertexRefVector trvertexes;
  TrackingParticleRefVector trdaughter;

  if (tp.isNonnull() && tp->decayVertices().isNonnull() ) {
    trvertexes = tp->decayVertices();
    //Loop on vector of TrackingVertex objects where the TrackingParticle decays. 
    for(TrackingVertexRefVector::iterator vi = trvertexes.begin(); vi != trvertexes.end(); vi++){
      //loop on all daughter tracks 
      for(TrackingParticleRefVector::iterator di = (**vi).daughterTracks().begin(); di != (**vi).daughterTracks().end(); di++){
	//Check if the daughter is the same as our mother tp particle
	if ( (*di) == tp  ) {continue;}
	trdaughter.push_back( (*di) );
      }//end on loop over daughter
    }//end on loop over vertices
    return trdaughter;
  } else {
    return TrackingParticleRefVector();
  }

}



// ------------ method called once each job just before starting event loop  ------------
void runPandora::beginJob()
{   
  std::cout << "I am beginning my job...LOOK!!!" << std::endl ; 
  firstEvent_ = true ;


// AP
  const char *pDisplay(::getenv("DISPLAY"));
  if (NULL == pDisplay) {
    std::cout << "DISPLAY environment not set" << std::endl;
  }  else {
    std::cout << "DISPLAY environment set to " << pDisplay << std::endl;
  }
  int argc = 0;
  char* argv = (char *)"";
  TApplication *m_pApplication;
  m_pApplication = gROOT->GetApplication();
  std::cout << "In runPandora::beginJob gVirtualX->GetDisplay()" << gVirtualX->GetDisplay() << std::endl;
  if(!m_pApplication){
    std::cout << "In if of m_pApplication in runPandora::beginJob " << std::endl;
    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
  } 
// END AP



  file = new TFile(_outputFileName.c_str(),"recreate");
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);
  // dbe_->book1D("TH1PFOEnergy","Energies of Pandora's PFOs",1000,0.,1000.);
  Epfos = new TH1F("TH1PFOEnergy","Energy of Reconstructed PFO",1000,0.,1000.);
  Egenpart = new TH1F("TH1GenParticlesEnergy","Energy of Generated Particles",1000,0.,1000.);
  Energy_res = new TH1F("TH1RES","#Delta E / E ",100,-2.,2.);

  h_sumCaloE = new TH1F("sumCaloE","sum hit E in Calos",1000,0,200);
  h_sumCaloEM = new TH1F("sumCaloEM","sum hit E in Calos",1000,0,200);
  h_sumCaloHad = new TH1F("sumCaloHad","sum hit E in Calos",1000,0,200);
  h_sumPfoE = new TH1F("hsumPfoE","sumPFOenergy",1000,0.,1000.);
  h_nbPFOs = new TH1F("hnbPfos","nb of rec PFOs",30,0.,30.);

  h2_Calo_EM_hcalEecalE = new TH2F("CalohcalEecalEem","",1000,0,200,1000,0,200);
  h2_Calo_Had_hcalEecalE = new TH2F("CalohcalEecalEhad","",1000,0,200,1000,0,200);


  h2_EM_hcalEecalE = new TH2F("hcalEecalEem","",400,0,40,400,0,40);
  h2_Had_hcalEecalE = new TH2F("hcalEecalEhad","",400,0,40,400,0,40);

  h_MIP_EE  = new TH1F("MIPEE" ,"Mip in EE ",1000,0,10);
  h_MIP_HEF = new TH1F("MIPHEF","Mip in HEF",1000,0,10);
  h_MIP_HEB = new TH1F("MIPHEB","Mip in HEB",1000,0,10);

  h_MIP_Corr_EE  = new TH1F("MIPCorrEE" ,"Mip corrected in EE ",1000,0,10);
  h_MIP_Corr_HEF = new TH1F("MIPCorrHEF","Mip corrected in HEF",1000,0,10);
  h_MIP_Corr_HEB = new TH1F("MIPCorrHEB","Mip corrected in HEB",1000,0,10);




  TH1::AddDirectory(oldAddDir);

  // read in calibration parameters
  initPandoraCalibrParameters();
  readCalibrParameterFile();


}


void runPandora::getLayerPropertiesEE (const HGCRecHit *eerh, int layer,
      float & nCellInteractionLengths, float & nCellRadiationLengths,
      float & absorberCorrectionEM, float & absorberCorrectionHAD
      )
{
   float absorberThickness_Pb[3] = { 1.6, 3.3, 5.6}; //mm
   float absorberThickness_Cu[3] = { 3. , 3. , 3. }; //mm

   float absorberRadLength_Pb   = 0.561 * 10; //mm
   float absorberRadLength_Cu   = 1.436 * 10; //mm

   float absorberIntLength_Pb = 17.59 * 10; //mm
   float absorberIntLength_Cu = 15.32 * 10; //mm

   float firstLayerCu = 3.0; //mm
   float firstLayerX0   = firstLayerCu / absorberRadLength_Cu;
   float firstLayerLdaI = firstLayerCu / absorberIntLength_Cu;

   int layerSet = 0;
   if ( 2 <= layer && layer <= 11 ) layerSet = 0;
   else if ( 12<= layer && layer <= 21 ) layerSet = 1;
   else if ( 22<= layer && layer <= 31 ) layerSet = 2;
   else {
      //std::cout << "OOUPSSSSSSSSSSSSSS, unknown layer!!!!" << std::endl;
      return;
   }

   float thisLayerX0 = absorberThickness_Pb[layerSet] / absorberRadLength_Pb
      + absorberThickness_Cu[layerSet] / absorberRadLength_Cu;
   float thisLayerLdaI = absorberThickness_Pb[layerSet] / absorberIntLength_Pb
      + absorberThickness_Cu[layerSet] / absorberIntLength_Cu;

   absorberCorrectionEM  = thisLayerX0   / firstLayerX0  ;
   absorberCorrectionHAD = thisLayerLdaI / firstLayerLdaI;

   nCellInteractionLengths = thisLayerLdaI;
   nCellRadiationLengths   = thisLayerX0  ;

   return;
}
void runPandora::getLayerPropertiesHEF (const HGCRecHit *eerh, int layer,
      float & nCellInteractionLengths, float & nCellRadiationLengths,
      float & absorberCorrectionEM, float & absorberCorrectionHAD
      )
{

   return;
}
void runPandora::getLayerPropertiesHEB (const HGCRecHit *eerh, int layer,
      float & nCellInteractionLengths, float & nCellRadiationLengths,
      float & absorberCorrectionEM, float & absorberCorrectionHAD
      )
{

   return;
}
// ------------ method called once each job just after ending the event loop  ------------
void 
runPandora::endJob() 
{
  file->Write();
  file->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
runPandora::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
runPandora::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
runPandora::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
runPandora::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
runPandora::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(runPandora);
