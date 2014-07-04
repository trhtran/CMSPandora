// -*- C++ -*-
//
// Package:    PFCal/runPandora
// Class:      runPandora
// 
/**\class runPandora runPandora.cc PFCal/runPandora/plugins/runPandora.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Andreas Psallidas
//         Created:  Mon, 11 Nov 2013 15:11:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"
#include "Api/PandoraApi.h"
#include "TLorentzVector.h"
#include "Objects/ParticleFlowObject.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/Vector3D.h"
//DQM services for histogram
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"

//#include "PFCal/runPandora/interface/steerManager.h"

#include <TH1.h>
#include <TFile.h>

//
// class declaration
//
namespace reco {class PFRecHit;}

// namespace pandora {class Pandora;}

class runPandora : public edm::EDAnalyzer {
public:
  explicit runPandora(const edm::ParameterSet&);
  ~runPandora();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static pandora::Pandora        *m_pPandora;

  void prepareTrack(math::XYZVector B_,const reco::RecoToSimCollection pRecoToSim,const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void prepareHits(edm::Handle<EcalRecHitCollection> ecalRecHitHandleEB,edm::Handle<HBHERecHitCollection> hcalRecHitHandleHBHE,edm::Handle<HGCeeRecHitCollection> HGCeeRecHitHandle,edm::Handle<HGChefRecHitCollection> HGChefRecHitHandle,edm::Handle<HGChebRecHitCollection> HGChebRecHitHandle, reco::Vertex& pv, const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void preparemcParticle(edm::Handle<std::vector<reco::GenParticle> > genpart);

  void preparePFO(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void prepareGeometry(const edm::EventSetup& iSetup);

  TrackingParticleRefVector getTpSiblings(TrackingParticleRef tp);
  TrackingParticleRefVector getTpDaughters(TrackingParticleRef tp);

  std::string _outputFileName;
  std::string     m_pandoraSettingsXmlFile;

  std::string     m_calibrationParameterFile;
  void initPandoraCalibrParameters();
  void readCalibrParameterFile();
  void getLayerPropertiesEE (const HGCRecHit *eerh, int layer,
        float & nCellInteractionLengths, float & nCellRadiationLengths,
        float & absorberCorrectionEM, float & absorberCorrectionHAD
        );
  void getLayerPropertiesHEF (const HGCRecHit *eerh, int layer,
        float & nCellInteractionLengths, float & nCellRadiationLengths,
        float & absorberCorrectionEM, float & absorberCorrectionHAD
        );
  void getLayerPropertiesHEB (const HGCRecHit *eerh, int layer,
        float & nCellInteractionLengths, float & nCellRadiationLengths,
        float & absorberCorrectionEM, float & absorberCorrectionHAD
        );




private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  // ----------access to event data
  edm::InputTag    inputTagEcalRecHitsEB_ ; 
  edm::InputTag    inputTagHcalRecHitsHBHE_ ; 
  edm::InputTag    inputTagHGCEErechit_;
  edm::InputTag    inputTagHGCHEFrechit_;
  edm::InputTag    inputTagHGCHEBrechit_;
  edm::InputTag    inputTagtPRecoTrackAsssociation_;
  edm::InputTag    inputTagGenParticles_;
  std::vector<edm::InputTag>  inputTagGeneralTracks_;
  TFile * file;
  TH1F * Epfos;
  TH1F * Egenpart;
  TH1F * Energy_res;

  TH1F * h_sumPfoE;
  TH1F * h_nbPFOs;

  TH1F * h_sumCaloE;
  TH1F * h_MIP_EE ;
  TH1F * h_MIP_HEF;
  TH1F * h_MIP_HEB;
  TH1F * h_MIP_Corr_EE ;
  TH1F * h_MIP_Corr_HEF;
  TH1F * h_MIP_Corr_HEB;


  TH2F * h2_hcalEecalE;



  float m_Calibr_ADC2GeV_EE     ;
  float m_Calibr_ADC2GeV_HEF    ;
  float m_Calibr_ADC2GeV_HEB    ;

  float m_hCalThresBarrel       ;
  float m_hCalThresEndCapHEF    ;
  float m_hCalThresEndCapHEB    ;
  float m_eCalThresBarrel       ;
  float m_eCalThresEndCap       ;

  float m_hCalMipThresBarrel    ;
  float m_hCalMipThresEndCapHEF ;
  float m_hCalMipThresEndCapHEB ;
  float m_eCalMipThresBarrel    ;
  float m_eCalMipThresEndCap    ;

  float m_eCalToMipEndCap       ;
  float m_eCalToMipBarrel       ;
  float m_hCalToMipEndCapHEF    ;
  float m_hCalToMipEndCapHEB    ;
  float m_hCalToMipBarrel       ;

  float m_eCalToEMGeVEndCap     ;
  float m_eCalToEMGeVBarrel     ;
  float m_hCalToEMGeVEndCapHEF  ;
  float m_hCalToEMGeVEndCapHEB  ;
  float m_hCalToEMGeVBarrel     ;

  float m_eCalToHadGeVEndCap    ;
  float m_eCalToHadGeVBarrel    ;
  float m_hCalToHadGeVEndCapHEF ;
  float m_hCalToHadGeVEndCapHEB ;
  float m_hCalToHadGeVBarrel    ;

  float m_muonToMip             ;

  bool firstEvent_ ; 
  short _debugLevel;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
