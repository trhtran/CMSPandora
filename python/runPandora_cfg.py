import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi")

# No data of type "HGCalGeometry" with label "HGCalEESensitive" in record "IdealGeometryRecord"
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')

#The below three lines were added to solve an error Exception Message:
#No "CaloGeometryRecord" record found in the EventSetup.
### process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
### process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
### process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

#Add the next two lines to solve an error Exception Message:
#No "IdealMagneticFieldRecord" record found in the EventSetup.
### process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#Add the next three lines to solve an error Exception Message:
#No "TransientTrackRecord" record found in the EventSetup.
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
# process.GlobalTag.globaltag = 'START70_V1::All'

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.EveService = cms.Service("EveService")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/t/trtran/slc5/testAbsCorr/CMSSW_6_2_0_SLHC14/sim/gamma/eta2/Pt10/step3.root'
#        'file:/afs/cern.ch/work/t/trtran/slc5/detidfix/CMSSW_6_2_0_SLHC14/SIM/first/step3.root'
#        'file:/afs/cern.ch/work/t/trtran/slc5/testAbsCorr/CMSSW_6_2_0_SLHC14/sim/muon/eta1.5-3/first/step3.root'
    )
)

process.pandorapfanew = cms.EDAnalyzer('runPandora',
    ecalRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    hcalRecHitsHBHE = cms.InputTag("reducedHcalRecHits","hbheUpgradeReco"),
    HGCEErechitCollection  = cms.InputTag('HGCalRecHit','HGCEERecHits'), 
    HGCHEFrechitCollection = cms.InputTag('HGCalRecHit','HGCHEFRecHits'), 
    HGCHEBrechitCollection = cms.InputTag('HGCalRecHit','HGCHEBRecHits'), 
    generaltracks = cms.VInputTag(cms.InputTag("generalTracks")),
    tPRecoTrackAsssociation= cms.InputTag("trackingParticleRecoTrackAsssociation"),
    genParticles= cms.InputTag("genParticles"),
    inputconfigfile = cms.string('PandoraSettingsDefault.test.xml'),
#    inputconfigfile = cms.string('PandoraSettings_test.xml'),
#    inputconfigfile = cms.string('PandoraSettingsBasic_WithoutMonitoring.xml'),
#    inputconfigfile = cms.string('PandoraSettingsMuon.xml'),

    calibrParFile = cms.string('pandoraCalibrPars.txt-test'),
    outputFile = cms.string('pandoraoutput.root')

)


process.p = cms.Path(process.pandorapfanew)
