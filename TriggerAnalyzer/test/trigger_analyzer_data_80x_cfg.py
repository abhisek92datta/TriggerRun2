import FWCore.ParameterSet.Config as cms
import sys
import os

from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#process.GlobalTag.globaltag = 'MCRUN2_74_V9'
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8' ##'74X_dataRun2_v2'#'74X_dataRun2_Express_v0'
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.ak4PFCHSL1Fastjet = cms.ESProducer(
    'L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK4PFchs'),
    srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )

process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsResidual  = ak4CaloResidual.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
	'ak4PFCHSL1Fastjet', 
        'ak4PFchsL2Relative', 
        'ak4PFchsL3Absolute',
	'ak4PFchsResidual'
	)
)

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None' ),  
  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
  runIVF=True,
  #btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
 isData=True,
)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring())

#process.source = cms.Source("PoolSource",
##        fileNames = cms.untracked.vstring(
            #'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/10590823-AA0C-E511-A3BC-00259073E388.root',
            #'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02B810EA-7214-E511-BDAB-0025905964C2.root',
            #'root://xrootd.unl.edu//store/data/Run2015D/SingleElectron/MINIAOD/16Dec2015-v1/20002/8C404194-47A7-E511-8390-0CC47A4D7654.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/0A7BD549-131A-E611-8287-02163E0134FC.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/0C65F947-0F1A-E611-A1E6-02163E0144EA.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/1CCC1100-0E1A-E611-98C7-02163E014332.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/2209DFC5-191A-E611-9809-02163E014272.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/361F045C-111A-E611-AEF0-02163E01457C.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/3A8562C4-2B1A-E611-96AC-02163E0141D2.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/429D71B2-1D1A-E611-A5A9-02163E013926.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/4CF12A3D-0D1A-E611-9C13-02163E011FB9.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/546F0CFC-131A-E611-906F-02163E013942.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/56DDF2F4-0F1A-E611-AD12-02163E012065.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/58B8145B-1E1A-E611-BBF1-02163E01262F.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/5A806FDD-121A-E611-AC73-02163E013929.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/6014DA64-0E1A-E611-BC7F-02163E011E1A.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/662BEDED-091A-E611-9450-02163E0125FE.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/689AF1F8-0B1A-E611-B9D7-02163E0139B9.root',
#	   '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/6EA788B9-0A1A-E611-9786-02163E011AB8.root',
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/7C4306F0-1A1A-E611-AD6A-02163E011CA6.root',
#	    )
#)

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_unprescaled.txt.txt').getVLuminosityBlockRange()


## override the L1 menu from an Xml file
#process.l1GtTriggerMenuXml = cms.ESProducer("L1GtTriggerMenuXmlProducer",
#  TriggerMenuLuminosity = cms.string('startup'),
#  DefXmlFile = cms.string('L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml'),
#  VmeXmlFile = cms.string('')
#)
#process.L1GtTriggerMenuRcdSource = cms.ESSource("EmptyESSource",
#  recordName = cms.string('L1GtTriggerMenuRcd'),
#  iovIsRunNotTime = cms.bool(True),
#  firstValid = cms.vuint32(1)
#)
#process.es_prefer_l1GtParameters = cms.ESPrefer('L1GtTriggerMenuXmlProducer','l1GtTriggerMenuXml')

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

process.triggeranalyzer = cms.EDAnalyzer('TriggerAnalyzer',
                                         HLTsource = cms.untracked.string("HLT"),
                                         PATsource = cms.untracked.string("RECO"),
                                         genTtbarId = cms.InputTag("categorizeGenTtbar", "genTtbarId"),
                                         isData = cms.bool(True),
                                         electronMVAvalues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                                         electronMVAcategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
    )

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('trigger_analyzer.root')
)

process.p = cms.Path(process.electronMVAValueMapProducer * process.fullPatMetSequence * process.triggeranalyzer)
