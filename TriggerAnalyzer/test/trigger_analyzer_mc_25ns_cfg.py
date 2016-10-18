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

#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v0' #'74X_mcRun2_asymptotic_v2' #'MCRUN2_74_V9'
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'#'74X_dataRun2_Express_v0'
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'

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

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
	'ak4PFCHSL1Fastjet', 
    'ak4PFchsL2Relative', 
    'ak4PFchsL3Absolute')
)

# HIP mitigted B-tagger

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None' ),  
  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
  runIVF=True,
  #btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
)

#process.source = cms.Source("PoolSource",
#        fileNames = cms.untracked.vstring())

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            #'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/10590823-AA0C-E511-A3BC-00259073E388.root',
            #'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/00759690-D16E-E511-B29E-00261894382D.root',
            #'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/06B7DC00-CA6D-E511-BBF3-002590DB925E.root',
            #'root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv1/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/30000/122A7134-C7A6-E511-8F6C-0CC47A78A436.root',
            #'/store/user/puigh/ttHTobb_M125_13TeV_powheg_pythia8_MINIAODSIM_PU25nsData2015v1_76X_mcRun2_asymptotic_v12_v1_30000_122A7134-C7A6-E511-8F6C-0CC47A78A436.root',
            #'/store/user/puigh/TTTo2L2Nu_13TeV-powheg_TTTo2L2Nu_13TeV-powheg_MiniAOD76_v2_160107_132443_0000_miniAOD-prod_PAT_1.root',
            #'root://xrootd.unl.edu//store/user/tarndt/TTTo2L2Nu_13TeV-powheg/TTTo2L2Nu_13TeV-powheg_MiniAOD76_v2/160107_132443/0000/miniAOD-prod_PAT_1.root',
#            'root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv1/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/30000/06984753-98A6-E511-8BC3-0025905A6070.root',
            #'root://xrootd.unl.edu//store/data/Run2015D/SingleElectron/MINIAOD/PromptReco-v3/000/256/630/00000/6E469C2A-165F-E511-9E77-02163E01414D.root',
           #'/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/049593B4-9E38-E611-832B-02163E014930.root' 
	   #'/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/0089CC67-6338-E611-947D-0025904C4E2A.root'
	   '/store/mc/RunIISpring16MiniAODv2/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/100000/086ED46A-1E76-E611-ABE5-180373FF8456.root'
	   )
)

### override the L1 menu from an Xml file
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
                                         PATsource = cms.untracked.string("PAT"),
                                         genTtbarId = cms.InputTag("categorizeGenTtbar", "genTtbarId"),
                                         isData = cms.bool(False),
                                         #electronMVAvalues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                                         #electronMVAcategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
                                         electronMVAvalues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                         electronMVAcategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
    )

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('trigger_analyzer.root')
)

process.p = cms.Path(process.electronMVAValueMapProducer * process.triggeranalyzer)
