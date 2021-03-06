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
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Combined B to H
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'

#for B to G
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'
#for H
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'

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

#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#updateJetCollection(
#  process,
#  jetSource = cms.InputTag('slimmedJets'),
#  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None' ),
#  btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
#  runIVF=True,
  #btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
#)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
 isData=True,
)

process.load("RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff")

process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
#myModulesInPathToRun += process.badGlobalMuonTaggerMAOD
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
#myModulesInPathToRun += process.cloneGlobalMuonTaggerMAOD

#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring())

#process.source = cms.Source("PoolSource",
##        fileNames = cms.untracked.vstring(
#	    '/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/00000/00099863-E799-E611-A876-141877343E6D.root',
#	    )
#)

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_unprescaled.txt').getVLuminosityBlockRange()

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

process.triggeranalyzer = cms.EDAnalyzer('TriggerAnalyzer',
                                         HLTsource = cms.untracked.string("HLT"),
                                         PATsource = cms.untracked.string("RECO"),
                                         genTtbarId = cms.InputTag("categorizeGenTtbar", "genTtbarId"),
                                         isData = cms.bool(True),
                                         electronMVAvalues = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
					 electronMVAcategories = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                                         MET_filter_names = cms.untracked.vstring("Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter", "Flag_globalTightHalo2016Filter"),
                                         #badchcandfilter = cms.InputTag("BadChargedCandidateFilter"),
                                         #badpfmufilter = cms.InputTag("BadPFMuonFilter"),
                                         badglobalmuontagger = cms.InputTag("badGlobalMuonTaggerMAOD"),
                                         cloneglobalmuontagger = cms.InputTag("cloneGlobalMuonTaggerMAOD"),
    )

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('trigger_analyzer.root')
)

process.p = cms.Path(process.electronMVAValueMapProducer * process.fullPatMetSequence
                     #* process.BadPFMuonFilter
                     #* process.BadChargedCandidateFilter
                     * process.badGlobalMuonTaggerMAOD
                     * process.cloneGlobalMuonTaggerMAOD
                     * process.triggeranalyzer)
