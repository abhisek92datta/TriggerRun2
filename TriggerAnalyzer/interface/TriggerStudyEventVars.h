#ifndef TriggerRun2_TriggerAnalyzer_TriggerStudyEventVars_h
#define TriggerRun2_TriggerAnalyzer_TriggerStudyEventVars_h

//
// Dependencies (#includes)
//
#include <iostream>
#include <vector>
#include "TLorentzVector.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

using namespace std;



typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<double> vdouble;
typedef std::vector<float> vfloat;
typedef std::vector<string> vstring;
typedef std::vector<bool> vbool;
typedef std::vector<int> vint;
typedef std::vector< TLorentzVector > vecTLorentzVector;

//
// Utility Class for Handling Event Variables
//

struct triggerStudyEventVars{


  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  explicit triggerStudyEventVars() { }

  int run_;
  int lumi_;
  long evt_;

  double  wgt_generator_;

  int numJets_;
  int numTags_;
  int numPVs_;
  int numTruePVs_;
  int numGenPVs_;

  int additionalJetEventId_;
  int ttbarDecayMode_;

  bool goodFirstVertex_;

  /////

  double rho_;
  double top_pt_;
  double antitop_pt_;

  double qscale_;
  double pthat_;
  double originalXWGTUP_;

  vdouble LHEEvent_weights_;

  double lheHT_;

  /////
  int pass_L1_SingleEG25_;
  int pass_L1_SingleMu16_;

  int pass_AlCa_SingleEle_WPVeryLoose_Gsf_v_;
  int pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v_;
  int pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v_;
  int pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v_;
  int pass_HLT_Ele25_WPTight_Gsf_v_;
  int pass_HLT_Ele25_eta2p1_WPTight_Gsf_v_;
  

  int pass_HLT_Ele27_eta2p1_WP75_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v_;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v_;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v_;

  int pass_HLT_Ele25WP60_SC4_Mass55_v_;
  int pass_HLT_Ele27_WP85_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_;

  int pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v_;

  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;

  int pass_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_;
  int pass_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_;
  int pass_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_;
  int pass_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_;
  int pass_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_;

  int pass_HLT_IsoMu20_v_;
  int pass_HLT_IsoMu18_v_;
  int pass_HLT_IsoMu17_eta2p1_v_;
  int pass_HLT_IsoTkMu20_v_;
  int pass_HLT_Ele23_WPLoose_Gsf_v_;

  int pass_HLT_Ele22_eta2p1_WP75_Gsf_v_;
  int pass_HLT_PFHT450_SixJet40_PFBTagCSV0p72_v_;
  int pass_HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v_;
  int pass_HLT_PFHT450_SixJet40_PFBTagCSV_v_;
  int pass_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v_;



  /////

  // pfMET
  double pfMET_pt_;
  double pfMET_phi_;

  double pfMET_pt_JESup_;
  double pfMET_phi_JESup_;
  double pfMET_pt_JESdown_;
  double pfMET_phi_JESdown_;

  double pfMET_pt_JERup_;
  double pfMET_phi_JERup_;
  double pfMET_pt_JERdown_;
  double pfMET_phi_JERdown_;

  // pfMETNoHF
  double pfMETNoHF_pt_;
  double pfMETNoHF_phi_;

  double pfMETNoHF_pt_JESup_;
  double pfMETNoHF_phi_JESup_;
  double pfMETNoHF_pt_JESdown_;
  double pfMETNoHF_phi_JESdown_;

  double pfMETNoHF_pt_JERup_;
  double pfMETNoHF_phi_JERup_;
  double pfMETNoHF_pt_JERdown_;
  double pfMETNoHF_phi_JERdown_;

  // puppiMET
  double puppiMET_pt_;
  double puppiMET_phi_;

  double puppiMET_pt_JESup_;
  double puppiMET_phi_JESup_;
  double puppiMET_pt_JESdown_;
  double puppiMET_phi_JESdown_;

  double puppiMET_pt_JERup_;
  double puppiMET_phi_JERup_;
  double puppiMET_pt_JERdown_;
  double puppiMET_phi_JERdown_;

  // puppiMET raw
  double puppiMET_Upt_;
  double puppiMET_Uphi_;


  double L1HTT_;
  double L1HTT_bxm2_;
  double L1HTT_bxm1_;
  double L1HTT_bxp1_;
  double L1HTT_bxp2_;

  double PV_x_;
  double PV_y_;
  double PV_z_;

  vdouble jet_pt_;
  vdouble jet_eta_;
  vdouble jet_phi_;
  vdouble jet_energy_;
  vdouble jet_csv_;
  vdouble jet_cmva_;
  vint    jet_partonFlavour_;
  vint    jet_hadronFlavour_;
  vdouble jet_pileupJetId_fullDiscriminant_;

  vdouble jet_nocc_pt_;
  vdouble jet_nocc_eta_;
  vdouble jet_nocc_phi_;
  vdouble jet_nocc_energy_;
  vdouble jet_nocc_csv_;
  vdouble jet_nocc_cmva_;
  vint    jet_nocc_partonFlavour_;
  vint    jet_nocc_hadronFlavour_;
  vdouble jet_nocc_pileupJetId_fullDiscriminant_;


  vdouble jet_JESup_pt_;
  vdouble jet_JESup_eta_;
  vdouble jet_JESup_phi_;
  vdouble jet_JESup_energy_;
  vdouble jet_JESup_csv_;
  vdouble jet_JESup_cmva_;
  vint    jet_JESup_partonFlavour_;
  vint    jet_JESup_hadronFlavour_;
  vdouble jet_JESup_pileupJetId_fullDiscriminant_;

  vdouble jet_JESdown_pt_;
  vdouble jet_JESdown_eta_;
  vdouble jet_JESdown_phi_;
  vdouble jet_JESdown_energy_;
  vdouble jet_JESdown_csv_;
  vdouble jet_JESdown_cmva_;
  vint    jet_JESdown_partonFlavour_;
  vint    jet_JESdown_hadronFlavour_;
  vdouble jet_JESdown_pileupJetId_fullDiscriminant_;

  vdouble jet_JERup_pt_;
  vdouble jet_JERup_eta_;
  vdouble jet_JERup_phi_;
  vdouble jet_JERup_energy_;
  vdouble jet_JERup_csv_;
  vdouble jet_JERup_cmva_;
  vint    jet_JERup_partonFlavour_;
  vint    jet_JERup_hadronFlavour_;
  vdouble jet_JERup_pileupJetId_fullDiscriminant_;

  vdouble jet_JERdown_pt_;
  vdouble jet_JERdown_eta_;
  vdouble jet_JERdown_phi_;
  vdouble jet_JERdown_energy_;
  vdouble jet_JERdown_csv_;
  vdouble jet_JERdown_cmva_;
  vint    jet_JERdown_partonFlavour_;
  vint    jet_JERdown_hadronFlavour_;
  vdouble jet_JERdown_pileupJetId_fullDiscriminant_;


  double mass_leplep_;
  double dR_leplep_;
  Int_t  oppositeLepCharge_;

  vint lepton_trkCharge_;
  vint lepton_charge_;
  vint lepton_isMuon_;
  vint lepton_isTight_;
  vint lepton_isLoose_;
  vint lepton_is_IDTight_;
  vint lepton_isPhys14L_;
  vint lepton_isPhys14M_;
  vint lepton_isPhys14T_;
  vint lepton_isSpring15L_;
  vint lepton_isSpring15M_;
  vint lepton_isSpring15T_;
  vint lepton_isTrigMVAM_;
  vint lepton_isNonTrigMVAM_;
  vint lepton_isTrigCutM_;
  vint lepton_genId_;
  vint lepton_genParentId_;
  vint lepton_genGrandParentId_;
  vdouble lepton_pt_;
  vdouble lepton_eta_;
  vdouble lepton_phi_;
  vdouble lepton_energy_;
  vdouble lepton_relIso_;
  vdouble lepton_relIsoR04_;
  vdouble lepton_iso_sumChargedHadronPt_;
  vdouble lepton_iso_sumNeutralHadronEt_;
  vdouble lepton_iso_sumPhotonEt_;
  vdouble lepton_iso_sumPUPt_;
  vdouble lepton_trigMVAOutput_;
  vint    lepton_trigMVACategory_;
  vint    lepton_passTrigPresel_;
  vdouble lepton_scSigmaIEtaIEta_;
  vdouble lepton_full5x5_scSigmaIEtaIEta_;
  vdouble lepton_hadronicOverEm_;
  vdouble lepton_hcalOverEcal_;
  vdouble lepton_ecalPFClusterIso_;
  vdouble lepton_hcalPFClusterIso_;
  vdouble lepton_dr03TkSumPt_;
  vdouble lepton_relEcalIso_;
  vdouble lepton_relHcalIso_;
  vdouble lepton_relTrackIso_;
  vdouble lepton_OneOESuperMinusOneOP_;
  vint lepton_numMissingHits_;
  vint lepton_isEB_;
  vint lepton_passHLTId_;
  vint lepton_passConversionVeto_;
  vint lepton_inCrack_;
  vdouble lepton_scEta_;
  vdouble lepton_dEtaSCTrackAtVtx_;
  vdouble lepton_dPhiSCTrackAtVtx_;
  vdouble lepton_d0_;
  vdouble lepton_dZ_;
  vint lepton_isGlobalMuon_;
  vint lepton_isTrackerMuon_;
  vint lepton_isPFMuon_;
  vdouble lepton_normalizedChi2_;
  vint lepton_numberOfValidMuonHits_;
  vint lepton_numberOfValidPixelHits_;
  vint lepton_trackerLayersWithMeasurement_;
  vint lepton_numberOfMatchedStations_;

  vint lepton_ele_matchHLT_hltL1sL1SingleEG25_;
  vint lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_;

  vint lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25_;
  vint lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_;

  vint lepton_ele_matchHLT_hltL1sL1EG25erHTT125_;
  vint lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter_;
  vint lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter_;

  vdouble hltL1SingleEG25_pt_;
  vdouble hltL1SingleEG25_eta_;
  vdouble hltL1SingleEG25_phi_;

  vdouble hltEle27WP85Gsf_pt_;
  vdouble hltEle27WP85Gsf_eta_;
  vdouble hltEle27WP85Gsf_phi_;
  vstring hltEle27WP85Gsf_filter_;

  vdouble hltJet30_pt_;
  vdouble hltJet30_eta_;
  vdouble hltJet30_phi_;

  vdouble hltBtagJet30_pt_;
  vdouble hltBtagJet30_eta_;
  vdouble hltBtagJet30_phi_;

  vdouble hltPFHT200Jet30_pt_;
  vdouble hltPFHT200Jet30_eta_;
  vdouble hltPFHT200Jet30_phi_;
  vint    hltPFHT200Jet30_id_;

  vint    flt_accept_;
  vstring flt_name_;

  void initialize();

};


typedef std::vector<triggerStudyEventVars> vtriggerStudyEventVars;


void triggerStudyEventVars::initialize(){

  run_  = -99;
  lumi_ = -99;
  evt_ = -99;

  wgt_generator_ = -99.9;

  numJets_ = -99;
  numTags_ = -99;
  numPVs_  = -99;
  numTruePVs_ = -99;
  numGenPVs_ = -99;

  additionalJetEventId_ = -99;
  ttbarDecayMode_ = -99;
  
  goodFirstVertex_ = false;
  
  rho_ = -99;
  top_pt_ = -99;
  antitop_pt_ = -99;
  
  qscale_ = -99;
  pthat_ = -99;
  originalXWGTUP_ = -99;

  LHEEvent_weights_.clear();

  lheHT_ = -99;

  pass_L1_SingleEG25_ = -99;
  pass_L1_SingleMu16_ = -99;
  
  pass_AlCa_SingleEle_WPVeryLoose_Gsf_v_ = -99;
  pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v_ = -99;
  pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v_ = -99;
  pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v_ = -99;
  pass_HLT_Ele25_WPTight_Gsf_v_ = -99;
  pass_HLT_Ele25_eta2p1_WPTight_Gsf_v_ = -99;
  
  pass_HLT_Ele27_eta2p1_WP75_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v_ = -99;

  pass_HLT_Ele25WP60_SC4_Mass55_v_ = -99;
  pass_HLT_Ele27_WP85_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_ = -99;


  pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v_ = -99;

  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_ = -99;


  pass_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = -99;
  pass_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = -99;
  pass_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_ = -99;
  pass_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = -99;
  pass_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_ = -99;

  pass_HLT_IsoMu20_v_ = -99;
  pass_HLT_IsoMu18_v_ = -99;
  pass_HLT_IsoMu17_eta2p1_v_ = -99;
  pass_HLT_IsoTkMu20_v_ = -99;
  pass_HLT_Ele23_WPLoose_Gsf_v_ = -99;

  pass_HLT_Ele22_eta2p1_WP75_Gsf_v_ = -99;
  pass_HLT_PFHT450_SixJet40_PFBTagCSV0p72_v_ = -99;
  pass_HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v_ = -99;
  pass_HLT_PFHT450_SixJet40_PFBTagCSV_v_ = -99;
  pass_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v_ = -99;


  pfMET_pt_ = -99;
  pfMET_phi_ = -99;

  pfMET_pt_JESup_ = -99;
  pfMET_phi_JESup_ = -99;
  pfMET_pt_JESdown_ = -99;
  pfMET_phi_JESdown_ = -99;

  pfMET_pt_JERup_ = -99;
  pfMET_phi_JERup_ = -99;
  pfMET_pt_JERdown_ = -99;
  pfMET_phi_JERdown_ = -99;


  pfMETNoHF_pt_ = -99;
  pfMETNoHF_phi_ = -99;

  pfMETNoHF_pt_JESup_ = -99;
  pfMETNoHF_phi_JESup_ = -99;
  pfMETNoHF_pt_JESdown_ = -99;
  pfMETNoHF_phi_JESdown_ = -99;

  pfMETNoHF_pt_JERup_ = -99;
  pfMETNoHF_phi_JERup_ = -99;
  pfMETNoHF_pt_JERdown_ = -99;
  pfMETNoHF_phi_JERdown_ = -99;


  puppiMET_pt_ = -99;
  puppiMET_phi_ = -99;

  puppiMET_pt_JESup_ = -99;
  puppiMET_phi_JESup_ = -99;
  puppiMET_pt_JESdown_ = -99;
  puppiMET_phi_JESdown_ = -99;

  puppiMET_pt_JERup_ = -99;
  puppiMET_phi_JERup_ = -99;
  puppiMET_pt_JERdown_ = -99;
  puppiMET_phi_JERdown_ = -99;


  puppiMET_Upt_ = -99;
  puppiMET_Uphi_ = -99;


  L1HTT_ = -1;
  L1HTT_bxm2_ = -1;
  L1HTT_bxm1_ = -1;
  L1HTT_bxp1_ = -1;
  L1HTT_bxp2_ = -1;

  PV_x_ = -99;
  PV_y_ = -99;
  PV_z_ = -99;

  jet_pt_.clear();
  jet_eta_.clear();
  jet_phi_.clear();
  jet_energy_.clear();
  jet_csv_.clear();
  jet_cmva_.clear();
  jet_partonFlavour_.clear();
  jet_hadronFlavour_.clear();
  jet_pileupJetId_fullDiscriminant_.clear();


  jet_nocc_pt_.clear();
  jet_nocc_eta_.clear();
  jet_nocc_phi_.clear();
  jet_nocc_energy_.clear();
  jet_nocc_csv_.clear();
  jet_nocc_cmva_.clear();
  jet_nocc_partonFlavour_.clear();
  jet_nocc_hadronFlavour_.clear();
  jet_nocc_pileupJetId_fullDiscriminant_.clear();


  jet_JESup_pt_.clear();
  jet_JESup_eta_.clear();
  jet_JESup_phi_.clear();
  jet_JESup_energy_.clear();
  jet_JESup_csv_.clear();
  jet_JESup_cmva_.clear();
  jet_JESup_partonFlavour_.clear();
  jet_JESup_hadronFlavour_.clear();
  jet_JESup_pileupJetId_fullDiscriminant_.clear();

  jet_JESdown_pt_.clear();
  jet_JESdown_eta_.clear();
  jet_JESdown_phi_.clear();
  jet_JESdown_energy_.clear();
  jet_JESdown_csv_.clear();
  jet_JESdown_cmva_.clear();
  jet_JESdown_partonFlavour_.clear();
  jet_JESdown_hadronFlavour_.clear();
  jet_JESdown_pileupJetId_fullDiscriminant_.clear();

  jet_JERup_pt_.clear();
  jet_JERup_eta_.clear();
  jet_JERup_phi_.clear();
  jet_JERup_energy_.clear();
  jet_JERup_csv_.clear();
  jet_JERup_cmva_.clear();
  jet_JERup_partonFlavour_.clear();
  jet_JERup_hadronFlavour_.clear();
  jet_JERup_pileupJetId_fullDiscriminant_.clear();

  jet_JERdown_pt_.clear();
  jet_JERdown_eta_.clear();
  jet_JERdown_phi_.clear();
  jet_JERdown_energy_.clear();
  jet_JERdown_csv_.clear();
  jet_JERdown_cmva_.clear();
  jet_JERdown_partonFlavour_.clear();
  jet_JERdown_hadronFlavour_.clear();
  jet_JERdown_pileupJetId_fullDiscriminant_.clear();


  mass_leplep_ = -99;
  dR_leplep_ = -99;
  oppositeLepCharge_ = -99;


  lepton_trkCharge_.clear();
  lepton_charge_.clear();
  lepton_isMuon_.clear();
  lepton_isTight_.clear();
  lepton_isLoose_.clear();
  lepton_is_IDTight_.clear();
  lepton_isPhys14L_.clear();
  lepton_isPhys14M_.clear();
  lepton_isPhys14T_.clear();
  lepton_isSpring15L_.clear();
  lepton_isSpring15M_.clear();
  lepton_isSpring15T_.clear();
  lepton_isTrigMVAM_.clear();
  lepton_isNonTrigMVAM_.clear();
  lepton_isTrigCutM_.clear();
  lepton_genId_.clear();
  lepton_genParentId_.clear();
  lepton_genGrandParentId_.clear();
  lepton_pt_.clear();
  lepton_eta_.clear();
  lepton_phi_.clear();
  lepton_energy_.clear();
  lepton_relIso_.clear();
  lepton_relIsoR04_.clear();
  lepton_iso_sumChargedHadronPt_.clear();
  lepton_iso_sumNeutralHadronEt_.clear();
  lepton_iso_sumPhotonEt_.clear();
  lepton_iso_sumPUPt_.clear();
  lepton_trigMVAOutput_.clear();
  lepton_trigMVACategory_.clear();
  lepton_passTrigPresel_.clear();
  lepton_scSigmaIEtaIEta_.clear();
  lepton_full5x5_scSigmaIEtaIEta_.clear();
  lepton_hadronicOverEm_.clear();
  lepton_hcalOverEcal_.clear();
  lepton_ecalPFClusterIso_.clear();
  lepton_hcalPFClusterIso_.clear();
  lepton_dr03TkSumPt_.clear();
  lepton_relEcalIso_.clear();
  lepton_relHcalIso_.clear();
  lepton_relTrackIso_.clear();
  lepton_OneOESuperMinusOneOP_.clear();
  lepton_numMissingHits_.clear();
  lepton_isEB_.clear();
  lepton_passHLTId_.clear();
  lepton_passConversionVeto_.clear();
  lepton_inCrack_.clear();
  lepton_scEta_.clear();
  lepton_dEtaSCTrackAtVtx_.clear();
  lepton_dPhiSCTrackAtVtx_.clear();
  lepton_d0_.clear();
  lepton_dZ_.clear();
  lepton_isGlobalMuon_.clear();
  lepton_isTrackerMuon_.clear();
  lepton_isPFMuon_.clear();
  lepton_normalizedChi2_.clear();
  lepton_numberOfValidMuonHits_.clear();
  lepton_numberOfValidPixelHits_.clear();
  lepton_trackerLayersWithMeasurement_.clear();
  lepton_numberOfMatchedStations_.clear();


  lepton_ele_matchHLT_hltL1sL1SingleEG25_.clear();
  lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_.clear();

  lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25_.clear();
  lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_.clear();

  lepton_ele_matchHLT_hltL1sL1EG25erHTT125_.clear();
  lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter_.clear();
  lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter_.clear();


  hltL1SingleEG25_pt_.clear();
  hltL1SingleEG25_eta_.clear();
  hltL1SingleEG25_phi_.clear();

  hltEle27WP85Gsf_pt_.clear();
  hltEle27WP85Gsf_eta_.clear();
  hltEle27WP85Gsf_phi_.clear();
  hltEle27WP85Gsf_filter_.clear();

  hltJet30_pt_.clear();
  hltJet30_eta_.clear();
  hltJet30_phi_.clear();

  hltBtagJet30_pt_.clear();
  hltBtagJet30_eta_.clear();
  hltBtagJet30_phi_.clear();

  hltPFHT200Jet30_pt_.clear();
  hltPFHT200Jet30_eta_.clear();
  hltPFHT200Jet30_phi_.clear();
  hltPFHT200Jet30_id_.clear();

  flt_accept_.clear();
  flt_name_.clear();


  return;
}

  

#endif
