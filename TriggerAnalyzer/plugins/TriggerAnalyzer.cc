// -*- C++ -*-
//
// Package:    TriggerRun2/TriggerAnalyzer
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc TriggerRun2/TriggerAnalyzer/plugins/TriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Darren Puigh
//         Created:  Wed, 02 Jul 2014 20:01:00 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerEvmReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom3.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "LHAPDF/LHAPDF.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"
#include "MiniAOD/MiniAODHelper/interface/LeptonSFHelper.h"

#include "TriggerRun2/TriggerAnalyzer/interface/TriggerStudyEventVars.h"


//
// class declaration
//

class TriggerAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TriggerAnalyzer(const edm::ParameterSet&);
      ~TriggerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const& iRun,edm::EventSetup const& iSetup) override;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> gtReadoutRecordToken;

  edm::EDGetTokenT <edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT <trigger::TriggerEvent> triggerEventToken;
  edm::EDGetTokenT <edm::TriggerResults> filterResultsToken;
  edm::EDGetTokenT <pat::TriggerObjectStandAloneCollection> triggerObjectsToken;

  edm::EDGetTokenT <reco::VertexCollection> vertexToken;
  edm::EDGetTokenT <reco::VertexCompositePtrCandidateCollection> secondaryVertexToken;
  edm::EDGetTokenT <edm::View<pat::Electron> > electronToken;
  edm::EDGetTokenT <pat::MuonCollection> muonToken;
  edm::EDGetTokenT <pat::JetCollection> jetToken;
  edm::EDGetTokenT <pat::METCollection> pfMetToken;
  edm::EDGetTokenT <pat::METCollection> pfMetNoHFToken;
  edm::EDGetTokenT <pat::METCollection> puppiMetToken;
  edm::EDGetTokenT<reco::GenJetCollection> token_genjets;

  edm::EDGetTokenT <pat::PackedCandidateCollection> packedpfToken;

  edm::EDGetTokenT <reco::BeamSpot> beamspotToken;
  edm::EDGetTokenT <reco::ConversionCollection> EDMConversionCollectionToken;
  edm::EDGetTokenT <double> rhoToken;
  edm::EDGetTokenT <reco::GenParticleCollection> mcparicleToken;
  //edm::EDGetTokenT <std::vector< PileupSummaryInfo > > puInfoToken;

  edm::EDGetTokenT<bool> BadChCandFilterToken;
  edm::EDGetTokenT<bool> BadPFMuonFilterToken;

  edm::EDGetTokenT <GenEventInfoProduct> genInfoProductToken;
  edm::EDGetTokenT <LHEEventProduct> lheEventProductToken;

  edm::EDGetTokenT <l1extra::L1EtMissParticleCollection> l1HTMissToken;

  edm::EDGetTokenT<edm::ValueMap<float> > eleMVAvaluesToken; // values of electron mva
  edm::EDGetTokenT<edm::ValueMap<int> >   eleMVAcategoriesToken;  // category of electron mva
  
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PU_info_Token;
  edm::EDGetTokenT<double> srcRho_Token;
  edm::EDGetTokenT<LHEEventProduct> lhep_Token;
  
  const edm::EDGetTokenT<int> genTtbarIdToken_;

  HLTConfigProvider hlt_config_;
  HLTConfigProvider filter_config_;

  std::string hltTag;
  std::string filterTag;

  std::map<std::string, int> l1talgo_cppath_;

  std::map<std::string, int> hlt_cppath_;
  std::map<std::string, int> flt_cppath_;

  bool verbose_;
  bool debug_;
  bool dumpHLT_;

  int numEvents_;


  int insample_;
  bool isData_;

  std::vector<std::string> hlt_triggerNames_;
  std::vector<std::string> flt_filterNames_;

  edm::Service<TFileService> fs_;

  TTree *m_ttree;
  triggerStudyEventVars *eve;


  double mySample_xSec_;
  double mySample_nGen_;
  double intLumi_;

  double ptmax;
  int NptBins;

  double jetptmax;
  int NjetptBins;

  int NjetMin;
  int NjetMax;
  int NjetBins;
  int NtagMin;
  int NtagMax;
  int NtagBins;
  int NpuMin;
  int NpuMax;
  int NpuBins;


  // Declare histograms
  
  TH1D *h_hlt;
  TH1D *h_flt;

  TH1D* h_numJet;
  TH1D* h_numTag;
  TH1D* h_numPV;
  TH1D* h_numTruePV;

  TH1D* h_numSecVtx;
  
  TRandom3 *r;

  // for JEC
  FactorizedJetCorrector *_jetCorrector_MC;
  FactorizedJetCorrector *_jetCorrector_BCD;
  FactorizedJetCorrector *_jetCorrector_EF;
  FactorizedJetCorrector *_jetCorrector_G;
  FactorizedJetCorrector *_jetCorrector_H;
  JetCorrectionUncertainty *_jetCorrectorUnc;

  // for JER
  std::vector<double>    JER_etaMin;
  std::vector<double>    JER_etaMax;
  std::vector<double>    JER_rhoMin;
  std::vector<double>    JER_rhoMax;
  std::vector<double>    JER_PtMin;
  std::vector<double>    JER_PtMax;
  std::vector<double>    JER_Par0;
  std::vector<double>    JER_Par1;
  std::vector<double>    JER_Par2;
  std::vector<double>    JER_Par3;

  MiniAODHelper miniAODhelper;
  LeptonSFHelper leptonSFhelper;
  
  // for PDF weight
  LHAPDF::PDFSet *PDF4LHC15_nlo_30_PDFSet;
  std::vector<LHAPDF::PDF *> _systPDFs;
  
  // for b-weights
  std::string inputFileHF;
  std::string inputFileLF;
  TFile *f_CSVwgt_HF;
  TFile *f_CSVwgt_LF;
  TH1D *h_csv_wgt_hf[9][6];
  TH1D *c_csv_wgt_hf[9][6];
  TH1D *h_csv_wgt_lf[9][4][3];
  
  double PU_x[100], PU_y[100];
  double PU_B_x[100], PU_B_y[100];
  double PU_C_x[100], PU_C_y[100];
  double PU_D_x[100], PU_D_y[100];
  double PU_E_x[100], PU_E_y[100];
  double PU_F_x[100], PU_F_y[100];
  double PU_G_x[100], PU_G_y[100];
  double PU_H_x[100], PU_H_y[100];
  double PU_BCDEF_x[100], PU_BCDEF_y[100];
  double PU_GH_x[100], PU_GH_y[100];


  std::vector<std::string> MET_filter_names_;
  inline bool Check_filters(edm::Handle<edm::TriggerResults>);

  inline void SetFactorizedJetCorrector(const sysType::sysType iSysType=sysType::NA);
  inline std::vector<pat::Jet> GetCorrectedJets(const std::vector<pat::Jet>&, const int &, const double &, const sysType::sysType iSysType=sysType::NA, const float& corrFactor = 1, const float& uncFactor = 1);
 
  inline double getPDFweight(const edm::Handle<GenEventInfoProduct> &, const int );
  inline double getPUweight(edm::Handle<std::vector<PileupSummaryInfo>>, const int );
  inline double getQ2weight(const edm::Handle<GenEventInfoProduct> &, const edm::Handle<LHEEventProduct> &, const string &);
  inline double getbweight( std::vector<pat::Jet> );
  inline void fillCSVHistos(TFile *, TFile *);
  inline double getCSVWeight(std::vector<double> , std::vector<double> , std::vector<double> , std::vector<int> , int , double &, double &, double &);
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig):
  //hltTag (iConfig.getParameter<string>("HLTsource")),
  //filterTag (iConfig.getParameter<string>("PATsource")),
  genTtbarIdToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarId"))),
  isData_(iConfig.getParameter<bool>("isData"))
  //MET_filter_names_(iConfig.getParameter<std::vector<string>>("MET_filter_names"))
{

  //now do what ever initialization is needed
  verbose_ = false;
  debug_ = false;
  dumpHLT_ = false;

  //hltTag = iConfig.getUntrackedParameter("HLTsource",std::string("HLT"));
  //filterTag = iConfig.getUntrackedParameter("PATsource",std::string("PAT"));
 
  hltTag = iConfig.getUntrackedParameter<string>("HLTsource");
  filterTag = iConfig.getUntrackedParameter<string>("PATsource");
  MET_filter_names_ = iConfig.getUntrackedParameter<std::vector<string>>("MET_filter_names");

  //for(int i=0; i<int(MET_filter_names_.size()); i++)
  //     std::cout<<MET_filter_names_[i]<<"\n";

  gtReadoutRecordToken = consumes <L1GlobalTriggerReadoutRecord> (edm::InputTag(std::string("gtDigis"), std::string(""), std::string("RECO")));

  triggerResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
  triggerEventToken = consumes <trigger::TriggerEvent> (edm::InputTag(std::string("hltTriggerSummaryAOD"), std::string(""), hltTag));
  filterResultsToken  = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), filterTag));
  triggerObjectsToken = consumes <pat::TriggerObjectStandAloneCollection> (edm::InputTag(std::string("selectedPatTrigger")));


  vertexToken = consumes <reco::VertexCollection> (edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
  secondaryVertexToken = consumes <reco::VertexCompositePtrCandidateCollection> (edm::InputTag(std::string("slimmedSecondaryVertices")));
  electronToken = consumes <edm::View<pat::Electron> > (edm::InputTag(std::string("slimmedElectrons")));
  muonToken = consumes <pat::MuonCollection> (edm::InputTag(std::string("slimmedMuons")));
  jetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJets")));
  //jetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("selectedUpdatedPatJets")));
  //pfMetToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETs")));
  pfMetToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETs"),std::string(""),std::string("MAOD")));
  //pfMetNoHFToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETsNoHF")));
  /// FIXME
  pfMetNoHFToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETs")));
  puppiMetToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETsPuppi")));
  token_genjets = consumes<reco::GenJetCollection>(edm::InputTag(std::string("slimmedGenJets")));

  BadChCandFilterToken = consumes<bool>(iConfig.getParameter<edm::InputTag>("badchcandfilter"));
  BadPFMuonFilterToken = consumes<bool>(iConfig.getParameter<edm::InputTag>("badpfmufilter"));

  packedpfToken = consumes <pat::PackedCandidateCollection> (edm::InputTag(std::string("packedPFCandidates")));

  beamspotToken = consumes <reco::BeamSpot> (edm::InputTag(std::string("offlineBeamSpot")));
  rhoToken = consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")));
  mcparicleToken = consumes <reco::GenParticleCollection> (edm::InputTag(std::string("prunedGenParticles")));
  //puInfoToken = consumes <std::vector< PileupSummaryInfo > > (edm::InputTag(std::string("slimmedAddPileupInfo")));
  EDMConversionCollectionToken = consumes <reco::ConversionCollection > (edm::InputTag("reducedEgamma","reducedConversions",""));

  genInfoProductToken = consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")));
  lheEventProductToken = consumes <LHEEventProduct> (edm::InputTag(std::string("externalLHEProducer")));


  l1HTMissToken = consumes <l1extra::L1EtMissParticleCollection> (edm::InputTag(std::string("l1extraParticles"), std::string("MHT")));

  eleMVAvaluesToken = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("electronMVAvalues"));
  eleMVAcategoriesToken = consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("electronMVAcategories"));
  
  //PU_info_Token = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup"));
  PU_info_Token = consumes <std::vector< PileupSummaryInfo > > (edm::InputTag(std::string("slimmedAddPileupInfo")));
  //srcRho_Token = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  //lhep_Token = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lhepprod"));

  ptmax = 500.;
  NptBins = int( ptmax/1. + 0.0001 );

  jetptmax = 1000.;
  NjetptBins = int( ptmax/1. + 0.0001 );

  NjetMin = 0;
  NjetMax = 15;
  NjetBins = NjetMax - NjetMin + 1;
  NtagMin = 0;
  NtagMax = 8;
  NtagBins = NtagMax - NtagMin + 1;
  NpuMin = 0;
  NpuMax = 50;
  NpuBins = NpuMax - NpuMin + 1;

  // book the TTree
  m_ttree = fs_->make<TTree>("triggerTree", "triggerTree");
  eve=0; 
  m_ttree->Branch("eve.", "triggerStudyEventVars", &eve, 8000, 1);


  h_numJet = fs_->make<TH1D>("h_numJet", ";Number of Jets", NjetBins, NjetMin-0.5, NjetMax+0.5 );
  h_numTag = fs_->make<TH1D>("h_numTag", ";Number of Tags", NtagBins, NtagMin-0.5, NtagMax+0.5 );
  h_numPV  = fs_->make<TH1D>("h_numPV", ";Number of Good Vertices", NpuBins, NpuMin-0.5, NpuMax+0.5 );
  h_numTruePV  = fs_->make<TH1D>("h_numTruePV", ";Number of True PU Vertices", NpuBins, NpuMin-0.5, NpuMax+0.5 );

  h_numSecVtx = fs_->make<TH1D>("h_numSecVtx", ";Number of Secondary Vertices", 6, 6-0.5, 6+0.5 );


  std::string era = "2012_53x";
  insample_ = 2500;

  mySample_xSec_ = 689.1;
  mySample_nGen_ = 25474122;
  intLumi_ = 20000;

  analysisType::analysisType iAnalysisType = analysisType::LJ;
  bool isData = isData_;

  miniAODhelper.SetUp(era, insample_, iAnalysisType, isData);

  r = new TRandom3(1);

  //miniAODhelper.SetJetCorrectorUncertainty();
  SetFactorizedJetCorrector();
  
  // Set up b-weights
  inputFileHF =
    	  "data/csv_weights/csv_rwt_fit_hf_v2_final_2017_1_10test.root";
    inputFileLF =
    	  "data/csv_weights/csv_rwt_fit_lf_v2_final_2017_1_10test.root";
    f_CSVwgt_HF = new TFile((inputFileHF).c_str());
    f_CSVwgt_LF = new TFile((inputFileLF).c_str());
    fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF);

  // Set up PU_weights
    ifstream fin;
    fin.open("data/PU_weight/PU_weight.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_x[i] >> PU_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_B.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_B_x[i] >> PU_B_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_C.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_C_x[i] >> PU_C_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_D.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_D_x[i] >> PU_D_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_E.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_E_x[i] >> PU_E_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_F.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_F_x[i] >> PU_F_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_G.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_G_x[i] >> PU_G_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_H.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_H_x[i] >> PU_H_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_BCDEF.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_BCDEF_x[i] >> PU_BCDEF_y[i];
    }
    fin.close();
    fin.open("data/PU_weight/PU_weight_GH.txt");
    for (int i = 0; i < 75; ++i) {
        fin >> PU_GH_x[i] >> PU_GH_y[i];
    }
    fin.close();

    // Set up PDF_weights
    PDF4LHC15_nlo_30_PDFSet = new LHAPDF::PDFSet("PDF4LHC15_nlo_30");
    _systPDFs = PDF4LHC15_nlo_30_PDFSet->mkPDFs();

}


TriggerAnalyzer::~TriggerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   r->SetSeed(0);
   delete PDF4LHC15_nlo_30_PDFSet;
}


//
// member functions
//

// for JEC

inline void TriggerAnalyzer::SetFactorizedJetCorrector(const sysType::sysType iSysType){


    std::vector<JetCorrectorParameters> corrParams_MC;
    std::vector<JetCorrectorParameters> corrParams_BCD;
    std::vector<JetCorrectorParameters> corrParams_EF;
    std::vector<JetCorrectorParameters> corrParams_G;
    std::vector<JetCorrectorParameters> corrParams_H;

    if (isData_) {

        JetCorrectorParameters *L3JetPar_BCD = new JetCorrectorParameters(
                                                                          "data/JEC/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar_BCD = new JetCorrectorParameters(
                                                                          "data/JEC/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar_BCD = new JetCorrectorParameters(
                                                                          "data/JEC/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt");
        JetCorrectorParameters *L2L3JetPar_BCD = new JetCorrectorParameters(
                                                                            "data/JEC/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt");

        JetCorrectorParameters *L3JetPar_EF = new JetCorrectorParameters(
                                                                         "data/JEC/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar_EF = new JetCorrectorParameters(
                                                                         "data/JEC/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar_EF = new JetCorrectorParameters(
                                                                         "data/JEC/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt");
        JetCorrectorParameters *L2L3JetPar_EF = new JetCorrectorParameters(
                                                                           "data/JEC/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt");

        JetCorrectorParameters *L3JetPar_G = new JetCorrectorParameters(
                                                                        "data/JEC/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar_G = new JetCorrectorParameters(
                                                                        "data/JEC/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar_G = new JetCorrectorParameters(
                                                                        "data/JEC/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt");
        JetCorrectorParameters *L2L3JetPar_G = new JetCorrectorParameters(
                                                                          "data/JEC/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt");

        JetCorrectorParameters *L3JetPar_H = new JetCorrectorParameters(
                                                                        "data/JEC/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar_H = new JetCorrectorParameters(
                                                                        "data/JEC/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar_H = new JetCorrectorParameters(
                                                                        "data/JEC/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt");
        JetCorrectorParameters *L2L3JetPar_H = new JetCorrectorParameters(
                                                                          "data/JEC/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt");

        corrParams_BCD.push_back(*L1JetPar_BCD);
        corrParams_BCD.push_back(*L2JetPar_BCD);
        corrParams_BCD.push_back(*L3JetPar_BCD);
        corrParams_BCD.push_back(*L2L3JetPar_BCD);
        _jetCorrector_BCD = new FactorizedJetCorrector(corrParams_BCD);
        corrParams_EF.push_back(*L1JetPar_EF);
        corrParams_EF.push_back(*L2JetPar_EF);
        corrParams_EF.push_back(*L3JetPar_EF);
        corrParams_EF.push_back(*L2L3JetPar_EF);
        _jetCorrector_EF = new FactorizedJetCorrector(corrParams_EF);
        corrParams_G.push_back(*L1JetPar_G);
        corrParams_G.push_back(*L2JetPar_G);
        corrParams_G.push_back(*L3JetPar_G);
        corrParams_G.push_back(*L2L3JetPar_G);
        _jetCorrector_G = new FactorizedJetCorrector(corrParams_G);
        corrParams_H.push_back(*L1JetPar_H);
        corrParams_H.push_back(*L2JetPar_H);
        corrParams_H.push_back(*L3JetPar_H);
        corrParams_H.push_back(*L2L3JetPar_H);
        _jetCorrector_H = new FactorizedJetCorrector(corrParams_H);

        std::string _JESUncFile =
        "data/JEC/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt";
        _jetCorrectorUnc = new JetCorrectionUncertainty(_JESUncFile);

        delete L3JetPar_BCD;
        delete L2JetPar_BCD;
        delete L1JetPar_BCD;
        delete L2L3JetPar_BCD;
        delete L3JetPar_EF;
        delete L2JetPar_EF;
        delete L1JetPar_EF;
        delete L2L3JetPar_EF;
        delete L3JetPar_G;
        delete L2JetPar_G;
        delete L1JetPar_G;
        delete L2L3JetPar_G;
        delete L3JetPar_H;
        delete L2JetPar_H;
        delete L1JetPar_H;
        delete L2L3JetPar_H;

    }
    
    else {
        JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(
                                                                      "data/JEC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(
                                                                      "data/JEC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(
                                                                      "data/JEC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt");
        
        corrParams_MC.push_back(*L1JetPar);
        corrParams_MC.push_back(*L2JetPar);
        corrParams_MC.push_back(*L3JetPar);
        _jetCorrector_MC = new FactorizedJetCorrector(corrParams_MC);
        
        std::string _JESUncFile =
        "data/JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt";
        _jetCorrectorUnc = new JetCorrectionUncertainty(_JESUncFile);
        
        delete L3JetPar;
        delete L2JetPar;
        delete L1JetPar;
    }
}

inline std::vector<pat::Jet> 
TriggerAnalyzer::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const int &run_nr, const double &rho, const sysType::sysType iSysType, const float& corrFactor , const float& uncFactor ){
	
  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    
    pat::Jet jet = (*it);
    double scale = 1.;

	// JEC

      if(!isData_) {   // MC
          _jetCorrector_MC->setJetPt(jet.pt());
          _jetCorrector_MC->setJetEta(jet.eta());
          _jetCorrector_MC->setJetA(jet.jetArea());
          _jetCorrector_MC->setRho(rho); //=fixedGridRhoFastjetAll
          scale = _jetCorrector_MC->getCorrection();
      }
      else {         // DATA
          if( run_nr>=272007 && run_nr<=276811 ) {
              _jetCorrector_BCD->setJetPt(jet.pt());
              _jetCorrector_BCD->setJetEta(jet.eta());
              _jetCorrector_BCD->setJetA(jet.jetArea());
              _jetCorrector_BCD->setRho(rho); //=fixedGridRhoFastjetAll
              scale = _jetCorrector_BCD->getCorrection();
          }
          else if( run_nr>=276831 && run_nr<=278801 ) {
              _jetCorrector_EF->setJetPt(jet.pt());
              _jetCorrector_EF->setJetEta(jet.eta());
              _jetCorrector_EF->setJetA(jet.jetArea());
              _jetCorrector_EF->setRho(rho); //=fixedGridRhoFastjetAll
              scale = _jetCorrector_EF->getCorrection();
          }
          else if( run_nr>=278802 && run_nr<=280385 ) {
              _jetCorrector_G->setJetPt(jet.pt());
              _jetCorrector_G->setJetEta(jet.eta());
              _jetCorrector_G->setJetA(jet.jetArea());
              _jetCorrector_G->setRho(rho); //=fixedGridRhoFastjetAll
              scale = _jetCorrector_G->getCorrection();
          }
          else if( run_nr>=280919 && run_nr<=284044 ) {
              _jetCorrector_H->setJetPt(jet.pt());
              _jetCorrector_H->setJetEta(jet.eta());
              _jetCorrector_H->setJetA(jet.jetArea());
              _jetCorrector_H->setRho(rho); //=fixedGridRhoFastjetAll
              scale = _jetCorrector_H->getCorrection();
          }
      }

      jet.scaleEnergy(scale);

    if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      _jetCorrectorUnc->setJetPt(jet.pt());
      _jetCorrectorUnc->setJetEta(jet.eta()); // here you must use the CORRECTED jet pt
      double unc = 1;
      double jes = 1;
      if( iSysType==sysType::JESup ){
	unc = _jetCorrectorUnc->getUncertainty(true);
	jes = 1 + unc;
      }
      else if( iSysType==sysType::JESdown ){
	unc = _jetCorrectorUnc->getUncertainty(false);
	jes = 1 - unc;
      }

      jet.scaleEnergy( jes );
    }
    outputJets.push_back(jet);
  }

  return outputJets;
}

inline bool TriggerAnalyzer::Check_filters(edm::Handle<edm::TriggerResults> filterResults)
{

    if (!filterResults.isValid()) {
        std::cerr << "Trigger results not valid for tag " << filterTag
        << std::endl;
        return 1;
    }


    bool pass = 1;
    for (std::vector<std::string>::const_iterator filter = MET_filter_names_.begin();
         filter != MET_filter_names_.end(); ++filter) {

        unsigned int filterIndex;
        std::string pathName = *filter;
        filterIndex = filter_config_.triggerIndex(pathName);
        if (filterIndex >= filterResults->size()){
            pass = pass*0;
            break;
        }
        if (filterResults->accept(filterIndex))
            pass=pass*1;
        else {
            pass = pass*0;
            break;
        }
    }

    return pass;
}

inline double TriggerAnalyzer::getPDFweight(const edm::Handle<GenEventInfoProduct> &genInfos, const int option)
{
    auto pdfInfos = genInfos->pdf();
    double pdfNominal = pdfInfos->xPDF.first * pdfInfos->xPDF.second;
   
    std::vector<double> pdfs;
    for (size_t j = 0; j < PDF4LHC15_nlo_30_PDFSet->size(); ++j) {
        double xpdf1 = _systPDFs[j]->xfxQ(pdfInfos->id.first, pdfInfos->x.first,
                                          pdfInfos->scalePDF);
        double xpdf2 = _systPDFs[j]->xfxQ(
            pdfInfos->id.second, pdfInfos->x.second, pdfInfos->scalePDF);
        pdfs.push_back(xpdf1 * xpdf2);
    }

    const LHAPDF::PDFUncertainty pdfUnc =
        PDF4LHC15_nlo_30_PDFSet->uncertainty(pdfs, 68.);

	double weight = 1.0;
    double weight_up = 1.0;
    double weight_down = 1.0;
    if (std::isfinite(1. / pdfNominal)) {
        weight = (pdfUnc.central) / pdfNominal;
        weight_up = (pdfUnc.central + pdfUnc.errplus) / pdfNominal;
        weight_down = (pdfUnc.central - pdfUnc.errminus) / pdfNominal;
    }
    if (option == 0)
    	return weight;
    else if (option == 1)
    	return weight_up;
    else
    	return weight_down;
}

inline double
TriggerAnalyzer::getPUweight(edm::Handle<std::vector<PileupSummaryInfo>> PupInfo, const int run)
{
    
    double pu_weight = -1;
    double numTruePV = -1;
    if ((PupInfo.isValid())) {
        for (std::vector<PileupSummaryInfo>::const_iterator PVI =
                 PupInfo->begin();
             PVI != PupInfo->end(); ++PVI) {
            int BX = PVI->getBunchCrossing();
            if (BX == 0) {
                numTruePV = PVI->getTrueNumInteractions();
            }
        }
    }

	if(run==0){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_x[i] + 1)) {
                pu_weight = PU_y[i];
                break;
            }
        }
    }
    else if(run==1){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_B_x[i] + 1)) {
                pu_weight = PU_B_y[i];
                break;
            }
        }
    }
    else if(run==2){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_C_x[i] + 1)) {
                pu_weight = PU_C_y[i];
                break;
            }
        }
    }
    else if(run==3){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_D_x[i] + 1)) {
                pu_weight = PU_D_y[i];
                break;
            }
        }
    }
    else if(run==4){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_E_x[i] + 1)) {
                pu_weight = PU_E_y[i];
                break;
            }
        }
    }
    else if(run==5){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_F_x[i] + 1)) {
                pu_weight = PU_F_y[i];
                break;
            }
        }
    }
    else if(run==6){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_G_x[i] + 1)) {
                pu_weight = PU_G_y[i];
                break;
            }
        }
    }
    else if(run==7){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_H_x[i] + 1)) {
                pu_weight = PU_H_y[i];
                break;
            }
        }
    }
    else if(run==8){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_BCDEF_x[i] + 1)) {
                pu_weight = PU_BCDEF_y[i];
                break;
            }
        }
    }
    else if(run==9){
        for (int i = 0; i < 75; ++i) {
            if (numTruePV < (PU_GH_x[i] + 1)) {
                pu_weight = PU_GH_y[i];
                break;
            }
        }
    }

    return pu_weight;
}

inline double
TriggerAnalyzer::getQ2weight(const edm::Handle<GenEventInfoProduct> &event_gen_info,
                        const edm::Handle<LHEEventProduct> &LHEEventProductHandle,
                        const string &ud)
{
    double theWeight;
    theWeight = event_gen_info->weight();
    unsigned int i;
    for (i = 0; i < LHEEventProductHandle->weights().size(); ++i) {
        if (!(ud.compare(LHEEventProductHandle->weights()[i].id)))
            theWeight *=
                LHEEventProductHandle->weights()[i].wgt / LHEEventProductHandle->originalXWGTUP();
    }
    return theWeight;
}

inline double TriggerAnalyzer::getbweight( std::vector<pat::Jet> selectedJets  )
{
    
    double csvWgtHF, csvWgtLF, csvWgtCF;
	double bweight = 1;

    std::vector<double> vec_jet_pt;
    std::vector<double> vec_jet_eta;
    std::vector<double> vec_jet_csv;
    std::vector<int> vec_jet_hadronFlavour;
    int iSys = 0;

        for (std::vector<pat::Jet>::const_iterator iJet =
                 selectedJets.begin();
             iJet != selectedJets.end(); ++iJet) {
            if( (iJet->pt()>30) && (fabs(iJet->eta())<2.4)  ){
            vec_jet_pt.push_back(iJet->pt());
            vec_jet_eta.push_back(iJet->eta());
            vec_jet_csv.push_back(miniAODhelper.GetJetCSV(
                *iJet, "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
            }
        }

	csvWgtHF = csvWgtLF = csvWgtCF = 0;
        bweight = getCSVWeight(vec_jet_pt, vec_jet_eta, vec_jet_csv,
                         vec_jet_hadronFlavour, iSys, csvWgtHF,
                         csvWgtLF, csvWgtCF);
	return bweight;
}

inline void TriggerAnalyzer::fillCSVHistos(TFile *fileHF, TFile *fileLF)
{
    for (int iSys = 0; iSys < 9; ++iSys) {
        for (int iPt = 0; iPt < 5; ++iPt)
            h_csv_wgt_hf[iSys][iPt] = NULL;
        for (int iPt = 0; iPt < 3; ++iPt) {
            for (int iEta = 0; iEta < 3; ++iEta)
                h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
        }
    }
    for (int iSys = 0; iSys < 5; ++iSys) {
        for (int iPt = 0; iPt < 5; ++iPt)
            c_csv_wgt_hf[iSys][iPt] = NULL;
    }

    // CSV reweighting /// only care about the nominal ones
    for (int iSys = 0; iSys < 9; ++iSys) {
        TString syst_csv_suffix_hf = "final";
        TString syst_csv_suffix_c = "final";
        TString syst_csv_suffix_lf = "final";

        switch (iSys) {
        case 0:
            // this is the nominal case
            break;
        case 1:
            // JESUp
            syst_csv_suffix_hf = "final_JESUp";
            syst_csv_suffix_lf = "final_JESUp";
            syst_csv_suffix_c = "final_cErr1Up";
            break;
        case 2:
            // JESDown
            syst_csv_suffix_hf = "final_JESDown";
            syst_csv_suffix_lf = "final_JESDown";
            syst_csv_suffix_c = "final_cErr1Down";
            break;
        case 3:
            // purity up
            syst_csv_suffix_hf = "final_LFUp";
            syst_csv_suffix_lf = "final_HFUp";
            syst_csv_suffix_c = "final_cErr2Up";
            break;
        case 4:
            // purity down
            syst_csv_suffix_hf = "final_LFDown";
            syst_csv_suffix_lf = "final_HFDown";
            syst_csv_suffix_c = "final_cErr2Down";
            break;
        case 5:
            // stats1 up
            syst_csv_suffix_hf = "final_Stats1Up";
            syst_csv_suffix_lf = "final_Stats1Up";
            break;
        case 6:
            // stats1 down
            syst_csv_suffix_hf = "final_Stats1Down";
            syst_csv_suffix_lf = "final_Stats1Down";
            break;
        case 7:
            // stats2 up
            syst_csv_suffix_hf = "final_Stats2Up";
            syst_csv_suffix_lf = "final_Stats2Up";
            break;
        case 8:
            // stats2 down
            syst_csv_suffix_hf = "final_Stats2Down";
            syst_csv_suffix_lf = "final_Stats2Down";
            break;
        }

        for (int iPt = 0; iPt < 5; ++iPt)
            h_csv_wgt_hf[iSys][iPt] = (TH1D *)fileHF->Get(
                Form("csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_hf.Data()));

        if (iSys < 5) {
            for (int iPt = 0; iPt < 5; ++iPt)
                c_csv_wgt_hf[iSys][iPt] = (TH1D *)fileHF->Get(Form(
                    "c_csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_c.Data()));
        }

        for (int iPt = 0; iPt < 4; ++iPt) {
            for (int iEta = 0; iEta < 3; ++iEta)
                h_csv_wgt_lf[iSys][iPt][iEta] =
                    (TH1D *)fileLF->Get(Form("csv_ratio_Pt%i_Eta%i_%s", iPt,
                                             iEta, syst_csv_suffix_lf.Data()));
        }
    }

    return;
}

inline double TriggerAnalyzer::getCSVWeight(std::vector<double> jetPts,
                                std::vector<double> jetEtas,
                                std::vector<double> jetCSVs,
                                std::vector<int> jetFlavors, int iSys,
                                double &csvWgtHF, double &csvWgtLF,
                                double &csvWgtCF)
{
    int iSysHF = 0;
    switch (iSys) {
    case 7:
        iSysHF = 1;
        break; // JESUp
    case 8:
        iSysHF = 2;
        break; // JESDown
    case 9:
        iSysHF = 3;
        break; // LFUp
    case 10:
        iSysHF = 4;
        break; // LFDown
    case 13:
        iSysHF = 5;
        break; // Stats1Up
    case 14:
        iSysHF = 6;
        break; // Stats1Down
    case 15:
        iSysHF = 7;
        break; // Stats2Up
    case 16:
        iSysHF = 8;
        break; // Stats2Down
    default:
        iSysHF = 0;
        break; // NoSys
    }

    int iSysC = 0;
    switch (iSys) {
    case 21:
        iSysC = 1;
        break;
    case 22:
        iSysC = 2;
        break;
    case 23:
        iSysC = 3;
        break;
    case 24:
        iSysC = 4;
        break;
    default:
        iSysC = 0;
        break;
    }

    int iSysLF = 0;
    switch (iSys) {
    case 7:
        iSysLF = 1;
        break; // JESUp
    case 8:
        iSysLF = 2;
        break; // JESDown
    case 11:
        iSysLF = 3;
        break; // HFUp
    case 12:
        iSysLF = 4;
        break; // HFDown
    case 17:
        iSysLF = 5;
        break; // Stats1Up
    case 18:
        iSysLF = 6;
        break; // Stats1Down
    case 19:
        iSysLF = 7;
        break; // Stats2Up
    case 20:
        iSysLF = 8;
        break; // Stats2Down
    default:
        iSysLF = 0;
        break; // NoSys
    }

    double csvWgthf = 1.;
    double csvWgtC = 1.;
    double csvWgtlf = 1.;

    for (int iJet = 0; iJet < int(jetPts.size()); ++iJet) {

        double csv = jetCSVs[iJet];
        double jetPt = jetPts[iJet];
        double jetAbsEta = fabs(jetEtas[iJet]);
        int flavor = jetFlavors[iJet];

        int iPt = -1;
        int iEta = -1;
        if (jetPt >= 19.99 && jetPt < 30)
            iPt = 0;
        else if (jetPt >= 30 && jetPt < 40)
            iPt = 1;
        else if (jetPt >= 40 && jetPt < 60)
            iPt = 2;
        else if (jetPt >= 60 && jetPt < 100)
            iPt = 3;
        else if (jetPt >= 100)
            iPt = 4;

        if (jetAbsEta >= 0 && jetAbsEta < 0.8)
            iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6)
            iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41)
            iEta = 2;

        if (iPt < 0 || iEta < 0)
            std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor "
                         "jet, jetPt = "
                      << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

        if (abs(flavor) == 5) {
            int useCSVBin =
                (csv >= 0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
            double iCSVWgtHF =
                h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
            if (iCSVWgtHF != 0)
                csvWgthf *= iCSVWgtHF;
        } else if (abs(flavor) == 4) {
            int useCSVBin =
                (csv >= 0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
            double iCSVWgtC =
                c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
            if (iCSVWgtC != 0)
                csvWgtC *= iCSVWgtC;
        } else {
            if (iPt >= 3)
                iPt =
                    3; /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
            int useCSVBin =
                (csv >= 0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
            double iCSVWgtLF =
                h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
            if (iCSVWgtLF != 0)
                csvWgtlf *= iCSVWgtLF;
        }
    }

    double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

    csvWgtHF = csvWgthf;
    csvWgtLF = csvWgtlf;
    csvWgtCF = csvWgtC;

    return csvWgtTotal;
}

// ------------ method called for each event  ------------
void
TriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  numEvents_++;

  int run  = iEvent.id().run();
  int lumi = iEvent.luminosityBlock();
  long evt = iEvent.id().event();

  eve->run_ = run;
  eve->lumi_ = lumi;
  eve->evt_ = evt;

  if( debug_ ) std::cout << " ====> test 0 " << std::endl;

  double minTightLeptonPt = 15.;
  double minLooseLeptonPt = 15.;

  double minLooseJetPt = 15.;

  if( debug_ ) std::cout << " ====> test 0.1 " << std::endl;

  
  h_hlt->Fill(0.,1);
  if( debug_ ) std::cout << " ====> test 0.2 " << std::endl;

  h_flt->Fill(0.,1);
  if( debug_ ) std::cout << " ====> test 0.3 " << std::endl;
  
  
  edm::Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;
  iEvent.getByToken(gtReadoutRecordToken, gtReadoutRecord);
  if( debug_ ) std::cout << " ====> test 0.4 " << std::endl;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults);
  if( debug_ ) std::cout << " ====> test 0.5 " << std::endl;

  edm::Handle<trigger::TriggerEvent> triggerEventAOD;
  iEvent.getByToken(triggerEventToken, triggerEventAOD);
  if( debug_ ) std::cout << " ====> test 0.6 " << std::endl;

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjectsToken, triggerObjects);

  if( debug_ ) std::cout << " ====> test 1 " << std::endl;

  // Get Trigger and Event Handles
 // edm::ESHandle<L1GtTriggerMenu> menuRcd;
  //iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  //const L1GtTriggerMenu* menu = menuRcd.product();

  // vint    l1t_accept;
  // vstring l1t_name;
  int pass_L1_SingleEG25 = -1;
  int pass_L1_SingleMu16 = -1;
  
  //cout<<"b";
  //bool passL1HTT100 = false;
  /*
  if( gtReadoutRecord.isValid() ){
    const DecisionWord& gtDecisionWord = gtReadoutRecord->decisionWord();
cout<<" c" ;
    // L1 algos
    for( CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
      int algoBitNumber = (algo->second).algoBitNumber();
      cout<<algoBitNumber<<" "<<gtDecisionWord.size();
cout<<"d";
      std::string name = (algo->second).algoName();
      cout<<" s";
      bool accept = gtDecisionWord.at(algoBitNumber);
      int prescaleFactor = -1;
      cout<<"e";

      //if( name=="L1_HTT100" && accept ) passL1HTT100 = true;

      if( name=="L1_SingleEG25" ) pass_L1_SingleEG25 = ( accept ) ? 1 : 0;
      if( name=="L1_SingleMu16" ) pass_L1_SingleMu16 = ( accept ) ? 1 : 0;

      // int is_accept = ( accept ) ? 1 : 0;
      // l1t_accept.push_back(is_accept);
      // l1t_name.push_back(name);
cout<<"f";
      if( accept ) l1talgo_cppath_[name]+=1;
      if( verbose_ ) std::cout << " =====>  L1T algo: path name = " << (algo->second).algoName() << ",\t prescale = " << prescaleFactor << ",\t pass = " << gtDecisionWord.at(algoBitNumber) << std::endl; 
    }
  }
  */
  //cout<<"22222 ";

  if( debug_ ) std::cout << " ====> test 2 " << std::endl;

  eve->pass_L1_SingleEG25_ = pass_L1_SingleEG25;
  eve->pass_L1_SingleMu16_ = pass_L1_SingleMu16;
  //eve->l1t_accept_ = l1t_accept;
  //eve->l1t_name_   = l1t_name;


  /*
  std::vector<TLorentzVector> triggerObjects_L1SingleEG25;
  std::vector<TLorentzVector> triggerObjects_hltEle;

  std::vector<TLorentzVector> triggerObjects_hltJet30;
  std::vector<TLorentzVector> triggerObjects_hltBtagJet30;
  std::vector<TLorentzVector> triggerObjects_hltPFHT200Jet30;

  std::vector<std::string> vec_hltEle27WP85Gsf_filter;

  std::vector<int> vec_hltPFHT200Jet30_id;
  */

  std::vector<TLorentzVector> triggerObjects_hlt_ele27wploose;
  std::vector<TLorentzVector> triggerObjects_hlt_ele27wptight;
  std::vector<TLorentzVector> triggerObjects_hlt_ele32wptight;
  std::vector<TLorentzVector> triggerObjects_hlt_ele27wplooseht200;


  if( triggerObjects.isValid() && triggerResults.isValid() ){
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
 
      // Print trigger object collection and type
      
      std::string collectionName = std::string( obj.collection() );

      TLorentzVector obj_TLV;
      obj_TLV.SetPxPyPzE( obj.px(), obj.py(), obj.pz(), obj.energy());

      //for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
        for (unsigned h = 0; h < obj.pathNames().size(); ++h){

            if( obj.pathNames()[h].find("HLT_Ele27_eta2p1_WPLoose_Gsf_v")!=std::string::npos )
                triggerObjects_hlt_ele27wploose.push_back(obj_TLV);
            if( obj.pathNames()[h].find("HLT_Ele27_eta2p1_WPTight_Gsf_v")!=std::string::npos )
                triggerObjects_hlt_ele27wptight.push_back(obj_TLV);
            if( obj.pathNames()[h].find("HLT_Ele32_eta2p1_WPTight_Gsf_v")!=std::string::npos )
                triggerObjects_hlt_ele32wptight.push_back(obj_TLV);
            if( obj.pathNames()[h].find("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v")!=std::string::npos )
                triggerObjects_hlt_ele27wplooseht200.push_back(obj_TLV);

    /*
	if( obj.filterLabels()[h]=="hltL1sL1SingleEG25" ) triggerObjects_L1SingleEG25.push_back(obj_TLV);

	if( obj.filterLabels()[h]=="hltL1sL1SingleEG25" ||
	    obj.filterLabels()[h]=="hltL1sL1SingleIsoEG22erOrSingleEG25" ||
	    obj.filterLabels()[h]=="hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20ORSingleIsoEG20erORSingleIsoEG22erORSingleEG25" ||
	    obj.filterLabels()[h]=="hltL1EG25Ele27WP85GsfTrackIsoFilter" ||
	    obj.filterLabels()[h]=="hltEle27WPLooseGsfTrackIsoFilter" ||
	    obj.filterLabels()[h]=="hltEle23WPLooseGsfTrackIsoFilter" ||
	    obj.filterLabels()[h]=="hltL1sL1SingleMu16" ||
	    obj.filterLabels()[h]=="hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09" ||
	    obj.filterLabels()[h].find("hltL1sL1EG25erHTT")!=std::string::npos ||
	    obj.filterLabels()[h]=="hltL1sL1EG25erHTT125" ||
	    obj.filterLabels()[h]=="hltL1EGHttEle27WPLooseGsfTrackIsoFilter" ||
	    obj.filterLabels()[h]=="hltL1EGHttEle27WP85GsfTrackIsoFilter" ){
	  triggerObjects_hltEle.push_back(obj_TLV);
	  vec_hltEle27WP85Gsf_filter.push_back(obj.filterLabels()[h]);
	}

	if( obj.filterLabels()[h]=="hltPFHT200Jet30" ){
	  triggerObjects_hltPFHT200Jet30.push_back(obj_TLV);
	  vec_hltPFHT200Jet30_id.push_back(obj.filterIds()[h]);
	}

	if( obj.filterLabels()[h]=="hltJetFilterSingleTopEle27" ) triggerObjects_hltJet30.push_back(obj_TLV);
	if( obj.filterLabels()[h]=="hltCSVFilterSingleTop" ) triggerObjects_hltBtagJet30.push_back(obj_TLV);

	if( (verbose_ && dumpHLT_) || (false && numEvents_<200) ) printf(" obj, filter %d: pt = %4.1f, filter = %s \n", h, obj.pt(), (obj.filterLabels()[h]).c_str());
    */
      }
    }
  }


    std::vector<double> vec_hlt_ele27wploose_pt;
    std::vector<double> vec_hlt_ele27wploose_eta;
    std::vector<double> vec_hlt_ele27wploose_phi;
    for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hlt_ele27wploose.size()); iHLTObj++ ){
        vec_hlt_ele27wploose_pt.push_back(triggerObjects_hlt_ele27wploose[iHLTObj].Pt());
        vec_hlt_ele27wploose_eta.push_back(triggerObjects_hlt_ele27wploose[iHLTObj].Eta());
        vec_hlt_ele27wploose_phi.push_back(triggerObjects_hlt_ele27wploose[iHLTObj].Phi());
    }
    eve->hltEle27WPLooseGsf_pt_  = vec_hlt_ele27wploose_pt;
    eve->hltEle27WPLooseGsf_eta_ = vec_hlt_ele27wploose_eta;
    eve->hltEle27WPLooseGsf_phi_ = vec_hlt_ele27wploose_phi;

    std::vector<double> vec_hlt_ele27wptight_pt;
    std::vector<double> vec_hlt_ele27wptight_eta;
    std::vector<double> vec_hlt_ele27wptight_phi;
    for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hlt_ele27wptight.size()); iHLTObj++ ){
        vec_hlt_ele27wptight_pt.push_back(triggerObjects_hlt_ele27wptight[iHLTObj].Pt());
        vec_hlt_ele27wptight_eta.push_back(triggerObjects_hlt_ele27wptight[iHLTObj].Eta());
        vec_hlt_ele27wptight_phi.push_back(triggerObjects_hlt_ele27wptight[iHLTObj].Phi());
    }
    eve->hltEle27WPTightGsf_pt_  = vec_hlt_ele27wptight_pt;
    eve->hltEle27WPTightGsf_eta_ = vec_hlt_ele27wptight_eta;
    eve->hltEle27WPTightGsf_phi_ = vec_hlt_ele27wptight_phi;

    std::vector<double> vec_hlt_ele32wptight_pt;
    std::vector<double> vec_hlt_ele32wptight_eta;
    std::vector<double> vec_hlt_ele32wptight_phi;
    for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hlt_ele32wptight.size()); iHLTObj++ ){
        vec_hlt_ele32wptight_pt.push_back(triggerObjects_hlt_ele32wptight[iHLTObj].Pt());
        vec_hlt_ele32wptight_eta.push_back(triggerObjects_hlt_ele32wptight[iHLTObj].Eta());
        vec_hlt_ele32wptight_phi.push_back(triggerObjects_hlt_ele32wptight[iHLTObj].Phi());
    }
    eve->hltEle32WPTightGsf_pt_  = vec_hlt_ele32wptight_pt;
    eve->hltEle32WPTightGsf_eta_ = vec_hlt_ele32wptight_eta;
    eve->hltEle32WPTightGsf_phi_ = vec_hlt_ele32wptight_phi;

    std::vector<double> vec_hlt_ele27wplooseht200_pt;
    std::vector<double> vec_hlt_ele27wplooseht200_eta;
    std::vector<double> vec_hlt_ele27wplooseht200_phi;
    for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hlt_ele27wplooseht200.size()); iHLTObj++ ){
        vec_hlt_ele27wplooseht200_pt.push_back(triggerObjects_hlt_ele27wplooseht200[iHLTObj].Pt());
        vec_hlt_ele27wplooseht200_eta.push_back(triggerObjects_hlt_ele27wplooseht200[iHLTObj].Eta());
        vec_hlt_ele27wplooseht200_phi.push_back(triggerObjects_hlt_ele27wplooseht200[iHLTObj].Phi());
    }
    eve->hltEle27WPLooseHT200Gsf_pt_  = vec_hlt_ele27wplooseht200_pt;
    eve->hltEle27WPLooseHT200Gsf_eta_ = vec_hlt_ele27wplooseht200_eta;
    eve->hltEle27WPLooseHT200Gsf_phi_ = vec_hlt_ele27wplooseht200_phi;

  /*
  std::vector<double> vec_hltL1SingleEG25_pt;
  std::vector<double> vec_hltL1SingleEG25_eta;
  std::vector<double> vec_hltL1SingleEG25_phi;
  for( int iL1TObj = 0; iL1TObj < int(triggerObjects_L1SingleEG25.size()); iL1TObj++ ){ 
    vec_hltL1SingleEG25_pt.push_back(triggerObjects_L1SingleEG25[iL1TObj].Pt());
    vec_hltL1SingleEG25_eta.push_back(triggerObjects_L1SingleEG25[iL1TObj].Eta());
    vec_hltL1SingleEG25_phi.push_back(triggerObjects_L1SingleEG25[iL1TObj].Phi());
  }
  eve->hltL1SingleEG25_pt_  = vec_hltL1SingleEG25_pt;
  eve->hltL1SingleEG25_eta_ = vec_hltL1SingleEG25_eta;
  eve->hltL1SingleEG25_phi_ = vec_hltL1SingleEG25_phi;

  std::vector<double> vec_hltEle27WP85Gsf_pt;
  std::vector<double> vec_hltEle27WP85Gsf_eta;
  std::vector<double> vec_hltEle27WP85Gsf_phi;
  for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hltEle.size()); iHLTObj++ ){ 
    vec_hltEle27WP85Gsf_pt.push_back(triggerObjects_hltEle[iHLTObj].Pt());
    vec_hltEle27WP85Gsf_eta.push_back(triggerObjects_hltEle[iHLTObj].Eta());
    vec_hltEle27WP85Gsf_phi.push_back(triggerObjects_hltEle[iHLTObj].Phi());
  }
  eve->hltEle27WP85Gsf_pt_  = vec_hltEle27WP85Gsf_pt;
  eve->hltEle27WP85Gsf_eta_ = vec_hltEle27WP85Gsf_eta;
  eve->hltEle27WP85Gsf_phi_ = vec_hltEle27WP85Gsf_phi;
  eve->hltEle27WP85Gsf_filter_ = vec_hltEle27WP85Gsf_filter;


  std::vector<double> vec_hltJet30_pt;
  std::vector<double> vec_hltJet30_eta;
  std::vector<double> vec_hltJet30_phi;
  for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hltJet30.size()); iHLTObj++ ){ 
    vec_hltJet30_pt.push_back(triggerObjects_hltJet30[iHLTObj].Pt());
    vec_hltJet30_eta.push_back(triggerObjects_hltJet30[iHLTObj].Eta());
    vec_hltJet30_phi.push_back(triggerObjects_hltJet30[iHLTObj].Phi());
  }
  eve->hltJet30_pt_  = vec_hltJet30_pt;
  eve->hltJet30_eta_ = vec_hltJet30_eta;
  eve->hltJet30_phi_ = vec_hltJet30_phi;


  std::vector<double> vec_hltBtagJet30_pt;
  std::vector<double> vec_hltBtagJet30_eta;
  std::vector<double> vec_hltBtagJet30_phi;
  for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hltBtagJet30.size()); iHLTObj++ ){ 
    vec_hltBtagJet30_pt.push_back(triggerObjects_hltBtagJet30[iHLTObj].Pt());
    vec_hltBtagJet30_eta.push_back(triggerObjects_hltBtagJet30[iHLTObj].Eta());
    vec_hltBtagJet30_phi.push_back(triggerObjects_hltBtagJet30[iHLTObj].Phi());
  }
  eve->hltBtagJet30_pt_  = vec_hltBtagJet30_pt;
  eve->hltBtagJet30_eta_ = vec_hltBtagJet30_eta;
  eve->hltBtagJet30_phi_ = vec_hltBtagJet30_phi;


  std::vector<double> vec_hltPFHT200Jet30_pt;
  std::vector<double> vec_hltPFHT200Jet30_eta;
  std::vector<double> vec_hltPFHT200Jet30_phi;
  for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hltPFHT200Jet30.size()); iHLTObj++ ){ 
    vec_hltPFHT200Jet30_pt.push_back(triggerObjects_hltPFHT200Jet30[iHLTObj].Pt());
    vec_hltPFHT200Jet30_eta.push_back(triggerObjects_hltPFHT200Jet30[iHLTObj].Eta());
    vec_hltPFHT200Jet30_phi.push_back(triggerObjects_hltPFHT200Jet30[iHLTObj].Phi());
  }
  eve->hltPFHT200Jet30_pt_  = vec_hltPFHT200Jet30_pt;
  eve->hltPFHT200Jet30_eta_ = vec_hltPFHT200Jet30_eta;
  eve->hltPFHT200Jet30_phi_ = vec_hltPFHT200Jet30_phi;
  eve->hltPFHT200Jet30_id_  = vec_hltPFHT200Jet30_id;
  */

  if( debug_ ) std::cout << " ====> test 3 " << std::endl;

  ////
  
  int pass_AlCa_SingleEle_WPVeryLoose_Gsf_v = -1;
  int pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v = -1;
  int pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v = -1;
  int pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v = -1;
  int pass_HLT_Ele25_WPTight_Gsf_v = -1;
  int pass_HLT_Ele25_eta2p1_WPTight_Gsf_v = -1;
  int pass_HLT_Ele32_eta2p1_WPTight_Gsf_v = -1;
  
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_v = -1;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v = -1;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v = -1;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v = -1;

  int pass_HLT_Ele25WP60_SC4_Mass55_v = -1;
  int pass_HLT_Ele27_WP85_Gsf_v = -1;
  int pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v = -1;


  int pass_HLT_Ele27_eta2p1_WPTight_Gsf_v = -1;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v = -1;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v = -1;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v = -1;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v = -1;

  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v = -1;


  int pass_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = -1;
  int pass_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = -1;
  
  int pass_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = -1;
  int pass_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v = -1;
  int pass_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = -1;
  int pass_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v = -1;

  int pass_HLT_IsoMu20_v = -1;
  int pass_HLT_IsoMu18_v = -1;
  int pass_HLT_IsoMu17_eta2p1_v = -1;
  int pass_HLT_IsoTkMu20_v = -1;
  int pass_HLT_Ele23_WPLoose_Gsf_v = -1;

  int pass_HLT_Ele22_eta2p1_WP75_Gsf_v = -1;
  int pass_HLT_PFHT450_SixJet40_PFBTagCSV0p72_v = -1;
  int pass_HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v = -1;
  int pass_HLT_PFHT450_SixJet40_PFBTagCSV_v = -1;
  int pass_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v = -1;


  //////

  if( triggerResults.isValid() ){
    std::vector<std::string> triggerNames = hlt_config_.triggerNames();

    for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
      std::string pathName = triggerNames[iPath];
      unsigned int hltIndex = hlt_config_.triggerIndex(pathName);

      if( hltIndex >= triggerResults->size() ) continue;

      int accept = triggerResults->accept(hltIndex);
      int prescale = -1;//hlt_config_.prescaleValue(iEvent, iSetup, pathName);

      if( verbose_ && dumpHLT_ ) std::cout << " =====>  HLT: path name = " << pathName << ",\t prescale = " << prescale << ",\t pass = " << accept << std::endl; 

      std::string pathNameNoVer = hlt_config_.removeVersion(pathName);

      if( accept ) hlt_cppath_[pathNameNoVer]+=1;

	 
      if( pathName.find("AlCa_SingleEle_WPVeryLoose_Gsf_v")!=std::string::npos )                          pass_AlCa_SingleEle_WPVeryLoose_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele22_eta2p1_WPLoose_Gsf_v")!=std::string::npos )          pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v = (accept) ? 1 : 0;
	  if( pathName.find("HLT_Ele24_eta2p1_WPLoose_Gsf_v")!=std::string::npos )                          pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele25_eta2p1_WPLoose_Gsf_v")!=std::string::npos )                          pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele25_WPTight_Gsf_v")!=std::string::npos )                          pass_HLT_Ele25_WPTight_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele25_eta2p1_WPTight_Gsf_v")!=std::string::npos )                          pass_HLT_Ele25_eta2p1_WPTight_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele32_eta2p1_WPTight_Gsf_v")!=std::string::npos )                          pass_HLT_Ele32_eta2p1_WPTight_Gsf_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_Ele27_eta2p1_WP75_Gsf_v")!=std::string::npos )                          pass_HLT_Ele27_eta2p1_WP75_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v")!=std::string::npos ) pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v")!=std::string::npos )        pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v")!=std::string::npos )  pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_Ele25WP60_SC4_Mass55_v")!=std::string::npos )        pass_HLT_Ele25WP60_SC4_Mass55_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_WP85_Gsf_v")!=std::string::npos )              pass_HLT_Ele27_WP85_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WP85_Gsf_HT200_v")!=std::string::npos ) pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_Ele27_eta2p1_WPTight_Gsf_v")!=std::string::npos )                          pass_HLT_Ele27_eta2p1_WPTight_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v")!=std::string::npos )                          pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v")!=std::string::npos ) pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v")!=std::string::npos )        pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v")!=std::string::npos )  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v")!=std::string::npos ) pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")!=std::string::npos ) pass_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")!=std::string::npos ) pass_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = (accept) ? 1 : 0;
      
      if( pathName.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")!=std::string::npos ) pass_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v")!=std::string::npos ) pass_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")!=std::string::npos ) pass_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")!=std::string::npos ) pass_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_IsoMu20_v")!=std::string::npos ) pass_HLT_IsoMu20_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_IsoMu18_v")!=std::string::npos ) pass_HLT_IsoMu18_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_IsoMu17_eta2p1_v")!=std::string::npos ) pass_HLT_IsoMu17_eta2p1_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_IsoTkMu20_v")!=std::string::npos ) pass_HLT_IsoTkMu20_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_Ele23_WPLoose_Gsf_v")!=std::string::npos ) pass_HLT_Ele23_WPLoose_Gsf_v = (accept) ? 1 : 0;

      if( pathName.find("HLT_Ele22_eta2p1_WP75_Gsf_v")!=std::string::npos ) pass_HLT_Ele22_eta2p1_WP75_Gsf_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_PFHT450_SixJet40_PFBTagCSV0p72_v")!=std::string::npos ) pass_HLT_PFHT450_SixJet40_PFBTagCSV0p72_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v")!=std::string::npos ) pass_HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_PFHT450_SixJet40_PFBTagCSV_v")!=std::string::npos ) pass_HLT_PFHT450_SixJet40_PFBTagCSV_v = (accept) ? 1 : 0;
      if( pathName.find("HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v")!=std::string::npos ) pass_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v = (accept) ? 1 : 0;

    
    if( accept ){
	TAxis * axis = h_hlt->GetXaxis();
	if( !axis ) continue;
	int bin_num = axis->FindBin(pathNameNoVer.c_str());
	int bn = bin_num - 1;
	h_hlt->Fill(bn, 1);
    }
     
    }
  }
  else{
    std::cout << "Trigger results not valid for tag " << hltTag << std::endl;
  }

  eve->pass_AlCa_SingleEle_WPVeryLoose_Gsf_v_ = pass_AlCa_SingleEle_WPVeryLoose_Gsf_v;
  eve->pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v_ = pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v;
  eve->pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v_ = pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v;
  eve->pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v_ = pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v;
  eve->pass_HLT_Ele25_WPTight_Gsf_v_ = pass_HLT_Ele25_WPTight_Gsf_v;
  eve->pass_HLT_Ele25_eta2p1_WPTight_Gsf_v_ = pass_HLT_Ele25_eta2p1_WPTight_Gsf_v;
  eve->pass_HLT_Ele32_eta2p1_WPTight_Gsf_v_ = pass_HLT_Ele32_eta2p1_WPTight_Gsf_v;

  eve->pass_HLT_Ele27_eta2p1_WP75_Gsf_v_                          = pass_HLT_Ele27_eta2p1_WP75_Gsf_v;
  eve->pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v_ = pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v;
  eve->pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v_        = pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v;
  eve->pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v_  = pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v;

  eve->pass_HLT_Ele25WP60_SC4_Mass55_v_        = pass_HLT_Ele25WP60_SC4_Mass55_v;
  eve->pass_HLT_Ele27_WP85_Gsf_v_              = pass_HLT_Ele27_WP85_Gsf_v;
  eve->pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_ = pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v;

  eve->pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_                          = pass_HLT_Ele27_eta2p1_WPTight_Gsf_v;
  eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_                          = pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v;
  eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v_ = pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v;
  eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v_        = pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v;
  eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v_  = pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v;

  eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_ = pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v;

  eve->pass_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = pass_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
  eve->pass_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = pass_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
  
  eve->pass_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = pass_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
  eve->pass_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_ = pass_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v;
  eve->pass_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = pass_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;
  eve->pass_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_ = pass_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;

  eve->pass_HLT_IsoMu20_v_ = pass_HLT_IsoMu20_v;
  eve->pass_HLT_IsoMu18_v_ = pass_HLT_IsoMu18_v;
  eve->pass_HLT_IsoMu17_eta2p1_v_ = pass_HLT_IsoMu17_eta2p1_v;
  eve->pass_HLT_IsoTkMu20_v_ = pass_HLT_IsoTkMu20_v;
  eve->pass_HLT_Ele23_WPLoose_Gsf_v_ = pass_HLT_Ele23_WPLoose_Gsf_v;

  eve->pass_HLT_Ele22_eta2p1_WP75_Gsf_v_ = pass_HLT_Ele22_eta2p1_WP75_Gsf_v;
  eve->pass_HLT_PFHT450_SixJet40_PFBTagCSV0p72_v_ = pass_HLT_PFHT450_SixJet40_PFBTagCSV0p72_v;
  eve->pass_HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v_ = pass_HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v;
  eve->pass_HLT_PFHT450_SixJet40_PFBTagCSV_v_ = pass_HLT_PFHT450_SixJet40_PFBTagCSV_v;
  eve->pass_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v_ = pass_HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v;


  if( debug_ ) std::cout << " ====> test 4 " << std::endl;



  if( debug_ ) std::cout << " ====> test 5 " << std::endl;



  if( debug_ ) std::cout << " ====> test 6 " << std::endl;

  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(filterResultsToken, filterResults);

  vint    flt_accept;
  vstring flt_name;
  if( filterResults.isValid() ){
    std::vector<std::string> triggerNames = filter_config_.triggerNames();

    for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
      std::string pathName = triggerNames[iPath];
      unsigned int hltIndex = filter_config_.triggerIndex(pathName);

      unsigned int filterResultsSize = filterResults->size();

      if( hltIndex >= filterResultsSize ) continue;

      int accept = filterResults->accept(hltIndex);

      int prescale = -1;//filter_config_.prescaleValue(iEvent, iSetup, pathName);

      if( verbose_ && dumpHLT_ ) std::cout << " =====>  Filter: path name = " << pathName << ",\t prescale = " << prescale << ",\t pass = " << accept << std::endl; 

      std::string pathNameNoVer = filter_config_.removeVersion(pathName);

      if( accept ) flt_cppath_[pathNameNoVer]+=1;

      flt_accept.push_back(accept);
      flt_name.push_back(pathName);

      if( accept ){
	TAxis * axis = h_flt->GetXaxis();
	if( !axis ) continue;
	int bin_num = axis->FindBin(pathNameNoVer.c_str());
	int bn = bin_num - 1;
	h_flt->Fill(bn, 1);
      }
    }
  }
  else{
    std::cout << "Trigger results not valid for tag " << filterTag << std::endl;
  }


  eve->flt_accept_ = flt_accept;
  eve->flt_name_   = flt_name;

  // MET Filters
  bool met_filters = Check_filters(filterResults);
  eve->met_filters = met_filters;

  bool filterbadPFMuon, filterbadChCandidate;
  filterbadChCandidate = filterbadPFMuon = 0;
  edm::Handle<bool> ifilterbadChCand;
  edm:: Handle<bool> ifilterbadPFMuon;
  iEvent.getByToken(BadChCandFilterToken, ifilterbadChCand);
  iEvent.getByToken(BadPFMuonFilterToken, ifilterbadPFMuon);
  filterbadChCandidate = *ifilterbadChCand;
  filterbadPFMuon = *ifilterbadPFMuon;
  eve->filterbadChCandidate = filterbadChCandidate;
  eve->filterbadPFMuon = filterbadPFMuon;

  if( debug_ ) std::cout << " ====> test 7 " << std::endl;


  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vertexToken,vtxHandle);
  reco::VertexCollection vtxs = *vtxHandle;

  edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVtxs;
  iEvent.getByToken(secondaryVertexToken,secondaryVtxs);

  //edm::Handle<pat::ElectronCollection> h_electrons;
  edm::Handle<edm::View<pat::Electron> > h_electrons;
  iEvent.getByToken(electronToken,h_electrons);

  edm::Handle<edm::ValueMap<float> > h_mvaValues; 
  iEvent.getByToken(eleMVAvaluesToken,h_mvaValues);

  edm::Handle<edm::ValueMap<int> > h_mvaCategories;
  iEvent.getByToken(eleMVAcategoriesToken,h_mvaCategories);
	
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken,muons);

  edm::Handle<pat::JetCollection> pfjets;
  iEvent.getByToken(jetToken,pfjets);

  edm::Handle<pat::METCollection> pfmet;
  iEvent.getByToken(pfMetToken,pfmet);

  edm::Handle<pat::METCollection> pfmetnohf;
  iEvent.getByToken(pfMetNoHFToken,pfmetnohf);

  edm::Handle<pat::METCollection> puppimet;
  iEvent.getByToken(puppiMetToken,puppimet);
 
  edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByToken(token_genjets, genjets);

  edm::Handle<pat::PackedCandidateCollection> packedPFcands;
  iEvent.getByToken(packedpfToken,packedPFcands);


  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(beamspotToken,bsHandle);

  edm::Handle<reco::ConversionCollection> h_conversioncollection;
  iEvent.getByToken( EDMConversionCollectionToken,h_conversioncollection );

  edm::Handle<reco::GenParticleCollection> mcparticles;
  iEvent.getByToken(mcparicleToken,mcparticles);

  /*
  double top_pt = -99;
  double antitop_pt = -99;
  if( mcparticles.isValid() ){
    for( size_t k = 0; k < mcparticles->size(); k++ ){
      const reco::Candidate & mcParticle = (*mcparticles)[k];
      int pdgId = mcParticle.pdgId();
      double pt = mcParticle.pt();
      if( pdgId==6 )       top_pt = pt;
      else if( pdgId==-6 ) antitop_pt = pt;
    }
  }

  eve->top_pt_ = top_pt;
  eve->antitop_pt_ = antitop_pt;
  */

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken,rhoHandle);
  ////------- set up rho for lepton effArea Isolation correction
  double rho_event = ( (rhoHandle.isValid()) ) ? *rhoHandle : -99;
  miniAODhelper.SetRho(rho_event);

  eve->rho_ = rho_event;

  //edm::Handle<std::vector< PileupSummaryInfo > > PupInfo;
  //iEvent.getByToken(puInfoToken,PupInfo);
 
  edm::Handle<std::vector< PileupSummaryInfo > > PupInfoHandle;
  iEvent.getByToken(PU_info_Token,PupInfoHandle);

  edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
  iEvent.getByToken(genInfoProductToken,GenEventInfoHandle);

  //edm::Handle<double> srcRhoHandle;
  //iEvent.getByToken(srcRho_Token, srcRhoHandle);

  double GenEventInfoWeight = 1.0;
  //double qScale=-99, pthat=-99;
  if( GenEventInfoHandle.isValid() ){
    GenEventInfoWeight = GenEventInfoHandle.product()->weight();
    //qScale = GenEventInfoHandle->qScale();
    //pthat = ( GenEventInfoHandle->hasBinningValues() ? (GenEventInfoHandle->binningValues())[0] : 0.0);
  }

  //eve->qscale_ = qScale;
  //eve->pthat_ = pthat;
  
  edm::Handle<LHEEventProduct> LHEEventProductHandle;
  iEvent.getByToken(lheEventProductToken,LHEEventProductHandle);
  
  //edm::Handle<LHEEventProduct> EvtHandle;
  //iEvent.getByToken(lhep_Token,EvtHandle);
  
  /*
  double originalXWGTUP=-99;
  vdouble LHEEvent_weights; 
  double lheHT = 0;
  if( LHEEventProductHandle.isValid() ){
    unsigned int NumEventWeights = LHEEventProductHandle->weights().size();
    originalXWGTUP = LHEEventProductHandle->originalXWGTUP();
    for( unsigned int iWgt=0; iWgt < NumEventWeights; iWgt++ ){
      double evtWgt = LHEEventProductHandle->weights()[iWgt].wgt;
      LHEEvent_weights.push_back(evtWgt);
    }

    for( int i = 0; i < LHEEventProductHandle->hepeup().NUP; ++i ){

      int id = LHEEventProductHandle->hepeup().IDUP[i];
      int status = LHEEventProductHandle->hepeup().ISTUP[i];
      int absId = abs(id);

      if( status==1 && ( absId==21 || (absId>0 && absId<7) ) ){
	double px = LHEEventProductHandle->hepeup().PUP[i][0];
	double py = LHEEventProductHandle->hepeup().PUP[i][1];

	lheHT += sqrt( px*px + py*py );
      }
    }
  }

  eve->originalXWGTUP_ = originalXWGTUP;
  eve->LHEEvent_weights_ = LHEEvent_weights;
  eve->lheHT_ = lheHT;
  */

  edm::Handle<int> genTtbarId;
  iEvent.getByToken(genTtbarIdToken_, genTtbarId);

  /*
  int additionalJetEventId = -99;
  if( genTtbarId.isValid() ) additionalJetEventId = *genTtbarId%100;
  eve->additionalJetEventId_ = additionalJetEventId;
  */

  //int ttbarDecayMode = -99;
  //MiniAODHelper::TTbarDecayMode ttbarDecayMode = ( isData_ ) ? MiniAODHelper::TTbarDecayMode::ChNotDefined : miniAODhelper.GetTTbarDecay(mcparticles);
  //eve->ttbarDecayMode_ = (int)ttbarDecayMode;

  edm::Handle<l1extra::L1EtMissParticleCollection> l1HTs;
  iEvent.getByToken(l1HTMissToken,l1HTs);


  if( debug_ ) std::cout << " ====> test 8 " << std::endl;


  math::XYZPoint beamSpotPosition;
  beamSpotPosition.SetCoordinates(0,0,0);
  double BSx=0,BSy=0,BSz=0;
  if( (bsHandle.isValid()) ){
    reco::BeamSpot bs = *bsHandle;
    BSx = bs.x0();
    BSy = bs.y0();
    BSz = bs.z0();
    beamSpotPosition = bsHandle->position();
  }


  if( verbose_ ) printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n", BSx, BSy, BSz );

  int numpv=0;
  reco::Vertex vertex;
  bool firstVertex = true;
  bool goodFirstVertex = false;
  if( vtxHandle.isValid() ){
    for( reco::VertexCollection::const_iterator vtx = vtxs.begin(); vtx!=vtxs.end(); ++vtx ){
      bool isGood = ( !(vtx->isFake()) &&
		      (vtx->ndof() >= 4.0) &&
		      (abs(vtx->z()) <= 24.0) &&
		      (abs(vtx->position().Rho()) <= 2.0) 
		      );

      if( firstVertex ){
	firstVertex = false;
	if( isGood ) goodFirstVertex = true;
      }

      if( !isGood ) continue;

      if( numpv==0 ) vertex = (*vtx);

      numpv++;
    }
  }

  eve->goodFirstVertex_ = ( goodFirstVertex );
  eve->PV_x_ = vertex.x();
  eve->PV_y_ = vertex.y();
  eve->PV_z_ = vertex.z();


  if( debug_ ) std::cout << " ====> test 9 " << std::endl;

  if( verbose_ ) printf("\t Event PV: x = %.3f,\t y = %.3f,\t z = %.3f \n", vertex.x(), vertex.y(), vertex.z());
  int numSecVtx = 0;
  for( reco::VertexCompositePtrCandidateCollection::const_iterator iVtx = secondaryVtxs->begin(); iVtx!=secondaryVtxs->end(); ++iVtx ){
    numSecVtx++;

    int numDgt = iVtx->numberOfDaughters();

    double dX = iVtx->vx() - vertex.x();
    double dY = iVtx->vy() - vertex.y();
    double dZ = iVtx->vz() - vertex.z();

    double dist = sqrt( dX*dX + dY*dY + dZ*dZ );

    if( verbose_ ) {
      printf("\t\t Secondary vtx %d: vx = %.3f,\t vy = %.3f,\t vz = %.3f,\t mass = %.1f,\t dist = %.2f,\t numDgt = %d \n", 
	     int(iVtx - secondaryVtxs->begin()), iVtx->vx(), iVtx->vy(), iVtx->vz(), iVtx->mass(), dist, numDgt);

      for( int iDau=0; iDau<numDgt; iDau++ ){
	printf("\t\t\t daughter %d: pdgId = %d,\t pt = %.1f,\t eta = %.2f,\t phi = %.2f \n", 
	       iDau, iVtx->daughter(iDau)->pdgId(), iVtx->daughter(iDau)->pt(), iVtx->daughter(iDau)->eta(), iVtx->daughter(iDau)->phi());
      }
    }
  }

  if( numpv>0 ) miniAODhelper.SetVertex(vertex);

  double numTruePV = -1;
  double numGenPV = -1;
  if( (PupInfoHandle.isValid()) ){
    for( std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfoHandle->begin(); PVI != PupInfoHandle->end(); ++PVI ) {
      int BX = PVI->getBunchCrossing();
      if( BX==0 ){
	numTruePV = PVI->getTrueNumInteractions();
	numGenPV  = PVI->getPU_NumInteractions();
      }
    }
  }



  if( debug_ ) std::cout << " ====> test 10 " << std::endl;

 
     //Use event setup for jet corrector 
  const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );   //Get the jet corrector from the event setup

  miniAODhelper.SetJetCorrector(corrector);
  miniAODhelper.SetJetCorrectorUncertainty(iSetup);

  if( debug_ ) std::cout << " ====> test 11 " << std::endl;


  double wgt = 1;

  double xSec = mySample_xSec_;
  double nGen = mySample_nGen_;

  wgt = intLumi_ * xSec / nGen;


  eve->wgt_generator_ = GenEventInfoWeight;

  eve->numPVs_ = numpv;
  eve->numTruePVs_ = numTruePV;
  eve->numGenPVs_ = numGenPV;

  /*
  double L1HTT = -999;
  double L1HTT_bxm2 = -999;
  double L1HTT_bxm1 = -999;
  double L1HTT_bxp1 = -999;
  double L1HTT_bxp2 = -999;
  */
  /*
  if( l1HTs.isValid() ){
    //L1HTT = l1HTs->front().etTotal();
    //std::cout << " use event L1HTT = " << L1HTT << "\t passL1HTT100 = " << passL1HTT100 << std::endl;
    for( l1extra::L1EtMissParticleCollection::const_iterator iL1HT = l1HTs->begin(); iL1HT!=l1HTs->end(); ++iL1HT ){
      double etTotal = iL1HT->etTotal();
      int bx = iL1HT->bx();

      if( bx==0 && L1HTT < -1 ) L1HTT = etTotal;

      if( bx==-2 )      L1HTT_bxm2 = etTotal;
      else if( bx==-1 ) L1HTT_bxm1 = etTotal;
      else if( bx==1  ) L1HTT_bxp1 = etTotal;
      else if( bx==2  ) L1HTT_bxp2 = etTotal;

      //printf(" iL1HT %d: etTotal = %4.1f \t bx = %d \n", int(iL1HT - l1HTs->begin()), etTotal, bx);
    }
  }
  else{ 
    std::cout << " l1HTs Handle not valid!! " << std::endl;
  }
  */
  /*
  eve->L1HTT_ = L1HTT;
  eve->L1HTT_bxm2_ = L1HTT_bxm2;
  eve->L1HTT_bxm1_ = L1HTT_bxm1;
  eve->L1HTT_bxp1_ = L1HTT_bxp1;
  eve->L1HTT_bxp2_ = L1HTT_bxp2;
  */

  if( debug_ ) std::cout << " ====> test 12 " << std::endl;

  /////////
  ///
  /// Electrons
  ///
  ////////

  //miniAODhelper.SetElectronMVAinfo(h_conversioncollection, bsHandle);

  std::vector<pat::Electron> electrons = miniAODhelper.GetElectronsWithMVAid( h_electrons, h_mvaValues, h_mvaCategories );
    
  // std::vector<pat::Electron> selectedElectrons_tight = miniAODhelper.GetSelectedElectrons( *electrons, minTightLeptonPt, electronID::electronSpring15M, 2.1 );
  //std::vector<pat::Electron> selectedElectrons_loose = miniAODhelper.GetSelectedElectrons( *electrons, minLooseLeptonPt, electronID::electronEndOf15MVAmedium, 2.4 );
  //std::vector<pat::Electron> selectedElectrons_loose = miniAODhelper.GetSelectedElectrons( electrons, minLooseLeptonPt, electronID::electronEndOf15MVA80, 2.4 );
  std::vector<pat::Electron> selectedElectrons_loose = miniAODhelper.GetSelectedElectrons( electrons, minLooseLeptonPt, electronID::electronNonTrigMVAid80, 2.4 );

  // int numTightElectrons = int(selectedElectrons_tight.size());
  // int numLooseElectrons = int(selectedElectrons_loose.size());// - numTightElectrons;

  if( debug_ ) std::cout << " ====> test 13 " << std::endl;


  /////////
  ///
  /// Muons
  ///
  ////////
  // std::vector<pat::Muon> selectedMuons_tight = miniAODhelper.GetSelectedMuons( *muons, minTightLeptonPt, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.1 );
  std::vector<pat::Muon> selectedMuons_loose = miniAODhelper.GetSelectedMuons( *muons, minLooseLeptonPt, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.4 );

  // int numTightMuons = int(selectedMuons_tight.size());
  // int numLooseMuons = int(selectedMuons_loose.size());// - numTightMuons;


  if( debug_ ) std::cout << " ====> test 14 " << std::endl;



  //////
  // Do jets stuff
  std::vector<pat::Jet> pfJets_ID = miniAODhelper.GetSelectedJets(*pfjets,0.,999,jetID::jetTight,'-');
  //std::vector<pat::Jet> pfJets_ID = miniAODhelper.GetSelectedJets(*pfjets,0.,999,jetID::jetLoose,'-');
  std::vector<pat::Jet> pfJets_ID_clean = miniAODhelper.GetDeltaRCleanedJets( pfJets_ID, selectedMuons_loose, selectedElectrons_loose, 0.4);
  std::vector<pat::Jet> rawJets = miniAODhelper.GetUncorrectedJets(pfJets_ID_clean);
  // Use JEC from GT
  bool doJER;
  if(isData_)
  	doJER = 0;
  else 
  	doJER = 1;
  // JEC
  std::vector<pat::Jet> correctedJets_noSys = GetCorrectedJets(rawJets, run, rho_event, sysType::NA);
  // JER
  correctedJets_noSys = miniAODhelper.GetCorrectedJets(correctedJets_noSys, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);

  std::vector<pat::Jet> selectedJets_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_noSys, minLooseJetPt, 3.0, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_tag_noSys_unsorted = miniAODhelper.GetSelectedJets( correctedJets_noSys, minLooseJetPt, 3.0, jetID::none, 'M' );


  //////
  if( debug_ ) std::cout << " ====> test 15 " << std::endl;

  // Sort jet collections by pT
  std::vector<pat::Jet> selectedJets       = miniAODhelper.GetSortedByPt( selectedJets_noSys_unsorted );
  std::vector<pat::Jet> selectedJets_tag   = miniAODhelper.GetSortedByPt( selectedJets_tag_noSys_unsorted );


  // Get numJets, numTags
  int numJet = int( selectedJets.size() );
  int numTag = int( selectedJets_tag.size() );

  eve->numJets_ = numJet;
  eve->numTags_ = numTag;

  /*
  /// jets without cc
  std::vector<pat::Jet> rawJets_nocc = miniAODhelper.GetUncorrectedJets(pfJets_ID);
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> correctedJets_nocc_noSys = GetCorrectedJets(rawJets_nocc, run, rho_event, sysType::NA);
  // JER
  correctedJets_nocc_noSys = miniAODhelper.GetCorrectedJets(correctedJets_nocc_noSys, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);

  std::vector<pat::Jet> selectedJets_nocc_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_nocc_noSys, minLooseJetPt, 3.0, jetID::none, '-' );

  //////
  if( debug_ ) std::cout << " ====> test 15 " << std::endl;

  // Sort jet collections by pT
  std::vector<pat::Jet> selectedJets_nocc = miniAODhelper.GetSortedByPt( selectedJets_nocc_noSys_unsorted );


  // Get Corrected MET (propagating JEC and JER)
  // pat::MET correctedMET = pfmet->front();//miniAODhelper.GetCorrectedMET( pfmets.at(0), pfJets_forMET, iSysType );
  std::vector<pat::Jet> oldJetsForMET = miniAODhelper.GetSelectedJets(*pfjets, 0., 999., jetID::jetMETcorrection, '-' );
  std::vector<pat::Jet> oldJetsForMET_uncorr = miniAODhelper.GetUncorrectedJets(oldJetsForMET);
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> newJetsForMET = GetCorrectedJets(oldJetsForMET_uncorr, run, rho_event, sysType::NA);
  // JER
  newJetsForMET = miniAODhelper.GetCorrectedJets(newJetsForMET, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);


  std::vector<pat::MET> newPfMETs = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET, *pfmet);
  std::vector<pat::MET> newPfMETsNoHF = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET, *pfmetnohf);
  std::vector<pat::MET> newPuppiMETs = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET, *puppimet);

  eve->pfMET_pt_  = newPfMETs.at(0).pt();
  eve->pfMET_phi_ = newPfMETs.at(0).phi();

  eve->pfMETNoHF_pt_  = newPfMETsNoHF.at(0).pt();
  eve->pfMETNoHF_phi_ = newPfMETsNoHF.at(0).phi();

  eve->puppiMET_pt_  = newPuppiMETs.at(0).pt();
  eve->puppiMET_phi_ = newPuppiMETs.at(0).phi();

  eve->puppiMET_Upt_  = puppimet->front().pt();
  eve->puppiMET_Uphi_ = puppimet->front().phi();


  // JESup
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> newJetsForMET_JESup = GetCorrectedJets(oldJetsForMET_uncorr, run, rho_event, sysType::JESup);
  // JER
  newJetsForMET_JESup = miniAODhelper.GetCorrectedJets(newJetsForMET_JESup, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);

  std::vector<pat::MET> newPfMETs_JESup = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JESup, *pfmet);
  std::vector<pat::MET> newPfMETsNoHF_JESup = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JESup, *pfmetnohf);
  std::vector<pat::MET> newPuppiMETs_JESup = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JESup, *puppimet);

  eve->pfMET_pt_JESup_  = newPfMETs_JESup.at(0).pt();
  eve->pfMET_phi_JESup_ = newPfMETs_JESup.at(0).phi();

  eve->pfMETNoHF_pt_JESup_  = newPfMETsNoHF_JESup.at(0).pt();
  eve->pfMETNoHF_phi_JESup_ = newPfMETsNoHF_JESup.at(0).phi();

  eve->puppiMET_pt_JESup_  = newPuppiMETs_JESup.at(0).pt();
  eve->puppiMET_phi_JESup_ = newPuppiMETs_JESup.at(0).phi();

  // JESdown
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> newJetsForMET_JESdown = GetCorrectedJets(oldJetsForMET_uncorr, run, rho_event, sysType::JESdown);
  // JER
  newJetsForMET_JESdown = miniAODhelper.GetCorrectedJets(newJetsForMET_JESdown, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);

  std::vector<pat::MET> newPfMETs_JESdown = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JESdown, *pfmet);
  std::vector<pat::MET> newPfMETsNoHF_JESdown = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JESdown, *pfmetnohf);
  std::vector<pat::MET> newPuppiMETs_JESdown = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JESdown, *puppimet);

  eve->pfMET_pt_JESdown_  = newPfMETs_JESdown.at(0).pt();
  eve->pfMET_phi_JESdown_ = newPfMETs_JESdown.at(0).phi();

  eve->pfMETNoHF_pt_JESdown_  = newPfMETsNoHF_JESdown.at(0).pt();
  eve->pfMETNoHF_phi_JESdown_ = newPfMETsNoHF_JESdown.at(0).phi();

  eve->puppiMET_pt_JESdown_  = newPuppiMETs_JESdown.at(0).pt();
  eve->puppiMET_phi_JESdown_ = newPuppiMETs_JESdown.at(0).phi();


  // JERup
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> newJetsForMET_JERup = GetCorrectedJets(oldJetsForMET_uncorr, run, rho_event, sysType::NA);
  // JER
  newJetsForMET_JERup = miniAODhelper.GetCorrectedJets(newJetsForMET_JERup, iEvent, iSetup, genjets, r, sysType::JERup, 0, doJER);

  std::vector<pat::MET> newPfMETs_JERup = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JERup, *pfmet);
  std::vector<pat::MET> newPfMETsNoHF_JERup = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JERup, *pfmetnohf);
  std::vector<pat::MET> newPuppiMETs_JERup = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JERup, *puppimet);

  eve->pfMET_pt_JERup_  = newPfMETs_JERup.at(0).pt();
  eve->pfMET_phi_JERup_ = newPfMETs_JERup.at(0).phi();

  eve->pfMETNoHF_pt_JERup_  = newPfMETsNoHF_JERup.at(0).pt();
  eve->pfMETNoHF_phi_JERup_ = newPfMETsNoHF_JERup.at(0).phi();

  eve->puppiMET_pt_JERup_  = newPuppiMETs_JERup.at(0).pt();
  eve->puppiMET_phi_JERup_ = newPuppiMETs_JERup.at(0).phi();

  // JERdown
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> newJetsForMET_JERdown = GetCorrectedJets(oldJetsForMET_uncorr, run, rho_event, sysType::NA);
  // JER
  newJetsForMET_JERdown = miniAODhelper.GetCorrectedJets(newJetsForMET_JERdown, iEvent, iSetup, genjets, r, sysType::JERdown, 0, doJER);
 
  std::vector<pat::MET> newPfMETs_JERdown = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JERdown, *pfmet);
  std::vector<pat::MET> newPfMETsNoHF_JERdown = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JERdown, *pfmetnohf);
  std::vector<pat::MET> newPuppiMETs_JERdown = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET_JERdown, *puppimet);

  eve->pfMET_pt_JERdown_  = newPfMETs_JERdown.at(0).pt();
  eve->pfMET_phi_JERdown_ = newPfMETs_JERdown.at(0).phi();

  eve->pfMETNoHF_pt_JERdown_  = newPfMETsNoHF_JERdown.at(0).pt();
  eve->pfMETNoHF_phi_JERdown_ = newPfMETsNoHF_JERdown.at(0).phi();

  eve->puppiMET_pt_JERdown_  = newPuppiMETs_JERdown.at(0).pt();
  eve->puppiMET_phi_JERdown_ = newPuppiMETs_JERdown.at(0).phi();

  */

  if( debug_ ) std::cout << " ====> test 16 " << std::endl;

  h_numJet->Fill(std::min(numJet,NjetMax),wgt);
  h_numTag->Fill(std::min(numTag,NtagMax),wgt);
  h_numPV->Fill(std::min(numpv,NpuMax),wgt);
  h_numTruePV->Fill(std::min(numTruePV,double(NpuMax)),wgt);

  h_numSecVtx->Fill(numSecVtx,wgt);


  // Loop over jets
  std::vector<double> vec_jet_pt;
  std::vector<double> vec_jet_eta;
  std::vector<double> vec_jet_phi;
  std::vector<double> vec_jet_energy;
  std::vector<double> vec_jet_csv;
  std::vector<double> vec_jet_cmva;
  //std::vector<int>    vec_jet_partonFlavour;
  //std::vector<int>    vec_jet_hadronFlavour;
  //std::vector<double> vec_jet_pileupJetId_fullDiscriminant;

  for( std::vector<pat::Jet>::const_iterator iJet = selectedJets.begin(); iJet != selectedJets.end(); iJet++ ){ 
    vec_jet_pt.push_back(iJet->pt());
    vec_jet_eta.push_back(iJet->eta());
    vec_jet_phi.push_back(iJet->phi());
    vec_jet_energy.push_back(iJet->energy());
    vec_jet_csv.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    vec_jet_cmva.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    //vec_jet_partonFlavour.push_back(iJet->partonFlavour());
    //vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
    //vec_jet_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
  } // end loop on jets

  eve->jet_pt_  = vec_jet_pt;
  eve->jet_eta_ = vec_jet_eta;
  eve->jet_phi_ = vec_jet_phi;
  eve->jet_energy_ = vec_jet_energy;
  eve->jet_csv_ = vec_jet_csv;
  eve->jet_cmva_ = vec_jet_cmva;
  //eve->jet_partonFlavour_ = vec_jet_partonFlavour;
  //eve->jet_hadronFlavour_ = vec_jet_hadronFlavour;
  //eve->jet_pileupJetId_fullDiscriminant_ = vec_jet_pileupJetId_fullDiscriminant;


  /////
  /*
  // Loop over jets
  std::vector<double> vec_jet_nocc_pt;
  std::vector<double> vec_jet_nocc_eta;
  std::vector<double> vec_jet_nocc_phi;
  std::vector<double> vec_jet_nocc_energy;
  std::vector<double> vec_jet_nocc_csv;
  std::vector<double> vec_jet_nocc_cmva;
  std::vector<int>    vec_jet_nocc_partonFlavour;
  std::vector<int>    vec_jet_nocc_hadronFlavour;
  std::vector<double> vec_jet_nocc_pileupJetId_fullDiscriminant;

  for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_nocc.begin(); iJet != selectedJets_nocc.end(); iJet++ ){ 
    vec_jet_nocc_pt.push_back(iJet->pt());
    vec_jet_nocc_eta.push_back(iJet->eta());
    vec_jet_nocc_phi.push_back(iJet->phi());
    vec_jet_nocc_energy.push_back(iJet->energy());
    vec_jet_nocc_csv.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    vec_jet_nocc_cmva.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    vec_jet_nocc_partonFlavour.push_back(iJet->partonFlavour());
    vec_jet_nocc_hadronFlavour.push_back(iJet->hadronFlavour());
    vec_jet_nocc_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
  } // end loop on jets

  eve->jet_nocc_pt_  = vec_jet_nocc_pt;
  eve->jet_nocc_eta_ = vec_jet_nocc_eta;
  eve->jet_nocc_phi_ = vec_jet_nocc_phi;
  eve->jet_nocc_energy_ = vec_jet_nocc_energy;
  eve->jet_nocc_csv_ = vec_jet_nocc_csv;
  eve->jet_nocc_cmva_ = vec_jet_nocc_cmva;
  eve->jet_nocc_partonFlavour_ = vec_jet_nocc_partonFlavour;
  eve->jet_nocc_hadronFlavour_ = vec_jet_nocc_hadronFlavour;
  eve->jet_nocc_pileupJetId_fullDiscriminant_ = vec_jet_nocc_pileupJetId_fullDiscriminant;



  ////


  // JESUp
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> correctedJets_JESup = GetCorrectedJets(rawJets, run, rho_event, sysType::JESup);
  // JER
  correctedJets_JESup = miniAODhelper.GetCorrectedJets(correctedJets_JESup, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);
 
  std::vector<pat::Jet> selectedJets_JESup_unsorted = miniAODhelper.GetSelectedJets(correctedJets_JESup, minLooseJetPt, 3.0, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_JESup = miniAODhelper.GetSortedByPt( selectedJets_JESup_unsorted );

  // Loop over jets
  std::vector<double> vec_jet_JESup_pt;
  std::vector<double> vec_jet_JESup_eta;
  std::vector<double> vec_jet_JESup_phi;
  std::vector<double> vec_jet_JESup_energy;
  std::vector<double> vec_jet_JESup_csv;
  std::vector<double> vec_jet_JESup_cmva;
  std::vector<int>    vec_jet_JESup_partonFlavour;
  std::vector<int>    vec_jet_JESup_hadronFlavour;
  std::vector<double> vec_jet_JESup_pileupJetId_fullDiscriminant;

  for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_JESup.begin(); iJet != selectedJets_JESup.end(); iJet++ ){ 
    vec_jet_JESup_pt.push_back(iJet->pt());
    vec_jet_JESup_eta.push_back(iJet->eta());
    vec_jet_JESup_phi.push_back(iJet->phi());
    vec_jet_JESup_energy.push_back(iJet->energy());
    vec_jet_JESup_csv.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    vec_jet_JESup_cmva.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    vec_jet_JESup_partonFlavour.push_back(iJet->partonFlavour());
    vec_jet_JESup_hadronFlavour.push_back(iJet->hadronFlavour());
    vec_jet_JESup_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
  } // end loop on jets

  eve->jet_JESup_pt_  = vec_jet_JESup_pt;
  eve->jet_JESup_eta_ = vec_jet_JESup_eta;
  eve->jet_JESup_phi_ = vec_jet_JESup_phi;
  eve->jet_JESup_energy_ = vec_jet_JESup_energy;
  eve->jet_JESup_csv_ = vec_jet_JESup_csv;
  eve->jet_JESup_cmva_ = vec_jet_JESup_cmva;
  eve->jet_JESup_partonFlavour_ = vec_jet_JESup_partonFlavour;
  eve->jet_JESup_hadronFlavour_ = vec_jet_JESup_hadronFlavour;
  eve->jet_JESup_pileupJetId_fullDiscriminant_ = vec_jet_JESup_pileupJetId_fullDiscriminant;


  //

  // JESDown
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> correctedJets_JESdown = GetCorrectedJets(rawJets, run, rho_event, sysType::JESdown);
  // JER
  correctedJets_JESdown = miniAODhelper.GetCorrectedJets(correctedJets_JESdown, iEvent, iSetup, genjets, r, sysType::NA, 0, doJER);
  
  std::vector<pat::Jet> selectedJets_JESdown_unsorted = miniAODhelper.GetSelectedJets(correctedJets_JESdown, minLooseJetPt, 3.0, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_JESdown = miniAODhelper.GetSortedByPt( selectedJets_JESdown_unsorted );

  // Loop over jets
  std::vector<double> vec_jet_JESdown_pt;
  std::vector<double> vec_jet_JESdown_eta;
  std::vector<double> vec_jet_JESdown_phi;
  std::vector<double> vec_jet_JESdown_energy;
  std::vector<double> vec_jet_JESdown_csv;
  std::vector<double> vec_jet_JESdown_cmva;
  std::vector<int>    vec_jet_JESdown_partonFlavour;
  std::vector<int>    vec_jet_JESdown_hadronFlavour;
  std::vector<double> vec_jet_JESdown_pileupJetId_fullDiscriminant;

  for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_JESdown.begin(); iJet != selectedJets_JESdown.end(); iJet++ ){ 
    vec_jet_JESdown_pt.push_back(iJet->pt());
    vec_jet_JESdown_eta.push_back(iJet->eta());
    vec_jet_JESdown_phi.push_back(iJet->phi());
    vec_jet_JESdown_energy.push_back(iJet->energy());
    vec_jet_JESdown_csv.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    vec_jet_JESdown_cmva.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    vec_jet_JESdown_partonFlavour.push_back(iJet->partonFlavour());
    vec_jet_JESdown_hadronFlavour.push_back(iJet->hadronFlavour());
    vec_jet_JESdown_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
  } // end loop on jets

  eve->jet_JESdown_pt_  = vec_jet_JESdown_pt;
  eve->jet_JESdown_eta_ = vec_jet_JESdown_eta;
  eve->jet_JESdown_phi_ = vec_jet_JESdown_phi;
  eve->jet_JESdown_energy_ = vec_jet_JESdown_energy;
  eve->jet_JESdown_csv_ = vec_jet_JESdown_csv;
  eve->jet_JESdown_cmva_ = vec_jet_JESdown_cmva;
  eve->jet_JESdown_partonFlavour_ = vec_jet_JESdown_partonFlavour;
  eve->jet_JESdown_hadronFlavour_ = vec_jet_JESdown_hadronFlavour;
  eve->jet_JESdown_pileupJetId_fullDiscriminant_ = vec_jet_JESdown_pileupJetId_fullDiscriminant;



  // JERUp
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> correctedJets_JERup = GetCorrectedJets(rawJets, run, rho_event, sysType::NA);
  // JER
  correctedJets_JERup = miniAODhelper.GetCorrectedJets(correctedJets_JERup, iEvent, iSetup, genjets, r, sysType::JERup, 0, doJER);
 
  std::vector<pat::Jet> selectedJets_JERup_unsorted = miniAODhelper.GetSelectedJets(correctedJets_JERup, minLooseJetPt, 3.0, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_JERup = miniAODhelper.GetSortedByPt( selectedJets_JERup_unsorted );

  // Loop over jets
  std::vector<double> vec_jet_JERup_pt;
  std::vector<double> vec_jet_JERup_eta;
  std::vector<double> vec_jet_JERup_phi;
  std::vector<double> vec_jet_JERup_energy;
  std::vector<double> vec_jet_JERup_csv;
  std::vector<double> vec_jet_JERup_cmva;
  std::vector<int>    vec_jet_JERup_partonFlavour;
  std::vector<int>    vec_jet_JERup_hadronFlavour;
  std::vector<double> vec_jet_JERup_pileupJetId_fullDiscriminant;

  for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_JERup.begin(); iJet != selectedJets_JERup.end(); iJet++ ){ 
    vec_jet_JERup_pt.push_back(iJet->pt());
    vec_jet_JERup_eta.push_back(iJet->eta());
    vec_jet_JERup_phi.push_back(iJet->phi());
    vec_jet_JERup_energy.push_back(iJet->energy());
    vec_jet_JERup_csv.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    vec_jet_JERup_cmva.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    vec_jet_JERup_partonFlavour.push_back(iJet->partonFlavour());
    vec_jet_JERup_hadronFlavour.push_back(iJet->hadronFlavour());
    vec_jet_JERup_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
  } // end loop on jets

  eve->jet_JERup_pt_  = vec_jet_JERup_pt;
  eve->jet_JERup_eta_ = vec_jet_JERup_eta;
  eve->jet_JERup_phi_ = vec_jet_JERup_phi;
  eve->jet_JERup_energy_ = vec_jet_JERup_energy;
  eve->jet_JERup_csv_ = vec_jet_JERup_csv;
  eve->jet_JERup_cmva_ = vec_jet_JERup_cmva;
  eve->jet_JERup_partonFlavour_ = vec_jet_JERup_partonFlavour;
  eve->jet_JERup_hadronFlavour_ = vec_jet_JERup_hadronFlavour;
  eve->jet_JERup_pileupJetId_fullDiscriminant_ = vec_jet_JERup_pileupJetId_fullDiscriminant;


  //

  // JERDown
  // Use JEC from GT
  // JEC
  std::vector<pat::Jet> correctedJets_JERdown = GetCorrectedJets(rawJets, run, rho_event, sysType::NA);
  // JER
  correctedJets_JERdown = miniAODhelper.GetCorrectedJets(correctedJets_JERdown, iEvent, iSetup, genjets, r, sysType::JERdown, 0, doJER);
  
  std::vector<pat::Jet> selectedJets_JERdown_unsorted = miniAODhelper.GetSelectedJets(correctedJets_JERdown, minLooseJetPt, 3.0, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_JERdown = miniAODhelper.GetSortedByPt( selectedJets_JERdown_unsorted );

  // Loop over jets
  std::vector<double> vec_jet_JERdown_pt;
  std::vector<double> vec_jet_JERdown_eta;
  std::vector<double> vec_jet_JERdown_phi;
  std::vector<double> vec_jet_JERdown_energy;
  std::vector<double> vec_jet_JERdown_csv;
  std::vector<double> vec_jet_JERdown_cmva;
  std::vector<int>    vec_jet_JERdown_partonFlavour;
  std::vector<int>    vec_jet_JERdown_hadronFlavour;
  std::vector<double> vec_jet_JERdown_pileupJetId_fullDiscriminant;

  for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_JERdown.begin(); iJet != selectedJets_JERdown.end(); iJet++ ){ 
    vec_jet_JERdown_pt.push_back(iJet->pt());
    vec_jet_JERdown_eta.push_back(iJet->eta());
    vec_jet_JERdown_phi.push_back(iJet->phi());
    vec_jet_JERdown_energy.push_back(iJet->energy());
    vec_jet_JERdown_csv.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    vec_jet_JERdown_cmva.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    vec_jet_JERdown_partonFlavour.push_back(iJet->partonFlavour());
    vec_jet_JERdown_hadronFlavour.push_back(iJet->hadronFlavour());
    vec_jet_JERdown_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
  } // end loop on jets

  eve->jet_JERdown_pt_  = vec_jet_JERdown_pt;
  eve->jet_JERdown_eta_ = vec_jet_JERdown_eta;
  eve->jet_JERdown_phi_ = vec_jet_JERdown_phi;
  eve->jet_JERdown_energy_ = vec_jet_JERdown_energy;
  eve->jet_JERdown_csv_ = vec_jet_JERdown_csv;
  eve->jet_JERdown_cmva_ = vec_jet_JERdown_cmva;
  eve->jet_JERdown_partonFlavour_ = vec_jet_JERdown_partonFlavour;
  eve->jet_JERdown_hadronFlavour_ = vec_jet_JERdown_hadronFlavour;
  eve->jet_JERdown_pileupJetId_fullDiscriminant_ = vec_jet_JERdown_pileupJetId_fullDiscriminant;

  */

  ////////////


  if( debug_ ) std::cout << " ====> test 17 " << std::endl;



  vint lepton_genId, lepton_charge, lepton_isMuon, lepton_isTight, lepton_is_IDTight, lepton_isLoose;
  //vint lepton_genParentId, lepton_genGrandParentId, lepton_trkCharge;
  //vint lepton_isPhys14L, lepton_isPhys14M, lepton_isPhys14T, lepton_isSpring15L, lepton_isSpring15M, lepton_isSpring15T;
  vint lepton_isTrigMVAM, lepton_isNonTrigMVAM, lepton_isTrigCutM;
  vdouble lepton_pt;
  vdouble lepton_px;
  vdouble lepton_py;
  vdouble lepton_pz;
  vdouble lepton_eta;
  vdouble lepton_sc_eta;
  vdouble lepton_phi;
  vdouble lepton_energy;
  vdouble lepton_relIso;
  vdouble lepton_relIsoR04;
  /*
  vdouble lepton_iso_sumChargedHadronPt;
  vdouble lepton_iso_sumNeutralHadronEt;
  vdouble lepton_iso_sumPhotonEt;
  vdouble lepton_iso_sumPUPt;
  vdouble lepton_trigMVAOutput;
  vint    lepton_trigMVACategory;
  vint    lepton_passTrigPresel;
  vint    lepton_numMissingHits;
  vdouble lepton_scSigmaIEtaIEta;
  vdouble lepton_full5x5_scSigmaIEtaIEta;
  vdouble lepton_hadronicOverEm;
  vdouble lepton_hcalOverEcal;
  vdouble lepton_ecalPFClusterIso;
  vdouble lepton_hcalPFClusterIso;
  vdouble lepton_dr03TkSumPt;
  vdouble lepton_relEcalIso;
  vdouble lepton_relHcalIso;
  vdouble lepton_relTrackIso;
  vdouble lepton_OneOESuperMinusOneOP;
  vint    lepton_isEB;
  vint    lepton_passHLTId;
  vint    lepton_passConversionVeto;
  vint    lepton_inCrack;
  vdouble lepton_scEta;
  vdouble lepton_dEtaSCTrackAtVtx;
  vdouble lepton_dPhiSCTrackAtVtx;
  vdouble lepton_d0;
  vdouble lepton_dZ;
  vint lepton_isGlobalMuon;
  vint lepton_isTrackerMuon;
  vint lepton_isPFMuon;
  vdouble lepton_normalizedChi2;
  vint lepton_numberOfValidMuonHits;
  vint lepton_numberOfValidPixelHits;
  vint lepton_trackerLayersWithMeasurement;
  vint lepton_numberOfMatchedStations;

  vint lepton_ele_matchHLT_hltL1sL1SingleEG25;
  vint lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter;

  vint lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25;
  vint lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter;

  vint lepton_ele_matchHLT_hltL1sL1EG25erHTT125;
  vint lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter;
  vint lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter;
  */

  std::vector<TLorentzVector> vec_TLV_lep;
  TLorentzVector sum_lepton_vect;
  sum_lepton_vect.SetPxPyPzE(0., 0., 0., 0.);

  // Loop over muons
  for( std::vector<pat::Muon>::const_iterator iMu = muons->begin(); iMu != muons->end(); iMu++ ){ 
 
    int genId=-99;
    //int genParentId=-99, genGrandParentId=-99;
    if( (iMu->genLepton()) ){
      genId = iMu->genLepton()->pdgId();
      //if( iMu->genLepton()->numberOfMothers()>=1 ){
	//genParentId = iMu->genLepton()->mother(0)->pdgId();
	//if( iMu->genLepton()->mother(0)->numberOfMothers()>=1 ) genGrandParentId = iMu->genLepton()->mother(0)->mother(0)->pdgId();
    //  }
    }

    //int trkCharge = -99;
    //if( iMu->muonBestTrack().isAvailable() ) trkCharge = iMu->muonBestTrack()->charge();

    int isTight = ( miniAODhelper.isGoodMuon(*iMu, minTightLeptonPt, 2.1, muonID::muonTight, coneSize::R04, corrType::deltaBeta) ) ? 1 : 0;
    int isLoose = ( miniAODhelper.isGoodMuon(*iMu, minLooseLeptonPt, 2.4, muonID::muonLoose, coneSize::R04, corrType::deltaBeta) ) ? 1 : 0;
	int is_IDTight = ( miniAODhelper.isGoodMuon(*iMu, minLooseLeptonPt, 2.4, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta) ) ? 1 : 0;

    /*
    int isPhys14L = false;
    int isPhys14M = false;
    int isPhys14T = false;

    int isSpring15L = false;
    int isSpring15M = false;
    int isSpring15T = false;
    */
    int isTrigMVAM = false;
    int isNonTrigMVAM = false;
    int isTrigCutM = false;

    /*
    double d0 = -999;
    double dZ = -999;
    if( iMu->muonBestTrack().isAvailable() ){
      d0 = iMu->muonBestTrack()->dxy(vertex.position());
      dZ = iMu->muonBestTrack()->dz(vertex.position());
    }
    
    double normalizedChi2 = -999;
    int numberOfValidMuonHits = -999;
    if( iMu->globalTrack().isAvailable() ){
      normalizedChi2 = iMu->globalTrack()->normalizedChi2();
      numberOfValidMuonHits = iMu->globalTrack()->hitPattern().numberOfValidMuonHits();
    }

    int numberOfValidPixelHits = -999;
    if( iMu->innerTrack().isAvailable() ){
      numberOfValidPixelHits = iMu->innerTrack()->hitPattern().numberOfValidPixelHits();
    }

    int trackerLayersWithMeasurement = -999;
    if( iMu->track().isAvailable() ){
      trackerLayersWithMeasurement = iMu->track()->hitPattern().trackerLayersWithMeasurement();
    }

    int numberOfMatchedStations = iMu->numberOfMatchedStations();
    */

    //lepton_trkCharge.push_back(trkCharge);
    lepton_charge.push_back(iMu->charge());
    lepton_isMuon.push_back(1);
    lepton_isTight.push_back(isTight);
    lepton_isLoose.push_back(isLoose);
    lepton_is_IDTight.push_back(is_IDTight);
    //lepton_isPhys14L.push_back(isPhys14L);
    //lepton_isPhys14M.push_back(isPhys14M);
    //lepton_isPhys14T.push_back(isPhys14T);
    //lepton_isSpring15L.push_back(isSpring15L);
    //lepton_isSpring15M.push_back(isSpring15M);
    //lepton_isSpring15T.push_back(isSpring15T);
    lepton_isTrigMVAM.push_back(isTrigMVAM);
    lepton_isNonTrigMVAM.push_back(isNonTrigMVAM);
    lepton_isTrigCutM.push_back(isTrigCutM);
    lepton_genId.push_back(genId);
    //lepton_genParentId.push_back(genParentId);
    //lepton_genGrandParentId.push_back(genGrandParentId);
    lepton_pt.push_back(iMu->pt());
    lepton_px.push_back(iMu->px());
    lepton_py.push_back(iMu->py()); 
    lepton_pz.push_back(iMu->pz());	  	  
    lepton_eta.push_back(iMu->eta());
    lepton_sc_eta.push_back(-999);
    lepton_phi.push_back(iMu->phi());
    lepton_energy.push_back(iMu->energy());
    lepton_relIso.push_back(miniAODhelper.GetMuonRelIso(*iMu));
    lepton_relIsoR04.push_back(miniAODhelper.GetMuonRelIso(*iMu, coneSize::R04, corrType::deltaBeta));
    /*   
    lepton_iso_sumChargedHadronPt.push_back(iMu->pfIsolationR03().sumChargedHadronPt);
    lepton_iso_sumNeutralHadronEt.push_back(iMu->pfIsolationR03().sumNeutralHadronEt);
    lepton_iso_sumPhotonEt.push_back(iMu->pfIsolationR03().sumPhotonEt);
    lepton_iso_sumPUPt.push_back(iMu->pfIsolationR03().sumPUPt);
    lepton_trigMVAOutput.push_back(-99);
    lepton_trigMVACategory.push_back(-99);
    lepton_passTrigPresel.push_back(-99);
    lepton_scSigmaIEtaIEta.push_back(-99);
    lepton_full5x5_scSigmaIEtaIEta.push_back(-99);
    lepton_hadronicOverEm.push_back(-99);
    lepton_hcalOverEcal.push_back(-99);
    lepton_ecalPFClusterIso.push_back(-99);
    lepton_hcalPFClusterIso.push_back(-99);
    lepton_dr03TkSumPt.push_back(-99);
    lepton_relEcalIso.push_back(-99);
    lepton_relHcalIso.push_back(-99);
    lepton_relTrackIso.push_back(-99);
    lepton_OneOESuperMinusOneOP.push_back(-99);
    lepton_numMissingHits.push_back(-99);
    lepton_isEB.push_back(-99);
    lepton_passHLTId.push_back(-99);
    lepton_passConversionVeto.push_back(-99);
    lepton_inCrack.push_back(-99);
    lepton_scEta.push_back(-99);
    lepton_dEtaSCTrackAtVtx.push_back(-99);
    lepton_dPhiSCTrackAtVtx.push_back(-99);

    lepton_d0.push_back(d0);
    lepton_dZ.push_back(dZ);
    lepton_isGlobalMuon.push_back(iMu->isGlobalMuon());
    lepton_isTrackerMuon.push_back(iMu->isTrackerMuon());
    lepton_isPFMuon.push_back(iMu->isPFMuon());
    lepton_normalizedChi2.push_back(normalizedChi2);
    lepton_numberOfValidMuonHits.push_back(numberOfValidMuonHits);
    lepton_numberOfValidPixelHits.push_back(numberOfValidPixelHits);
    lepton_trackerLayersWithMeasurement.push_back(trackerLayersWithMeasurement);
    lepton_numberOfMatchedStations.push_back(numberOfMatchedStations);

    lepton_ele_matchHLT_hltL1sL1SingleEG25.push_back(-99);
    lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter.push_back(-99);
    lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25.push_back(-99);
    lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter.push_back(-99);
    lepton_ele_matchHLT_hltL1sL1EG25erHTT125.push_back(-99);
    lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter.push_back(-99);
    lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter.push_back(-99);
    */

    // Get muon 4Vector and add to vecTLorentzVector for muons
    TLorentzVector leptonP4;	  
    leptonP4.SetPxPyPzE(iMu->px(),iMu->py(),iMu->pz(),iMu->energy());
    vec_TLV_lep.push_back(leptonP4);

    sum_lepton_vect += leptonP4;
  }

  // Loop over electrons
  for( std::vector<pat::Electron>::const_iterator iEle = electrons.begin(); iEle != electrons.end(); iEle++ ){ 

    int genId=-99;
    //int genParentId=-99, genGrandParentId=-99;
    if( (iEle->genLepton()) ){
      genId = iEle->genLepton()->pdgId();
    //  if( iEle->genLepton()->numberOfMothers()>=1 ){
	//genParentId = iEle->genLepton()->mother(0)->pdgId();
	//if( iEle->genLepton()->mother(0)->numberOfMothers()>=1 ) genGrandParentId = iEle->genLepton()->mother(0)->mother(0)->pdgId();
    //  }
    }

    //int trkCharge = -99;
    //if( iEle->gsfTrack().isAvailable() ) trkCharge = iEle->gsfTrack()->charge();

    int isTight = ( miniAODhelper.isGoodElectron(*iEle, minTightLeptonPt, 2.1, electronID::electronTight) ) ? 1 : 0;
    int isLoose = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronLoose) ) ? 1 : 0;
    int is_IDTight = false;

    //int isPhys14L = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronPhys14L) ) ? 1 : 0;
    //int isPhys14M = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronPhys14M) ) ? 1 : 0;
    //int isPhys14T = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronPhys14T) ) ? 1 : 0;

    //int isSpring15L = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronSpring15L) ) ? 1 : 0;
    //int isSpring15M = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronSpring15M) ) ? 1 : 0;
    //int isSpring15T = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronSpring15T) ) ? 1 : 0;
    int isTrigMVAM = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronEndOf15MVA80) ) ? 1 : 0;
    int isNonTrigMVAM = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electronNonTrigMVAid80) ) ? 1 : 0;
	int isTrigCutM = ( miniAODhelper.isGoodElectron(*iEle, minLooseLeptonPt, 2.4, electronID::electron80XCutBasedM) ) ? 1 : 0;


    //double trigMVAOutput = iEle->userFloat("mvaValue");//miniAODhelper.GetElectronMVAIDValue(*iEle);//myMVATrig->mvaValue(*iEle,false);
    //int category = iEle->userInt("mvaCategory");

    /*
    bool myTrigPresel = ( iEle->pt()>15 && 
			  ( ( abs(iEle->superCluster()->position().eta()) < 1.4442 && 
			      iEle->full5x5_sigmaIetaIeta() < 0.012 && 
			      iEle->hcalOverEcal() < 0.09 && 
			      (iEle->ecalPFClusterIso() / iEle->pt()) < 0.37 && 
			      (iEle->hcalPFClusterIso() / iEle->pt()) < 0.25 && 
			      (iEle->dr03TkSumPt() / iEle->pt()) < 0.18 && 
			      abs(iEle->deltaEtaSuperClusterTrackAtVtx()) < 0.0095 && 
			      abs(iEle->deltaPhiSuperClusterTrackAtVtx()) < 0.065 ) || 
			    ( abs(iEle->superCluster()->position().eta()) > 1.5660 && 
			      iEle->full5x5_sigmaIetaIeta() < 0.033 && 
			      iEle->hcalOverEcal() <0.09 && 
			      (iEle->ecalPFClusterIso() / iEle->pt()) < 0.45 && 
			      (iEle->hcalPFClusterIso() / iEle->pt()) < 0.28 && 
			      (iEle->dr03TkSumPt() / iEle->pt()) < 0.18 
			      ))
			  );

    bool inCrack = false;
    double scEta = -99;
    if( iEle->superCluster().isAvailable() ){
      inCrack = ( fabs(iEle->superCluster()->position().eta())>1.4442 && fabs(iEle->superCluster()->position().eta())<1.5660 );
      scEta = iEle->superCluster()->position().eta();
    }

    double ooEmooP = -999;
    if( iEle->ecalEnergy() == 0 ) ooEmooP = 1e30;
    else if( !std::isfinite(iEle->ecalEnergy()) ) ooEmooP = 1e30;
    else ooEmooP = fabs(1.0/iEle->ecalEnergy() - iEle->eSuperClusterOverP()/iEle->ecalEnergy() );

    bool passHLTId = false;
    double OneOESuperMinusOneOP = ooEmooP;
    int numMissingHits = 99;
    int isEB = -1;
    double relEcalIso = 99, relHcalIso = 99, relTrackIso = 99;
    double d0 = -999;
    double dZ = -999;
    if( (iEle->superCluster().isAvailable()) && (iEle->gsfTrack().isAvailable()) ){
      double SCenergy = iEle->superCluster()->energy();
      double SCeta = iEle->superCluster()->position().eta();
      double absSCeta = fabs( SCeta );
      double SCet = SCenergy * sin (2*atan(exp(-SCeta))); 

      numMissingHits = iEle->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

      d0 = fabs(iEle->gsfTrack()->dxy(vertex.position()));
      dZ = fabs(iEle->gsfTrack()->dz(vertex.position()));

      relEcalIso = iEle->ecalIso()/SCet;
      relHcalIso = iEle->hcalIso()/SCet;
      relTrackIso = iEle->trackIso()/SCet;

      if( absSCeta < 1.479 ){
	isEB = 1;
	passHLTId = ( iEle->scSigmaIEtaIEta() <= 0.011 &&
		      iEle->hadronicOverEm() <= 0.15 &&
		      iEle->ecalIso()/SCet <= 0.16 &&
		      iEle->hcalIso()/SCet <= 0.20 &&
		      iEle->trackIso()/SCet <= 0.05 &&
		      OneOESuperMinusOneOP <= 0.012 &&
		      numMissingHits <= 999
		      );
      }
      else{
	isEB = 0;
	passHLTId = ( iEle->scSigmaIEtaIEta() <= 0.033 &&
		      iEle->hadronicOverEm() <= 0.20 &&
		      iEle->ecalIso()/SCet <= 0.12 &&
		      iEle->hcalIso()/SCet <= 0.30 &&
		      iEle->trackIso()/SCet <= 0.05 &&
		      OneOESuperMinusOneOP <= 0.009 &&
		      numMissingHits <= 1
		      );
      }
    }
    */

    /*
    bool matchHLT_hltL1sL1SingleEG25 = false;
    bool matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter = false;
    bool matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25 = false;
    bool matchHLT_hltEle27WPLooseGsfTrackIsoFilter = false;
    bool matchHLT_hltL1sL1EG25erHTT125 = false;
    bool matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter = false;
    bool matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter = false;

    for( int iHLTObj = 0; iHLTObj < int(triggerObjects_hltEle.size()); iHLTObj++ ){ 
      std::string name = vec_hltEle27WP85Gsf_filter[iHLTObj];
      double deltaR = reco::deltaR( iEle->eta(), iEle->phi(), triggerObjects_hltEle[iHLTObj].Eta(), triggerObjects_hltEle[iHLTObj].Phi() );
      if( deltaR<0.5 ){
	if( name=="hltL1sL1SingleEG25" ) matchHLT_hltL1sL1SingleEG25 = true;
	if( name=="hltL1EG25Ele27WP85GsfTrackIsoFilter" ) matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter = true;
	if( name=="hltL1sL1SingleIsoEG22erOrSingleEG25" ) matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25 = true;
	if( name=="hltEle27WPLooseGsfTrackIsoFilter" ) matchHLT_hltEle27WPLooseGsfTrackIsoFilter = true;
	//if( name=="hltL1sL1EG25erHTT125" ) matchHLT_hltL1sL1EG25erHTT125 = true;
	if( name.find("hltL1sL1EG25erHTT")!=std::string::npos ) matchHLT_hltL1sL1EG25erHTT125 = true;
	if( name=="hltL1sL1EG25erHTT125" ) matchHLT_hltL1sL1EG25erHTT125 = true;
	if( name=="hltL1EGHttEle27WP85GsfTrackIsoFilter" ) matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter = true;
	if( name=="hltL1EGHttEle27WPLooseGsfTrackIsoFilter" ) matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter = true;
      }
    }
    */
    //lepton_trkCharge.push_back(trkCharge);
    lepton_charge.push_back(iEle->charge());
    lepton_isMuon.push_back(0);
    lepton_isTight.push_back(isTight);
    lepton_isLoose.push_back(isLoose);
    lepton_is_IDTight.push_back(is_IDTight);
    //lepton_isPhys14L.push_back(isPhys14L);
    //lepton_isPhys14M.push_back(isPhys14M);
    //lepton_isPhys14T.push_back(isPhys14T);
    //lepton_isSpring15L.push_back(isSpring15L);
    //lepton_isSpring15M.push_back(isSpring15M);
    //lepton_isSpring15T.push_back(isSpring15T);
    lepton_isTrigMVAM.push_back(isTrigMVAM);
    lepton_isNonTrigMVAM.push_back(isNonTrigMVAM);
    lepton_isTrigCutM.push_back(isTrigCutM);
    lepton_genId.push_back(genId);
    //lepton_genParentId.push_back(genParentId);
    //lepton_genGrandParentId.push_back(genGrandParentId);
    lepton_pt.push_back(iEle->pt());
    lepton_px.push_back(iEle->px());
    lepton_py.push_back(iEle->py());
    lepton_pz.push_back(iEle->pz());
    lepton_eta.push_back(iEle->eta());
    lepton_sc_eta.push_back(iEle->superCluster()->position().eta());
    lepton_phi.push_back(iEle->phi());
    lepton_energy.push_back(iEle->energy());
    lepton_relIso.push_back(miniAODhelper.GetElectronRelIso(*iEle,coneSize::R03,corrType::rhoEA,effAreaType::spring15));
    lepton_relIsoR04.push_back(miniAODhelper.GetElectronRelIso(*iEle));
    /*
    lepton_iso_sumChargedHadronPt.push_back(iEle->pfIsolationVariables().sumChargedHadronPt);
    lepton_iso_sumNeutralHadronEt.push_back(iEle->pfIsolationVariables().sumNeutralHadronEt);
    lepton_iso_sumPhotonEt.push_back(iEle->pfIsolationVariables().sumPhotonEt);
    lepton_iso_sumPUPt.push_back(iEle->pfIsolationVariables().sumPUPt);
    lepton_trigMVAOutput.push_back(trigMVAOutput);
    lepton_trigMVACategory.push_back(category);
    lepton_passTrigPresel.push_back( (myTrigPresel) ? 1 : 0);
    lepton_scSigmaIEtaIEta.push_back(iEle->scSigmaIEtaIEta());
    lepton_full5x5_scSigmaIEtaIEta.push_back(iEle->full5x5_sigmaIetaIeta());
    lepton_hadronicOverEm.push_back(iEle->hadronicOverEm());
    lepton_hcalOverEcal.push_back(iEle->hcalOverEcal());
    lepton_ecalPFClusterIso.push_back(iEle->ecalPFClusterIso());
    lepton_hcalPFClusterIso.push_back(iEle->hcalPFClusterIso());
    lepton_dr03TkSumPt.push_back(iEle->dr03TkSumPt());
    lepton_relEcalIso.push_back(relEcalIso);
    lepton_relHcalIso.push_back(relHcalIso);
    lepton_relTrackIso.push_back(relTrackIso);
    lepton_OneOESuperMinusOneOP.push_back(OneOESuperMinusOneOP);
    lepton_numMissingHits.push_back(numMissingHits);
    lepton_isEB.push_back(isEB);
    lepton_passHLTId.push_back(passHLTId);
    lepton_passConversionVeto.push_back(iEle->passConversionVeto());
    lepton_inCrack.push_back(inCrack);
    lepton_scEta.push_back(scEta);
    lepton_dEtaSCTrackAtVtx.push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
    lepton_dPhiSCTrackAtVtx.push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
    lepton_d0.push_back(d0);
    lepton_dZ.push_back(dZ);
    lepton_isGlobalMuon.push_back(-99);
    lepton_isTrackerMuon.push_back(-99);
    lepton_isPFMuon.push_back(-99);
    lepton_normalizedChi2.push_back(-99);
    lepton_numberOfValidMuonHits.push_back(-99);
    lepton_numberOfValidPixelHits.push_back(-99);
    lepton_trackerLayersWithMeasurement.push_back(-99);
    lepton_numberOfMatchedStations.push_back(-99);


    lepton_ele_matchHLT_hltL1sL1SingleEG25.push_back(matchHLT_hltL1sL1SingleEG25);
    lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter.push_back(matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter);
    lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25.push_back(matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25);
    lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter.push_back(matchHLT_hltEle27WPLooseGsfTrackIsoFilter);
    lepton_ele_matchHLT_hltL1sL1EG25erHTT125.push_back(matchHLT_hltL1sL1EG25erHTT125);
    lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter.push_back(matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter);
    lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter.push_back(matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter);
    */

    // Get electron 4Vector and add to vecTLorentzVector for electrons
    TLorentzVector leptonP4;	  
    leptonP4.SetPxPyPzE(iEle->px(),iEle->py(),iEle->pz(),iEle->energy());
    vec_TLV_lep.push_back(leptonP4);

    sum_lepton_vect += leptonP4;
  }

  //eve->lepton_trkCharge_        = lepton_trkCharge;
  eve->lepton_charge_           = lepton_charge;
  eve->lepton_isMuon_           = lepton_isMuon;
  eve->lepton_isTight_          = lepton_isTight;
  eve->lepton_isLoose_          = lepton_isLoose;
  eve->lepton_is_IDTight_		= lepton_is_IDTight;
  //eve->lepton_isPhys14L_        = lepton_isPhys14L;
  //eve->lepton_isPhys14M_        = lepton_isPhys14M;
  //eve->lepton_isPhys14T_        = lepton_isPhys14T;
  //eve->lepton_isSpring15L_      = lepton_isSpring15L;
  //eve->lepton_isSpring15M_      = lepton_isSpring15M;
  //eve->lepton_isSpring15T_      = lepton_isSpring15T;
  eve->lepton_isTrigMVAM_       = lepton_isTrigMVAM;
  eve->lepton_isNonTrigMVAM_    = lepton_isNonTrigMVAM;
  eve->lepton_isTrigCutM_       = lepton_isTrigCutM;
  eve->lepton_genId_            = lepton_genId;
  //eve->lepton_genParentId_      = lepton_genParentId;
  //eve->lepton_genGrandParentId_ = lepton_genGrandParentId;
  eve->lepton_pt_               = lepton_pt;
  eve->lepton_px_               = lepton_px;
  eve->lepton_py_               = lepton_py;
  eve->lepton_pz_               = lepton_pz;
  eve->lepton_eta_              = lepton_eta;
  eve->lepton_sc_eta_           = lepton_sc_eta;
  eve->lepton_phi_              = lepton_phi;
  eve->lepton_energy_           = lepton_energy;
  eve->lepton_relIso_           = lepton_relIso;
  eve->lepton_relIsoR04_           = lepton_relIsoR04;
  /*
  eve->lepton_iso_sumChargedHadronPt_ = lepton_iso_sumChargedHadronPt;
  eve->lepton_iso_sumNeutralHadronEt_ = lepton_iso_sumNeutralHadronEt;
  eve->lepton_iso_sumPhotonEt_        = lepton_iso_sumPhotonEt;
  eve->lepton_iso_sumPUPt_            = lepton_iso_sumPUPt;
  eve->lepton_trigMVAOutput_     = lepton_trigMVAOutput;
  eve->lepton_trigMVACategory_   = lepton_trigMVACategory;
  eve->lepton_passTrigPresel_    = lepton_passTrigPresel;
  eve->lepton_scSigmaIEtaIEta_  = lepton_scSigmaIEtaIEta;
  eve->lepton_full5x5_scSigmaIEtaIEta_ = lepton_full5x5_scSigmaIEtaIEta;
  eve->lepton_hadronicOverEm_   = lepton_hadronicOverEm;
  eve->lepton_hcalOverEcal_     = lepton_hcalOverEcal;
  eve->lepton_ecalPFClusterIso_ = lepton_ecalPFClusterIso;
  eve->lepton_hcalPFClusterIso_ = lepton_hcalPFClusterIso;
  eve->lepton_dr03TkSumPt_      = lepton_dr03TkSumPt;
  eve->lepton_relEcalIso_       = lepton_relEcalIso;
  eve->lepton_relHcalIso_       = lepton_relHcalIso;
  eve->lepton_relTrackIso_      = lepton_relTrackIso;
  eve->lepton_OneOESuperMinusOneOP_ = lepton_OneOESuperMinusOneOP;
  eve->lepton_numMissingHits_   = lepton_numMissingHits;
  eve->lepton_isEB_             = lepton_isEB;
  eve->lepton_passHLTId_        = lepton_passHLTId;
  eve->lepton_passConversionVeto_ = lepton_passConversionVeto;
  eve->lepton_inCrack_          = lepton_inCrack;
  eve->lepton_scEta_            = lepton_scEta;
  eve->lepton_dEtaSCTrackAtVtx_ = lepton_dEtaSCTrackAtVtx;
  eve->lepton_dPhiSCTrackAtVtx_ = lepton_dPhiSCTrackAtVtx;
  eve->lepton_d0_ = lepton_d0;
  eve->lepton_dZ_ = lepton_dZ;
  eve->lepton_isGlobalMuon_ = lepton_isGlobalMuon;
  eve->lepton_isTrackerMuon_ = lepton_isTrackerMuon;
  eve->lepton_isPFMuon_ = lepton_isPFMuon;
  eve->lepton_normalizedChi2_ = lepton_normalizedChi2;
  eve->lepton_numberOfValidMuonHits_ = lepton_numberOfValidMuonHits;
  eve->lepton_numberOfValidPixelHits_ = lepton_numberOfValidPixelHits;
  eve->lepton_trackerLayersWithMeasurement_ = lepton_trackerLayersWithMeasurement;
  eve->lepton_numberOfMatchedStations_ = lepton_numberOfMatchedStations;


  eve->lepton_ele_matchHLT_hltL1sL1SingleEG25_ = lepton_ele_matchHLT_hltL1sL1SingleEG25;
  eve->lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_ = lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter;
  eve->lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25_ = lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25;
  eve->lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_ = lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter;
  eve->lepton_ele_matchHLT_hltL1sL1EG25erHTT125_ = lepton_ele_matchHLT_hltL1sL1EG25erHTT125;
  eve->lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter_ = lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter;
  eve->lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter_ = lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter;
  */

  /*
  int oppositeLepCharge = -9;
  if( lepton_charge.size()==2 ){
    int chg0 = lepton_charge[0];
    int chg1 = lepton_charge[1];

    if( (chg0 * chg1)==-1 )     oppositeLepCharge = 1;
    else if( (chg0 * chg1)==1 ) oppositeLepCharge = 0;
    else if( chg0==-99 )        oppositeLepCharge = -1;
    else if( chg1==-99 )        oppositeLepCharge = -2;
    else                        oppositeLepCharge = -3;
  }
  eve->oppositeLepCharge_ = oppositeLepCharge;

  double mass_leplep = -99;
  if( vec_TLV_lep.size()==2 ){
    mass_leplep = sum_lepton_vect.M();
    eve->mass_leplep_ = mass_leplep;
    eve->dR_leplep_ = vec_TLV_lep[0].DeltaR(vec_TLV_lep[1]);
  }
  */

  if( debug_ ) std::cout << " ====> test 18 " << std::endl;
  
  // Event Weights
  
  double gen_weight = 1;
  double csv_weight = 1;
  double PU_weight = 1;
  double PU_weight_B = 1;
  double PU_weight_C = 1;
  double PU_weight_D = 1;
  double PU_weight_E = 1;
  double PU_weight_F = 1;
  double PU_weight_G = 1;
  double PU_weight_H = 1;
  double PU_weight_BCDEF = 1;
  double PU_weight_GH = 1;
  double PDF_weight = 1;
  double PDF_weight_up = 1;
  double PDF_weight_down = 1;
  double Q2_weight = 1;
  double Q2_weight_up = 1;
  double Q2_weight_down = 1;
  vdouble lepton_id_sf;
  vdouble lepton_iso_sf;
  vdouble lepton_gsf_sf;
  vdouble lepton_trig_sf;
  vdouble lepton_hip_sf;

  eve->isData_ = isData_;
  
  if(!isData_) {
		
	//Generator weight
	gen_weight = GenEventInfoHandle->weight();
	
	//CSV weight
	csv_weight = getbweight(selectedJets);
  
    //PDF weight
    PDF_weight = getPDFweight(GenEventInfoHandle, 0);
    PDF_weight_up = getPDFweight(GenEventInfoHandle, 1);
    PDF_weight_down = getPDFweight(GenEventInfoHandle, -1);

    //PU_weight
   	PU_weight = getPUweight(PupInfoHandle, 0);
   	PU_weight_B = getPUweight(PupInfoHandle, 1);
   	PU_weight_C = getPUweight(PupInfoHandle, 2);
	PU_weight_D = getPUweight(PupInfoHandle, 3);
	PU_weight_E = getPUweight(PupInfoHandle, 4);
	PU_weight_F = getPUweight(PupInfoHandle, 5);
	PU_weight_G = getPUweight(PupInfoHandle, 6);
	PU_weight_H = getPUweight(PupInfoHandle, 7);
    PU_weight_BCDEF = getPUweight(PupInfoHandle, 8);
    PU_weight_GH = getPUweight(PupInfoHandle, 9);
   		
    //Q2_weight
	Q2_weight_up = getQ2weight(GenEventInfoHandle, LHEEventProductHandle , "1005");
	Q2_weight_down = getQ2weight(GenEventInfoHandle, LHEEventProductHandle , "1009");

 }
 
 
  int N_lep = lepton_pt.size();
  
  for(int i=0; i<N_lep; i++) {

    if(!isData_) {

    //ID
    if(lepton_isMuon[i] == 1)
   		lepton_id_sf.push_back(leptonSFhelper.GetMuonSF(
                lepton_pt[i], lepton_eta[i], 0, "ID"));
    else
    	lepton_id_sf.push_back(leptonSFhelper.GetElectronSF(
                lepton_pt[i], lepton_sc_eta[i], 0, "ID")); 
    //ISO
    if(lepton_isMuon[i] == 1)
   		lepton_iso_sf.push_back(leptonSFhelper.GetMuonSF(
                lepton_pt[i], lepton_eta[i], 0, "Iso"));
    else
    	lepton_iso_sf.push_back(leptonSFhelper.GetElectronSF(
                lepton_pt[i], lepton_sc_eta[i], 0, "Iso")); 
    //GSF
    if(lepton_isMuon[i] == 1)
   		lepton_gsf_sf.push_back(1);
    else
    	lepton_gsf_sf.push_back(leptonSFhelper.GetElectronSF(
                lepton_pt[i], lepton_sc_eta[i], 0, "Gsf")); 
    //TRIGGER
    if(lepton_isMuon[i] == 1)
   		lepton_trig_sf.push_back(leptonSFhelper.GetMuonSF(
                lepton_pt[i], lepton_eta[i], 0, "Trigger"));
    else
    	lepton_trig_sf.push_back(leptonSFhelper.GetElectronSF(
                lepton_pt[i], lepton_sc_eta[i], 0, "Trigger"));
    //HIP
    if(lepton_isMuon[i] == 1)
   		lepton_hip_sf.push_back(leptonSFhelper.GetMuonSF(
                lepton_pt[i], lepton_eta[i], 0, "HIP"));
    else
    	lepton_hip_sf.push_back(1);
    }

    else {
        lepton_id_sf.push_back(1);
        lepton_iso_sf.push_back(1);
        lepton_gsf_sf.push_back(1);
        lepton_trig_sf.push_back(1);
        lepton_hip_sf.push_back(1);
    }

    }
   
  eve->gen_weight_ = gen_weight;
  eve->csv_weight_ = csv_weight;
  eve->PU_weight_ = PU_weight;
  eve->PU_weight_B_ = PU_weight_B;
  eve->PU_weight_C_ = PU_weight_C;
  eve->PU_weight_D_ = PU_weight_D;
  eve->PU_weight_E_ = PU_weight_E;
  eve->PU_weight_F_ = PU_weight_F;
  eve->PU_weight_G_ = PU_weight_G;
  eve->PU_weight_H_ = PU_weight_H;
  eve->PU_weight_BCDEF_ = PU_weight_BCDEF;
  eve->PU_weight_GH_ = PU_weight_GH;
  eve->PDF_weight_ = PDF_weight;
  eve->PDF_weight_up_ = PDF_weight_up;
  eve->PDF_weight_down_ = PDF_weight_down;
  eve->Q2_weight_ = Q2_weight;
  eve->Q2_weight_up_ = Q2_weight_up;
  eve->Q2_weight_down_ = Q2_weight_down;
  eve->lepton_id_sf_ = lepton_id_sf;
  eve->lepton_iso_sf_ = lepton_iso_sf;
  eve->lepton_gsf_sf_ = lepton_gsf_sf;
  eve->lepton_trig_sf_ = lepton_trig_sf;
  eve->lepton_hip_sf_ = lepton_hip_sf;

  m_ttree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzer::beginJob()
{

  TH1::SetDefaultSumw2(true);

  numEvents_         = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzer::endJob() 
{
  // report on triggers fired
  if( dumpHLT_ || true ){
    std::cout << "***************************************************************" << std::endl;
    std::cout << "  Summary for L1T Algo: Total number of events = " << numEvents_ << std::endl;
    int bit=0;
    for( std::map<std::string, int>::const_iterator iter = l1talgo_cppath_.begin(); iter != l1talgo_cppath_.end(); iter++){
      std::string name = iter->first;
      double eff = 100*double(iter->second)/double(numEvents_);
      printf("\t %d,\t %55s,\t %5d,\t %.1f\n",bit++, name.c_str(), iter->second, eff);
    }
    std::cout << "***************************************************************" << std::endl;
    std::cout << "      Trigger menu = " << hlt_config_.tableName() << std::endl;
    std::cout << "  Summary for HLT: Total number of events = " << numEvents_ << std::endl;
    for( std::map<std::string, int>::const_iterator iter = hlt_cppath_.begin(); iter != hlt_cppath_.end(); iter++){
      std::string name = iter->first;
      double eff = 100*double(iter->second)/double(numEvents_);
      printf("\t %90s,\t %5d,\t %.1f \n", name.c_str(), iter->second, eff);
    }
    std::cout << "***************************************************************" << std::endl;
    std::cout << "  Summary for Filters: Total number of events = " << numEvents_ << std::endl;
    for( std::map<std::string, int>::const_iterator iter = flt_cppath_.begin(); iter != flt_cppath_.end(); iter++){
      std::string name = iter->first;
      double eff = 100*double(iter->second)/double(numEvents_);
      printf("\t %55s,\t %5d,\t %.1f \n",name.c_str(), iter->second, eff);
    }
    std::cout << "***************************************************************" << std::endl;
  }

  std::cout << "***************************************************************" << std::endl;
  std::cout << "  Total number of events = " << numEvents_ << std::endl;
  std::cout << "***************************************************************" << std::endl;


}

// ------------ method called when starting to processes a run  ------------
/*
void 
TriggerAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
void 
TriggerAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  bool hltchanged = true;
  if (!hlt_config_.init(iRun, iSetup, hltTag, hltchanged)) {
    std::cout << "Warning, didn't find trigger process HLT with input tag:" << hltTag << std::endl;
    return;
  }
  bool filterchanged = true;
  if (!filter_config_.init(iRun, iSetup, filterTag, filterchanged)) {
    std::cout << "Warning, didn't find filter process HLT with input tag:" << filterTag << std::endl;
    return;
  }


  // Zero out map
  std::vector<std::string> triggerNames = hlt_config_.triggerNames();
  std::vector<std::string> filterNames  = filter_config_.triggerNames();

  hlt_triggerNames_.clear();
  flt_filterNames_.clear();

  hlt_triggerNames_.push_back("All");
  std::string prefix = "HLT_";
  for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
    std::string name = triggerNames[iPath];
    std::string pathNameNoVer = hlt_config_.removeVersion(name);
    hlt_cppath_[pathNameNoVer] = 0;
    if( name.compare(0, prefix.length(), prefix) == 0 ) hlt_triggerNames_.push_back(pathNameNoVer);
  }

  flt_filterNames_.push_back("All");
  for( unsigned int iPath=0; iPath<filterNames.size(); iPath++ ){
    std::string name = filterNames[iPath];
    std::string pathNameNoVer = filter_config_.removeVersion(name);
    flt_cppath_[pathNameNoVer] = 0;
    flt_filterNames_.push_back(pathNameNoVer);
  }

  unsigned int numHLT = hlt_triggerNames_.size();
  unsigned int numFLT = flt_filterNames_.size();

  h_hlt = fs_->make<TH1D>("h_hlt",";HLT path", numHLT , 0 , numHLT );
  h_flt = fs_->make<TH1D>("h_flt",";Filter path", numFLT , 0 , numFLT );
  
  for( unsigned int iPath=0; iPath<numHLT; iPath++ ){
    std::string pathNameNoVer = hlt_triggerNames_[iPath];
    int jPath = iPath+1;
    if( h_hlt ){
      TAxis * axis = h_hlt->GetXaxis();
      if( axis ) axis->SetBinLabel(jPath, pathNameNoVer.c_str());
    }
  }

  for( unsigned int iPath=0; iPath<numFLT; iPath++ ){
    std::string pathNameNoVer = flt_filterNames_[iPath];
    int jPath = iPath+1;
    if( h_flt ){
      TAxis * axis = h_flt->GetXaxis();
      if( axis ) axis->SetBinLabel(jPath, pathNameNoVer.c_str());
    }
  }

  // Get Trigger and Event Handles
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();

  // L1 algos
  for( CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
    std::string name = (algo->second).algoName();
    l1talgo_cppath_[name] = 0;
  }


  /*
  edm::Handle<LHERunInfoProduct> run; 
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
 
  iRun.getByLabel( "externalLHEProducer", run );
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());
 
  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }
  */

}

// ------------ method called when ending the processing of a run  ------------
/*
void 
TriggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TriggerAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TriggerAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzer);
