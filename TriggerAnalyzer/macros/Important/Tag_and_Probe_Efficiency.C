#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"
#include "TH1D.h"
#include <time.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "TriggerStudyEventVars.h"

#endif

typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;
typedef std::vector<std::vector<int> >         vvint;

void Control_Region_Sel( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

  clock_t t;
  t = clock();
  
  int N_total = 0;
  int N_eve = 0;
  int N_eve_loose_27 = 0;
  int N_eve_tight_27 = 0;
  int N_eve_tight_32 = 0;
  int N_eve_loose_27_ht200 = 0;
  int N_eve_tight_27_OR_loose_27_ht200 = 0;

  double effi_num_loose_27 = 0;
  double effi_den_loose_27 = 0;
  double effi_loose_27;
  double effi_num_loose_27_BCDEF = 0;
  double effi_den_loose_27_BCDEF = 0;
  double effi_loose_27_BCDEF;
  double effi_num_loose_27_GH = 0;
  double effi_den_loose_27_GH = 0;
  double effi_loose_27_GH;
  double effi_num_tight_27 = 0;
  double effi_den_tight_27 = 0;
  double effi_tight_27;
  double effi_num_tight_27_BCDEF = 0;
  double effi_den_tight_27_BCDEF = 0;
  double effi_tight_27_BCDEF;
  double effi_num_tight_27_GH = 0;
  double effi_den_tight_27_GH = 0;
  double effi_tight_27_GH;
  double effi_num_tight_32 = 0;
  double effi_den_tight_32 = 0;
  double effi_tight_32;
  double effi_num_tight_32_BCDEF = 0;
  double effi_den_tight_32_BCDEF = 0;
  double effi_tight_32_BCDEF;
  double effi_num_tight_32_GH = 0;
  double effi_den_tight_32_GH = 0;
  double effi_tight_32_GH;
  double effi_num_loose_27_ht200 = 0;
  double effi_den_loose_27_ht200 = 0;
  double effi_loose_27_ht200;
  double effi_num_loose_27_ht200_BCDEF = 0;
  double effi_den_loose_27_ht200_BCDEF = 0;
  double effi_loose_27_ht200_BCDEF;
  double effi_num_loose_27_ht200_GH = 0;
  double effi_den_loose_27_ht200_GH = 0;
  double effi_loose_27_ht200_GH;
  double effi_num_tight_27_or_loose_27_ht200 = 0;
  double effi_den_tight_27_or_loose_27_ht200 = 0;
  double effi_tight_27_or_loose_27_ht200;
  double effi_num_tight_27_or_loose_27_ht200_BCDEF = 0;
  double effi_den_tight_27_or_loose_27_ht200_BCDEF = 0;
  double effi_tight_27_or_loose_27_ht200_BCDEF;
  double effi_num_tight_27_or_loose_27_ht200_GH = 0;
  double effi_den_tight_27_or_loose_27_ht200_GH = 0;
  double effi_tight_27_or_loose_27_ht200_GH;

  double sum_gen_weight = 0;

  std::cout << "   ===> load the root files! " << std::endl;
  
  double wgt=1;
  
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string treefilename = "/eos/uscms/store/user/adatta/Trigger_Analysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/triggerTree_MC_80X_drell_yan/161104_212140/0000/trigger_analyzer_*.root";
  //std::string treefilename2 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0001/trigger_analyzer_*.root";
  //std::string treefilename3 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0002/trigger_analyzer_*.root";
  
  std::string s_end = "_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = ".root";

  std::string histofilename = "Tag_Probe_Efficiency_histos" + s_end;
 
  std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

  TChain *chain = new TChain("triggeranalyzer/triggerTree");
 
  chain->Add(treefilename.c_str());
  //chain->Add(treefilename2.c_str());
  //chain->Add(treefilename3.c_str());
  
  ////////////////////////////////////
   
  triggerStudyEventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();

  //TEfficiency
	 
  TH1::SetDefaultSumw2();

    int n_ptBins = 10;
    double x_ptBins[10] = {30, 40, 50, 60, 70, 80, 90 , 100, 150, 200, 500 };

    int n_scetaBins = 8;
    double x_scetaBins[8]
    = {  -2.2, // range -2.2_-2.1 is to contain electron with |eta|<2.1 while |SC_eta|>2.1.
        -2.1 + 0.6 * 1,
        -2.1 + 0.6 * 2,
        -2.1 + 0.6 * 3,
        -2.1 + 0.6 * 4,
        -2.1 + 0.6 * 5,
        -2.1 + 0.6 * 6,
        -2.1 + 0.6 * 7 + 0.1
    };

    int n_htBins = 18;
    double x_htBins[18] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90 , 100, 150, 200, 300, 400, 500, 750, 1000 };

    TEfficiency* Eff_pt_WPLoose_27 = new TEfficiency("Eff_pt_WPLoose_27","Efficiency vs pT for WPLoose_27;pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPLoose_27 = new TEfficiency("Eff_sceta_WPLoose_27","Efficiency vs SuperCluster eta for WPLoose_27;SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPLoose_27 = new TEfficiency("Eff_phi_WPLoose_27","Efficiency vs phi for WPLoose_27;#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPLoose_27 = new TEfficiency("Eff_2d_WPLoose_27","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPLoose_27 = new TEfficiency("Eff_HT_WPLoose_27","Efficiency vs HT for WPLoose_27;HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPLoose_27 = new TEfficiency("Eff_numPV_WPLoose_27","Efficiency vs numPV for WPLoose_27;numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPLoose_27 = new TEfficiency("Eff_nJets_WPLoose_27","Efficiency vs nJets for WPLoose_27;nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPLoose_27_BCDEF = new TEfficiency("Eff_pt_WPLoose_27_BCDEF","Efficiency vs pT for WPLoose_27 (for B to F);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPLoose_27_BCDEF = new TEfficiency("Eff_sceta_WPLoose_27_BCDEF","Efficiency vs SuperCluster eta for WPLoose_27 (for B to F);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPLoose_27_BCDEF = new TEfficiency("Eff_phi_WPLoose_27_BCDEF","Efficiency vs phi for WPLoose_27 (for B to F);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPLoose_27_BCDEF = new TEfficiency("Eff_2d_WPLoose_27_BCDEF","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPLoose_27_BCDEF = new TEfficiency("Eff_HT_WPLoose_27_BCDEF","Efficiency vs HT for WPLoose_27 (for B to F);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPLoose_27_BCDEF = new TEfficiency("Eff_numPV_WPLoose_27_BCDEF","Efficiency vs numPV for WPLoose_27 (for B to F);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPLoose_27_BCDEF = new TEfficiency("Eff_nJets_WPLoose_27_BCDEF","Efficiency vs nJets for WPLoose_27 (for B to F);nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPLoose_27_GH = new TEfficiency("Eff_pt_WPLoose_27_GH","Efficiency vs pT for WPLoose_27 (for G and H);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPLoose_27_GH = new TEfficiency("Eff_sceta_WPLoose_27_GH","Efficiency vs SuperCluster eta for WPLoose_27 (for G and H);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPLoose_27_GH = new TEfficiency("Eff_phi_WPLoose_27_GH","Efficiency vs phi for WPLoose_27 (for G and H);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPLoose_27_GH = new TEfficiency("Eff_2d_WPLoose_27_GH","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPLoose_27_GH = new TEfficiency("Eff_HT_WPLoose_27_GH","Efficiency vs HT for WPLoose_27 (for G and H);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPLoose_27_GH = new TEfficiency("Eff_numPV_WPLoose_27_GH","Efficiency vs numPV for WPLoose_27 (for G and H);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPLoose_27_GH = new TEfficiency("Eff_nJets_WPLoose_27_GH","Efficiency vs nJets for WPLoose_27 (for G and H);nJets;Efficiency",10,0,10);


    TEfficiency* Eff_pt_WPTight_27 = new TEfficiency("Eff_pt_WPTight_27","Efficiency vs pT for WPTight_27;pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_27 = new TEfficiency("Eff_sceta_WPTight_27","Efficiency vs SuperCluster eta for WPTight_27;SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_27 = new TEfficiency("Eff_phi_WPTight_27","Efficiency vs phi for WPTight_27;#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_27 = new TEfficiency("Eff_2d_WPTight_27","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_27 = new TEfficiency("Eff_HT_WPTight_27","Efficiency vs HT for WPTight_27;HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_27 = new TEfficiency("Eff_numPV_WPTight_27","Efficiency vs numPV for WPTight_27;numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_27 = new TEfficiency("Eff_nJets_WPTight_27","Efficiency vs nJets for WPTight_27;nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPTight_27_BCDEF = new TEfficiency("Eff_pt_WPTight_27_BCDEF","Efficiency vs pT for WPTight_27 (for B to F);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_27_BCDEF = new TEfficiency("Eff_sceta_WPTight_27_BCDEF","Efficiency vs SuperCluster eta for WPTight_27 (for B to F);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_27_BCDEF = new TEfficiency("Eff_phi_WPTight_27_BCDEF","Efficiency vs phi for WPTight_27 (for B to F);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_27_BCDEF = new TEfficiency("Eff_2d_WPTight_27_BCDEF","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_27_BCDEF = new TEfficiency("Eff_HT_WPTight_27_BCDEF","Efficiency vs HT for WPTight_27 (for B to F);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_27_BCDEF = new TEfficiency("Eff_numPV_WPTight_27_BCDEF","Efficiency vs numPV for WPTight_27 (for B to F);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_27_BCDEF = new TEfficiency("Eff_nJets_WPTight_27_BCDEF","Efficiency vs nJets for WPTight_27 (for B to F);nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPTight_27_GH = new TEfficiency("Eff_pt_WPTight_27_GH","Efficiency vs pT for WPTight_27 (for G and H);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_27_GH = new TEfficiency("Eff_sceta_WPTight_27_GH","Efficiency vs SuperCluster eta for WPTight_27 (for G and H);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_27_GH = new TEfficiency("Eff_phi_WPTight_27_GH","Efficiency vs phi for WPTight_27 (for G and H);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_27_GH = new TEfficiency("Eff_2d_WPTight_27_GH","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_27_GH = new TEfficiency("Eff_HT_WPTight_27_GH","Efficiency vs HT for WPTight_27 (for G and H);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_27_GH = new TEfficiency("Eff_numPV_WPTight_27_GH","Efficiency vs numPV for WPTight_27 (for G and H);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_27_GH = new TEfficiency("Eff_nJets_WPTight_27_GH","Efficiency vs nJets for WPTight_27 (for G and H);nJets;Efficiency",10,0,10);


    TEfficiency* Eff_pt_WPTight_32 = new TEfficiency("Eff_pt_WPTight_32","Efficiency vs pT for WPTight_32;pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_32 = new TEfficiency("Eff_sceta_WPTight_32","Efficiency vs SuperCluster eta for WPTight_32;SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_32 = new TEfficiency("Eff_phi_WPTight_32","Efficiency vs phi for WPTight_32;#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_32 = new TEfficiency("Eff_2d_WPTight_32","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_32 = new TEfficiency("Eff_HT_WPTight_32","Efficiency vs HT for WPTight_32;HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_32 = new TEfficiency("Eff_numPV_WPTight_32","Efficiency vs numPV for WPTight_32;numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_32 = new TEfficiency("Eff_nJets_WPTight_32","Efficiency vs nJets for WPTight_32;nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPTight_32_BCDEF = new TEfficiency("Eff_pt_WPTight_32_BCDEF","Efficiency vs pT for WPTight_32 (for B to F);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_32_BCDEF = new TEfficiency("Eff_sceta_WPTight_32_BCDEF","Efficiency vs SuperCluster eta for WPTight_32 (for B to F);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_32_BCDEF = new TEfficiency("Eff_phi_WPTight_32_BCDEF","Efficiency vs phi for WPTight_32 (for B to F);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_32_BCDEF = new TEfficiency("Eff_2d_WPTight_32_BCDEF","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_32_BCDEF = new TEfficiency("Eff_HT_WPTight_32_BCDEF","Efficiency vs HT for WPTight_32 (for B to F);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_32_BCDEF = new TEfficiency("Eff_numPV_WPTight_32_BCDEF","Efficiency vs numPV for WPTight_32 (for B to F);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_32_BCDEF = new TEfficiency("Eff_nJets_WPTight_32_BCDEF","Efficiency vs nJets for WPTight_32 (for B to F);nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPTight_32_GH = new TEfficiency("Eff_pt_WPTight_32_GH","Efficiency vs pT for WPTight_32 (for G and H);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_32_GH = new TEfficiency("Eff_sceta_WPTight_32_GH","Efficiency vs SuperCluster eta for WPTight_32 (for G and H);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_32_GH = new TEfficiency("Eff_phi_WPTight_32_GH","Efficiency vs phi for WPTight_32 (for G and H);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_32_GH = new TEfficiency("Eff_2d_WPTight_32_GH","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_32_GH = new TEfficiency("Eff_HT_WPTight_32_GH","Efficiency vs HT for WPTight_32 (for G and H);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_32_GH = new TEfficiency("Eff_numPV_WPTight_32_GH","Efficiency vs numPV for WPTight_32 (for G and H);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_32_GH = new TEfficiency("Eff_nJets_WPTight_32_GH","Efficiency vs nJets for WPTight_32 (for G and H);nJets;Efficiency",10,0,10);


    TEfficiency* Eff_pt_WPLoose_27_HT200 = new TEfficiency("Eff_pt_WPLoose_27_HT200","Efficiency vs pT for WPLoose_27_HT200;pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPLoose_27_HT200 = new TEfficiency("Eff_sceta_WPLoose_27_HT200","Efficiency vs SuperCluster eta for WPLoose_27_HT200;SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPLoose_27_HT200 = new TEfficiency("Eff_phi_WPLoose_27_HT200","Efficiency vs phi for WPLoose_27_HT200;#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPLoose_27_HT200 = new TEfficiency("Eff_2d_WPLoose_27_HT200","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPLoose_27_HT200 = new TEfficiency("Eff_HT_WPLoose_27_HT200","Efficiency vs HT for WPLoose_27_HT200;HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPLoose_27_HT200 = new TEfficiency("Eff_numPV_WPLoose_27_HT200","Efficiency vs numPV for WPLoose_27_HT200;numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPLoose_27_HT200 = new TEfficiency("Eff_nJets_WPLoose_27_HT200","Efficiency vs nJets for WPLoose_27_HT200;nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_pt_WPLoose_27_HT200_BCDEF","Efficiency vs pT for WPLoose_27_HT200 (for B to F);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_sceta_WPLoose_27_HT200_BCDEF","Efficiency vs SuperCluster eta for WPLoose_27_HT200 (for B to F);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_phi_WPLoose_27_HT200_BCDEF","Efficiency vs phi for WPLoose_27_HT200 (for B to F);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_2d_WPLoose_27_HT200_BCDEF","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_HT_WPLoose_27_HT200_BCDEF","Efficiency vs HT for WPLoose_27_HT200 (for B to F);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_numPV_WPLoose_27_HT200_BCDEF","Efficiency vs numPV for WPLoose_27_HT200 (for B to F);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_nJets_WPLoose_27_HT200_BCDEF","Efficiency vs nJets for WPLoose_27_HT200 (for B to F);nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPLoose_27_HT200_GH = new TEfficiency("Eff_pt_WPLoose_27_HT200_GH","Efficiency vs pT for WPLoose_27_HT200 (for G and H);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPLoose_27_HT200_GH = new TEfficiency("Eff_sceta_WPLoose_27_HT200_GH","Efficiency vs SuperCluster eta for WPLoose_27_HT200 (for G and H);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPLoose_27_HT200_GH = new TEfficiency("Eff_phi_WPLoose_27_HT200_GH","Efficiency vs phi for WPLoose_27_HT200 (for G and H);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPLoose_27_HT200_GH = new TEfficiency("Eff_2d_WPLoose_27_HT200_GH","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPLoose_27_HT200_GH = new TEfficiency("Eff_HT_WPLoose_27_HT200_GH","Efficiency vs HT for WPLoose_27_HT200 (for G and H);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPLoose_27_HT200_GH = new TEfficiency("Eff_numPV_WPLoose_27_HT200_GH","Efficiency vs numPV for WPLoose_27_HT200 (for G and H);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPLoose_27_HT200_GH = new TEfficiency("Eff_nJets_WPLoose_27_HT200_GH","Efficiency vs nJets for WPLoose_27_HT200 (for G and H);nJets;Efficiency",10,0,10);


    TEfficiency* Eff_pt_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_pt_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_sceta_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs SuperCluster eta for WPTight_27_OR_WPLoose_27_HT200;SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_phi_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200;#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_2d_WPTight_27_OR_WPLoose_27_HT200","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_HT_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200;HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_numPV_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs numPV for WPTight_27_OR_WPLoose_27_HT200;numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_nJets_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs nJets for WPTight_27_OR_WPLoose_27_HT200;nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_pt_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200 (for B to F);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_sceta_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency vs SuperCluster eta for WPTight_27_OR_WPLoose_27_HT200 (for B to F);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_phi_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200 (for B to F);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_2d_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_HT_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200 (for B to F);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_numPV_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency vs numPV for WPTight_27_OR_WPLoose_27_HT200 (for B to F);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_27_OR_WPLoose_27_HT200_BCDEF = new TEfficiency("Eff_nJets_WPTight_27_OR_WPLoose_27_HT200_BCDEF","Efficiency vs nJets for WPTight_27_OR_WPLoose_27_HT200 (for B to F);nJets;Efficiency",10,0,10);

    TEfficiency* Eff_pt_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_pt_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200 (for G and H);pT (GeV);Efficiency",n_ptBins,x_ptBins);
    TEfficiency* Eff_sceta_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_sceta_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency vs SuperCluster eta for WPTight_27_OR_WPLoose_27_HT200 (for G and H);SuperCluster #eta;Efficiency",n_scetaBins,x_scetaBins);
    TEfficiency* Eff_phi_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_phi_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200 (for G and H);#phi;x_scetaBins",30,-3,3);
    TEfficiency* Eff_2d_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_2d_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency;pT;SuperCluster #eta",n_ptBins,x_ptBins,n_scetaBins,x_scetaBins);
    TEfficiency* Eff_HT_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_HT_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200 (for G and H);HT (GeV);Efficiency",n_htBins,x_htBins);
    TEfficiency* Eff_numPV_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_numPV_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency vs numPV for WPTight_27_OR_WPLoose_27_HT200 (for G and H);numPV;Efficiency",15,0,75);
    TEfficiency* Eff_nJets_WPTight_27_OR_WPLoose_27_HT200_GH = new TEfficiency("Eff_nJets_WPTight_27_OR_WPLoose_27_HT200_GH","Efficiency vs nJets for WPTight_27_OR_WPLoose_27_HT200 (for G and H);nJets;Efficiency",10,0,10);

   //////////////////////////////////////////////////////////////////////////
  ///  Weights
  //////////////////////////////////////////////////////////////////////////
 
  int nentries = chain->GetEntries();
  std::cout << "\n\t Number of entries = " << nentries << std::endl;
  std::cout << "\t Max number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;
  if( jobN==1 ) firstEvent = 0;
 
  std::cout << "========  Starting Event Loop  ========" << std::endl;
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t
  
  		N_total++;
  		
	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;
	  
	  
	  	chain->GetEntry(ievt);
		
		//check Primary Vertices 
		
		bool GoodFirstPV = eve->goodFirstVertex_;
		if(!GoodFirstPV) continue;
		
		
		//Grab Specific Lepton information from trees
		
		int event_nr = eve->evt_;
		vint lepton_charge = eve->lepton_charge_;
        vdouble lepton_pt = eve->lepton_pt_;
        vdouble lepton_px = eve->lepton_px_;
        vdouble lepton_py = eve->lepton_py_;
        vdouble lepton_pz = eve->lepton_pz_;
		vdouble lepton_eta = eve->lepton_eta_;
        vdouble lepton_sc_eta = eve->lepton_sc_eta_;
		vdouble lepton_phi = eve->lepton_phi_;
		vdouble lepton_energy = eve->lepton_energy_;
		vdouble lepton_rel_Iso = eve->lepton_relIso_;
		vdouble lepton_rel_IsoR04 = eve->lepton_relIsoR04_;
		
		vint is_IDTight = eve->lepton_is_IDTight_;
		vint isTight = eve->lepton_isTight_;
		vint isLoose = eve->lepton_isLoose_;
		vint isMuon = eve->lepton_isMuon_; 
		vint isTrigMVAM = eve->lepton_isTrigMVAM_;
		vint isNonTrigMVAM = eve->lepton_isNonTrigMVAM_;
		vint isTrigCutM = eve->lepton_isTrigCutM_;

        int opp_charge= 0;
        double mass_leplep = 0;
	   
		vdouble jet_pt = eve->jet_pt_;
		vdouble jet_eta = eve->jet_eta_;
		vdouble jet_csv = eve->jet_csv_;

        // FOR MC
        int isData = eve->isData_;
        // FOR DATA
        //int isData = 1;

		double gen_weight = eve->gen_weight_;
  		double csv_weight = eve->csv_weight_;
  		double PU_weight = eve->PU_weight_ ;
        if(!isData){
            double PU_weight_B = eve->PU_weight_B_;
            double PU_weight_C = eve->PU_weight_C_;
            double PU_weight_D = eve->PU_weight_D_;
            double PU_weight_E = eve->PU_weight_E_;
            double PU_weight_F = eve->PU_weight_F_;
            double PU_weight_G = eve->PU_weight_G_;
            double PU_weight_H = eve->PU_weight_H_;
            double PU_weight_BCDEF = eve->PU_weight_BCDEF_;
            double PU_weight_GH = eve->PU_weight_GH_;
        }
        else{
            double PU_weight_B = 1;
            double PU_weight_C = 1;
            double PU_weight_D = 1;
            double PU_weight_E = 1;
            double PU_weight_F = 1;
            double PU_weight_G = 1;
            double PU_weight_H = 1;
            double PU_weight_BCDEF = 1;
            double PU_weight_GH = 1;
        }
  		double PDF_weight = eve->PDF_weight_;
  		double Q2_weight = eve->Q2_weight_;
  		vdouble lepton_id_sf = eve->lepton_id_sf_;
  		vdouble lepton_iso_sf = eve->lepton_iso_sf_;
  		vdouble lepton_gsf_sf = eve->lepton_gsf_sf_;
  		vdouble lepton_trig_sf = eve->lepton_trig_sf_;
  		vdouble lepton_hip_sf = eve->lepton_hip_sf_;

        vdouble hltEle27WPLooseGsf_pt = eve->hltEle27WPLooseGsf_pt_;
        vdouble hltEle27WPLooseGsf_eta = eve->hltEle27WPLooseGsf_eta_;
        vdouble hltEle27WPLooseGsf_phi = eve->hltEle27WPLooseGsf_phi_;

        vdouble hltEle27WPTightGsf_pt = eve->hltEle27WPTightGsf_pt_;
        vdouble hltEle27WPTightGsf_eta = eve->hltEle27WPTightGsf_eta_;
        vdouble hltEle27WPTightGsf_phi = eve->hltEle27WPTightGsf_phi_;

        vdouble hltEle32WPTightGsf_pt = eve->hltEle32WPTightGsf_pt_;
        vdouble hltEle32WPTightGsf_eta = eve->hltEle32WPTightGsf_eta_;
        vdouble hltEle32WPTightGsf_phi = eve->hltEle32WPTightGsf_phi_;

        vdouble hltEle27WPLooseHT200Gsf_pt = eve->hltEle27WPLooseHT200Gsf_pt_;
        vdouble hltEle27WPLooseHT200Gsf_eta = eve->hltEle27WPLooseHT200Gsf_eta_;
        vdouble hltEle27WPLooseHT200Gsf_phi = eve->hltEle27WPLooseHT200Gsf_phi_;
  		
  		double tot_weight = 1;          // Data and MC (all combined)
        double tot_weight_BCDEF = 1;    // MC (B to F)
        double tot_weight_GH = 1;       // MC (G and H)

        sum_gen_weight = sum_gen_weight + gen_weight;
		
		double HT = 0;
		int numPV = eve->numPVs_;
		
		int numJets = 0;
		int numTags = 0;
		
		//Organize the lepton information into a LorentzVector and Count Tight and Loose Leptons for selection
		
		vvdouble vvLEPTON;
		int numLooseEle=0;
		int numMu = 0;
		int numTightEle=0;
		for(int i=0;i<int(lepton_pt.size());i++){
		
			 if(is_IDTight[i]==1) {
				numMu++;			
			 }
			
			else {
			    if (isNonTrigMVAM[i]==1 && lepton_rel_Iso[i]<0.15 ) {  			 // for Non-Triggering MVA electron ID
				//if (isTrigCutM[i]==1 && lepton_rel_Iso[i]<0.15 ) {  			 // for Cut based electron ID
				//if (isTrigMVAM[i]==1 && lepton_rel_Iso[i]<0.15 ) {               // for Triggering MVA electron ID
					numLooseEle++;
					vdouble vlepton;
					vlepton.push_back(lepton_pt[i]);
					vlepton.push_back(lepton_eta[i]);
					vlepton.push_back(lepton_phi[i]);
					vlepton.push_back(lepton_energy[i]);
                    vlepton.push_back(lepton_charge[i]);
                    vlepton.push_back(lepton_px[i]);
                    vlepton.push_back(lepton_py[i]);
                    vlepton.push_back(lepton_pz[i]);
				    //Trigger SF not applied
					vlepton.push_back(lepton_id_sf[i]*lepton_iso_sf[i]*lepton_gsf_sf[i]*lepton_hip_sf[i]);
                    vlepton.push_back(lepton_sc_eta[i]);
					vvLEPTON.push_back(vlepton);
					if ( lepton_pt[i]>30 && fabs(lepton_eta[i])<2.1  ) {
						numTightEle++;
					}
				}
			}	
		}
		
		double jet1_pt, jet1_csv;
		jet1_pt = jet1_csv = -9999;
		
		for(int i=0;i<int(jet_pt.size());i++) {
			if ( jet_pt[i]>30 && fabs(jet_eta[i])<2.4 ) {
				numJets++;
				HT = HT + jet_pt[i];
				if(jet_pt[i]>=jet1_pt) {
					jet1_pt = jet_pt[i];
					jet1_csv = jet_csv[i];
				}						
				if ( jet_csv[i] > 0.8 )
					numTags++;
			}
		}
		
		//Selection
		
		
		if(numMu!=0)continue;
		if(numLooseEle!=2)continue;
		if(numTightEle!=2)continue;
		
        if((vvLEPTON[0][4]*vvLEPTON[1][4]) == -1)
            opp_charge = 1;
        double E = vvLEPTON[0][3] + vvLEPTON[1][3];
        double px = vvLEPTON[0][5] + vvLEPTON[1][5];
        double py = vvLEPTON[0][6] + vvLEPTON[1][6];
        double pz = vvLEPTON[0][7] + vvLEPTON[1][7];
        mass_leplep = sqrt( E*E - px*px - py*py - pz*pz );
	  
		// Opposite charge and dilepton mass window
		if(opp_charge!=1) continue;
		//if( mass_leplep < 88.7 || mass_leplep > 93.7) continue;   // 5 GeV window
		if( mass_leplep < 86.2 || mass_leplep > 96.2) continue;   // 10 GeV window
        if( mass_leplep < 81.2 || mass_leplep > 101.2) continue;   // 20 GeV window
		
		//Additionally we require 4Jets and 2Btag Jets
			
		//if( numJets<2 ) continue;
		//if( numJets<4 || numTags<2 ) continue;	
		
		N_eve++;
		
		int pass_WPLoose_27 = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_;
		int pass_WPLoose_27_HT200 = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;
		int pass_WPTight_27 = eve->pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_;
		int pass_WPTight_32 = eve->pass_HLT_Ele32_eta2p1_WPTight_Gsf_v_;
		int pass_WPTight_27_OR_WPLoose_27_HT200 = 0;
		if (pass_WPLoose_27_HT200 ==1 || pass_WPTight_27 ==1 ) 
			pass_WPTight_27_OR_WPLoose_27_HT200 = 1;

        int tag, probe;
        double deta, dphi, dR;
        bool trigger_pass;

		// PDF Weight, Q2 Weight not applied
	    //tot_weight = vvLEPTON[0][8]*vvLEPTON[1][8]*csv_weight*PU_weight*gen_weight;
        tot_weight = vvLEPTON[0][8]*vvLEPTON[1][8]*PU_weight;
        tot_weight_BCDEF = vvLEPTON[0][8]*vvLEPTON[1][8]*PU_weight_BCDEF;
        tot_weight_GH = vvLEPTON[0][8]*vvLEPTON[1][8]*PU_weight_GH;

        // Fill Histograms

        if(pass_WPLoose_27==1) {

            tag = probe = -9999;
            for(int i=0; i<int(hltEle27WPLooseGsf_eta.size());i++){
                deta = hltEle27WPLooseGsf_eta[i] - vvLEPTON[0][1];
                dphi = hltEle27WPLooseGsf_phi[i] - vvLEPTON[0][2];
                dR = sqrt( deta*deta + dphi*dphi );
                if(dR<0.1){
                    tag = 0;
                    probe = 1;
                    break;
                }
            }
            if(tag!=0){
                for(int i=0; i<int(hltEle27WPLooseGsf_eta.size());i++){
                    deta = hltEle27WPLooseGsf_eta[i] - vvLEPTON[1][1];
                    dphi = hltEle27WPLooseGsf_phi[i] - vvLEPTON[1][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 1;
                        probe = 0;
                        break;
                    }
                }
            }

            if(tag==0 || tag==1){
                trigger_pass = 0;
                for(int i=0; i<int(hltEle27WPLooseGsf_eta.size());i++){
                    deta = hltEle27WPLooseGsf_eta[i] - vvLEPTON[probe][1];
                    dphi = hltEle27WPLooseGsf_phi[i] - vvLEPTON[probe][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        trigger_pass = 1;
                        break;
                    }
                }

                Eff_pt_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0]);
                Eff_sceta_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][9]);
                Eff_phi_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][2]);
                Eff_2d_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,HT);
                Eff_numPV_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,numPV);
                Eff_nJets_WPLoose_27-> FillWeighted(trigger_pass,tot_weight,nJets);

                Eff_pt_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0]);
                Eff_sceta_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][9]);
                Eff_phi_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][2]);
                Eff_2d_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,HT);
                Eff_numPV_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,numPV);
                Eff_nJets_WPLoose_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,nJets);

                Eff_pt_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0]);
                Eff_sceta_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][9]);
                Eff_phi_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][2]);
                Eff_2d_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,HT);
                Eff_numPV_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,numPV);
                Eff_nJets_WPLoose_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,nJets);

                effi_num_loose_27 = effi_num_loose_27 + trigger_pass*tot_weight;
                effi_den_loose_27 = effi_den_loose_27 + tot_weight;
                effi_num_loose_27_BCDEF = effi_num_loose_27_BCDEF + trigger_pass*tot_weight_BCDEF;
                effi_den_loose_27_BCDEF = effi_den_loose_27_BCDEF + tot_weight_BCDEF;
                effi_num_loose_27_GH = effi_num_loose_27_GH + trigger_pass*tot_weight_GH;
                effi_den_loose_27_GH = effi_den_loose_27_GH + tot_weight_GH;

            }

			N_eve_loose_27++;
		}
			
		if(pass_WPTight_27==1) {

            tag = probe = -9999;
            for(int i=0; i<int(hltEle27WPTightGsf_eta.size());i++){
                deta = hltEle27WPTightGsf_eta[i] - vvLEPTON[0][1];
                dphi = hltEle27WPTightGsf_phi[i] - vvLEPTON[0][2];
                dR = sqrt( deta*deta + dphi*dphi );
                if(dR<0.1){
                    tag = 0;
                    probe = 1;
                    break;
                }
            }
            if(tag!=0){
                for(int i=0; i<int(hltEle27WPTightGsf_eta.size());i++){
                    deta = hltEle27WPTightGsf_eta[i] - vvLEPTON[1][1];
                    dphi = hltEle27WPTightGsf_phi[i] - vvLEPTON[1][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 1;
                        probe = 0;
                        break;
                    }
                }
            }

            if(tag==0 || tag==1){
                trigger_pass = 0;
                for(int i=0; i<int(hltEle27WPTightGsf_eta.size());i++){
                    deta = hltEle27WPTightGsf_eta[i] - vvLEPTON[probe][1];
                    dphi = hltEle27WPTightGsf_phi[i] - vvLEPTON[probe][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        trigger_pass = 1;
                        break;
                    }
                }

                Eff_pt_WPTight_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][9]);
                Eff_phi_WPTight_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][2]);
                Eff_2d_WPTight_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_27-> FillWeighted(trigger_pass,tot_weight,HT);
                Eff_numPV_WPTight_27-> FillWeighted(trigger_pass,tot_weight,numPV);
                Eff_nJets_WPTight_27-> FillWeighted(trigger_pass,tot_weight,nJets);

                Eff_pt_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][9]);
                Eff_phi_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][2]);
                Eff_2d_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,HT);
                Eff_numPV_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,numPV);
                Eff_nJets_WPTight_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,nJets);

                Eff_pt_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][9]);
                Eff_phi_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][2]);
                Eff_2d_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,HT);
                Eff_numPV_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,numPV);
                Eff_nJets_WPTight_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,nJets);

                effi_num_tight_27 = effi_num_tight_27 + trigger_pass*tot_weight;
                effi_den_tight_27 = effi_den_tight_27 + tot_weight;
                effi_num_tight_27_BCDEF = effi_num_tight_27_BCDEF + trigger_pass*tot_weight_BCDEF;
                effi_den_tight_27_BCDEF = effi_den_tight_27_BCDEF + tot_weight_BCDEF;
                effi_num_tight_27_GH = effi_num_tight_27_GH + trigger_pass*tot_weight_GH;
                effi_den_tight_27_GH = effi_den_tight_27_GH + tot_weight_GH;
                
            }

			N_eve_tight_27++;
		}
		
		if(pass_WPTight_32==1) {

            tag = probe = -9999;
            for(int i=0; i<int(hltEle32WPTightGsf_eta.size());i++){
                deta = hltEle32WPTightGsf_eta[i] - vvLEPTON[0][1];
                dphi = hltEle32WPTightGsf_phi[i] - vvLEPTON[0][2];
                dR = sqrt( deta*deta + dphi*dphi );
                if(dR<0.1){
                    tag = 0;
                    probe = 1;
                    break;
                }
            }
            if(tag!=0){
                for(int i=0; i<int(hltEle32WPTightGsf_eta.size());i++){
                    deta = hltEle32WPTightGsf_eta[i] - vvLEPTON[1][1];
                    dphi = hltEle32WPTightGsf_phi[i] - vvLEPTON[1][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 1;
                        probe = 0;
                        break;
                    }
                }
            }

            if(tag==0 || tag==1){
                trigger_pass = 0;
                for(int i=0; i<int(hltEle32WPTightGsf_eta.size());i++){
                    deta = hltEle32WPTightGsf_eta[i] - vvLEPTON[probe][1];
                    dphi = hltEle32WPTightGsf_phi[i] - vvLEPTON[probe][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        trigger_pass = 1;
                        break;
                    }
                }

                Eff_pt_WPTight_32-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_32-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][9]);
                Eff_phi_WPTight_32-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][2]);
                Eff_2d_WPTight_32-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_32-> FillWeighted(trigger_pass,tot_weight,HT);
                Eff_numPV_WPTight_32-> FillWeighted(trigger_pass,tot_weight,numPV);
                Eff_nJets_WPTight_32-> FillWeighted(trigger_pass,tot_weight,nJets);

                Eff_pt_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][9]);
                Eff_phi_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][2]);
                Eff_2d_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,HT);
                Eff_numPV_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,numPV);
                Eff_nJets_WPTight_32_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,nJets);

                Eff_pt_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][9]);
                Eff_phi_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][2]);
                Eff_2d_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,HT);
                Eff_numPV_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,numPV);
                Eff_nJets_WPTight_32_GH-> FillWeighted(trigger_pass,tot_weight_GH,nJets);

                effi_num_tight_32 = effi_num_tight_32 + trigger_pass*tot_weight;
                effi_den_tight_32 = effi_den_tight_32 + tot_weight;
                effi_num_tight_32_BCDEF = effi_num_tight_32_BCDEF + trigger_pass*tot_weight_BCDEF;
                effi_den_tight_32_BCDEF = effi_den_tight_32_BCDEF + tot_weight_BCDEF;
                effi_num_tight_32_GH = effi_num_tight_32_GH + trigger_pass*tot_weight_GH;
                effi_den_tight_32_GH = effi_den_tight_32_GH + tot_weight_GH;

            }

			N_eve_tight_32++;
		}
		
		if(pass_WPLoose_27_HT200==1) {
		
            tag = probe = -9999;
            for(int i=0; i<int(hltEle27WPLooseHT200Gsf_eta.size());i++){
                deta = hltEle27WPLooseHT200Gsf_eta[i] - vvLEPTON[0][1];
                dphi = hltEle27WPLooseHT200Gsf_phi[i] - vvLEPTON[0][2];
                dR = sqrt( deta*deta + dphi*dphi );
                if(dR<0.1){
                    tag = 0;
                    probe = 1;
                    break;
                }
            }
            if(tag!=0){
                for(int i=0; i<int(hltEle27WPLooseHT200Gsf_eta.size());i++){
                    deta = hltEle27WPLooseHT200Gsf_eta[i] - vvLEPTON[1][1];
                    dphi = hltEle27WPLooseHT200Gsf_phi[i] - vvLEPTON[1][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 1;
                        probe = 0;
                        break;
                    }
                }
            }

            if(tag==0 || tag==1){
                trigger_pass = 0;
                for(int i=0; i<int(hltEle27WPLooseHT200Gsf_eta.size());i++){
                    deta = hltEle27WPLooseHT200Gsf_eta[i] - vvLEPTON[probe][1];
                    dphi = hltEle27WPLooseHT200Gsf_phi[i] - vvLEPTON[probe][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        trigger_pass = 1;
                        break;
                    }
                }

                Eff_pt_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0]);
                Eff_sceta_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][9]);
                Eff_phi_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][2]);
                Eff_2d_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,HT);
                Eff_numPV_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,numPV);
                Eff_nJets_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,nJets);

                Eff_pt_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0]);
                Eff_sceta_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][9]);
                Eff_phi_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][2]);
                Eff_2d_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,HT);
                Eff_numPV_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,numPV);
                Eff_nJets_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,nJets);

                Eff_pt_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0]);
                Eff_sceta_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][9]);
                Eff_phi_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][2]);
                Eff_2d_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,HT);
                Eff_numPV_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,numPV);
                Eff_nJets_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,nJets);

                effi_num_loose_ht200_27 = effi_num_loose_ht200_27 + trigger_pass*tot_weight;
                effi_den_loose_ht200_27 = effi_den_loose_ht200_27 + tot_weight;
                effi_num_loose_ht200_27_BCDEF = effi_num_loose_ht200_27_BCDEF + trigger_pass*tot_weight_BCDEF;
                effi_den_loose_ht200_27_BCDEF = effi_den_loose_ht200_27_BCDEF + tot_weight_BCDEF;
                effi_num_loose_ht200_27_GH = effi_num_loose_ht200_27_GH + trigger_pass*tot_weight_GH;
                effi_den_loose_ht200_27_GH = effi_den_loose_ht200_27_GH + tot_weight_GH;
                
            }
			
			N_eve_loose_27_ht200++;
		}
		
		if(pass_WPTight_27_OR_WPLoose_27_HT200==1) {
		
            tag = probe = -9999;
            for(int i=0; i<int(hltEle27WPLooseHT200Gsf_eta.size());i++){
                deta = hltEle27WPLooseHT200Gsf_eta[i] - vvLEPTON[0][1];
                dphi = hltEle27WPLooseHT200Gsf_phi[i] - vvLEPTON[0][2];
                dR = sqrt( deta*deta + dphi*dphi );
                if(dR<0.1){
                    tag = 0;
                    probe = 1;
                    break;
                }
            }
            if(tag!=0){
                for(int i=0; i<int(hltEle27WPTightGsf_eta.size());i++){
                    deta = hltEle27WPTightGsf_eta[i] - vvLEPTON[0][1];
                    dphi = hltEle27WPTightGsf_phi[i] - vvLEPTON[0][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 0;
                        probe = 1;
                        break;
                    }
                }
            }

            if(tag!=0){
                for(int i=0; i<int(hltEle27WPLooseHT200Gsf_eta.size());i++){
                    deta = hltEle27WPLooseHT200Gsf_eta[i] - vvLEPTON[1][1];
                    dphi = hltEle27WPLooseHT200Gsf_phi[i] - vvLEPTON[1][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 1;
                        probe = 0;
                        break;
                    }
                }
            }
            if(tag!=0){
                for(int i=0; i<int(hltEle27WPTightGsf_eta.size());i++){
                    deta = hltEle27WPTightGsf_eta[i] - vvLEPTON[1][1];
                    dphi = hltEle27WPTightGsf_phi[i] - vvLEPTON[1][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        tag = 1;
                        probe = 0;
                        break;
                    }
                }
            }

            if(tag==0 || tag==1){
                trigger_pass = 0;
                for(int i=0; i<int(hltEle27WPLooseHT200Gsf_eta.size());i++){
                    deta = hltEle27WPLooseHT200Gsf_eta[i] - vvLEPTON[probe][1];
                    dphi = hltEle27WPLooseHT200Gsf_phi[i] - vvLEPTON[probe][2];
                    dR = sqrt( deta*deta + dphi*dphi );
                    if(dR<0.1){
                        trigger_pass = 1;
                        break;
                    }
                }

                Eff_pt_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][9]);
                Eff_phi_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][2]);
                Eff_2d_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,HT);
                Eff_numPV_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,numPV);
                Eff_nJets_WPTight_27_OR_WPLoose_HT200_27-> FillWeighted(trigger_pass,tot_weight,nJets);

                Eff_pt_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][9]);
                Eff_phi_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][2]);
                Eff_2d_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,HT);
                Eff_numPV_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,numPV);
                Eff_nJets_WPTight_27_OR_WPLoose_HT200_27_BCDEF-> FillWeighted(trigger_pass,tot_weight_BCDEF,nJets);

                Eff_pt_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0]);
                Eff_sceta_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][9]);
                Eff_phi_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][2]);
                Eff_2d_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,vvLEPTON[probe][0],vvLEPTON[probe][9]);
                Eff_HT_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,HT);
                Eff_numPV_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,numPV);
                Eff_nJets_WPTight_27_OR_WPLoose_HT200_27_GH-> FillWeighted(trigger_pass,tot_weight_GH,nJets);

                effi_num_tight_27_or_loose_ht200_27 = effi_num_tight_27_or_loose_ht200_27 + trigger_pass*tot_weight;
                effi_den_tight_27_or_loose_ht200_27 = effi_den_tight_27_or_loose_ht200_27 + tot_weight;
                effi_num_tight_27_or_loose_ht200_27_BCDEF = effi_num_tight_27_or_loose_ht200_27_BCDEF + trigger_pass*tot_weight_BCDEF;
                effi_den_tight_27_or_loose_ht200_27_BCDEF = effi_den_tight_27_or_loose_ht200_27_BCDEF + tot_weight_BCDEF;
                effi_num_tight_27_or_loose_ht200_27_GH = effi_num_tight_27_or_loose_ht200_27_GH + trigger_pass*tot_weight_GH;
                effi_den_tight_27_or_loose_ht200_27_GH = effi_den_tight_27_or_loose_ht200_27_GH + tot_weight_GH;
                
            }
			
			N_eve_tight_27_OR_loose_27_ht200++;
		}
				
   }
   
  t=clock()-t;



  effi_loose_27 = effi_num_loose_27 / effi_den_loose_27;
  effi_loose_27_BCDEF = effi_num_loose_27_BCDEF / effi_den_loose_27_BCDEF;
  effi_loose_27_GH = effi_num_loose_27_GH / effi_den_loose_27_GH;
  effi_tight_27 = effi_num_tight_27 / effi_den_tight_27;
  effi_tight_27_BCDEF = effi_num_tight_27_BCDEF / effi_den_tight_27_BCDEF;
  effi_tight_27_GH = effi_num_tight_27_GH / effi_den_tight_27_GH;
  effi_tight_32 = effi_num_tight_32 / effi_den_tight_32;
  effi_tight_32_BCDEF = effi_num_tight_32_BCDEF / effi_den_tight_32_BCDEF;
  effi_tight_32_GH = effi_num_tight_32_GH / effi_den_tight_32_GH;
  effi_loose_27_ht200 = effi_num_loose_27_ht200 / effi_den_loose_27_ht200;
  effi_loose_27_ht200_BCDEF = effi_num_loose_27_ht200_BCDEF / effi_den_loose_27_ht200_BCDEF;
  effi_loose_27_ht200_GH = effi_num_loose_27_ht200_GH / effi_den_loose_27_ht200_GH;
  effi_tight_27_or_loose_27_ht200 = effi_num_tight_27_or_loose_27_ht200 / effi_den_tight_27_or_loose_27_ht200;
  effi_tight_27_or_loose_27_ht200_BCDEF = effi_num_tight_27_or_loose_27_ht200_BCDEF / effi_den_tight_27_or_loose_27_ht200_BCDEF;
  effi_tight_27_or_loose_27_ht200_GH = effi_num_tight_27_or_loose_27_ht200_GH / effi_den_tight_27_or_loose_27_ht200_GH;

  std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
  std::cout<<"**********************************************************************************************\n";
  std::cout<<"Total No. of events : "<<N_total<<"\n";
  std::cout<<"No. of events passing event selection only : "<<N_eve<<"\n\n";
  std::cout<<"No. of events passing event selection plus WP_Loose_Ele27 Trigger : "<<N_eve_loose_27<<"\n";
  std:cout<<"Combined Efficiency B to H :"<<effi_loose_27<<"\n";
  std:cout<<"Efficiency B to F :"<<effi_loose_27_BCDEF<<"\n";
  std:cout<<"Efficiency G and H :"<<effi_loose_27_GH<<"\n\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele27 Trigger : "<<N_eve_tight_27<<"\n";
  std:cout<<"Combined Efficiency B to H :"<<effi_tight_27<<"\n";
  std:cout<<"Efficiency B to F :"<<effi_tight_27_BCDEF<<"\n";
  std:cout<<"Efficiency G and H :"<<effi_tight_27_GH<<"\n\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele32 Trigger : "<<N_eve_tight_32<<"\n";
  std:cout<<"Combined Efficiency B to H :"<<effi_tight_32<<"\n";
  std:cout<<"Efficiency B to F :"<<effi_tight_32_BCDEF<<"\n";
  std:cout<<"Efficiency G and H :"<<effi_tight_32_GH<<"\n\n";
  std::cout<<"No. of events passing event selection plus WP_Loose_Ele27_HT200 Trigger : "<<N_eve_loose_27_ht200<<"\n";
  std:cout<<"Combined Efficiency B to H :"<<effi_loose_27_ht200<<"\n";
  std:cout<<"Efficiency B to F :"<<effi_loose_27_ht200_BCDEF<<"\n";
  std:cout<<"Efficiency G and H :"<<effi_loose_27_ht200_GH<<"\n\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele27 OR WP_Loose_Ele27_HT200 Trigger : "<<N_eve_tight_27_OR_loose_27_ht200<<"\n";
  std:cout<<"Combined Efficiency B to H :"<<effi_tight_27_or_loose_27_ht200<<"\n";
  std:cout<<"Efficiency B to F :"<<effi_tight_27_or_loose_27_ht200_BCDEF<<"\n";
  std:cout<<"Efficiency G and H :"<<effi_tight_27_or_loose_27_ht200_GH<<"\n\n";
  std::cout<<"Sum of Generator Weights for sample : "<<sum_gen_weight<<"\n";
  std::cout<<"**********************************************************************************************\n";
  histofile.Write();
  histofile.Close();
  
}


