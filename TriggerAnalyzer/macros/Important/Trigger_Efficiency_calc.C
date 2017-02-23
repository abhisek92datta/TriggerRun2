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
#include "TEfficiency.h"

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

void Trigger_Efficiency_calc( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

  clock_t t;
  t = clock();
  
  int N_total = 0;
  int N_eve = 0;
  int N_eve_control = 0;
  int N_eve_tight_27 = 0;
  int N_eve_tight_32 = 0;
  int N_eve_loose_27_ht200 = 0;
  int N_eve_tight_27_OR_loose_27_ht200 = 0;
  
  std::cout << "   ===> load the root files! " << std::endl;
  
  double wgt=1;
  
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string treefilename = "/eos/uscms/store/user/adatta/Trigger_Analysis/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/triggerTree_MC_80X_tthbb/161104_195458/0000/trigger_analyzer_*.root";
  //std::string treefilename2 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0001/trigger_analyzer_*.root";
  //std::string treefilename3 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0002/trigger_analyzer_*.root";
  
  std::string s_end = "_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = ".root";

  std::string histofilename = "Efficiency_Trigger_histos" + s_end;
 
  std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

  TChain *chain = new TChain("triggeranalyzer/triggerTree");
 
  chain->Add(treefilename.c_str());
  //chain->Add(treefilename2.c_str());
  //chain->Add(treefilename3.c_str());
  
  ////////////////////////////////////
   
  std::vector<std::string> cat_labels;
  cat_labels.push_back("4j2t"); //0
  cat_labels.push_back("5j2t"); //1
  cat_labels.push_back("6j2t"); //2
  cat_labels.push_back("4j3t"); //3
  cat_labels.push_back("5j3t"); //4
  cat_labels.push_back("6j3t"); //5
  cat_labels.push_back("4j4t"); //6
  cat_labels.push_back("5j4t"); //7
  cat_labels.push_back("6j4t"); //8
  
  int NumCat = int(cat_labels.size());

  /////////////
   
  triggerStudyEventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();

  //Histograms
	 
  TH1::SetDefaultSumw2();

  //Basic Category Yields//
	 
  TH1D* h_catyield;
  TH1D* h_catyield_CrossEleTrig;
  TH1D* h_catyield_SingleEleTrig;	
	 
  h_catyield = new TH1D((std::string("h_catyield")).c_str(), "Category Yield", NumCat, 0, NumCat );
  h_catyield_SingleEleTrig = new TH1D((std::string("h_catyield_SingleEleTrig")).c_str(), "Category Yield Single", NumCat, 0, NumCat );
  h_catyield_CrossEleTrig = new TH1D((std::string("h_catyield_CrossEleTrig")).c_str(), "Category Yield Cross", NumCat, 0, NumCat );
	
  for(int b=0; b<NumCat; b++){
		h_catyield->GetXaxis()->SetBinLabel(b+1,cat_labels[b].c_str());
		h_catyield_SingleEleTrig->GetXaxis()->SetBinLabel(b+1,cat_labels[b].c_str());
		h_catyield_CrossEleTrig->GetXaxis()->SetBinLabel(b+1,cat_labels[b].c_str());
  }
	
  //TEfficiency things!
  //So the idea for TEfficiency is you are creating a histogram that looks bin by bin like (Passed Events)/(All Events) 
	
  int n_ptBins;
  double x_ptBins[200];	
  int s=0;
  int n_bin=0;
  while(s<100){
  	x_ptBins[n_bin++]=s;
  	s=s+2;
  }
  while(s<150){
  	x_ptBins[n_bin++]=s;
  	s=s+5;
  }
  while(s<200){
  	x_ptBins[n_bin++]=s;
  	s=s+10;
  }
  while(s<=300){
  	x_ptBins[n_bin++]=s;
  	s=s+25;
  }
  n_ptBins = n_bin-1;
  
  s=n_bin=0;
  int n_htBins;
  double x_htBins[200];	
  while(s<200){
  	x_htBins[n_bin++]=s;
  	s=s+20;
  }
  while(s<600){
  	x_htBins[n_bin++]=s;
  	s=s+10;
  }
  while(s<700){
  	x_htBins[n_bin++]=s;
  	s=s+20;
  }
  while(s<=1000){
  	x_htBins[n_bin++]=s;
  	s=s+25;
  }
  n_htBins = n_bin-1;
	
  //TEfficiency* Eff_pt_WPTight_27 = new TEfficiency("Eff_pt_WPTight_27","Efficiency vs pT for WPTight_27;pT (GeV);(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",150,0,300);
  TEfficiency* Eff_pt_WPTight_27 = new TEfficiency("Eff_pt_WPTight_27","Efficiency vs pT for WPTight_27;pT (GeV);(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_ptBins,x_ptBins);
  TEfficiency* Eff_eta_WPTight_27 = new TEfficiency("Eff_eta_WPTight_27","Efficiency vs eta for WPTight_27;#eta;(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPTight_27 = new TEfficiency("Eff_phi_WPTight_27","Efficiency vs phi for WPTight_27;#phi;(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPTight_27 = new TEfficiency("Eff_2d_WPTight_27","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  //TEfficiency* Eff_HT_WPTight_27 = new TEfficiency("Eff_HT_WPTight_27","Efficiency vs HT for WPTight_27;HT (GeV);(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",200,0,1000);
  TEfficiency* Eff_HT_WPTight_27 = new TEfficiency("Eff_HT_WPTight_27","Efficiency vs HT for WPTight_27;HT (GeV);(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_htBins,x_htBins);	
  TEfficiency* Eff_numPV_WPTight_27 = new TEfficiency("Eff_numPV_WPTight_27","Efficiency vs numPV for WPTight_27;numPV;(WPTight_27 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);		
	
  TEfficiency* Eff_pt_WPTight_32 = new TEfficiency("Eff_pt_WPTight_32","Efficiency vs pT for WPTight_32;pT (GeV);(WPTight_32 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_ptBins,x_ptBins);
  TEfficiency* Eff_eta_WPTight_32 = new TEfficiency("Eff_eta_WPTight_32","Efficiency vs eta for WPTight_32;#eta;(WPTight_32 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPTight_32 = new TEfficiency("Eff_phi_WPTight_32","Efficiency vs phi for WPTight_32;#phi;(WPTight_32 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPTight_32 = new TEfficiency("Eff_2d_WPTight_32","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  TEfficiency* Eff_HT_WPTight_32 = new TEfficiency("Eff_HT_WPTight_32","Efficiency vs HT for WPTight_32;HT (GeV);(WPTight_32 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_htBins,x_htBins);	
  TEfficiency* Eff_numPV_WPTight_32 = new TEfficiency("Eff_numPV_WPTight_32","Efficiency vs numPV for WPTight_32;numPV;(WPTight_32 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);	
	
  TEfficiency* Eff_pt_WPLoose_27_HT200 = new TEfficiency("Eff_pt_WPLoose_27_HT200","Efficiency vs pT for WPLoose_27_HT200;pT (GeV);(WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_ptBins,x_ptBins);
  TEfficiency* Eff_eta_WPLoose_27_HT200 = new TEfficiency("Eff_eta_WPLoose_27_HT200","Efficiency vs eta for WPLoose_27_HT200;#eta;(WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPLoose_27_HT200 = new TEfficiency("Eff_phi_WPLoose_27_HT200","Efficiency vs phi for WPLoose_27_HT200;#phi;(WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPLoose_27_HT200 = new TEfficiency("Eff_2d_WPLoose_27_HT200","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  TEfficiency* Eff_HT_WPLoose_27_HT200 = new TEfficiency("Eff_HT_WPLoose_27_HT200","Efficiency vs HT for WPLoose_27_HT200;HT (GeV);(WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_htBins,x_htBins);	
  TEfficiency* Eff_numPV_WPLoose_27_HT200 = new TEfficiency("Eff_numPV_WPLoose_27_HT200","Efficiency vs numPV for WPLoose_27_HT200;numPV;(WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);	
  
  TEfficiency* Eff_pt_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_pt_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);(WPTight_27_OR_WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_ptBins,x_ptBins);
  TEfficiency* Eff_eta_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_eta_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs eta for WPTight_27_OR_WPLoose_27_HT200;#eta;(WPTight_27_OR_WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_phi_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200;#phi;(WPTight_27_OR_WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_2d_WPTight_27_OR_WPLoose_27_HT200","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  TEfficiency* Eff_HT_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_HT_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200;HT (GeV);(WPTight_27_OR_WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",n_htBins,x_htBins);	
  TEfficiency* Eff_numPV_WPTight_27_OR_WPLoose_27_HT200 = new TEfficiency("Eff_numPV_WPTight_27_OR_WPLoose_27_HT200","Efficiency vs numPV for WPTight_27_OR_WPLoose_27_HT200;numPV;(WPTight_27_OR_WPLoose_27_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);		
  	
  //Sometimes you are asked to do specific binning for your efficiency plots (in the case of making Trigger Scale Factors)
	
  //int n_ptBins = 5;
  //double x_ptBins[6] = {15.0, 25.0, 35.0, 45.0, 55.0, 200.0};
	
  //int n_etaBins = 10;
  //double x_etaBins[11] = {-2.5,-2.1, -1.5660, -1.4442, -0.8, 0, 0.8, 1.4442, 1.5660, 2.1, 2.5};
	
  //TEfficiency* Eff_bin_pt = new TEfficiency("Eff_bin_pt","my efficiency;x;#epsilon",n_ptBins,x_ptBins);
  //TEfficiency* Eff_bin_eta = new TEfficiency("Eff_bin_eta","my efficiency;x;#epsilon",n_etaBins,x_etaBins);
  //TEfficiency* Eff_bin_2d = new TEfficiency("Eff_bin_2d","my efficiency;x;#epsilon",n_ptBins,x_ptBins,n_etaBins,x_etaBins);
	
	
  //Something interesting you can do when viewing the root file is if you want to investigate the indvidual numerator and denominator
  //TEfficiency actually has these saved so you can just type in the command line
  //TH1* numerator_name = Eff_pt->GetCopyPassedHistogram() or TH1* denominator_name = Eff_pt->GetCopyTotalHistogram()

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

        // check MET Filters
        bool met_filters_ = eve->met_filters;
        // check MET Filters
        bool met_filters_ = eve->met_filters;
        //bool filterbadChCandidate_ = eve->filterbadChCandidate;
        //bool filterbadPFMuon_ = eve->filterbadPFMuon;
        //if(!met_filters_ || !filterbadChCandidate_ || !filterbadPFMuon_)
        //    continue;
        bool badGlobalMuonTagger_ = eve->badGlobalMuonTagger;
        bool cloneGlobalMuonTagger_ = eve->cloneGlobalMuonTagger;
        if(!met_filters_ || !badGlobalMuonTagger_ || !cloneGlobalMuonTagger_)
            continue;


		//Grab Specific Lepton information from trees
		
		int event_nr = eve->evt_;
		vdouble lepton_pt = eve->lepton_pt_;
		vdouble lepton_eta = eve->lepton_eta_;
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
		
		vdouble jet_pt = eve->jet_pt_;
		vdouble jet_eta = eve->jet_eta_;
		vdouble jet_csv = eve->jet_csv_;
		
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
			    //if (isNonTrigMVAM[i]==1 && lepton_rel_Iso[i]<0.15 ) {  			 // for Non-Triggering MVA electron ID
				if (isTrigCutM[i]==1 && lepton_rel_Iso[i]<0.15 ) {  			 // for Cut based electron ID
				//if (isTrigMVAM[i]==1 && lepton_rel_Iso[i]<0.15 ) {               // for Triggering MVA electron ID
					numLooseEle++;
					if ( lepton_pt[i]>30 && fabs(lepton_eta[i])<2.1  ) {
						numTightEle++;
						vdouble vlepton;
						vlepton.push_back(lepton_pt[i]);
						vlepton.push_back(lepton_eta[i]);
						vlepton.push_back(lepton_phi[i]);
						vlepton.push_back(lepton_energy[i]);
						vvLEPTON.push_back(vlepton);
					}
				}
			}	
		}
		
		for(int i=0;i<int(jet_pt.size());i++) {
			if ( jet_pt[i]>30 && fabs(jet_eta[i])<2.4 ) {
				numJets++;
				HT = HT + jet_pt[i];
				if ( jet_csv[i] >= 0.8484 )
					numTags++;
			}
		}
		
		//Selection
		
		//Generally speaking for ttH in the lepton+jets H->bb channel we require exactly 1 Tight Lepton + 0 additional Loose Leptons
		//Since we are only going to look at Electron triggers we require exactly 1 Tight Electron and 0 additional Loose Leptons 
		
		if(numMu!=0)continue;
		if(numTightEle!=1)continue;
		if(numLooseEle!=1)continue;
		
		//Additionally we require 4Jets and 2Btag Jets
			
		if( numJets<4 || numTags<2 ) continue;	
			
		int this_category = -1;
  	    if( numJets==4 && numTags==2) this_category=0;
  	  	if( numJets==5 && numTags==2) this_category=1;
  	    if( numJets>=6 && numTags==2) this_category=2;	
  	    if( numJets==4 && numTags==3) this_category=3;
  	    if( numJets==5 && numTags==3) this_category=4;
  	    if( numJets>=6 && numTags==3) this_category=5;
  	    if( numJets==4 && numTags>=4) this_category=6;
  	    if( numJets==5 && numTags>=4) this_category=7;
  	    if( numJets>=6 && numTags>=4) this_category=8;
		
		//if(this_category==-1)continue;
		
		N_eve++;
		
		// Checking control trigger selection
		
		if (eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_ != 1) continue;		
		
		h_catyield->Fill(this_category);	
		
		int pass_WPLoose_27_HT200 = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;
		int pass_WPTight_27 = eve->pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_;
		int pass_WPTight_32 = eve->pass_HLT_Ele32_eta2p1_WPTight_Gsf_v_;
		int pass_WPTight_27_OR_WPLoose_27_HT200 = 0;
		if (pass_WPLoose_27_HT200 ==1 || pass_WPTight_27 ==1 ) 
			pass_WPTight_27_OR_WPLoose_27_HT200 = 1;
		
		
		if(pass_WPLoose_27_HT200==1)h_catyield_CrossEleTrig->Fill(this_category);
		if(pass_WPTight_27==1)h_catyield_SingleEleTrig->Fill(this_category);
		
		N_eve_control++;
		if(pass_WPTight_27==1) N_eve_tight_27++;
		if(pass_WPTight_32==1) N_eve_tight_32++;
		if(pass_WPLoose_27_HT200==1) N_eve_loose_27_ht200++;
		if(pass_WPTight_27_OR_WPLoose_27_HT200==1) N_eve_tight_27_OR_loose_27_ht200++;
	
		
		//Using TEfficiency
		
		Eff_pt_WPTight_27->Fill(pass_WPTight_27,vvLEPTON[0][0]);
		Eff_eta_WPTight_27->Fill(pass_WPTight_27,vvLEPTON[0][1]);
		Eff_phi_WPTight_27->Fill(pass_WPTight_27,vvLEPTON[0][2]);
		Eff_2d_WPTight_27->Fill(pass_WPTight_27,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPTight_27->Fill(pass_WPTight_27,HT);
		Eff_numPV_WPTight_27->Fill(pass_WPTight_27,numPV);
		
		Eff_pt_WPTight_32->Fill(pass_WPTight_32,vvLEPTON[0][0]);
		Eff_eta_WPTight_32->Fill(pass_WPTight_32,vvLEPTON[0][1]);
		Eff_phi_WPTight_32->Fill(pass_WPTight_32,vvLEPTON[0][2]);
		Eff_2d_WPTight_32->Fill(pass_WPTight_32,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPTight_32->Fill(pass_WPTight_32,HT);
		Eff_numPV_WPTight_32->Fill(pass_WPTight_32,numPV);
		
		Eff_pt_WPLoose_27_HT200->Fill(pass_WPLoose_27_HT200,vvLEPTON[0][0]);
		Eff_eta_WPLoose_27_HT200->Fill(pass_WPLoose_27_HT200,vvLEPTON[0][1]);
		Eff_phi_WPLoose_27_HT200->Fill(pass_WPLoose_27_HT200,vvLEPTON[0][2]);
		Eff_2d_WPLoose_27_HT200->Fill(pass_WPLoose_27_HT200,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPLoose_27_HT200->Fill(pass_WPLoose_27_HT200,HT);
		Eff_numPV_WPLoose_27_HT200->Fill(pass_WPLoose_27_HT200,numPV);
		
		Eff_pt_WPTight_27_OR_WPLoose_27_HT200->Fill(pass_WPTight_27_OR_WPLoose_27_HT200,vvLEPTON[0][0]);
		Eff_eta_WPTight_27_OR_WPLoose_27_HT200->Fill(pass_WPTight_27_OR_WPLoose_27_HT200,vvLEPTON[0][1]);
		Eff_phi_WPTight_27_OR_WPLoose_27_HT200->Fill(pass_WPTight_27_OR_WPLoose_27_HT200,vvLEPTON[0][2]);
		Eff_2d_WPTight_27_OR_WPLoose_27_HT200->Fill(pass_WPTight_27_OR_WPLoose_27_HT200,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPTight_27_OR_WPLoose_27_HT200->Fill(pass_WPTight_27_OR_WPLoose_27_HT200,HT);
		Eff_numPV_WPTight_27_OR_WPLoose_27_HT200->Fill(pass_WPTight_27_OR_WPLoose_27_HT200,numPV);
		
		//Eff_bin_pt->Fill(passsing,vvLEPTON[0][0]);
		//Eff_bin_eta->Fill(passsing,vvLEPTON[0][1]);
		//Eff_bin_2d->Fill(passsing,vvLEPTON[0][0],vvLEPTON[0][1]);
				
   }
   
  t=clock()-t;
	
  std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
  std::cout<<"**********************************************************************************************\n";
  std::cout<<"Total No. of events : "<<N_total<<"\n";
  std::cout<<"No. of events passing event selection only : "<<N_eve<<"\n";
  std::cout<<"No. of events passing event selection plus Control Trigger  : "<<N_eve_control<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele27 Trigger : "<<N_eve_tight_27<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele32 Trigger : "<<N_eve_tight_32<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Loose_Ele27_HT200 Trigger : "<<N_eve_loose_27_ht200<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele27 OR WP_Loose_Ele27_HT200 Trigger : "<<N_eve_tight_27_OR_loose_27_ht200<<"\n";
  std::cout<<"**********************************************************************************************\n";
  histofile.Write();
  histofile.Close();
  
}


