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
  int N_eve_control = 0;
  int N_eve_loose_27 = 0;
  int N_eve_tight_27 = 0;
  int N_eve_tight_32 = 0;
  int N_eve_loose_27_ht200 = 0;
  int N_eve_tight_27_OR_loose_27_ht200 = 0;

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

  std::string histofilename = "Control_Region_histos" + s_end;
 
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

  //Histograms
	 
  TH1::SetDefaultSumw2();

  TH1D* pt_WPLoose_27 = new TH1D("Pt_WPLoose_27","Tag Electron Pt Distribution for WPLoose_27;pT (GeV);Nr. of Events",150,0,300);
  TH1D* pt2_WPLoose_27 = new TH1D("Pt2_WPLoose_27","Probe Electron Pt Distribution for WPLoose_27;pT (GeV);Nr. of Events",150,0,300);
  TH1D* eta_WPLoose_27 = new TH1D("Eta_WPLoose_27","Eta Distribution for WPLoose_27;#eta;Nr. of Events",30,-3,3);
  TH1D* phi_WPLoose_27 = new TH1D("Phi_WPLoose_27","Phi Distribution for WPLoose_27;#phi;Nr. of Events",30,-3,3);
  TH1D* HT_WPLoose_27 = new TH1D("HT_WPLoose_27","HT Distribution for WPLoose_27;HT (GeV);Nr. of Events",200,0,1000);	
  TH1D* numPV_WPLoose_27 = new TH1D("NumPV_WPLoose_27","NumPV Distribution for WPLoose_27;numPV;Nr. of Events",50,0,50);		
  TH1D* jet1_pt_WPLoose_27 = new TH1D("Jet1_pt_WPLoose_27","Jet1 pT Distribution for WPLoose_27;pT (GeV);Nr. of Events",250,0,500);
  TH1D* jet1_csv_WPLoose_27 = new TH1D("Jet1_csv_WPLoose_27","Jet1 csv Distribution for WPLoose_27;csv;Nr. of Events",50,0,1);
  TH1D* njets_WPLoose_27 = new TH1D("Njets_WPLoose_27","N_jets Distribution for WPLoose_27;Nr. of jets;Nr. of Events",10,0,10);
  TH1D* nbtags_WPLoose_27 = new TH1D("Nbtags_WPLoose_27","N_btags Distribution for WPLoose_27;Nr. of btags;Nr. of Events",5,0,5);	
	
  TH1D* pt_WPTight_27 = new TH1D("Pt_WPTight_27","Tag Electron Pt Distribution for WPTight_27;pT (GeV);Nr. of Events",150,0,300);
  TH1D* pt2_WPTight_27 = new TH1D("Pt2_WPTight_27","Probe Electron Pt Distribution for WPTight_27;pT (GeV);Nr. of Events",150,0,300);
  TH1D* eta_WPTight_27 = new TH1D("Eta_WPTight_27","Eta Distribution for WPTight_27;#eta;Nr. of Events",30,-3,3);
  TH1D* phi_WPTight_27 = new TH1D("Phi_WPTight_27","Phi Distribution for WPTight_27;#phi;Nr. of Events",30,-3,3);
  TH1D* HT_WPTight_27 = new TH1D("HT_WPTight_27","HT Distribution for WPTight_27;HT (GeV);Nr. of Events",200,0,1000);	
  TH1D* numPV_WPTight_27 = new TH1D("NumPV_WPTight_27","NumPV Distribution for WPTight_27;numPV;Nr. of Events",50,0,50);		
  TH1D* jet1_pt_WPTight_27 = new TH1D("Jet1_pt_WPTight_27","Jet1 pT Distribution for WPTight_27;pT (GeV);Nr. of Events",250,0,500);
  TH1D* jet1_csv_WPTight_27 = new TH1D("Jet1_csv_WPTight_27","Jet1 csv Distribution for WPTight_27;csv;Nr. of Events",50,0,1);
  TH1D* njets_WPTight_27 = new TH1D("Njets_WPTight_27","N_jets Distribution for WPTight_27;Nr. of jets;Nr. of Events",10,0,10);
  TH1D* nbtags_WPTight_27 = new TH1D("Nbtags_WPTight_27","N_btags Distribution for WPTight_27;Nr. of btags;Nr. of Events",5,0,5);	
	
  TH1D* pt_WPTight_32 = new TH1D("Pt_WPTight_32","Tag Electron Pt Distribution for WPTight_32;pT (GeV);Nr. of Events",150,0,300);
  TH1D* pt2_WPTight_32 = new TH1D("Pt2_WPTight_32","Probe Electron Pt Distribution for WPTight_32;pT (GeV);Nr. of Events",150,0,300);
  TH1D* eta_WPTight_32 = new TH1D("Eta_WPTight_32","Eta Distribution for WPTight_32;#eta;Nr. of Events",30,-3,3);
  TH1D* phi_WPTight_32 = new TH1D("Phi_WPTight_32","Phi Distribution for WPTight_32;#phi;Nr. of Events",30,-3,3);
  TH1D* HT_WPTight_32 = new TH1D("HT_WPTight_32","HT Distribution for WPTight_32;HT (GeV);Nr. of Events",200,0,1000);	
  TH1D* numPV_WPTight_32 = new TH1D("NumPV_WPTight_32","NumPV Distribution for WPTight_32;numPV;Nr. of Events",50,0,50);
  TH1D* jet1_pt_WPTight_32 = new TH1D("Jet1_pt_WPTight_32","Jet1 pT Distribution for WPTight_32;pT (GeV);Nr. of Events",250,0,500);
  TH1D* jet1_csv_WPTight_32 = new TH1D("Jet1_csv_WPTight_32","Jet1 csv Distribution for WPTight_32;csv;Nr. of Events",50,0,1);	
  TH1D* njets_WPTight_32 = new TH1D("Njets_WPTight_32","N_jets Distribution for WPTight_32;Nr. of jets;Nr. of Events",10,0,10);
  TH1D* nbtags_WPTight_32 = new TH1D("Nbtags_WPTight_32","N_btags Distribution for WPTight_32;Nr. of btags;Nr. of Events",5,0,5);	
  
  TH1D* pt_WPLoose_27_HT200 = new TH1D("Pt_WPLoose_27_HT200","Tag Electron Pt Distribution for WPLoose_27_HT200;pT (GeV);Nr. of Events",150,0,300);
  TH1D* pt2_WPLoose_27_HT200 = new TH1D("Pt2_WPLoose_27_HT200","Probe Electron Pt Distribution for WPLoose_27_HT200;pT (GeV);Nr. of Events",150,0,300);
  TH1D* eta_WPLoose_27_HT200 = new TH1D("Eta_WPLoose_27_HT200","Eta Distribution for WPLoose_27_HT200;#eta;Nr. of Events",30,-3,3);
  TH1D* phi_WPLoose_27_HT200 = new TH1D("Phi_WPLoose_27_HT200","Phi Distribution for WPLoose_27_HT200;#phi;Nr. of Events",30,-3,3);
  TH1D* HT_WPLoose_27_HT200 = new TH1D("HT_WPLoose_27_HT200","HT Distribution for WPLoose_27_HT200;HT (GeV);Nr. of Events",200,0,1000);	
  TH1D* numPV_WPLoose_27_HT200 = new TH1D("NumPV_WPLoose_27_HT200","NumPV Distribution for WPLoose_27_HT200;numPV;Nr. of Events",50,0,50);
  TH1D* jet1_pt_WPLoose_27_HT200 = new TH1D("Jet1_pt_WPLoose_27_HT200","Jet1 pT Distribution for WPLoose_27_HT200;pT (GeV);Nr. of Events",250,0,500);
  TH1D* jet1_csv_WPLoose_27_HT200 = new TH1D("Jet1_csv_WPLoose_27_HT200","Jet1 csv Distribution for WPLoose_27_HT200;csv;Nr. of Events",50,0,1);		
  TH1D* njets_WPLoose_27_HT200 = new TH1D("Njets_WPLoose_27_HT200","N_jets Distribution for WPLoose_27_HT200;Nr. of jets;Nr. of Events",10,0,10);
  TH1D* nbtags_WPLoose_27_HT200 = new TH1D("Nbtags_WPLoose_27_HT200","N_btags Distribution for WPLoose_27_HT200;Nr. of btags;Nr. of Events",5,0,5);	
  
  TH1D* pt_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Pt_WPTight_27_OR_WPLoose_27_HT200","Tag Electron Pt Distribution for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);Nr. of Events",150,0,300);
  TH1D* pt2_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Pt2_WPTight_27_OR_WPLoose_27_HT200","Probe Electron Pt Distribution for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);Nr. of Events",150,0,300);
  TH1D* eta_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Eta_WPTight_27_OR_WPLoose_27_HT200","Eta Distribution for WPTight_27_OR_WPLoose_27_HT200;#eta;Nr. of Events",30,-3,3);
  TH1D* phi_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Phi_WPTight_27_OR_WPLoose_27_HT200","Phi Distribution for WPTight_27_OR_WPLoose_27_HT200;#phi;Nr. of Events",30,-3,3);
  TH1D* HT_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("HT_WPTight_27_OR_WPLoose_27_HT200","HT Distribution for WPTight_27_OR_WPLoose_27_HT200;HT (GeV);Nr. of Events",200,0,1000);	
  TH1D* numPV_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("NumPV_WPTight_27_OR_WPLoose_27_HT200","NumPV Distribution for WPTight_27_OR_WPLoose_27_HT200;numPV;Nr. of Events",50,0,50);
  TH1D* jet1_pt_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200","Jet1 pT Distribution for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);Nr. of Events",250,0,500);
  TH1D* jet1_csv_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200","Jet1 csv Distribution for WPTight_27_OR_WPLoose_27_HT200;csv;Nr. of Events",50,0,1);			
  TH1D* njets_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Njets_WPTight_27_OR_WPLoose_27_HT200","N_jets Distribution for WPTight_27_OR_WPLoose_27_HT200;Nr. of jets;Nr. of Events",10,0,10);
  TH1D* nbtags_WPTight_27_OR_WPLoose_27_HT200 = new TH1D("Nbtags_WPTight_27_OR_WPLoose_27_HT200","N_btags Distribution for WPTight_27_OR_WPLoose_27_HT200;Nr. of btags;Nr. of Events",5,0,5);

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
		
		//int opp_charge = eve->oppositeLepCharge_;
		//double mass_leplep = eve->mass_leplep_;
	        int opp_charge= 0;	
	        double mass_leplep = 0;     
	   
		vdouble jet_pt = eve->jet_pt_;
		vdouble jet_eta = eve->jet_eta_;
		vdouble jet_csv = eve->jet_csv_;
		
		double gen_weight = eve->gen_weight_;
  		double csv_weight = eve->csv_weight_;
  		double PU_weight = eve->PU_weight_ ;
  		double PDF_weight = eve->PDF_weight_;
  		double Q2_weight = eve->Q2_weight_;
  		vdouble lepton_id_sf = eve->lepton_id_sf_;
  		vdouble lepton_iso_sf = eve->lepton_iso_sf_;
  		vdouble lepton_gsf_sf = eve->lepton_gsf_sf_;
  		vdouble lepton_trig_sf = eve->lepton_trig_sf_;
  		vdouble lepton_hip_sf = eve->lepton_hip_sf_;
  		
  		double tot_weight = 1;

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
		if(numTightEle<1)continue;
		
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
		
		//Additionally we require 4Jets and 2Btag Jets
			
		//if( numJets<2 ) continue;
		//if( numJets<4 || numTags<2 ) continue;	
		
		N_eve++;
		
		// Checking control trigger selection
		
		if (eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_ == 1) N_eve_control++;		
		
		int pass_WPLoose_27 = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_;
		int pass_WPLoose_27_HT200 = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;
		int pass_WPTight_27 = eve->pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_;
		int pass_WPTight_32 = eve->pass_HLT_Ele32_eta2p1_WPTight_Gsf_v_;
		int pass_WPTight_27_OR_WPLoose_27_HT200 = 0;
		if (pass_WPLoose_27_HT200 ==1 || pass_WPTight_27 ==1 ) 
			pass_WPTight_27_OR_WPLoose_27_HT200 = 1;


        //vvLEPTON[0][8] = 1; // for Data
        //vvLEPTON[1][8] = 1; // for Data

		// PDF Weight, Q2 Weight not applied
	    tot_weight = vvLEPTON[0][8]*vvLEPTON[1][8]*csv_weight*PU_weight*gen_weight;
	
		// Fill Histograms
		
		int lead, sub_lead;
		lead = sub_lead = 0;
		if(vvLEPTON[0][0] >= vvLEPTON[1][0])
			sub_lead = 1;
		else 
			lead = 1;
		
		
		if(pass_WPLoose_27==1) {
		
			pt_WPLoose_27->Fill(vvLEPTON[lead][0], tot_weight );
			pt2_WPLoose_27->Fill(vvLEPTON[sub_lead][0], tot_weight );
			eta_WPLoose_27->Fill(vvLEPTON[lead][1], tot_weight );
			phi_WPLoose_27->Fill(vvLEPTON[lead][2], tot_weight );
			HT_WPLoose_27->Fill(HT, tot_weight );
			numPV_WPLoose_27->Fill(numPV, tot_weight );
			jet1_pt_WPLoose_27->Fill(jet1_pt, tot_weight );
			jet1_csv_WPLoose_27->Fill(jet1_csv, tot_weight );
			njets_WPLoose_27->Fill(numJets, tot_weight );
			nbtags_WPLoose_27->Fill(numTags, tot_weight );
			
			N_eve_loose_27++;
		}
			
		if(pass_WPTight_27==1) {
		
			pt_WPTight_27->Fill(vvLEPTON[lead][0], tot_weight );
			pt2_WPTight_27->Fill(vvLEPTON[sub_lead][0], tot_weight );
			eta_WPTight_27->Fill(vvLEPTON[lead][1], tot_weight );
			phi_WPTight_27->Fill(vvLEPTON[lead][2], tot_weight );
			HT_WPTight_27->Fill(HT, tot_weight );
			numPV_WPTight_27->Fill(numPV, tot_weight );
			jet1_pt_WPTight_27->Fill(jet1_pt, tot_weight );
			jet1_csv_WPTight_27->Fill(jet1_csv, tot_weight );
			njets_WPTight_27->Fill(numJets, tot_weight );
			nbtags_WPTight_27->Fill(numTags, tot_weight );
			
			N_eve_tight_27++;
		}
		
		if(pass_WPTight_32==1) {
		
			pt_WPTight_32->Fill(vvLEPTON[lead][0], tot_weight );
			pt2_WPTight_32->Fill(vvLEPTON[sub_lead][0], tot_weight );
			eta_WPTight_32->Fill(vvLEPTON[lead][1], tot_weight );
			phi_WPTight_32->Fill(vvLEPTON[lead][2], tot_weight );
			HT_WPTight_32->Fill(HT, tot_weight );
			numPV_WPTight_32->Fill(numPV, tot_weight );
			jet1_pt_WPTight_32->Fill(jet1_pt, tot_weight );
			jet1_csv_WPTight_32->Fill(jet1_csv, tot_weight );
			njets_WPTight_32->Fill(numJets, tot_weight );
			nbtags_WPTight_32->Fill(numTags, tot_weight );
			
			N_eve_tight_32++;
		}
		
		if(pass_WPLoose_27_HT200==1) {
		
			pt_WPLoose_27_HT200->Fill(vvLEPTON[lead][0], tot_weight );
			pt2_WPLoose_27_HT200->Fill(vvLEPTON[sub_lead][0], tot_weight );
			eta_WPLoose_27_HT200->Fill(vvLEPTON[lead][1], tot_weight );
			phi_WPLoose_27_HT200->Fill(vvLEPTON[lead][2], tot_weight );
			HT_WPLoose_27_HT200->Fill(HT, tot_weight );
			numPV_WPLoose_27_HT200->Fill(numPV, tot_weight );
			jet1_pt_WPLoose_27_HT200->Fill(jet1_pt, tot_weight );
			jet1_csv_WPLoose_27_HT200->Fill(jet1_csv, tot_weight );
			njets_WPLoose_27_HT200->Fill(numJets, tot_weight );
			nbtags_WPLoose_27_HT200->Fill(numTags, tot_weight );
			
			N_eve_loose_27_ht200++;
		}
		
		if(pass_WPTight_27_OR_WPLoose_27_HT200==1) {
		
			pt_WPTight_27_OR_WPLoose_27_HT200->Fill(vvLEPTON[lead][0], tot_weight );
			pt2_WPTight_27_OR_WPLoose_27_HT200->Fill(vvLEPTON[sub_lead][0], tot_weight );
			eta_WPTight_27_OR_WPLoose_27_HT200->Fill(vvLEPTON[lead][1], tot_weight );
			phi_WPTight_27_OR_WPLoose_27_HT200->Fill(vvLEPTON[lead][2], tot_weight );
			HT_WPTight_27_OR_WPLoose_27_HT200->Fill(HT, tot_weight );
			numPV_WPTight_27_OR_WPLoose_27_HT200->Fill(numPV, tot_weight );
			jet1_pt_WPTight_27_OR_WPLoose_27_HT200->Fill(jet1_pt, tot_weight );
			jet1_csv_WPTight_27_OR_WPLoose_27_HT200->Fill(jet1_csv, tot_weight );
			njets_WPTight_27_OR_WPLoose_27_HT200->Fill(numJets, tot_weight );
			nbtags_WPTight_27_OR_WPLoose_27_HT200->Fill(numTags, tot_weight );
			
			N_eve_tight_27_OR_loose_27_ht200++;
		}
				
   }
   
  t=clock()-t;
	
  std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
  std::cout<<"**********************************************************************************************\n";
  std::cout<<"Total No. of events : "<<N_total<<"\n";
  std::cout<<"No. of events passing event selection only : "<<N_eve<<"\n";
  std::cout<<"No. of events passing event selection plus Control Trigger  : "<<N_eve_control<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Loose_Ele27 Trigger : "<<N_eve_loose_27<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele27 Trigger : "<<N_eve_tight_27<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele32 Trigger : "<<N_eve_tight_32<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Loose_Ele27_HT200 Trigger : "<<N_eve_loose_27_ht200<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight_Ele27 OR WP_Loose_Ele27_HT200 Trigger : "<<N_eve_tight_27_OR_loose_27_ht200<<"\n";
  std::cout<<"Sum of Generator Weights for sample : "<<sum_gen_weight<<"\n";
  std::cout<<"**********************************************************************************************\n";
  histofile.Write();
  histofile.Close();
  
}


