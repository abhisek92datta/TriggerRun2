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

/*
bool passesMuonPOGIdTight(const double normalizedChi2, const int numberOfValidMuonHits, const double d0, const double dz, const int numberOfValidPixelHits , const int trackerLayersWithMeasurement, const int numberOfMatchedStations, const int isPFMuon ){

    bool passesGlobalTrackID = ( ( normalizedChi2 < 10.) 
				 && (numberOfValidMuonHits > 0)
				 );
    if(!passesGlobalTrackID) return false;

    bool passesMuonBestTrackID = ( ( fabs(d0) < 0.2)
				   && (fabs(dz) < 0.5)
				   );
    if(!passesMuonBestTrackID) return false;

    bool passesInnerTrackID = ( numberOfValidPixelHits > 0);
    if(!passesInnerTrackID) return false;

    bool passesTrackID = (trackerLayersWithMeasurement > 5);
    if(!passesTrackID) return false;

    if(numberOfMatchedStations <= 1) return false;

    if(!isPFMuon) return false;

    return true;

}
*/

void Trigger_Express_June( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

  clock_t t;
  t = clock();
  
  int N_total = 0;
  int N_eve = 0;
  int N_eve_control = 0;
  int N_eve_loose = 0;
  int N_eve_tight = 0;
  int N_eve_loose_ht200 = 0;
  
  std::cout << "   ===> load the root files! " << std::endl;
  
  double wgt=1;
  
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string treefilename = "/eos/uscms/store/user/adatta/Trigger_Analysis/ttHTobb_M125_13TeV_powheg_pythia8/triggerTree_MC_80X_tth/160906_030939/0000/trigger_analyzer_*.root";
  //std::string treefilename2 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0001/trigger_analyzer_*.root";
  //std::string treefilename3 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0002/trigger_analyzer_*.root";
  
  std::string s_end = "_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = ".root";

  std::string histofilename = "Express_Trigger_June_Data" + s_end;
 
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
	
  TEfficiency* Eff_pt_WPLoose = new TEfficiency("Eff_pt_WPLoose","Efficiency vs pT for WPLoose;pT (GeV);(WPLoose + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",150,0,300);
  TEfficiency* Eff_eta_WPLoose = new TEfficiency("Eff_eta_WPLoose","Efficiency vs eta for WPLoose;#eta;(WPLoose + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPLoose = new TEfficiency("Eff_phi_WPLoose","Efficiency vs phi for WPLoose;#phi;(WPLoose + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPLoose = new TEfficiency("Eff_2d_WPLoose","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  TEfficiency* Eff_HT_WPLoose = new TEfficiency("Eff_HT_WPLoose","Efficiency vs HT for WPLoose;HT (GeV);(WPLoose + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",200,0,1000);	
  TEfficiency* Eff_numPV_WPLoose = new TEfficiency("Eff_numPV_WPLoose","Efficiency vs numPV for WPLoose;numPV;(WPLoose + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);		
	
  TEfficiency* Eff_pt_WPTight = new TEfficiency("Eff_pt_WPTight","Efficiency vs pT for WPTight;pT (GeV);(WPTight + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",150,0,300);
  TEfficiency* Eff_eta_WPTight = new TEfficiency("Eff_eta_WPTight","Efficiency vs eta for WPTight;#eta;(WPTight + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPTight = new TEfficiency("Eff_phi_WPTight","Efficiency vs phi for WPTight;#phi;(WPTight + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPTight = new TEfficiency("Eff_2d_WPTight","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  TEfficiency* Eff_HT_WPTight = new TEfficiency("Eff_HT_WPTight","Efficiency vs HT for WPTight;HT (GeV);(WPTight + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",200,0,1000);	
  TEfficiency* Eff_numPV_WPTight = new TEfficiency("Eff_numPV_WPTight","Efficiency vs numPV for WPTight;numPV;(WPTight + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);	
	
	
  TEfficiency* Eff_pt_WPLoose_HT200 = new TEfficiency("Eff_pt_WPLoose_HT200","Efficiency vs pT for WPLoose_HT200;pT (GeV);(WPLoose_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",150,0,300);
  TEfficiency* Eff_eta_WPLoose_HT200 = new TEfficiency("Eff_eta_WPLoose_HT200","Efficiency vs eta for WPLoose_HT200;#eta;(WPLoose_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_phi_WPLoose_HT200 = new TEfficiency("Eff_phi_WPLoose_HT200","Efficiency vs phi for WPLoose_HT200;#phi;(WPLoose_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",30,-3,3);
  TEfficiency* Eff_2d_WPLoose_HT200 = new TEfficiency("Eff_2d_WPLoose_HT200","my efficiency;x;#epsilon",150,0,300,30,-3,3);
  TEfficiency* Eff_HT_WPLoose_HT200 = new TEfficiency("Eff_HT_WPLoose_HT200","Efficiency vs HT for WPLoose_HT200;HT (GeV);(WPLoose_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",200,0,1000);	
  TEfficiency* Eff_numPV_WPLoose_HT200 = new TEfficiency("Eff_numPV_WPLoose_HT200","Efficiency vs numPV for WPLoose_HT200;numPV;(WPLoose_HT200 + Control Trigger + Event_sel)/(Control Trigger + Event_sel)",50,0,50);		
  	
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
		
		
		//Grab Specific Lepton information from trees
		
		int event_nr = eve->evt_;
		vdouble lepton_pt = eve->lepton_pt_;
		vdouble lepton_eta = eve->lepton_eta_;
		vdouble lepton_phi = eve->lepton_phi_;
		vdouble lepton_energy = eve->lepton_energy_;
		vdouble lepton_rel_Iso = eve->lepton_relIso_;
		vdouble lepton_rel_IsoR04 = eve->lepton_relIsoR04_;
		
		/*
		vdouble lepton_d0 = eve->lepton_d0_;
  		vdouble lepton_dZ = eve->lepton_dZ_;
  		vint lepton_isPFMuon = eve->lepton_isPFMuon_;
  		vdouble lepton_normalizedChi2 = eve->lepton_normalizedChi2_;
  		vint lepton_numberOfValidMuonHits = eve->lepton_numberOfValidMuonHits_;
  		vint lepton_numberOfValidPixelHits = eve->lepton_numberOfValidPixelHits_;
  		vint lepton_trackerLayersWithMeasurement = eve->lepton_trackerLayersWithMeasurement_;
		vint lepton_numberOfMatchedStations = eve->lepton_numberOfMatchedStations_;
		*/
		
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
		
		//int numJets = eve->numJets_;
		//int numTags = eve->numTags_;
		int numJets = 0;
		int numTags = 0;
		
		//Organize the lepton information into a LorentzVector and Count Tight and Loose Leptons for selection
		
		vvdouble vvLEPTON;
		int numLooseEle=0;
		int numMu = 0;
		//int numLooseMu=0;
		int numTightEle=0;
		//int numTightMu=0;
		for(int i=0;i<int(lepton_pt.size());i++){
			//vdouble vlepton;
			//vlepton.push_back(lepton_pt[i]);
			//vlepton.push_back(lepton_eta[i]);
			//vlepton.push_back(lepton_phi[i]);
			//vlepton.push_back(lepton_energy[i]);
			//vvLEPTON.push_back(vlepton);
			
		/*	
			if(isMuon[i]==1) {
				if (lepton_pt[i]>15 && fabs(lepton_eta[i])<2.4) {
					if (lepton_rel_IsoR04[i]<0.25) {
						if (passesMuonPOGIdTight(lepton_normalizedChi2[i], lepton_numberOfValidMuonHits[i], lepton_d0[i], lepton_dZ[i], lepton_numberOfValidPixelHits[i] , lepton_trackerLayersWithMeasurement[i], lepton_numberOfMatchedStations[i], lepton_isPFMuon[i] ) == true)
							numLooseMu++;
					}
				}
			}
		*/		
			 if(is_IDTight[i]==1) {
				numMu++;			
			 }
			
			else {
			    if (isNonTrigMVAM[i]==1 && lepton_rel_Iso[i]<0.15 ) {  			 // for Non-Triggering MVA electron ID
				//if (isTrigCutM[i]==1 && lepton_rel_Iso[i]<0.15 ) {  			 // for Cut based electron ID
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
			
			//if(isTight[i]==1 && isMuon[i]==1)numTightMu++;
			//if(isTight[i]==1 && isMuon[i]==0)numTightEle++;
			//if(isTight[i]==0 && isMuon[i]==1)numLooseMu++;
			//if(isTight[i]==0 && isMuon[i]==0)numLooseEle++;		
		}
		
		for(int i=0;i<int(jet_pt.size());i++) {
			if ( jet_pt[i]>30 && fabs(jet_eta[i])<2.4 ) {
				numJets++;
				HT = HT + jet_pt[i];
				if ( jet_csv[i] > 0.8 )
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
		
		//if (eve->pass_AlCa_SingleEle_WPVeryLoose_Gsf_v_ != 1) continue;
		if (eve->pass_HLT_Ele22_eta2p1_WPLoose_Gsf_v_ != 1) continue;
		//if (eve->pass_HLT_Ele23_WPLoose_Gsf_v_ != 1) continue;
  		//if (eve->pass_HLT_Ele24_eta2p1_WPLoose_Gsf_v_ != 1) continue;
  		//if (eve->pass_HLT_Ele25_eta2p1_WPLoose_Gsf_v_ != 1) continue;
  		//if (eve->pass_HLT_Ele25_WPTight_Gsf_v_ != 1) continue;
  		//if (eve->pass_HLT_Ele25_eta2p1_WPTight_Gsf_v_ != 1) continue;
		
		
		h_catyield->Fill(this_category);	
		
		int pass_WPLoose_HT200 = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;
		int pass_WPTight = eve->pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_;
		int pass_WPLoose = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_;
		
		if(pass_WPLoose_HT200==1)h_catyield_CrossEleTrig->Fill(this_category);
		if(pass_WPTight==1)h_catyield_SingleEleTrig->Fill(this_category);
		
		N_eve_control++;
		if(pass_WPLoose==1) N_eve_loose++;
		if(pass_WPTight==1) N_eve_tight++;
		if(pass_WPLoose_HT200==1) N_eve_loose_ht200++;
	
		
		//Using TEfficiency
		
		//Eff_var->Fill(bool_pass,variable);
		//Eff_2d->Fill(bool_pass,variable1,variable2);
		
		Eff_pt_WPTight->Fill(pass_WPTight,vvLEPTON[0][0]);
		Eff_eta_WPTight->Fill(pass_WPTight,vvLEPTON[0][1]);
		Eff_phi_WPTight->Fill(pass_WPTight,vvLEPTON[0][2]);
		Eff_2d_WPTight->Fill(pass_WPTight,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPTight->Fill(pass_WPTight,HT);
		Eff_numPV_WPTight->Fill(pass_WPTight,numPV);
		
		Eff_pt_WPLoose->Fill(pass_WPLoose,vvLEPTON[0][0]);
		Eff_eta_WPLoose->Fill(pass_WPLoose,vvLEPTON[0][1]);
		Eff_phi_WPLoose->Fill(pass_WPLoose,vvLEPTON[0][2]);
		Eff_2d_WPLoose->Fill(pass_WPLoose,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPLoose->Fill(pass_WPLoose,HT);
		Eff_numPV_WPLoose->Fill(pass_WPLoose,numPV);
		
		Eff_pt_WPLoose_HT200->Fill(pass_WPLoose_HT200,vvLEPTON[0][0]);
		Eff_eta_WPLoose_HT200->Fill(pass_WPLoose_HT200,vvLEPTON[0][1]);
		Eff_phi_WPLoose_HT200->Fill(pass_WPLoose_HT200,vvLEPTON[0][2]);
		Eff_2d_WPLoose_HT200->Fill(pass_WPLoose_HT200,vvLEPTON[0][0],vvLEPTON[0][1]);
		Eff_HT_WPLoose_HT200->Fill(pass_WPLoose_HT200,HT);
		Eff_numPV_WPLoose_HT200->Fill(pass_WPLoose_HT200,numPV);
		
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
  std::cout<<"No. of events passing event selection plus WP_Loose Trigger : "<<N_eve_loose<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Tight Trigger : "<<N_eve_tight<<"\n";
  std::cout<<"No. of events passing event selection plus WP_Loose_HT200 Trigger : "<<N_eve_loose_ht200<<"\n";
  std::cout<<"**********************************************************************************************\n";
  histofile.Write();
  histofile.Close();
  
}


