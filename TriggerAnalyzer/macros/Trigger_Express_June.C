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


#include "TriggerRun2/TriggerAnalyzer/interface/TriggerStudyEventVars.h"

#endif




typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;
typedef std::vector<std::vector<int> >         vvint;




void Trigger_Express_June( int maxNentries=-1, int Njobs=1, int jobN=1 ) {



clock_t t;
t = clock();


  
  std::cout << "   ===> load the root files! " << std::endl;
  
  double wgt=1;
  

 
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string treefilename = "/eos/uscms/store/user/sflowers/treeMaker/Summer16_June13_Trigger/SingleElectron/Run2015D-16Run2016B-v1_SilverJSON_2016_01_31_triggerTree_v1_211/160615_131259/0000/trigger_analyzer_*.root";
  std::string treefilename2 = "/eos/uscms/store/user/sflowers/treeMaker/Summer16_June13_Trigger/SingleElectron/Run2015D-16Run2016B-v1_SilverJSON_2016_01_31_triggerTree_v1_211/160615_131259/0001/trigger_analyzer_*.root";
  
 
  std::string s_end = "_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = ".root";

   std::string histofilename = "Express_Trigger_June_Data" + s_end;
 
  
  std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

 TChain *chain = new TChain("triggeranalyzer/triggerTree");
 
  chain->Add(treefilename.c_str());
  chain->Add(treefilename2.c_str());
 
  
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
	
	TEfficiency* Eff_pt = new TEfficiency("Eff_pt","my efficiency;x;#epsilon",100,0,200);
	TEfficiency* Eff_eta = new TEfficiency("Eff_eta","my efficiency;x;#epsilon",25,-2.5,2.5);
	TEfficiency* Eff_2d = new TEfficiency("Eff_2d","my efficiency;x;#epsilon",100,0,200,25,-2.5,2.5);
	
	//Sometimes you are asked to do specific binning for your efficiency plots (in the case of making Trigger Scale Factors)
	
	int n_ptBins = 5;
	double x_ptBins[6] = {15.0, 25.0, 35.0, 45.0, 55.0, 200.0};
	
	int n_etaBins = 10;
	double x_etaBins[11] = {-2.5,-2.1, -1.5660, -1.4442, -0.8, 0, 0.8, 1.4442, 1.5660, 2.1, 2.5};
	
	TEfficiency* Eff_bin_pt = new TEfficiency("Eff_bin_pt","my efficiency;x;#epsilon",n_ptBins,x_ptBins);
	TEfficiency* Eff_bin_eta = new TEfficiency("Eff_bin_eta","my efficiency;x;#epsilon",n_etaBins,x_etaBins);
	TEfficiency* Eff_bin_2d = new TEfficiency("Eff_bin_2d","my efficiency;x;#epsilon",n_ptBins,x_ptBins,n_etaBins,x_etaBins);
	
	
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
		
		vdouble lepton_pt = eve->lepton_pt_;
		vdouble lepton_eta = eve->lepton_eta_;
		vdouble lepton_phi = eve->lepton_phi_;
		vdouble lepton_energy = eve->lepton_energy_;
		vint isTight = eve->lepton_isTight_;
		vint isLoose = eve->lepton_isLoose_;
		vint isMuon = eve->lepton_isMuon_; 
		int numJets = eve->numJets_;
		int numTags = eve->numTags_;
		
		//Organize the lepton information into a LorentzVector and Count Tight and Loose Leptons for selection
		
		vvdouble vvLEPTON;
		int numLooseEle=0;
		int numLooseMu=0;
		int numTightEle=0;
		int numTightMu=0;
		for(int i=0;i<int(lepton_pt.size());i++){
			vdouble vlepton;
			vlepton.push_back(lepton_pt[i]);
			vlepton.push_back(lepton_eta[i]);
			vlepton.push_back(lepton_phi[i]);
			vlepton.push_back(lepton_energy[i]);
			vvLEPTON.push_back(vlepton);
			if(isTight[i]==1 && isMuon[i]==1)numTightMu++;
			if(isTight[i]==1 && isMuon[i]==0)numTightEle++;
			if(isTight[i]==0 && isMuon[i]==1)numLooseMu++;
			if(isTight[i]==0 && isMuon[i]==0)numLooseEle++;
		}
		
		//Selection
		
		//Generally speaking for ttH in the lepton+jets H->bb channel we require exactly 1 Tight Lepton + 0 additional Loose Leptons
		//Since we are only going to look at Electron triggers we require exactly 1 Tight Electron and 0 additional Loose Leptons 
		
		if(numLooseMu!=0)continue;
		if(numTightEle!=1)continue;
		if(numLooseEle!=1)continue;
		
		//Additionally we require 4Jets and 2Btag Jets
		
			
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
		
		if(this_category==-1)continue;
		
		h_catyield->Fill(this_category);	
		
		int passcross = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;
		int passsing = eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_;
		
		if(passcross==1)h_catyield_CrossEleTrig->Fill(this_category);
		if(passsing==1)h_catyield_SingleEleTrig->Fill(this_category);
		
		//Using TEfficiency
		
		//Eff_var->Fill(bool_pass,variable);
		//Eff_2d->Fill(bool_pass,variable1,variable2);
		
		Eff_pt->Fill(passcross,vvLEPTON[0][0]);
		Eff_eta->Fill(passcross,vvLEPTON[0][1]);
		Eff_2d->Fill(passcross,vvLEPTON[0][0],vvLEPTON[0][1]);
		
		Eff_bin_pt->Fill(passcross,vvLEPTON[0][0]);
		Eff_bin_eta->Fill(passcross,vvLEPTON[0][1]);
		Eff_bin_2d->Fill(passcross,vvLEPTON[0][0],vvLEPTON[0][1]);
		
		
		
		
		
		
   }
   
t=clock()-t;
	
	std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;

    histofile.Write();
    histofile.Close();

	

	 

}


