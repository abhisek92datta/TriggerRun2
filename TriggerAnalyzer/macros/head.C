{

  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy1 = makeshared.ReplaceAll("-W ", "");
  gSystem->SetMakeSharedLib(makeshared);
  TString dummy2 = makeshared.ReplaceAll("-Wshadow ", "-std=c++0x ");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();

  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libMiniAODMiniAODHelper.so");
  gSystem->Load("libTriggerRun2TriggerAnalyzer.so");

  TString path = gSystem->GetIncludePath();
  path.Append(" -I$ROOTSYS/include/root -I./include  ");
  gSystem->SetIncludePath(path);


  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  gStyle->SetPalette(1);

}

