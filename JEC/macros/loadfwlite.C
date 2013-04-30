void loadfwlite(){


// Set up FW Lite for automatic loading of CMS libraries
// and data formats.   As you may have other user-defined setup
// in your rootlogon.C, the CMS setup is executed only if the CMS
// environment is set up.
//


  TString cmsswbase = getenv("CMSSW_BASE");
  if (cmsswbase.Length() > 0) {
    //
    // The CMSSW environment is defined (this is true even for FW Lite)
    // so set up the rest.
    //
    cout << "Loading FW Lite setup." << endl;

    gSystem->Load("libFWCoreFWLite");
    gSystem->Load("libDataFormatsFWLite");
    gSystem->Load("libDataFormatsCommon");
    gSystem->Load("libDataFormatsCaloTowers");
    gSystem->Load("libDataFormatsHeavyIonEvent");
    gSystem->Load("libDataFormatsPatCandidates");
    gSystem->Load("libDataFormatsRecoCandidate");
    
    gSystem->Load("libSimDataFormatsHiGenData");
    gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
    gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
    AutoLibraryLoader::enable();

    cout<<"Loading FW lite done : "<<endl;

}





}
