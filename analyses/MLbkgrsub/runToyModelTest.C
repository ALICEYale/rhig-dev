/// \file runToyModelTest.C
/// \brief  macro to run a test of toy model task
///
/// 
/// This macros illustrates how to run a jet analysis. It can run locally or on
/// grid (test/full/terminate modes).
/// The script runEMCalJetSampleTask.sh in the same folder allows to easily set
/// the input values
///
/// \author Hannah Bossi <hannah.bossi@cern.ch>, Yale University
/// \date Feb 4th, 2019 modified from runEMCalJetSampleTask.C created by Salvatore Aiola


class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalCorrectionTask;
class AliEmcalJetTask;
class AliAnalysisTaskRho;
class AliAnalysisTaskRhoMass;
class AliAnalysisTaskEmcalJetPerformance;
class AliAnalysisTaskJetExtractor;
class AliAnalysisGrid;
class AliAnalysisAlien;
class AliJetResponseMaker;
class AliAnalysisTaskChargedJetsHadronToy;

// Include AddTask macros for ROOT 6 compatibility
#ifdef __CLING__
// Tell ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "OADB/macros/AddTaskPhysicsSelection.C"
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"
#include "PWGPP/PilotTrain/AddTaskCDBconnect.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskRhoMass.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetPerformance.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskChargedJetsHadronToy.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskJetExtractor.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C"
#endif

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runToyModelTest(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC15o",                               // set the run period
    const char   *cLocalFiles    = "files_LHC15o_AOD.txt",          // set the local list file
    const UInt_t  iNumEvents     = 500,                                     // number of events to be analyzed
    const UInt_t  kComPhysSel       = AliVEvent::kAnyINT |
    AliVEvent::kCentral | AliVEvent::kSemiCentral,                          // physics selection
    const char   *cTaskName      = "EMCalJetAna",                           // sets name of analysis manager
    const Bool_t  bDoChargedJets = kTRUE,
    const Bool_t  bDoFullJets    = kTRUE,
    const char   *obsolete       = "",                                      // Previous handled the ocdb settings, but obsolete due to CDBconnect task
    // 0 = only prepare the analysis manager but do not start the analysis
    // 1 = prepare the analysis manager and start the analysis
    // 2 = launch a grid analysis
    Int_t         iStartAnalysis = 1,
    const UInt_t  iNumFiles      = 100,                                     // number of files analyzed locally
    const char   *cGridMode      = "test"
)
{
  TString sRunPeriod(cRunPeriod);
  sRunPeriod.ToLower();

  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp;

  Bool_t bIsRun2 = kFALSE;
  if (sRunPeriod.Length() == 6 && sRunPeriod.BeginsWith("lhc15")) bIsRun2 = kTRUE;

  if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
    iBeamType = AliAnalysisTaskEmcal::kAA;
  }
  else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f" ||
      sRunPeriod == "lhc16q" || sRunPeriod == "lhc16r" || sRunPeriod == "lhc16s" ||
      sRunPeriod == "lhc16t") {
    iBeamType = AliAnalysisTaskEmcal::kpA;
  }

  Double_t kGhostArea = 0.01;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) kGhostArea = 0.005;

  //AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);

  //Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// Global Train Configuration file
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*Latest changes:  
  31.8.15 TH2 >THn for comparisons
  28.8.15 
  -Updated configuration based on EMC_JE pp MC train to ensure comparability of published JE results and results from HFCJ
  Please note:
  i. minimum jet p_T, jet finder (tracks) has been lowerd from 1 GeV/c -> 0 GeV/c
  ii.  Response maker matching distance increased for
  THN based output from 0.25 > 1.2; is still 0.25 for TH2
  boost::v1_53_0 cgal::v4.4 fastjet::v3.0.6_1.012
  */
  
  const Bool_t   bDoEmcalCorrections  = kTRUE;
  Bool_t fSpecialPtHard =kTRUE;
  Bool_t fGlobalParam_VetoMB= kFALSE;
  Bool_t fParamDoSmearing=kTRUE;
  Bool_t fParamRunCompositionCorrection=kTRUE;
  /*************************************************************************/
  //Jet Response Maker Configuration
  /*************************************************************************/


  Bool_t doZaxis = kFALSE;
  Bool_t doNEFaxis = kFALSE;
  Bool_t kDoWeighting = kFALSE;
  Int_t kHistoType = 0;  //0=TH2, 1=THnSparse
  UInt_t kMatching = 1; //1=geometrical, 2=MClabel
  Int_t kPtHard = -999;
  Double_t kMaxDistance04 = 1.2;
  Double_t kAreaCut = 0.;  
  Double_t kJetLeadingTrackBias = 0;
  Double_t kClusPtCut = 0.30;
  Double_t kTrackPtCut = 0.15;
  Double_t kPartLevPtCut = 0.15;
  //Double_t kGhostArea = 0.01;
  Double_t kMaxTrackPt = 1000;
  Double_t kJetPtCut = 1;

  Int_t kType = 0; // 0 = charged, 1 = neutral, 2 = both
   
  if (kHistoType == 0) {
    kAreaCut = 0.557;
    if (kMatching == 1) {
      kMaxDistance04 = 0.25;

    } 
    else if (kMatching == 2) {
      kMaxDistance04 = 0.99;
    }
  }
  /*************************************************************************/
  //Global Names
  /*************************************************************************/

  char pTString[200];
  if (kClusPtCut == 0)
    sprintf(pTString,"0000");
  else if (kClusPtCut < 1.0)
    sprintf(pTString,"0%3.0f",kTrackPtCut*1000.0);
  else
    sprintf(pTString,"%4.0f",kTrackPtCut*1000.0);

  char ETString[200];
  if (kClusPtCut == 0)
    sprintf(ETString,"0000");
  else if (kTrackPtCut < 1.0)
    sprintf(ETString,"0%3.0f",kClusPtCut*1000.0);
  else
    sprintf(ETString,"%4.0f",kClusPtCut*1000.0);

  char PartLevString[200];
  if (kPartLevPtCut == 0)
    sprintf(PartLevString,"0000");
  else if (kPartLevPtCut < 1.0)
    sprintf(PartLevString,"0%3.0f",kPartLevPtCut*1000.0);
  else
    sprintf(PartLevString,"%4.0f",kPartLevPtCut*1000.0);

  ///TRACKS  PARTICLES
  TString kMCTracksName("MCParticlesSelected");
  TString kTracksName("PicoTracks");
  //////////////////////////////////////////////////////////////////////////////////////////
  TString kJets04ChargedName(Form("Jet_AKTChargedR040_%s_pT%s_pt_scheme",kTracksName.Data(),pTString));

  TString kJetAkTName(Form("Jet_AKTChargedR040_%s_pT%s_pt_scheme",kTracksName.Data(),pTString));
  TString kJetAkTNameMC(Form("Jet_AKTChargedR040_%s_pT%s_pt_scheme",kMCTracksName.Data(),PartLevString));

  TString kJetRhokTName(Form("Jet_KTChargedR040_%s_pT%s_pt_scheme",kTracksName.Data(),pTString));
  TString kJetRhokTNameMC(Form("Jet_KTChargedR040_%s_pT%s_pt_scheme",kMCTracksName.Data(),PartLevString));

  TString kRhoTaskName       ="ExternalRhoTask";
  TString kRhoTaskNameMC   ="ExternalRhoTaskMC";

  //////////////// OTHER EMCAL VARIABLE //////
  TString kMatchingChainStr="16:1:includeNoITS=kTRUE doProp=kFALSE doAttemptProp=kTRUE isMC=kTRUE";
  TString kClusterName="EmcCaloClusters";
  // Variables for Ruediger
  const Double_t kRHJetRadius = 0.4;
  const Bool_t kRHBuiltinEventSelection = kFALSE;
  const Bool_t kRHEmbedParticleLevel   = kFALSE;

  enum eDataType { kAod, kEsd };

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD")) {
    iDataType = kEsd;
  }
  else if (!strcmp(cDataType, "AOD")) {
    iDataType = kAod;
  }
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("datatype = AOD or ESD");
    return 0;
  }

  Printf("%s analysis chosen.", cDataType);

  TString sLocalFiles(cLocalFiles);
  if (iStartAnalysis == 1) {
    if (sLocalFiles == "") {
      Printf("You need to provide the list of local files!");
      return 0;
    }
    Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);
  }

  #ifndef __CLING__
  LoadMacros();
  #endif

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);
  pMgr->SetDebugLevel(3);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AliAnalysisTaskEmcal::AddAODHandler();
  }
  else {
    AliESDInputHandler* pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
  }

  // Physics selection task
  if (iDataType == kEsd) {
    AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();
  }

  // Centrality task
  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp) {
    AliMultSelectionTask *pMultiplicityTask = AddTaskMultSelection(kFALSE);
    pMultiplicityTask->SelectCollisionCandidates(kComPhysSel);
  }

  //  CDBconnect task
  if (bDoFullJets || iDataType == kEsd) {
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    taskCDB->SetFallBackToRaw(kTRUE);
  }
  


  // ============================================= JetExtractor_R04_Toymodel =======================================================
  AliAnalysisTaskChargedJetsHadronToy *pJetExtractor_R04_Toymodel = 0;
  pJetExtractor_R04_Toymodel = AddTaskChargedJetsHadronToy("tracks_toy");
  pJetExtractor_R04_Toymodel->SetVzRange(-10,10);
  pJetExtractor_R04_Toymodel->SetUseNewCentralityEstimation(kFALSE);
  pJetExtractor_R04_Toymodel->SetNeedEmcalGeom(kFALSE);
  if(kRHEmbedParticleLevel)
    pJetExtractor_R04_Toymodel->AddTracksFromInputEvent("mcparticles");
  else
    pJetExtractor_R04_Toymodel->AddTracksFromInputEvent("tracks");
    pJetExtractor_R04_Toymodel->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
    pJetExtractor_R04_Toymodel->AddCellsFromToy("alien:///alice/cern.ch/user/h/hbossi/Distributions_Clusters.root");
  pJetExtractor_R04_Toymodel->AddTracksFromToy( "alien:///alice/cern.ch/user/r/rhaake/Toy/Distributions_PbPb_FlatMultiplicity.root");
  pJetExtractor_R04_Toymodel->SetOutputToyCellsName("toy_cells");
  // ============================================= EMCAL CORRECTIONS ======================================================
  if (bDoEmcalCorrections) {
    // Configuration of the Correction Task is handled via a YAML file, which is setup below
    // NOTE: Calling this function is equivalent to the AddTask macro, just the function is defined in the class.
    //       The AddTask macro still exists for use on the LEGO train
    // AliEmcalCorrectionTask * correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
    // correctionTask->SelectCollisionCandidates(kComPhysSel);
    // correctionTask->SetUseNewCentralityEstimation(bIsRun2);
    // correctionTask->SetNCentBins(5);
    // correctionTask->SetForceBeamType(static_cast<AliEmcalCorrectionTask::BeamType>(iBeamType));
    // 
    // // Configure and initialize
    // correctionTask->SetUserConfigurationFilename("ToyModelConfig_Combine.yaml");
    // correctionTask->Initialize();
    Printf(" ****** Starting Corrections *********** ");
    TObjArray correctionTasks;
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("data"));
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("toy"));
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("combined"));

    // Loop over all of the correction tasks to configure them
    AliEmcalCorrectionTask * tempCorrectionTask = 0;
    TIter next(&correctionTasks);
    while (( tempCorrectionTask = static_cast<AliEmcalCorrectionTask *>(next())))
    {
      Printf(" ****** Moving on to next correction task *********** ");
      tempCorrectionTask->SelectCollisionCandidates(kComPhysSel);
      tempCorrectionTask->SetUseNewCentralityEstimation(bIsRun2);
      tempCorrectionTask->SetNCentBins(5);
      tempCorrectionTask->SetForceBeamType(static_cast<AliEmcalCorrectionTask::BeamType>(iBeamType));
      // Local configuration
      tempCorrectionTask->SetUserConfigurationFilename("ToyModelConfig_Combine.yaml");
      tempCorrectionTask->Initialize();
    }
  }

  // ============================================= JetExtractor_R04_Rho_JFKT =======================================================
  // AliEmcalJetTask *pJetExtractor_R04_Rho_JFKT = AliEmcalJetTask::AddTaskEmcalJet("","", AliJetContainer::kt_algorithm, 0.2, AliJetContainer::kChargedJet, 0.15, 0.3, 0.005, AliJetContainer::E_scheme, "Jet", 0., kTRUE, kFALSE);
  // pJetExtractor_R04_Rho_JFKT->SetNeedEmcalGeom(kFALSE);
  // AliTrackContainer* trackCont03 = new AliTrackContainer("tracks_toy");
  // trackCont03->SetFilterHybridTracks(kTRUE);
  // trackCont03->SetParticlePtCut(kTrackPtCut);
  // pJetExtractor_R04_Rho_JFKT->AdoptParticleContainer(trackCont03);
  // pJetExtractor_R04_Rho_JFKT->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
  // //
  // // // ============================================= JetExtractor_R04_Rho =======================================================
  // AliAnalysisTaskRho *pJetExtractor_R04_Rho = 0;
  // pJetExtractor_R04_Rho = AddTaskRhoNew("tracks_toy", "", "RhoR020", 0.2, AliJetContainer::kTPCfid, AliJetContainer::kChargedJet, kTRUE, AliJetContainer::E_scheme, "Rho_ExLJ"); 
  // pJetExtractor_R04_Rho->SetVzRange(-10,10);
  // pJetExtractor_R04_Rho->SetNeedEmcalGeom(kFALSE);
  // pJetExtractor_R04_Rho->SetExcludeLeadJets(2); // default
  // //pJetExtractor_R04_Rho->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
  // //
  // // // ============================================= JetExtractor_R04_Rho_Mass ======================================================
  // AliAnalysisTaskRhoMass *pJetExtractor_R04_RhoMass = 0;
  // pJetExtractor_R04_RhoMass = AddTaskRhoMass("Jet_KTChargedR020_tracks_toy_pT0150_E_scheme", "tracks_toy", "", "RhoR020_mass", 0.2, "TPC", 0.01,0,0,2,kTRUE,"RhoMass"); 
  // pJetExtractor_R04_RhoMass->SetVzRange(-10,10);
  // pJetExtractor_R04_RhoMass->SetNeedEmcalGeom(kFALSE);
  // pJetExtractor_R04_RhoMass->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
  // 
  // // // ============================================= JetExtractor_R04_JF_Toy =============================================
  // AliEmcalJetTask *pJetExtractor_R04_JF_Toy = AliEmcalJetTask::AddTaskEmcalJet( "","",AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.005, AliJetContainer::E_scheme, "Jet", 0., kTRUE, kFALSE);
  // pJetExtractor_R04_JF_Toy->SetNeedEmcalGeom(kFALSE);
  // AliTrackContainer* trackCont02 = new AliTrackContainer("tracks_toy");
  // //trackCont02->SetFilterHybridTracks(kTRUE);
  // trackCont02->SetParticlePtCut(kTrackPtCut);
  // pJetExtractor_R04_JF_Toy->AdoptParticleContainer(trackCont02);
  // 
  // AliEmcalJetUtilityGenSubtractor* genSub = (AliEmcalJetUtilityGenSubtractor*)pJetExtractor_R04_JF_Toy->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
  // genSub->SetGenericSubtractionJetMass(kTRUE);
  // genSub->SetGenericSubtractionExtraJetShapes(kTRUE);
  // genSub->SetUseExternalBkg(kTRUE);
  // genSub->SetRhoName("RhoR020");
  // genSub->SetRhomName("RhoR020_mass");
  // pJetExtractor_R04_JF_Toy->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
  // 
  // // 
  // // // ============================================= JetExtractor_R04_allJets_Toy ========================================
  // AliAnalysisTaskJetExtractor *pJetExtractor_R04_allJets_Toy = 0;
  // pJetExtractor_R04_allJets_Toy = AddTaskJetExtractor("tracks_toy", "", "Jet_AKTChargedR040_tracks_toy_pT0150_E_scheme", "RhoR020", 0.4, "", "allJets");
  // pJetExtractor_R04_allJets_Toy->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
  // pJetExtractor_R04_allJets_Toy->SetVzRange(-10,10);
  // pJetExtractor_R04_allJets_Toy->SetNeedEmcalGeom(kFALSE);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->SetSaveConstituentPID(0);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->SetSaveConstituentsIP(0);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->SetSaveMCInformation(kTRUE);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(-20,10, 0.001);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(10,20, 0.01);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(20,30, 0.1);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(30,40, 0.5);
  // pJetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(40,200, 1.0);
  // pJetExtractor_R04_allJets_Toy->ActivateTrueJetMatching(Form("JetPartLevel_AKTChargedR%03d_mcparticles_pT0150_E_scheme", (Int_t)(kRHJetRadius*100)), "");
  // pJetExtractor_R04_allJets_Toy->SetTrueJetMatchingRadius(0.25);
  // pJetExtractor_R04_allJets_Toy->SetTrueJetRhoMassName("");
  // pJetExtractor_R04_allJets_Toy->SetTrueParticleArrayName("mcparticles");
  // pJetExtractor_R04_allJets_Toy->GetJetContainer(0)->SetRhoMassName("RhoR020_mass");
  // pJetExtractor_R04_allJets_Toy->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
  // //
  // // 
  // 






  // ============================================= Other stuff, don't touch! =======================================================
  TObjArray *pTopTasks = pMgr->GetTasks();
  for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
    AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
    if (!pTask) continue;
    if (pTask->InheritsFrom("AliAnalysisTaskEmcal")) {
      AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
      Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcal->GetName());
      pTaskEmcal->SetForceBeamType(iBeamType);
    }
  }

  if (!pMgr->InitAnalysis()) return 0;
  pMgr->PrintStatus();

  // pMgr->SetUseProgressBar(kTRUE, 250);

  TFile *pOutFile = new TFile("train.root","RECREATE");
  pOutFile->cd();
  pMgr->Write();
  pOutFile->Close();
  delete pOutFile;

  if (iStartAnalysis == 1) { // start local analysis
    TChain* pChain = 0;
    if (iDataType == kAod) {
      #ifdef __CLING__
      std::stringstream aodChain;
      aodChain << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/CreateAODChain.C(";
      aodChain << "\"" << sLocalFiles.Data() << "\", ";
      aodChain << iNumEvents << ", ";
      aodChain << 0 << ", ";
      aodChain << std::boolalpha << kFALSE << ");";
      pChain = reinterpret_cast<TChain *>(gROOT->ProcessLine(aodChain.str().c_str()));
      #else
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
      #endif
    }
    else {
      #ifdef __CLING__
      std::stringstream esdChain;
      esdChain << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/CreateESDChain.C(";
      esdChain << "\"" << sLocalFiles.Data() << "\", ";
      esdChain << iNumEvents << ", ";
      esdChain << 0 << ", ";
      esdChain << std::boolalpha << kFALSE << ");";
      pChain = reinterpret_cast<TChain *>(gROOT->ProcessLine(esdChain.str().c_str()));
      #else
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
      #endif
    }

    // start analysis
    Printf("Starting Analysis...");
    pMgr->StartAnalysis("local", pChain, iNumEvents);
  }
  else if (iStartAnalysis == 2) {  // start grid analysis
    StartGridAnalysis(pMgr, cTaskName, cGridMode);
  }

  return pMgr;
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetPerformance.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskChargedJetsHadronToy.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetExtractor.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");

}

void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode)
{
  Int_t maxFilesPerWorker = 4;
  Int_t workerTTL = 7200;
  const char* runNumbers = "246945";
  const char* pattern = "pass1/AOD194/*/AliAOD.root";
  const char* gridDir = "/alice/data/2015/LHC15o";
  const char* additionalCXXs = "";
  const char* additionalHs = "";

  AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, cGridMode, runNumbers, pattern, additionalCXXs, additionalHs, maxFilesPerWorker, workerTTL, kFALSE);
  pMgr->SetGridHandler(plugin);

  // start analysis
   Printf("Starting GRID Analysis...");
   pMgr->SetDebugLevel(0);
   pMgr->StartAnalysis("grid");
}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC)
{
  TDatime currentTime;
  TString tmpName(uniqueName);

  // Only add current date and time when not in terminate mode! In this case the exact name has to be supplied by the user
  if (strcmp(gridMode, "terminate")) {
    tmpName += "_";
    tmpName += currentTime.GetDate();
    tmpName += "_";
    tmpName += currentTime.GetTime();
  }

  TString macroName("");
  TString execName("");
  TString jdlName("");
  macroName = Form("%s.C", tmpName.Data());
  execName = Form("%s.sh", tmpName.Data());
  jdlName = Form("%s.jdl", tmpName.Data());

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetRunMode(gridMode);

  // Here you can set the (Ali)PHYSICS version you want to use
  plugin->SetAliPhysicsVersion("vAN-20181128-1");

  plugin->SetGridDataDir(gridDir); // e.g. "/alice/sim/LHC10a6"
  plugin->SetDataPattern(pattern); //dir structure in run directory

  if (!isMC) plugin->SetRunPrefix("000");

  plugin->AddRunList(runNumbers);

  plugin->SetGridWorkingDir(Form("work/%s",tmpName.Data()));
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

  plugin->SetAnalysisSource(additionalCode.Data());

  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetAnalysisMacro(macroName.Data());
  plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
  plugin->SetExecutable(execName.Data());
  plugin->SetTTL(workerTTL);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(jdlName.Data());
  plugin->SetPrice(1);
  plugin->SetSplitMode("se");

  // merging via jdl
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);

  return plugin;
}
