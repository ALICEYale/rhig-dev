/// \file runEMCalJetEmbedDataTask.C
/// \embedding MB data into detector level pythia and using the jet extractor toy to extract info for ML
///
/// \ingroup EMCALJETFW
/// The script runEMCalJetEmbedDataTask.sh in the same folder allows to easily set
/// the input values
///
/// \author Laura Havener <laura.brittany.havener@cern.ch>, Yale University and Hannah Bossi <hannah.bossi@yale.edu>, Yale University
/// \date Mar 12, 2019

/** this comes from $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/runEmbeddingAnalysis.C
  * Includes the ability to match jets from the hybrid level to embedded particle level.                       
  *                                                                                                            
  * The various defined jet collections include:                                                               
  * - "hybridLevelJets": PbPb + embedded detector level jets                                                   
  * - "PbPbLevelJets": PbPb only level jets                                                                    
  * - "detLevelJets": embedded detector level jets                                                             
  * - "partLevelJets": embeddded particle level jets                                                           
  *                                                                                                            
  * Such collections should be available in the case of embedding MC productions. Embedding                    
  * a real pp data set (which would not have embedded particle level information available)                    
  * must be setup with more care.                                                                              
  *                                                                                                            
  * Some options can be modified in the section under the comment "Configuration options".
  */

class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisGrid;
class AliAnalysisManager;
class AliAnalysisAlien;
class AliPhysicsSelectionTask;
class AliMultSelectionTask;
class AliCentralitySelectionTask;
class AliTaskCDBconnect;

class AliClusterContainer;
class AliParticleContainer;
class AliJetContainer;

class AliAnalysisTaskEmcalEmbeddingHelper;
class AliEmcalCorrectionTask;
class AliEmcalJetTask;
class AliAnalysisTaskEmcalJetSample;
class AliJetResponseMaker;

namespace PWGJE {
  namespace EMCALJetTasks {
    class AliAnalysisTaskEmcalJetPerfTree;
  }
}

namespace PWGJE {
  namespace EMCALJetTasks {
    class AliEmcalJetTaggerTaskFast;
  }
}

// Include AddTask macros for ROOT 6 compatibility
#ifdef __CLING__
// Tell ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "OADB/macros/AddTaskPhysicsSelection.C"
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"
#include "OADB/macros/AddTaskCentrality.C"
#include "PWGPP/PilotTrain/AddTaskCDBconnect.C"
#include "PWG/EMCAL/macros/AddTaskEmcalSample.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskJetExtractor.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetPerfTree.C"
#include "PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C"
#endif

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runEMCalJetEmbedDataTask(
				       const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
				       const char   *cRunPeriod     = "LHC11h",                                // set the run period
				       const char   *cLocalFiles    = "files_test.txt",   // set the local list file
				       const UInt_t  iNumEvents     = 1000,                                    // number of events to be analyzed
				       const UInt_t  kPhysSel       = AliVEvent::kAnyINT |
				           AliVEvent::kCentral | AliVEvent::kSemiCentral,                          // physics selection
				       // const UInt_t kPhysSel        = 0,                                      // for MC!
				       const char   *cTaskName      = "EMCalEmbeddingDataExtractor",                           // sets name of analysis manager
				       //				       const Bool_t  bDoChargedJets = kTRUE,
				       //				       const Bool_t  bDoFullJets    = kFALSE,
				       const char   *obsolete       = "",                                      // Previous handled the ocdb settings, but obsolete due to CDBconnect task
				       // 0 = only prepare the analysis manager but do not start the analysis
				       // 1 = prepare the analysis manager and start the analysis
				       // 2 = launch a grid analysis
				       Int_t         iStartAnalysis = 1,
				       const UInt_t  iNumFiles      = 5,                                     // number of files analyzed locally
				       const char   *cGridMode      = "test"
				       )
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Global Train Configuration file//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
/*  Bool_t fSpecialPtHard =kTRUE;
  Bool_t fGlobalParam_VetoMB= kFALSE;
  Bool_t fParamDoSmearing=kTRUE;
  Bool_t fParamRunCompositionCorrection=kTRUE;*/
  /*************************************************************************/
  //Jet Response Maker Configuration
  /*************************************************************************/
  
  
  /*  Bool_t doZaxis = kFALSE;
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
  Double_t kGhostArea = 0.01;
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
    }*/
  /*************************************************************************/
  //Global Names
  /*************************************************************************/
  
  /*  char pTString[200];
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

  // Variable for Ruediger
  Double_t kRHJetRadius = 0.4;
  Bool_t kRHBuiltinEventSelection = kFALSE;
  Bool_t kRHUseHFTracks = kTRUE;
  Bool_t kRHEmbedParticleLevel   = kTRUE;*/

  ////////end of the stuff from config file from Ruediger and Hannah/////
  
  //setup period
  TString sRunPeriod(cRunPeriod);
  sRunPeriod.ToLower();

  // set Run 2
  Bool_t bIsRun2 = kFALSE;
  if (sRunPeriod.Length() == 6 && sRunPeriod.BeginsWith("lhc15")) bIsRun2 = kTRUE;

  cout << "run2: " << bIsRun2 << endl;

  //set beam type
  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp;
  if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
    iBeamType = AliAnalysisTaskEmcal::kAA;
  }
  else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f" ||
      sRunPeriod == "lhc16q" || sRunPeriod == "lhc16r" || sRunPeriod == "lhc16s" ||
      sRunPeriod == "lhc16t") {
    iBeamType = AliAnalysisTaskEmcal::kpA;
    }

  //ghost area
  Double_t kGhostArea = 0.01;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) kGhostArea = 0.005;

  //setup track container
  AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);
  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());

  //set data file type
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

  //load macros needed for the analysis
  #ifndef __CLING__
  LoadMacros();
  #endif

  ////////////////////////                                                                                                                              // Configuration options                                                                                                                            
  ////////////////////////                                                                                                                             
  // Use full or charged jets                                                                                                                         
  const bool fullJets =true;
  // Embedding files list                                                                                                                             
  const std::string embeddedFilesList = "files_mctest.txt";
  // If true, events that are not selected in the PbPb will not be used for embedding.                                                                  // This ensures that good embedded events are not wasted on bad PbPb events.                                                                      
  const bool internalEventSelection = false;
  // Background subtraction                                                                                                                         
  const bool enableBackgroundSubtraction = false;
  // Do jet matching                                                                                                                                 
  const bool useJetTagger = true;
  // General track and cluster cuts (used particularly for jet finding)                                                                             
  const Double_t minTrackPt = 0.15;
  const Double_t minClusterPt = 0.3;

  // Debug options                                                                                                                                       //AliLog::SetClassDebugLevel("AliAnalysisTaskEmcalEmbeddingHelper", AliLog::kDebug+0); 

  // Determine track, cluster, and cell names                                                                                                        
  const bool IsEsd = (iDataType == kEsd);
  // These names correspond to the _uncombined_ input objects that we are interestd in the external (embedded) event                                
  TString externalEventParticlesName = "mcparticles";
  // Empty because there are no particle level clusters                                                                                              
  TString externalEventClustersName = "";
  
  // General input object names                                                                                                                      
  TString tracksName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack, IsEsd);
  TString emcalCellsName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCaloCells, IsEsd);
  TString emcalCellsNameCombined = emcalCellsName + "Combined";
  TString clustersName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster, IsEsd);
  // Combined (PbPb + embedded det level) clusters                                                                                                  
  TString clustersNameCombined = clustersName + "Combined";

  // Handle charged jets                                                                                                                             
  if (fullJets == false) {
    emcalCellsName = "";
    emcalCellsNameCombined = "";
    clustersName = "";
    clustersNameCombined = "";
  }

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);
  
  // Create Input Handler                                                                                                                            
  if (iDataType == kAod) {
    AliAODInputHandler * pESDHandler = AliAnalysisTaskEmcal::AddAODHandler();
  }
  else {
    AliESDInputHandler * pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
  }
  
  // Physics selection task
  if (iDataType == kEsd) {
    AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();
  }

  // Centrality task
  // The Run 2 condition is too restrictive, but until the switch to MultSelection is complete, it is the best we can do                              
  //  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp && bIsRun2 == kFALSE) {
  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  // AliMultSelection                                                                                                                                    // Works for pp, pPb, and PbPb for the periods that it is calibrated                                                                             
  if (bIsRun2 == kTRUE) {
    AliMultSelectionTask * pMultSelectionTask = AddTaskMultSelection(kFALSE);
    pMultSelectionTask->SelectCollisionCandidates(kPhysSel);
  }

  // CDBconnect task                                                                                                                               
  //  if (fullJets) {
    AliTaskCDBconnect * taskCDB = AddTaskCDBconnect();
    taskCDB->SetFallBackToRaw(kTRUE);
    //  }
  /* if (bDoFullJets || iDataType == kEsd) {
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    taskCDB->SetFallBackToRaw(kTRUE);
    }*/

  // Create the Embedding Helper
  AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
  // Set the file pattern. This example uses ptHardBin 4 of LHC12a15e_fix.
  // The pT hard bin and anchor run can also be set by adding a printf() wild card to the string (ie %d)
  // See the documentation of AliAnalysisTaskEmcalEmbeddingHelper::GetFilenames()
  // For the grid. Include "alien://" and don't use this locally!! It will be very slow and cause strain on the grid!
  //  embeddingHelper->SetFilePattern("~/alice/sim/2018/LHC18b8_cent_woSDD/20/282306/AliAOD.root");
  // For local use. File should be formatted the same as the normal list of input files (one filename per line).
  // Again, don't use AliEn locally!! It will be very slow and cause strain on the grid!
  // The pt hard bin should be set via the filenames in this file                                                                                        // If using a file pattern, it could be configured via embeddingHelper->SetPtHardBin(ptHardBin);                                                 
  embeddingHelper->SetFileListFilename(embeddedFilesList.c_str());
    // If the embedded file is an ESD, then set:
  //embeddingHelper->SetESD();
  // Add additional configure as desired.
  // For some examples...
  // ... randomly select which file to start from:
  //  embeddingHelper->SetRandomFileAccess(kTRUE);
  // ... Start from a random event within each file
  //  embeddingHelper->SetRandomEventNumberAccess(kTRUE);
  // ... Set pt hard bin properties
  //  embeddingHelper->SetPtHardBin(4);
  // Some example settings for LHC12a15e_fix (anchored to LHC11h)                                                                                     
  embeddingHelper->SelectCollisionCandidates(kPhysSel);
  embeddingHelper->SetNPtHardBins(11);
  embeddingHelper->SetMCRejectOutliers();
  // Setup internal event selection and additional configuration options                                                                             
  embeddingHelper->SetConfigurationPath("EmbeddingConfigurationExample.yaml");

  // etc..
  // As your last step, always initialize the helper!
  embeddingHelper->Initialize();

  if (fullJets) {
    // EMCal corrections                                                                                                                            
    TObjArray correctionTasks;
    // Create the Correction Tasks                                                                                                                        // "data" corresponds to the PbPb level                                                                                                              // "embed" corresponds to the embedded detector level                                                                                            
    // "combined" corresponds to the hybrid (PbPb + embedded detector) level                                                                          
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("data"));
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("embed"));
    // It is important that combined is last!                                                                                                        
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("combined"));
    // Loop over all of the correction tasks to configure them                                                                                        
    AliEmcalCorrectionTask * tempCorrectionTask = 0;
    TIter next(&correctionTasks);
    while (( tempCorrectionTask = static_cast<AliEmcalCorrectionTask *>(next())))
      {
	tempCorrectionTask->SelectCollisionCandidates(kPhysSel);
	// Configure centrality                                                                                                                      
	tempCorrectionTask->SetNCentBins(5);
	if (bIsRun2) {
	  tempCorrectionTask->SetUseNewCentralityEstimation(kTRUE);
	}
	tempCorrectionTask->SetUserConfigurationFilename("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/EMCalCorrectionTaskEmbeddingExample.yaml");

	tempCorrectionTask->Initialize();
      }
  }

  // Background                                                                                                                                        
  std::string sRhoChargedName = "";
  std::string sRhoFullName = "";
  if (iBeamType != AliAnalysisTaskEmcal::kpp && enableBackgroundSubtraction == true) {
    const AliJetContainer::EJetAlgo_t rhoJetAlgorithm = AliJetContainer::kt_algorithm;
    const AliJetContainer::EJetType_t rhoJetType = AliJetContainer::kChargedJet;
    const AliJetContainer::ERecoScheme_t rhoRecoScheme = AliJetContainer::pt_scheme;
    const double rhoJetRadius = 0.4;
    sRhoChargedName = "Rho";
    sRhoFullName = "Rho_Scaled";

    AliEmcalJetTask *pKtChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", rhoJetAlgorithm, rhoJetRadius, rhoJetType, minTrackPt, 0, kGhostArea, rhoRecoScheme, "Jet", 0., kFALSE, kFALSE);
    pKtChJetTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskRho * pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChargedName.c_str(), rhoJetRadius);
    pRhoTask->SetExcludeLeadJets(2);
    pRhoTask->SelectCollisionCandidates(kPhysSel);
    pRhoTask->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
    pRhoTask->SetUseNewCentralityEstimation(kTRUE);
    
    if (fullJets)
      {
	TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
	TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
	pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
      }
  }


    /*  // Background
  TString sRhoChName;
  TString sRhoFuName;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) {
    sRhoChName = "Rho";
    sRhoFuName = "Rho_Scaled";

    AliEmcalJetTask *pKtChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
    pKtChJetTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskRho* pRhoTask =AddTaskRhoNew("usedefault", "usedefault", sRhoChName, 0.4);
    pRhoTask->SetExcludeLeadJets(2);
    pRhoTask->SelectCollisionCandidates(kPhysSel);

    if (bDoFullJets) {
      TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
      TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
      pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
    }
    }*/

   // Jet finding                                                                                                                                   
  const AliJetContainer::EJetAlgo_t jetAlgorithm = AliJetContainer::antikt_algorithm;
  const Double_t jetRadius = 0.2;
  AliJetContainer::EJetType_t jetType = AliJetContainer::kFullJet;
  const AliJetContainer::ERecoScheme_t recoScheme = AliJetContainer::pt_scheme;
  const char * label = "Jet";
  const Double_t minJetPt = 1;
  const Bool_t lockTask = kTRUE;
  const Bool_t fillGhosts = kFALSE;

  // Do not pass clusters if we are only looking at charged jets                                                                                     
  if (fullJets == false) {
    jetType = AliJetContainer::kChargedJet;
  }

  ///////                                                                                                                                               // Particle level PYTHIA jet finding                                                                                                                 ///////                                                                                                                                           
  AliEmcalJetTask * pFullJetTaskPartLevel = AliEmcalJetTask::AddTaskEmcalJet("mcparticles", "", jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "partLevelJets", minJetPt, lockTask, fillGhosts);
  pFullJetTaskPartLevel->SelectCollisionCandidates(kPhysSel);
  pFullJetTaskPartLevel->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  pFullJetTaskPartLevel->SetUseNewCentralityEstimation(kTRUE);

  ///////                                                                                                                                                // External event (called embedding) settings for particle level PYTHIA jet finding                                                               
  ///////                                                                                                                                           
  // Setup the tracks properly to be retrieved from the external event                                                                                // It does not matter here if it's a Particle Container or MCParticleContainer                                                                       
  AliParticleContainer * partLevelTracks = pFullJetTaskPartLevel->GetParticleContainer(0);
  // Called Embedded, but really just means get from an external event!                                                                             
  partLevelTracks->SetIsEmbedding(kTRUE);

  ///////                                                                                                                                               // Detector level PYTHIA jet finding                                                                                                                 ///////                                                                                                                                              
  AliEmcalJetTask * pFullJetTaskDetLevel = AliEmcalJetTask::AddTaskEmcalJet(tracksName.Data(), clustersName.Data(),jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "detLevelJets", minJetPt, lockTask, fillGhosts);
  pFullJetTaskDetLevel->SelectCollisionCandidates(kPhysSel);
  pFullJetTaskDetLevel->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  pFullJetTaskDetLevel->SetUseNewCentralityEstimation(kTRUE);

  ///////                                                                                                                                                // External event (embedding) settings for det level PYTHIA jet finding                                                                                ///////                                                                                                                                               // Tracks                                                                                                                                             // Uses the name of the container passed into AliEmcalJetTask                                                                                   
  AliTrackContainer * tracksDetLevel = pFullJetTaskDetLevel->GetTrackContainer(0);
  // Get the det level tracks from the external event!                                                                                               
  tracksDetLevel->SetIsEmbedding(kTRUE);


  // Clusters                                                                                                                                       
  if (fullJets) {
    // Uses the name of the container passed into AliEmcalJetTask                                                                                   
    AliClusterContainer * clustersDetLevel = pFullJetTaskDetLevel->GetClusterContainer(0);
    // Get the det level clusters from the external event!                                                                                           
    clustersDetLevel->SetIsEmbedding(kTRUE);

    // Additional configuration                                                                                                                      
    clustersDetLevel->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  ///////                                                                                                                                              
  // PbPb jet finding                                                                                                                                  
  ///////                                                                                                                                              
  AliEmcalJetTask * pFullJetTask = AliEmcalJetTask::AddTaskEmcalJet(tracksName.Data(), clustersName.Data(),
								    jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "PbPbLevelJets", minJetPt, lockTask, fillGhosts);
  pFullJetTask->SelectCollisionCandidates(kPhysSel);
  pFullJetTask->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  pFullJetTask->SetUseNewCentralityEstimation(kTRUE);
  if (fullJets) {
    pFullJetTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  ///////                                                                                                                                              
  // Hybrid (PbPb + Detector) level PYTHIA jet finding                                                                                                 
  ///////                                                                                                                                              
  // Sets up PbPb tracks and clusters                                                                                                                  
  // NOTE: The clusters name is different here since we output to a different branch!                                                                  
  AliEmcalJetTask * pFullJetTaskHybrid = AliEmcalJetTask::AddTaskEmcalJet(tracksName.Data(), clustersNameCombined.Data(),
									  jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "hybridLevelJets", minJetPt, lockTask, fillGhosts);
  pFullJetTaskHybrid->SelectCollisionCandidates(kPhysSel);
  pFullJetTaskHybrid->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  pFullJetTaskHybrid->SetUseNewCentralityEstimation(kTRUE);
  if (fullJets) {
    pFullJetTaskHybrid->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }
  ///////                                                                                                                                              
  // External event (ie embedding) settings for PbPb jet finding (adds detector level PYTHIA)                                                          
  ///////                                                                                                                                              
  // Add embedded tracks and clusters to jet finder                                                                                                    
  // Tracks                                                                                                                                            
  AliTrackContainer * tracksEmbedDetLevel = new AliTrackContainer(tracksName.Data());
  // Get the det level tracks from the external event!                                                                                                 
  tracksEmbedDetLevel->SetIsEmbedding(kTRUE);

  tracksEmbedDetLevel->SetParticlePtCut(minTrackPt);
  pFullJetTaskHybrid->AdoptTrackContainer(tracksEmbedDetLevel);
  // Clusters                                                                                                                                          
  // Already combined in clusterizer, so we shouldn't add an additional cluster container here  

  //////////////////////////////////////////                                                                                                           
  // Sample Jet Tasks                                                                                                                                  
  //////////////////////////////////////////                                                                                                           

  // Use Sample task to show how the embedding did                                                                                                     
  // Need both the type and the string for various classes...                                                                                          
  AliEmcalJet::JetAcceptanceType acceptanceType = AliEmcalJet::kEMCALfid;
  const std::string acceptanceTypeStr = "EMCALfid";

  ///////                                                                                                                                              
  // Particle level PYTHIA sample task                                                                                                                 
  ///////                                                                                                                                              
  AliAnalysisTaskEmcalJetSample * sampleTaskPartLevel = AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample("mcparticles", "", "", "partLevelJets");

  // Set embedding                                                                                                                                     
  AliParticleContainer * partCont = sampleTaskPartLevel->GetParticleContainer(0);
  partCont->SetIsEmbedding(kTRUE);

  partCont->SetParticlePtCut(0.15);
  sampleTaskPartLevel->SetHistoBins(600, 0, 300);
  sampleTaskPartLevel->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  sampleTaskPartLevel->SelectCollisionCandidates(kPhysSel);
  sampleTaskPartLevel->SetUseNewCentralityEstimation(kTRUE);

  AliJetContainer* jetContPartLevel02 = sampleTaskPartLevel->AddJetContainer(pFullJetTaskPartLevel->GetName(), acceptanceType, jetRadius);

  ///////                                                                                                                                              
  // Detector level PYTHIA sample task                                                                                                                 
  ///////                                                                                                                                              
  // Cells are left empty because special care is needed to retrieve the embedded cells. Such code is not in the                                       
  // sample task, but is available in the EMCal Correction Framework                                                                                   
  AliAnalysisTaskEmcalJetSample * sampleTaskDetLevel = AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample(tracksName.Data(), clustersName.Data(), "", "detLevelJets");
  sampleTaskDetLevel->SetUseNewCentralityEstimation(kTRUE);

  // Tracks                                                                                                                                            
  // Set embedding                                                                                                                                     
  AliParticleContainer * detLevelPartCont = sampleTaskDetLevel->GetParticleContainer(0);
  detLevelPartCont->SetIsEmbedding(kTRUE);

  // Settings                                                                                                                                          
  detLevelPartCont->SetParticlePtCut(0.15);

  if (fullJets) {
    // Clusters                                                                                                                                        
    // Set embedding                                                                                                                                   
    sampleTaskDetLevel->GetClusterContainer(0)->SetIsEmbedding(kTRUE);
    // Settings                                                                                                                                        
    sampleTaskDetLevel->GetClusterContainer(0)->SetClusECut(0.);
    sampleTaskDetLevel->GetClusterContainer(0)->SetClusPtCut(0.);
    sampleTaskDetLevel->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
    sampleTaskDetLevel->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
    sampleTaskDetLevel->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  sampleTaskDetLevel->SetHistoBins(600, 0, 300);
  sampleTaskDetLevel->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  sampleTaskDetLevel->SelectCollisionCandidates(kPhysSel);

  AliJetContainer* detLevelJetCont02 = sampleTaskDetLevel->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, "detLevelJets");

  ///////                                                                                                                                              
  // PbPb sample task                                                                                                                                  
  ///////                                                                                                                                              
  AliAnalysisTaskEmcalJetSample * sampleTask = AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample(tracksName.Data(), clustersName.Data(), emcalCellsName.Data(), "PbPbLevelJets");
  if (fullJets) {
    sampleTask->GetClusterContainer(0)->SetClusECut(0.);
    sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
    sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
    sampleTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
    sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }
  sampleTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  sampleTask->SetHistoBins(600, 0, 300);
  sampleTask->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  sampleTask->SelectCollisionCandidates(kPhysSel);
  sampleTask->SetUseNewCentralityEstimation(kTRUE);


  AliJetContainer* jetContPbPb02 = sampleTask->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, "PbPbLevelJets");

  ///////                                                                                                                                              
  // PbPb + Detector level PYTHIA sample task                                                                                                          
  ///////                                                                                                                                              
  // NOTE: The clusters name is different here since we output to a different branch!                                                                  
  AliAnalysisTaskEmcalJetSample * sampleTaskHybrid = AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample(tracksName.Data(), clustersNameCombined.Data\
													  (), emcalCellsNameCombined.Data(), "hybridLevelJets");

  // PbPb tracks settings                                                                                                                              
  sampleTaskHybrid->GetParticleContainer(0)->SetParticlePtCut(0.15);
  sampleTaskHybrid->SetUseNewCentralityEstimation(kTRUE);


  // Embed tracks                                                                                                                                      
  tracksDetLevel = new AliTrackContainer(tracksName.Data());
  // Get the det level tracks from the external event!                                                                                                 
  tracksDetLevel->SetIsEmbedding(kTRUE);
  // Settings                                                                                                                                          
  tracksDetLevel->SetParticlePtCut(0.15);
  sampleTaskHybrid->AdoptTrackContainer(tracksDetLevel);

  if (fullJets) {
    // PbPb + Detector Level clusters settings                                                                                                         
    sampleTaskHybrid->GetClusterContainer(0)->SetClusECut(0.);
    sampleTaskHybrid->GetClusterContainer(0)->SetClusPtCut(0.);
    sampleTaskHybrid->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
    sampleTaskHybrid->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
    sampleTaskHybrid->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  sampleTaskHybrid->SetHistoBins(600, 0, 300);
  sampleTaskHybrid->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
  sampleTaskHybrid->SelectCollisionCandidates(kPhysSel);

  AliJetContainer * jetContHybrid02 = sampleTaskHybrid->AddJetContainer(pFullJetTaskHybrid->GetName(), acceptanceType, jetRadius);

  // Jet Matching Tasks                                                                                                                                
  // 0.3 is the default max matching distance from the EMCal Jet Tagger                                                                                
  double maxGeoMatchingDistance = 0.3;
  double fractionSharedMomentum = 0.5;

  if (useJetTagger) {
    /**                                                                                                                                                
     * Jet taggers:                                                                                                                                    
     *                                                                                                                                                 
     * For hybrid -> part level for a response matrix:                                                                                                 
     * strategy is to tag:           hybrid level <-> embedded detector level                                                                          
     * and:                    embedded det level <-> embedded particle level                                                                          
     *                                                                                                                                                 
     * For hybrid -> det level for pp embedding:                                                                                                       
     * Strategy is to only tag       hybrid level <-> embedded detector level                                                                          
     * because that is sufficient for those purposes                                                                                                   
     *                                                                                                                                                 
     * NOTE: The shared momentum fraction is set here, but it only affects the                                                                         
     * output histograms (not the matchig itself). For more customization,                                                                             
     * the user may want to handle this at the user task level.                                                                                        
     *                                                                                                                                                 
     * NOTE: The PbPb level jets are not necessary for this matching. Only the                                                                         
     * other three jet collections.                                                                                                                    
     */

    // Hybrid (PbPb + embed) jets are the "base" jet collection                                                                                        
    const std::string hybridLevelJetsName = pFullJetTaskHybrid->GetName();
    // Embed det level jets are the "tag" jet collection                                                                                               
    const std::string detLevelJetsName = pFullJetTaskDetLevel->GetName();
    // Centrality estimotor                                                                                                                            
    const std::string centralityEstimator = "V0M";
    // NOTE: The AddTask macro is "AddTaskEmcalJetTaggerFast" ("Emcal" is removed for the static definition...)                                        
    PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast * jetTaggerDetLevel = PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::AddTaskJetTaggerFast(                  
     hybridLevelJetsName.c_str(),         // "Base" jet collection which will be tagged											 
     detLevelJetsName.c_str(),       // "Tag" jet collection which will be used to tag (and will be tagged)                                         
     jetRadius,                      // Jet radius                                                                                                  
     "",                             // Hybrid ("base") rho name                                                                                    
     "",                             // Det level ("tag") rho name                                                                                  
     "",                             // tracks to attach to the jet containers. Not meaningful here, so left empty                                  
     "",                             // clusters to attach to the jet conatiners. Not meaingful here, so left empty (plus, it's not the same for the two jet collections!)
     acceptanceTypeStr.c_str(),      // Jet acceptance type for the "base" collection                                                               
     centralityEstimator.c_str(),    // Centrality estimator                                                                                        
     kPhysSel,                       // Physics selection                                                                                           
     ""                              // Trigger class. We can just leave blank, as it's only used in the task name                                  
     );
    // Task level settings                                                                                                                             
    jetTaggerDetLevel->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
    // Tag via geometrical matching                     

    jetTaggerDetLevel->SetUseNewCentralityEstimation(kTRUE);
                                                                                               
    jetTaggerDetLevel->SetJetTaggingMethod(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::kGeo);
    // Tag the closest jet                                                                                                                             
    jetTaggerDetLevel->SetJetTaggingType(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::kClosest);
    // Don't impose any additional acceptance cuts beyond the jet containers                                                                           
    jetTaggerDetLevel->SetTypeAcceptance(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::kNoLimit);
    // Use default matching distance                                                                                                                   
    jetTaggerDetLevel->SetMaxDistance(maxGeoMatchingDistance);
    // Redundant, but done for completeness                                                                                                            
    jetTaggerDetLevel->SelectCollisionCandidates(kPhysSel);
    // Set fraction shared momentum for output hists                                                                                                   
    jetTaggerDetLevel->SetMinFractionShared(fractionSharedMomentum);

    // Reapply the max track pt cut off to maintain energy resolution and avoid fake tracks                                                            
    AliJetContainer * hybridJetCont = jetTaggerDetLevel->GetJetContainer(0);
    hybridJetCont->SetMaxTrackPt(100);
    AliJetContainer * detLevelJetCont = jetTaggerDetLevel->GetJetContainer(1);
    detLevelJetCont->SetMaxTrackPt(100);

    // Embed det level jets are the "base" jet collection                                                                                              
    // Embed part level jets are the "tag" jet collection                                                                                              
    const std::string partLevelJetsName = pFullJetTaskPartLevel->GetName();

    PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast * jetTaggerPartLevel = PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::AddTaskJetTaggerFast(
    detLevelJetsName.c_str(),       // "Base" jet collection which will be tagged               									   
    partLevelJetsName.c_str(),      // "Tag" jet collection which will be used to tag (and will be tagged)								
    jetRadius,                      // Jet radius                                                                                                    
    "",                             // Det level ("base") rho name                                                                                   
    "",                             // Part level ("tag") rho name                                                                                   
    "",                             // tracks to attach to the jet containers. Not meaningful here, so left empty            						 
    "",                             // clusters to attach to the jet conatiners. Not meaingful here, so left empty (plus, it's not the same for the two jet collections!)
    acceptanceTypeStr.c_str(),      // Jet acceptance type for the "base" collection                                                                 
    centralityEstimator.c_str(),    // Centrality estimator                                                                                          
    kPhysSel,                       // Physics selection                                                                                             
    ""                              // Trigger class. We can just leave blank, as it's only used in the task name                                    
    );
    // Task level settings                                                                                                                             
    jetTaggerPartLevel->SetRecycleUnusedEmbeddedEventsMode(internalEventSelection);
    // Tag via geometrical matching                                                                                                                    
    jetTaggerPartLevel->SetUseNewCentralityEstimation(kTRUE);

    jetTaggerPartLevel->SetJetTaggingMethod(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::kGeo);
    // Tag the closest jet                                                                                                                             
    jetTaggerPartLevel->SetJetTaggingType(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::kClosest);
    // Don't impose any additional acceptance cuts beyond the jet containers                                                                           
    jetTaggerPartLevel->SetTypeAcceptance(PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::kNoLimit);
    // Use default matching distance                                                                                                                   
    jetTaggerPartLevel->SetMaxDistance(maxGeoMatchingDistance);
    // Redundant, but done for completeness                                                                                                            
    jetTaggerPartLevel->SelectCollisionCandidates(kPhysSel);
    // We don't want to apply a shared momentum fraction cut here, as it's not meaningful                                                              

    // Reapply the max track pt cut off to maintain energy resolution and avoid fake tracks                                                            
    // However, don't apply to the particle level jets which don't suffer this effect                                                                  
    detLevelJetCont = jetTaggerPartLevel->GetJetContainer(0);
    detLevelJetCont->SetMaxTrackPt(100);
  }

  /*//jets from pythia at particle level
  AliEmcalJetTask *JetExtractor_R04_JF_PYTHIA_PartLevel = AliEmcalJetTask::AddTaskEmcalJet("","",AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.005, AliJetContainer::E_scheme, "JetPartLevel", 0., kTRUE, kFALSE);
  JetExtractor_R04_JF_PYTHIA_PartLevel->SetNeedEmcalGeom(kFALSE);
  AliMCParticleContainer* trackCont = new AliMCParticleContainer("mcparticles");
  //trackCont->SetFilterHybridTracks(kTRUE);
  trackCont->SetParticlePtCut(kTrackPtCut);
  JetExtractor_R04_JF_PYTHIA_PartLevel->AdoptParticleContainer(trackCont);
  JetExtractor_R04_JF_PYTHIA_PartLevel->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);
  
  //jets from pythia at detector level                                                                                                                               
  AliEmcalJetTask *JetExtractor_R04_JF_PYTHIA_DetLevel = AliEmcalJetTask::AddTaskEmcalJet("","",AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.005, AliJetContainer::E_scheme, "JetPY", 0., kTRUE, kFALSE);
  JetExtractor_R04_JF_PYTHIA_DetLevel->SetNeedEmcalGeom(kFALSE);
  AliTrackContainer* trackCont_det = new AliTrackContainer("tracks");
  trackCont_det->SetFilterHybridTracks(kTRUE);
  trackCont_det->SetParticlePtCut(kTrackPtCut);
  JetExtractor_R04_JF_PYTHIA_DetLevel->AdoptParticleContainer(trackCont_det);
  JetExtractor_R04_JF_PYTHIA_DetLevel->GetTrackContainer(0)->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
  JetExtractor_R04_JF_PYTHIA_DetLevel->GetTrackContainer(0)->SetAODFilterBits((1<<4)|(1<<9));
    
    //background 

    //add jets for background estimation
    AliEmcalJetTask *JetExtractor_R04_Rho_JFKT = AliEmcalJetTask::AddTaskEmcalJet("","", AliJetContainer::kt_algorithm, 0.2, AliJetContainer::kChargedJet, 0.15, 0.3, 0.005, AliJetContainer::E_scheme, "Jet", 0., kTRUE, kFALSE);
    JetExtractor_R04_Rho_JFKT->SetNeedEmcalGeom(kFALSE);
    AliTrackContainer* trackCont_rho = new AliTrackContainer("tracks_toy");
    trackCont_rho->SetFilterHybridTracks(kTRUE);
    trackCont_rho->SetParticlePtCut(kTrackPtCut);
    JetExtractor_R04_Rho_JFKT->AdoptParticleContainer(trackCont_rho);
    JetExtractor_R04_Rho_JFKT->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);

    //estimate rho mass
    AliAnalysisTaskRhoMass* JetExtractor_R04_Rho_Mass =AddTaskRhoMass("Jet_KTChargedR020_tracks_toy_pT0150_E_scheme", "tracks_toy", "", "RhoR020_mass", 0.2, "TPC", 0.01,0,0,2,kTRUE,"RhoMass");
    JetExtractor_R04_Rho_Mass->SetVzRange(-10,10);
    JetExtractor_R04_Rho_Mass->SetNeedEmcalGeom(kFALSE);
    //  __R_ADDTASK__->SetExcludeLeadJets(2); // default
    JetExtractor_R04_Rho_Mass->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);

    //estimate rho                                                                                                                                              
    AliAnalysisTaskRho* JetExtractor_R04_Rho =AddTaskRhoNew("tracks_toy", "", "RhoR020", 0.2, AliJetContainer::kTPCfid, AliJetContainer::kChargedJet, kTRUE, AliJetContainer::E_scheme, "Rho_ExLJ");
    JetExtractor_R04_Rho->SetVzRange(-10,10);
    JetExtractor_R04_Rho->SetNeedEmcalGeom(kFALSE);
    JetExtractor_R04_Rho->SetExcludeLeadJets(2); // default
    JetExtractor_R04_Rho->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);

    // jet from toy 
    AliEmcalJetTask *JetExtractor_R04_JF_Toy = AliEmcalJetTask::AddTaskEmcalJet("","",AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.005, AliJetContainer::E_scheme, "Jet", 0., kTRUE, kFALSE);
    JetExtractor_R04_JF_Toy->SetNeedEmcalGeom(kFALSE);
    AliTrackContainer* trackCont_toy = new AliTrackContainer("tracks_toy");
    trackCont_toy->SetFilterHybridTracks(kTRUE);
    trackCont_toy->SetParticlePtCut(kTrackPtCut);
    JetExtractor_R04_JF_Toy->AdoptParticleContainer(trackCont_toy);
    AliEmcalJetUtilityGenSubtractor* genSub = (AliEmcalJetUtilityGenSubtractor*)JetExtractor_R04_JF_Toy->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
    genSub->SetGenericSubtractionJetMass(kTRUE);
    genSub->SetGenericSubtractionExtraJetShapes(kTRUE);
    genSub->SetUseExternalBkg(kTRUE);
    genSub->SetRhoName("RhoR020");
    genSub->SetRhomName("RhoR020_mass");
    JetExtractor_R04_JF_Toy->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);

    //make jet response
    AliJetResponseMaker* JetExtractor_R04_ResponseMaker = AddTaskJetResponseMaker("tracks_toy", "", Form("Jet_AKTChargedR%03d_tracks_toy_pT0150_E_scheme", (Int_t)(kRHJetRadius*100)), "RhoR020", kRHJetRadius, "mcparticles", "", Form("JetPartLevel_AKTChargedR%03d_mcparticles_pT0150_E_scheme",(Int_t)(kRHJetRadius*100)), "", kRHJetRadius, 1, 0.707, 0, 0, AliJetResponseMaker::kGeometrical, 0.6*kRHJetRadius, 0.6*kRHJetRadius, "TPC", -999, -999, -999, "AliJetResponseMakerMB", kFALSE, 0, -10, 10);
    //  __R_ADDTASK__->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
    JetExtractor_R04_ResponseMaker->SetNumberOfPtHardBins(21);
    Int_t arr[22] = {0,5,7,9,12,16,21,28,36,45,57,70,85,99,115,132,150,169,190,212,235,1000000};
    TArrayI bins(22, arr);
    JetExtractor_R04_ResponseMaker->SetUserPtHardBinning(bins);

    //jet extractor for true?
    AliAnalysisTaskJetExtractor* JetExtractor_R04_allJets_Toy = AddTaskJetExtractor("tracks_toy", "", "Jet_AKTChargedR040_tracks_toy_pT0150_E_scheme", "RhoR020", 0.4,"allJets");
    JetExtractor_R04_allJets_Toy->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
    JetExtractor_R04_allJets_Toy->SetVzRange(-10,10);
    JetExtractor_R04_allJets_Toy->SetNeedEmcalGeom(kFALSE);
    JetExtractor_R04_allJets_Toy->SetSaveConstituentPID(0);
    JetExtractor_R04_allJets_Toy->SetSaveConstituentsIP(0);
    JetExtractor_R04_allJets_Toy->SetSaveMCInformation(kTRUE);
    JetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(-20,10, 0.001);
    JetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(10,20, 0.01);
    JetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(20,30, 0.1);
    JetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(30,40, 0.5);
    JetExtractor_R04_allJets_Toy->GetJetTree()->AddExtractionPercentage(40,200, 1.0);
    JetExtractor_R04_allJets_Toy->ActivateJetMatching(Form("JetPartLevel_AKTChargedR%03d_mcparticles_pT0150_E_scheme", (Int_t)(kRHJetRadius*100)), "");
    JetExtractor_R04_allJets_Toy->SetHadronMatchingRadius(0.25);
    //    JetExtractor_R04_allJets_Toy->SetJetRhoMassName("");
    JetExtractor_R04_allJets_Toy->SetMCParticleArrayName("mcparticles");
    JetExtractor_R04_allJets_Toy->GetJetContainer(0)->SetRhoMassName("RhoR020_mass");
    JetExtractor_R04_allJets_Toy->SetUseBuiltinEventSelection(kRHBuiltinEventSelection);

    AliAnalysisTaskJetExtractor* JetExtractor_R04_allJets_PYTHIA_DetLevel = AddTaskJetExtractor("tracks", "", "JetPY_AKTChargedR040_tracks_pT0150_E_scheme", "", 0.4, "allJets");
    JetExtractor_R04_allJets_PYTHIA_DetLevel->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->SetVzRange(-10,10);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->SetNeedEmcalGeom(kFALSE);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->SetSaveConstituentPID(0);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->SetSaveConstituentsIP(0);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->SetSaveMCInformation(kTRUE);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->GetJetTree()->AddExtractionPercentage(0,15, 0.05);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->GetJetTree()->AddExtractionPercentage(15,30, 0.1);
    JetExtractor_R04_allJets_PYTHIA_DetLevel->GetJetTree()->AddExtractionPercentage(30,200, 0.4);*/

  TObjArray *pTopTasks = pMgr->GetTasks();
  for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
    AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
    if (!pTask) continue;
    if (pTask->InheritsFrom("AliAnalysisTaskEmcal")) {
      AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
      Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcal->GetName());
      pTaskEmcal->SetForceBeamType(iBeamType);
    }
    if (pTask->InheritsFrom("AliEmcalCorrectionTask")) {
      AliEmcalCorrectionTask * pTaskEmcalCorrection = static_cast<AliEmcalCorrectionTask*>(pTask);
      Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcalCorrection->GetName());
      pTaskEmcalCorrection->SetForceBeamType(static_cast<AliEmcalCorrectionTask::BeamType>(iBeamType));
    }
  }

  if (!pMgr->InitAnalysis()) return 0;
  pMgr->PrintStatus();
    
  pMgr->SetUseProgressBar(kTRUE, 250);
  
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
  //  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetPerfTree.C");
}

void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode)
{
  Int_t maxFilesPerWorker = 4;
  Int_t workerTTL = 7200;
  const char* runNumbers = "180720";
  const char* pattern = "pass2/AOD/*/AliAOD.root";
  const char* gridDir = "/alice/data/2012/LHC12c";
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
  plugin->SetAliPhysicsVersion("vAN-20170628-1");

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
