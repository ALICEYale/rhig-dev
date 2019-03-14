/// \file runJetHAnalysis.C
/// \brief Embedding run macro
///
/// \ingroup EMCALJETFW
/// Embedding run macro to create a response matrix
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
/// \date Jul 27, 2016


class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliAnalysisGrid;

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runNewEmbeddingAnalysis(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    const char   *cLocalFiles    = "aodFiles.txt",                          // set the local list file
    const UInt_t  iNumEvents     = 1000,                                    // number of events to be analyzed
    const UInt_t  kPhysSel       = AliVEvent::kEMCEGA | AliVEvent::kMB |
                    AliVEvent::kCentral | AliVEvent::kSemiCentral, //AliVEvent::kAny,                         // physics selection
    const char   *cTaskName      = "EMCalEmbeddingAnalysis",                     // sets name of analysis manager
    // 0 = only prepare the analysis manager but do not start the analysis
    // 1 = prepare the analysis manager and start the analysis
    // 2 = launch a grid analysis
    Int_t         iStartAnalysis = 1,
    const UInt_t  iNumFiles      = 5,                                     // number of files analyzed locally
    const char   *cGridMode      = "test"
)
{
  // Setup period
  TString sRunPeriod(cRunPeriod);
  sRunPeriod.ToLower();

  // Set Run 2
  Bool_t bIsRun2 = kFALSE;
  if (sRunPeriod.Length() == 6 && (sRunPeriod.BeginsWith("lhc15") || sRunPeriod.BeginsWith("lhc16"))) bIsRun2 = kTRUE;

  // Set beam type
  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp;
  if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
    iBeamType = AliAnalysisTaskEmcal::kAA;
  }
  else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f") {
    iBeamType = AliAnalysisTaskEmcal::kpA;
  }

  // Ghost area
  Double_t kGhostArea = 0.01;
  if (iBeamType != AliAnalysisTaskEmcal::kpp) kGhostArea = 0.005;

  // Setup track container
  AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);
  Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());

  // Set data file type
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

  // Load macros needed for the analysis
  #ifndef __CLING__
  LoadMacros();
  #endif

  ////////////////////
  // Configure options
  ////////////////////
  // Select which embedding framework to use
  const bool newFramework = true;
  const bool oldFramework = false;
  const bool fullJets = true;

  // Configure corrections
  const Bool_t newCorrectionFramework = kTRUE;
  const Bool_t oldCorrectionFramework = kFALSE;
  Bool_t comparisonBetweenTasks = kFALSE;

  // Will enable a comparison between the two if both are enabled
  if (newCorrectionFramework == kTRUE && oldCorrectionFramework == kTRUE) {
    comparisonBetweenTasks = kTRUE;
  }

  TString newFrameworkTracks = "tracks";
  TString newFrameworkClusters = "caloClusters";
  TString newFrameworkCells = "emcalCells";
  if (comparisonBetweenTasks == kTRUE) {
    newFrameworkTracks += "New";
    newFrameworkClusters += "New";
    newFrameworkCells += "New";
  }

  // Control background subtraction
  Bool_t bEnableBackgroundSubtraction = kFALSE;

  // Setup pt cuts
  // TODO: Implement these throughout
  const Double_t kClusECut            = 0.30;
  const Double_t kClusPtCut           = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kJetPtCut            = 1.;
  const Double_t kJetAreaCut          = 0;

  // Define relevant variables
  const bool IsEsd = (iDataType == kEsd);
  // Note what tag we are using
  const TString tag = "new";
  TString truthParticlesName = "";
  TString emcalCellsName = "emcalCells";
  TString clustersName = "";
  TString clustersNameCombined = "";
  TString tracksName = "";

  ///////////////////////////////
  // Setup and Configure Analysis
  ///////////////////////////////

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
    #ifdef __CLING__
    std::stringstream physSel;
    physSel << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/OADB/macros/AddTaskPhysicsSelection.C()";
    AliPhysicsSelectionTask * pPhysSelTask = reinterpret_cast<AliPhysicsSelectionTask *>(gROOT->ProcessLine(physSel.str().c_str()));
    #else
    AliPhysicsSelectionTask * pPhysSelTask = AddTaskPhysicsSelection();
    #endif
  }

  // Centrality task
  // The Run 2 condition is too restrictive, but until the switch to MultSelection is complete, it is the best we can do
  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp && bIsRun2 == kFALSE) {
    #ifdef __CLING__
    std::stringstream centralityTask;
    centralityTask << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/OADB/macros/AddTaskCentrality.C()";
    AliCentralitySelectionTask * pCentralityTask = reinterpret_cast<AliCentralitySelectionTask *>(gROOT->ProcessLine(centralityTask.str().c_str()));
    #else
    AliCentralitySelectionTask * pCentralityTask = AddTaskCentrality(kTRUE);
    #endif
    pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
  // AliMultSelection
  // Works for both pp and PbPb for the periods that it is calibrated
  if (bIsRun2 == kTRUE) {
    #ifdef __CLING__
    std::stringstream multSelection;
    multSelection << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C(kFALSE)";
    AliMultSelectionTask * pMultSelectionTask = reinterpret_cast<AliMultSelectionTask *>(gROOT->ProcessLine(multSelection.str().c_str()));
    #else
    AliMultSelectionTask * pMultSelectionTask = AddTaskMultSelection(kFALSE);
    #endif
    pMultSelectionTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  // CDBconnect task
  #ifdef __CLING__
  std::stringstream cdbConnect;
  cdbConnect << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWGPP/PilotTrain/AddTaskCDBconnect.C()";
  AliTaskCDBconnect * taskCDB = reinterpret_cast<AliTaskCDBconnect *>(gROOT->ProcessLine(cdbConnect.str().c_str()));
  #else
  AliTaskCDBconnect * taskCDB = AddTaskCDBconnect();
  #endif
  taskCDB->SetFallBackToRaw(kTRUE);

  if (fullJets) {
    emcalCellsName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCaloCells, IsEsd);
    clustersName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster, IsEsd);
    clustersNameCombined = clustersName;
  }

  if (newFramework) {
    // Setup relevant values
    truthParticlesName = "mcparticles";
    tracksName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack, IsEsd);
    if (fullJets) {
      clustersNameCombined += "Combined";
    }

    // Debugging
    //AliLog::SetClassDebugLevel("AliAnalysisTaskEmcalEmbeddingHelper", AliLog::kDebug+0);

    // Setup embedding task
    AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
    embeddingHelper->SetFileListFilename("aodFilesEmbed.txt");
    embeddingHelper->SetTriggerMask(AliVEvent::kAny); // Equivalent to 0
    embeddingHelper->SelectCollisionCandidates(kPhysSel);
    //embeddingHelper->SetRandomAccess(kTRUE);
    //embeddingHelper->SetMaxVertexDistance(3);
    //embeddingHelper->SetZVertexCut(3);
    // Runs successfully starting with 15, but not if starting at 1...
    //embeddingHelper->SetStartingFileIndex(15);
    embeddingHelper->Initialize();
  }

  // EMCal corrections
  if (newCorrectionFramework == kTRUE && newFramework == true)
  {
    TObjArray correctionTasks;

    // Debugging
    // NOTE: kDebug = 4. Can go up or down from there!
    //AliLog::SetClassDebugLevel("AliEmcalCorrectionTask", AliLog::kDebug-2);

    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("data"));
    correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("embed"));
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
        tempCorrectionTask->SetUseNewCentralityEstimator(kTRUE);
      }

      // The configuration file for all three is the same! They take advantage of component specialization
      // Local configuration
      tempCorrectionTask->SetUserConfigurationFilename("userConfiguration.yaml");
      // Grid configuration
      //tempCorrectionTask->SetUserConfigurationFilename("alien:///alice/cern.ch/user/r/rehlersi/embedding/userConfiguration.yaml");
      tempCorrectionTask->Initialize();    
    }
  }

  // Background
  TString sRhoChName;
  TString sRhoFuName;
  if (iBeamType != AliAnalysisTaskEmcal::kpp && bEnableBackgroundSubtraction == kTRUE) {
    sRhoChName = "Rho";
    sRhoFuName = "Rho_Scaled";

    AliEmcalJetTask *pKtChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
    pKtChJetTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChName, 0.4);
    pRhoTask->SetExcludeLeadJets(2);
    pRhoTask->SelectCollisionCandidates(kPhysSel);

    TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
    TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
    pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
  }

  // Jet finding
  const AliJetContainer::EJetAlgo_t jetAlgorithm = AliJetContainer::antikt_algorithm;
  const Double_t jetRadius = 0.2;
  const AliJetContainer::EJetType_t jetType = AliJetContainer::kFullJet;
  const Double_t minTrackPt = 0.15;//3.0;
  const Double_t minClusterPt = 0.30;//3.0;
  const AliJetContainer::ERecoScheme_t recoScheme = AliJetContainer::pt_scheme;
  const char * label = "Jet";
  const Double_t minJetPt = 1; // TEMP: Consider at 1! Usually: 0
  const Bool_t lockTask = kTRUE;
  const Bool_t fillGhosts = kFALSE;

  // Do not pass clusters if we are only looking at charged jets
  if (jetType == AliJetContainer::kChargedJet) {
    newFrameworkClusters = "";
  }

  AliEmcalJetTask * pFullJetTaskTruthLevelNew = 0;
  AliEmcalJetTask * pFullJetTaskDetLevelNew = 0;
  AliEmcalJetTask * pFullJetTaskNew = 0;
  AliEmcalJetTask * pFullJetTaskHybridNew = 0;
  if (newCorrectionFramework)
  {
    ///////
    // Truth level PYTHIA jet finding
    ///////
    pFullJetTaskTruthLevelNew = AliEmcalJetTask::AddTaskEmcalJet("mcparticles", "",
            jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "truthLevel", minJetPt, lockTask, fillGhosts);
    pFullJetTaskTruthLevelNew->SelectCollisionCandidates(kPhysSel);

    ///////
    // External event (called embedding) settings for truth level PYTHIA jet finding
    ///////
    // Setup the tracks properly to be retrieved from the external event
    AliMCParticleContainer * truthTracks = pFullJetTaskTruthLevelNew->GetMCParticleContainer(0);
    // Called Embedded, but really just means get from an external event!
    truthTracks->SetIsEmbedding(kTRUE);
    
    ///////
    // Detector level PYTHIA jet finding
    ///////
    /*if (comparisonBetweenTasks) {
      Printf("WARNING: More care is needed in configuring pFullJetTaskDetLevelNew whe comparing tasks! In particular, the objects are not being copied, so corrections are being applied twice. This needs to be fixed!");
      std::exit(1);
    }
    pFullJetTaskDetLevelNew = AliEmcalJetTask::AddTaskEmcalJet(newFrameworkTracks.Data(), newFrameworkClusters.Data(),
            jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "detLevel", minJetPt, lockTask, fillGhosts);
    pFullJetTaskDetLevelNew->SelectCollisionCandidates(kPhysSel);
    pFullJetTaskDetLevelNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

    ///////
    // External event (called embedding) settings for det level PYTHIA jet finding
    ///////
    // Tracks
    // Uses the name of the container passed into AliEmcalJetTask
    AliTrackContainer * tracksDetLevel = pFullJetTaskDetLevelNew->GetTrackContainer(0);
    // Get the det level tracks from the external event!
    tracksDetLevel->SetIsEmbedding(kTRUE);
    // Clusters
    if (jetType != AliJetContainer::kChargedJet) {
      // Uses the name of the container passed into AliEmcalJetTask
      AliClusterContainer * clustersDetLevel = pFullJetTaskDetLevelNew->GetClusterContainer(0);
      // Get the det level clusters from the external event!
      clustersDetLevel->SetIsEmbedding(kTRUE);
    }*/

    ///////
    // PbPb jet finding
    ///////
    /*pFullJetTaskNew = AliEmcalJetTask::AddTaskEmcalJet(newFrameworkTracks.Data(), newFrameworkClusters.Data(),
            jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "PbPbJets", minJetPt, lockTask, fillGhosts);
    pFullJetTaskNew->SelectCollisionCandidates(kPhysSel);
    if (jetType != AliJetContainer::kChargedJet) {
      pFullJetTaskNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }
    // Handles proper track container when comparing for new framework.
    if (comparisonBetweenTasks)
    {
      // Remove wrong particle container
      pFullJetTaskNew->RemoveParticleContainer(0);
      // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
      AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
      newTracks->SetParticlePtCut(minTrackPt);
      pFullJetTaskNew->AdoptTrackContainer(newTracks);
    }*/

    ///////
    // PbPb + Detector level PYTHIA jet finding
    ///////
    // Sets up PbPb tracks and clusters
    // NOTE: The clusters name is different here since we output to a different branch!
    pFullJetTaskHybridNew = AliEmcalJetTask::AddTaskEmcalJet(newFrameworkTracks.Data(), "caloClustersCombined",
    //pFullJetTaskHybridNew = AliEmcalJetTask::AddTaskEmcalJet(newFrameworkTracks.Data(), newFrameworkClusters.Data(),
            jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, "hybridJets", minJetPt, lockTask, fillGhosts);
    pFullJetTaskHybridNew->SelectCollisionCandidates(kPhysSel);
    if (jetType != AliJetContainer::kChargedJet) {
      pFullJetTaskHybridNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

    // Handles proper track container when comparing for new framework.
    if (comparisonBetweenTasks)
    {
      // Remove wrong particle container
      pFullJetTaskHybridNew->RemoveParticleContainer(0);
      // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
      AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
      newTracks->SetParticlePtCut(minTrackPt);
      pFullJetTaskHybridNew->AdoptTrackContainer(newTracks);
    }

    ///////
    // External event (ie embedding) settings for PbPb jet finding (adds detector level PYTHIA)
    ///////
    // NOTE: These will break when comparing frameworks!
    // Add embedded tracks and clusters to jet finder
    // Tracks
    tracksDetLevel = new AliTrackContainer(newFrameworkTracks.Data());
    // Get the det level tracks from the external event!
    tracksDetLevel->SetIsEmbedding(kTRUE);
    tracksDetLevel->SetParticlePtCut(minTrackPt);
    pFullJetTaskHybridNew->AdoptTrackContainer(tracksDetLevel);
    // Clusters
    /*if (jetType != AliJetContainer::kChargedJet) {
      AliClusterContainer * clustersDetLevel = new AliClusterContainer(newFrameworkClusters.Data());
      // Get the det level clusters from the external event!
      clustersDetLevel->SetIsEmbedding(kTRUE);
      clustersDetLevel->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      pFullJetTaskHybridNew->AdoptClusterContainer(clustersDetLevel);
    }*/
  }
  /*if (oldCorrectionFramework)
  {
    AliEmcalJetTask *pFullJet02Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault",
            jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, label, minJetPt, lockTask, fillGhosts);
    pFullJet02Task->SelectCollisionCandidates(kPhysSel);
    pFullJet02Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }*/

  // TEMP
  //pMgr->AddClassDebug("AliEmcalJetTask", AliLog::kDebug+2);
  //pMgr->AddClassDebug("AliAnalysisTaskEmcalEmbeddingHelper", AliLog::kDebug+2);
  
  //////////////////////////////////////////
  // Jet Tasks
  //////////////////////////////////////////

  // Use Sample task to access how the embedding did
  AliEmcalJet::JetAcceptanceType acceptanceType = AliEmcalJet::kEMCALfid;
  TString jetTag = "Jet";
  AliAnalysisTaskEmcalJetSample * sampleTaskTruthLevelNew = 0;
  AliAnalysisTaskEmcalJetSample * sampleTaskDetLevelNew = 0;
  AliAnalysisTaskEmcalJetSample * sampleTaskNew = 0;
  AliAnalysisTaskEmcalJetSample * sampleTaskHybridNew = 0;
  if (newCorrectionFramework) {
    ///////
    // Truth level PYTHIA sample task
    ///////
    sampleTaskTruthLevelNew = AddTaskEmcalJetSample("mcparticles", "", "", "truthLevel");
    // Not necessary, since there is explicit support for "mcparticles"
    /*if (comparisonBetweenTasks)
    {
      // Remove wrong particle container
      sampleTaskTruthLevelNew->RemoveParticleContainer(0);
      // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
      AliMCParticleContainer * newMCParticles = new AliMCParticleContainer(newFrameworkTracks.Data());
      sampleTaskTruthLevelNew->AdoptMCParticleContainer(newMCParticles);
    }*/

    // Set embedding
    AliParticleContainer * partCont = sampleTaskTruthLevelNew->GetParticleContainer(0);
    partCont->SetIsEmbedding(kTRUE);
    // Set name to ensure no clashes
    //partCont->SetName(partCont->GetName() + "_Emb");

    partCont->SetParticlePtCut(0.15);
    sampleTaskTruthLevelNew->SetHistoBins(600, 0, 300);
    sampleTaskTruthLevelNew->SelectCollisionCandidates(kPhysSel);

    AliJetContainer* jetCont02 = sampleTaskTruthLevelNew->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, "truthLevel");

    ///////
    // Detector level PYTHIA sample task
    ///////
    /*sampleTaskDetLevelNew = AddTaskEmcalJetSample(newFrameworkTracks.Data(), newFrameworkClusters.Data(), newFrameworkCells.Data(), "detLevel");
    if (comparisonBetweenTasks)
    {
      // Remove wrong particle container
      sampleTaskDetLevelNew->RemoveParticleContainer(0);
      // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
      AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
      sampleTaskDetLevelNew->AdoptTrackContainer(newTracks);
    }

    // Tracks
    // Set embedding
    partCont = sampleTaskDetLevelNew->GetParticleContainer(0);
    partCont->SetIsEmbedding(kTRUE);
    // Set name to ensure no clashes
    //partCont->SetName(partCont->GetName() + "_Emb");
    // Settings
    partCont->SetParticlePtCut(0.15);
    if (jetType != AliJetContainer::kChargedJet) {
      // Clusters
      // Set embedding
      sampleTaskDetLevelNew->GetClusterContainer(0)->SetIsEmbedding(kTRUE);
      // Set name to ensure no clashes
      //sampleTaskDetLevelNew->GetClusterContainer(0)->SetName(sampleTaskDetLevelNew->GetClusterContainer(0)->GetName() + "_Emb");
      // Settings
      sampleTaskDetLevelNew->GetClusterContainer(0)->SetClusECut(0.);
      sampleTaskDetLevelNew->GetClusterContainer(0)->SetClusPtCut(0.);
      sampleTaskDetLevelNew->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
      sampleTaskDetLevelNew->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
      sampleTaskDetLevelNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

    sampleTaskDetLevelNew->SetHistoBins(600, 0, 300);
    sampleTaskDetLevelNew->SelectCollisionCandidates(kPhysSel);

    AliJetContainer* jetCont02 = sampleTaskDetLevelNew->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, "detLevel");*/

    ///////
    // PbPb sample task
    ///////
    /*sampleTaskNew = AddTaskEmcalJetSample(newFrameworkTracks.Data(), newFrameworkClusters.Data(), newFrameworkCells.Data(), "PbPbJets");
    if (comparisonBetweenTasks)
    {
      // Remove wrong particle container
      sampleTaskNew->RemoveParticleContainer(0);
      // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
      AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
      sampleTaskNew->AdoptTrackContainer(newTracks);
    }

    if (jetType != AliJetContainer::kChargedJet) {
      sampleTaskNew->GetClusterContainer(0)->SetClusECut(0.);
      sampleTaskNew->GetClusterContainer(0)->SetClusPtCut(0.);
      sampleTaskNew->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
      sampleTaskNew->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
      sampleTaskNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }
    sampleTaskNew->GetParticleContainer(0)->SetParticlePtCut(0.15);
    sampleTaskNew->SetHistoBins(600, 0, 300);
    sampleTaskNew->SelectCollisionCandidates(kPhysSel);

    AliJetContainer* jetCont02 = sampleTaskNew->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, "PbPbJets");*/

    ///////
    // PbPb + Detector level PYTHIA sample task
    ///////
    // NOTE: The clusters name is different here since we output to a different branch!
    sampleTaskHybridNew = AddTaskEmcalJetSample(newFrameworkTracks.Data(), "caloClustersCombined", "emcalCellsCombined", "hybridJets");
    //sampleTaskHybridNew = AddTaskEmcalJetSample(newFrameworkTracks.Data(), newFrameworkClusters.Data(), "emcalCellsCombined", "hybridJets");
    if (comparisonBetweenTasks)
    {
      // Remove wrong particle container
      sampleTaskHybridNew->RemoveParticleContainer(0);
      // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
      AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
      sampleTaskHybridNew->AdoptTrackContainer(newTracks);
    }

    // PbPb tracks settings
    sampleTaskHybridNew->GetParticleContainer(0)->SetParticlePtCut(0.15);

    // Embed tracks
    tracksDetLevel = new AliTrackContainer(newFrameworkTracks.Data());
    // Get the det level tracks from the external event!
    tracksDetLevel->SetIsEmbedding(kTRUE);
    // Set name to ensure no clashes
    TString previousName = tracksDetLevel->GetName();
    previousName += "_Emb";
    tracksDetLevel->SetName(previousName.Data());
    // Settings
    tracksDetLevel->SetParticlePtCut(0.15);
    sampleTaskHybridNew->AdoptTrackContainer(tracksDetLevel);

    if (jetType != AliJetContainer::kChargedJet) {
      // PbPb clusters settings
      sampleTaskHybridNew->GetClusterContainer(0)->SetClusECut(0.);
      sampleTaskHybridNew->GetClusterContainer(0)->SetClusPtCut(0.);
      sampleTaskHybridNew->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
      sampleTaskHybridNew->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
      sampleTaskHybridNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

      // Embed clusters
      /*AliClusterContainer * clustersDetLevel = new AliClusterContainer(newFrameworkClusters.Data());
      // Get the det level clusters from the external event!
      clustersDetLevel->SetIsEmbedding(kTRUE);
      // Set name to ensure no clashes
      previousName = clustersDetLevel->GetName();
      previousName += "_Emb";
      clustersDetLevel->SetName(previousName.Data());
      // Settings
      clustersDetLevel->SetClusECut(0.);
      clustersDetLevel->SetClusPtCut(0.);
      clustersDetLevel->SetClusNonLinCorrEnergyCut(0.);
      clustersDetLevel->SetClusHadCorrEnergyCut(0.30);
      clustersDetLevel->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      sampleTaskHybridNew->AdoptClusterContainer(clustersDetLevel);*/
    }

    sampleTaskHybridNew->SetHistoBins(600, 0, 300);
    sampleTaskHybridNew->SelectCollisionCandidates(kPhysSel);

    AliJetContainer* jetCont02 = sampleTaskHybridNew->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, "hybridJets");
  }

  // Jet Reponse Maker
  /*AliTrackContainer * trackCont = new AliTrackContainer(newFrameworkTracks.Data());
  //AliClusterContainer * clusCont = new AliClusterContainer(newFrameworkClusters.Data());
  AliClusterContainer * clusCont = new AliClusterContainer("caloClustersCombined");
  clusCont->SetClusECut(minClusterPt);
  TString detLevelJetName = AliJetContainer::GenerateJetName(jetType, jetAlgorithm, recoScheme, jetRadius, trackCont, clusCont, "hybridJets");*/
  TString detLevelJetName = pFullJetTaskHybridNew->GetName();
  //std::cout << "detLevelJetName:  " << detLevelJetName << std::endl << "detLevelJetName2: " << detLevelJetName2 << std::endl;
  TString mcTracksName = "mcparticles";
  /*AliMCParticleContainer * mcPartCont = new AliMCParticleContainer(mcTracksName.Data());
  TString partLevelJetName = AliJetContainer::GenerateJetName(jetType, jetAlgorithm, recoScheme, jetRadius, mcPartCont, 0, "truthLevel");*/
  TString partLevelJetName = pFullJetTaskTruthLevelNew->GetName();
  //std::cout << "partLevelJetName:  " << partLevelJetName << std::endl << "partLevelJetName2: " << partLevelJetName2 << std::endl;
  AliJetResponseMaker::MatchingType matchingType = AliJetResponseMaker::kGeometrical;
  // Since this applies equally to both Jet Containers and we are comparing against MC truth,
  // it is probably best to disable this and just use any bias applied at the jet finder level!
  Double_t kJetLeadingTrackBias = .3;      // Previously was 0! Should consider what we need here, but perhaps not huge since the jets are already biased.
  Int_t kLeadHadType = 1;                 // 0 = charged, 1 = neutral, 2 = both
  Double_t kMaxGeoDistance02 = 1.2;       // 02 corresponds to R=0.2 jets
  Int_t kNcent = 1;                       // TODO: Check this value!
  Int_t kHistoType = 1;                   // 1 = THnSparse, 0 = TH1/TH2/TH3
  //AliJetResponseMaker * jetResponseMatrix = AddTaskJetRespPtHard(newFrameworkTracks, newFrameworkClusters, detLevelJetName, "", 0.2,
  AliJetResponseMaker * jetResponseMatrix = AddTaskJetResponseMaker(newFrameworkTracks, "caloClustersCombined", detLevelJetName, "", 0.2,
                                  mcTracksName, "", partLevelJetName, "", 0.2,
                                  kJetPtCut, kJetAreaCut, kJetLeadingTrackBias, kLeadHadType,
                                  matchingType, kMaxGeoDistance02, kMaxGeoDistance02, "EMCAL", // Cut type
                                  -999,-999,-999);
                                  //-999,-999,kNcent);
  // Set mc particles to be from the embedded event!
  //for (Int_t i = 0; i < 3; i++)
  //{
    //pMgr->AddClassDebug("AliEmcalJetTask", AliLog::kDebug+2);
    //pMgr->AddClassDebug("AliJetContainer", AliLog::kDebug+7);
    //pMgr->AddClassDebug("AliJetResponseMaker", AliLog::kDebug+2);
    jetResponseMatrix->GetJetContainer(1)->GetParticleContainer()->SetIsEmbedding(kTRUE);
    jetResponseMatrix->GetJetContainer(1)->SetJetAcceptanceType(acceptanceType);
    jetResponseMatrix->GetJetContainer(0)->SetJetAcceptanceType(acceptanceType);
    jetResponseMatrix->SelectCollisionCandidates(kPhysSel);
    // TODO: Set embedding on containers
    jetResponseMatrix->SetHistoType(kHistoType);
    jetResponseMatrix->SetJetRelativeEPAngleAxis(kTRUE);
  //}

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
      pTaskEmcalCorrection->SetForceBeamType(iBeamType);
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
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
    }
    else {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
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
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetRespPtHard.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
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
  plugin->SetAliPhysicsVersion("vAN-20160203-1");

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
