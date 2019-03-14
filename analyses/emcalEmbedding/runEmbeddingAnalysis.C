/// \file runJetHAnalysis.C
/// \brief Jet-Hadron away-side analysis
///
/// \ingroup EMCALJETFW
/// Embedding correction macro for jet energy scale in the Jet-H Analysis
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
/// \date Jul 27, 2016


class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalSetupTask;
class AliAnalysisGrid;

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runEmbeddingAnalysis(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    const char   *cLocalFiles    = "aodFiles.txt",                          // set the local list file
    const UInt_t  iNumEvents     = 5000,                                    // number of events to be analyzed
    const UInt_t  kPhysSel       = AliVEvent::kAny,                         // physics selection
    const char   *cTaskName      = "embeddingAnalysis",                     // sets name of analysis manager
    const char   *cOCDBpath      = "uselocal", //"raw://",                                // change to "raw://" if running on the grid
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

  // Configure corrections
  const Bool_t   bDoTender            = kTRUE;
  const Bool_t   bDoHadCorr           = kTRUE;
  const Double_t kHadCorrF            = 2.;

  // Control background subtraction
  Bool_t bEnableBackgroundSubtraction = kFALSE;

  // Setup pt cuts
  // TODO: Implement these
  const Double_t kClusECut            = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kJetPtCut            = 1.;

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
  LoadMacros();

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AddAODHandler();
  }
  else {  
    AliESDInputHandler* pESDHandler = AddESDHandler();
  }

  // Physics selection task
  if (iDataType == kEsd) {
    AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();
  }

  // Centrality task
  // The Run 2 condition is too restrictive, but until the switch to MultSelection is complete, it is the best we can do
  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp && bIsRun2 == kFALSE) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
  // AliMultSelection
  // Works for both pp and PbPb for the periods that it is calibrated
  if (bIsRun2 == kTRUE) {
    AliMultSelectionTask * pMultSelectionTask = AddTaskMultSelection(kFALSE);
    pMultSelectionTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  // CDBconnect task
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  taskCDB->SetFallBackToRaw(kTRUE);

  // AODTrackMaker
  // Create AOD tracks. This is using real data
  AliEmcalAodTrackFilterTask *aodTrackMaker = AddTaskEmcalAodTrackFilter("AODFilterTracks", "tracks", sRunPeriod.Data());
  aodTrackMaker->SelectCollisionCandidates(kPhysSel);
  aodTrackMaker->SetAttemptProp(kTRUE);

  // PicoTrackMaker
  // Create pico tracks according to the AOD track filter. This is using PbPb data.
  TString tracksName = "PicoTracks";
  TString clustersName = "EmcCaloClusters";
  AliEmcalPicoTrackMaker * picoTracks = AddTaskEmcalPicoTrackMaker(tracksName, "AODFilterTracks");
  picoTracks->SelectCollisionCandidates(kPhysSel);

  // EmbeddingFromPythia
  // Extracts information from PYTHIA and injects the tracks into the pico tracks, and cells into a AliVCelss
  // TODO: Extract name from the correction task
  TString emcalCellsName = "emcalCells";
  TString mcTracksName = "MCSelectedParticles";
  TString sPYTHIAPath = "alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root";
  Int_t kNpTHardBins = 11;
  Double_t kPtHardBinsScaling[11] = {0.000000E+00, 5.206101E-05, 5.859497E-06, 4.444755E-07, 4.344664E-08, 5.154750E-09, 6.956634E-10, 1.149828E-10, 2.520137E-11, 6.222240E-12, 2.255832E-12};
  Double_t kMinPythiaJetPt = 0;
  AliJetEmbeddingFromPYTHIATask * pythiaEmbedding = AddTaskJetEmbeddingFromPYTHIA(tracksName, "",
                  emcalCellsName,                         // EMCal cells name
                  mcTracksName,                           // Name of the MC at particle level
                  sPYTHIAPath,                            // Path to PYTHIA files. Seach string
                  kNpTHardBins, kPtHardBinsScaling,       // Pt hard bin settings
                  "aodTree", "tracks", "", "emcalCells", "mcparticles",   // AOD standard branch information. The combined tracks and cells will be injected into these branches
                  "lhc12a15e",                            // Run period for MC particles to embed
                  kFALSE,                                 // False includes the ITS
                  -1, -1,                                 // Centrality settings
                  0,                                      // Trigger Mask
                  kMinPythiaJetPt,
                  kFALSE, kTRUE);                         // Copy array and make QA

  
  if (bDoTender) {
    // Only cell energy/time recalibration (and bad channel) is switched on
    Bool_t   bCalibEnergy    = kTRUE;
    Bool_t   bCalibTime      = kTRUE;
    Bool_t   bRemBC          = kTRUE;
    Bool_t   bUpdateCellOnly = kTRUE;

    // All these parameters are irrelevant for the tender
    const char *cPass        = 0;
    Bool_t   bDistBC         = kFALSE;
    Bool_t   bRecalibClus    = kFALSE;
    Bool_t   bRecalcClusPos  = kFALSE;
    Bool_t   bNonLinearCorr  = kFALSE;
    Bool_t   bRemExoticCell  = kFALSE;
    Bool_t   bRemExoticClus  = kFALSE;
    Bool_t   bFidRegion      = kFALSE;
    UInt_t   iNonLinFunct    = AliEMCALRecoUtils::kNoCorrection;
    Bool_t   bReclusterize   = kFALSE;
    Float_t  fSeedEThresh    = 0.1;      // 100 MeV
    Float_t  fCellEThresh    = 0.05;     // 50 MeV
    UInt_t   iClusterizer    = AliEMCALRecParam::kClusterizerv2;
    Bool_t   bTrackMatch     = kFALSE;
    AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
        bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedEThresh,
        fCellEThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, 0, 1e6, 1e6, cPass);
    pTenderTask->SelectCollisionCandidates(kPhysSel);

    // Configure clusterizer
    if (iBeamType != AliAnalysisTaskEmcal::kpp) {
      iClusterizer         = AliEMCALRecParam::kClusterizerv2;
    }
    else {
      iClusterizer         = AliEMCALRecParam::kClusterizerv1;
    }
    fSeedEThresh             =  0.1;      // 100 MeV
    fCellEThresh             =  0.05;     // 50 MeV
    // Time cuts are switched off at cell level. However, they are enabled for clusters
    // Nominal values
    /*
    Float_t  fEMCtimeMin     = -50e-6;
    Float_t  fEMCtimeMax     =  50e-6;
    Float_t  fEMCtimeCut     =  1e6;
    */
    // Values on the PbPb Train
    Float_t  fEMCtimeMin     = -50e-6;
    Float_t  fEMCtimeMax     =  100e-6;
    Float_t  fEMCtimeCut     =  75e6;

    AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", iClusterizer,
        fCellEThresh, fSeedEThresh, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut,
        kFALSE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
    pClusterizerTask->SelectCollisionCandidates(kPhysSel);

    // Configure cluster maker
    bRemExoticClus  = kTRUE;
    iNonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrectedv3;

    // NOTE: ClusECut is 0.15 in original embedding macro
    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(iNonLinFunct, bRemExoticClus, "usedefault", clustersName, 0., kTRUE);
    pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  }

  //////////////////////////////////////////////////////////
  // EmcalParticleMaker is now obsolete!!
  //
  // The logic is included in the cluster-track matcher now!
  //////////////////////////////////////////////////////////

  // Cluster-track matcher task
  AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher(tracksName, clustersName, 0.1, kFALSE, kTRUE, kTRUE, kTRUE);
  pMatcherTask->SelectCollisionCandidates(kPhysSel);
  pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  pMatcherTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
  pMatcherTask->GetClusterContainer(0)->SetClusECut(0.);
  pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);

  // emcalTracksName and emcalClustersName are the names of the desired objects after the cluster-track matcher
  // Track output name for cluster-track matcher
  TString emcalTracksName = "EmcalTracks_";
  emcalTracksName += tracksName;
  // Clusters output name for cluster-track matcher
  TString emcalClustersName = "EmcalClusters_";
  emcalClustersName += clustersName;

  if (iDataType == kEsd) {
    pMatcherTask->SetDoPropagation(kTRUE);
  }

  TString hadCorrClusName = "";
  if (bDoHadCorr) {
    // Configure hadCorr
    hadCorrClusName = "CaloClustersCorr";
    // Hadronic correction task
    Double_t minPt = 0.15; // Applies to clusters and tracks, but overridden below
    AliHadCorrTask *pHadCorrTask = AddTaskHadCorr(emcalTracksName, emcalClustersName, hadCorrClusName,
        kHadCorrF, minPt, 0.030, 0.015, 0, kTRUE, kTRUE);
    pHadCorrTask->SelectCollisionCandidates(kPhysSel);
    pHadCorrTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(minPt);
    pHadCorrTask->GetClusterContainer(0)->SetClusECut(0);
    pHadCorrTask->GetClusterContainer(0)->SetClusPtCut(0.);
    // Configure histograms
    // NOTE: Different than in original run macro
    pHadCorrTask->SetHistoBins(200, 0, 30);
    // Enable embedding mode
    pHadCorrTask->SetIsEmbedded(kTRUE);
  }

  // Background
  TString sRhoChName;
  TString sRhoFuName;
  // NOTE: The tracks and clusters names are probably not correctly configured here
  if (iBeamType != AliAnalysisTaskEmcal::kpp && bEnableBackgroundSubtraction == kTRUE) {
    Printf("Running background subtraction tasks!");
    sRhoChName = "Rho";
    sRhoFuName = "Rho_Scaled";

    AliEmcalJetTask *pKtChJetTask = AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
    pKtChJetTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChName, 0.4);
    pRhoTask->SetExcludeLeadJets(2);
    pRhoTask->SelectCollisionCandidates(kPhysSel);

    TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
    TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
    pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
  }

  // Full jet finder
  // TODO: Should tracksName -> emcalTracksName???
  AliEmcalJetTask *pFullJet02Task = AddTaskEmcalJet(tracksName, hadCorrClusName, AliJetContainer::antikt_algorithm, 0.2, AliJetContainer::kFullJet, 3, 3, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kTRUE, kFALSE);
  pFullJet02Task->SelectCollisionCandidates(kPhysSel);
  pFullJet02Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

  // MC particle level jet finder
  // We want full jets here despite not passing clusters, since the particle level knowns where it is charged or neutral.
  // If we only took charged jets, then all particles with a neutral charge would be thrown out.
  AliEmcalJetTask *pMCJet02Task = AddTaskEmcalJet(mcTracksName, "", AliJetContainer::antikt_algorithm, 0.2, AliJetContainer::kFullJet, 3, 3, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kTRUE, kFALSE);

  // Response matrix
  // Configure the response matrix
  // TODO: Should tracksName -> emcalTracksName??? No, it should be fine, since everything in the later cluster corrections is stored in the cluster
  // TODO: Attempt to generate jet name
  //TString detLevelJetName = AliJetContainer::GenerateJetName();
  //TString partLevelJetName = AliJetContainer::GenerateJetName();
  // Needed to change ET -> E
  TString detLevelJetName(Form("Jet_AKTFullR020_%s_pT3000_%s_E3000_pt_scheme", tracksName.Data(), hadCorrClusName.Data()));
  TString partLevelJetName(Form("Jet_AKTFullR020_%s_pT3000_pt_scheme",mcTracksName.Data()));
  // GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag);
  // NOTE: These are from the original embedding macro. These could and probably should be changed
  Double_t kJetAreaCut = 0;
  Double_t kJetLeadingTrackBias = 0;
  Int_t kLeadHadType = 1; // 0 = charged, 1 = neutral, 2 = both
  Double_t kMaxGeoDistance02 = 1.2; // 02 corresponds to R=0.2 jets
  Int_t kNcent = 1;
  Int_t kHistoType = 1; // 1 = THnSparse, 0 = TH1/TH2/TH3
  AliJetResponseMaker * jetResponseMatrix = AddTaskJetRespPtHard(tracksName, hadCorrClusName, detLevelJetName, "", 0.2, mcTracksName, "", partLevelJetName, "", 0.2, kJetPtCut, kJetAreaCut, kJetLeadingTrackBias, kLeadHadType, AliJetResponseMaker::kGeometrical, kMaxGeoDistance02, kMaxGeoDistance02, "EMCAL", // Cut type
            -999,-999,kNcent);
  jetResponseMatrix->SelectCollisionCandidates(kPhysSel);
  // Related to handling MC handling in AliAnalysisTaskEmcal. Not relevant here.
  jetResponseMatrix->SetIsPythia(kFALSE);
  jetResponseMatrix->SetIsEmbedded(kTRUE);
  jetResponseMatrix->SetHistoType(kHistoType);
  jetResponseMatrix->SetDebugLevel(3);

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
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromPYTHIA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetRespPtHard.C");
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
