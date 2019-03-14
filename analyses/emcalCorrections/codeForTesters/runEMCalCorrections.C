// runEMCalCorrections.C

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
AliAnalysisManager* runEMCalCorrections(
        const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
        const char   *cLocalFiles    = "aodFiles.txt",   // set the local list file
        UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
        UInt_t        iNumEvents     = 1000,                                    // number of events to be analyzed
        const char   *cRunPeriod     = "LHC11h",                                // set the run period
        UInt_t        kPhysSel       = AliVEvent::kAnyINT,//kAnyINT,                      // physics selection - Set to 0 for MC productions!
        const char   *cTaskName      = "EMCalCorrections",                      // sets name of analysis manager
        // 0 = only prepare the analysis manager but do not start the analysis
        // 1 = prepare the analysis manager and start the analysis
        // 2 = launch a grid analysis
        Int_t         iStartAnalysis = 1,
        const char   *cGridMode      = "test"
        )
{
    TString sRunPeriod(cRunPeriod);
    sRunPeriod.ToLower();

    // Configure corrections
    const Bool_t newCorrectionFramework = kTRUE;
    const Bool_t oldCorrectionFramework = kTRUE;
    Bool_t comparisonBetweenTasks = kFALSE;
    // Will enable a comparison between the two if both are enabled
    if (newCorrectionFramework == kTRUE && oldCorrectionFramework == kTRUE) {
      comparisonBetweenTasks = kTRUE;
    }

    // Setup input object names
    TString newFrameworkTracks = "tracks";
    TString newFrameworkClusters = "caloClusters";
    TString newFrameworkCells = "emcalCells";
    if (comparisonBetweenTasks == kTRUE) {
      newFrameworkTracks += "New";
      newFrameworkClusters += "New";
      newFrameworkCells += "New";
    }
    // Change the names properly for ESDs
    if (cDataType == "ESD") {
      newFrameworkCells.Replace(0, 5, "EMCAL");
      newFrameworkClusters.Replace(0, 1, "C");
      newFrameworkTracks.Replace(0, 1, "T");
    }

    AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp;

    Bool_t bIsRun2 = kFALSE;
    if (sRunPeriod.Length() == 6 && sRunPeriod.BeginsWith("lhc15")) bIsRun2 = kTRUE;

    if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
        iBeamType = AliAnalysisTaskEmcal::kAA;
    }
    else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
            sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f") {
        iBeamType = AliAnalysisTaskEmcal::kpA;
    }

    Bool_t bSeparateEMCalDCal = bIsRun2;

    AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);
    Printf("Default track cut period set to: %s", AliTrackContainer::GetDefTrackCutsPeriod().Data());

    Bool_t   bDoTender            = kFALSE;
    Bool_t   bDoHadCorr           = kFALSE;
    const Double_t kClusPtCut           = 0.30;
    const Double_t kTrackPtCut          = 0.15;
    const Double_t kHadCorrF            = 2.;

    // Enable tender and hadCorr
    if (oldCorrectionFramework == kTRUE)
    {
      bDoTender = kTRUE;
      bDoHadCorr = kTRUE;
    }

    const AliAnalysisTaskEmcalJetSpectraQA::EHistoType_t kHistoType = AliAnalysisTaskEmcalJetSpectraQA::kTHnSparse;

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
    if (sLocalFiles == "") {
        Printf("You need to provide the list of local files!");
        return 0;
    }
    Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);

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

    // Copy collections
    if (comparisonBetweenTasks == kTRUE)
    {
      // Cells
      AliEmcalCorrectionTask::InputObject_t inputObject = AliEmcalCorrectionTask::kCaloCells;
      bool IsEsd = (iDataType == kEsd);
      TString inputObjectBranchName = AliEmcalCorrectionTask::DetermineUseDefaultName(inputObject, IsEsd);
      AliEmcalCopyCollection * copyTaskCells = AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newFrameworkCells.Data());

      // Clusters
      // We don't need to copy clusters since we are reclusterizing
      // Tracks
      inputObject = AliEmcalCorrectionTask::kTrack;
      inputObjectBranchName = AliEmcalCorrectionTask::DetermineUseDefaultName(inputObject, IsEsd);
      AliEmcalCopyCollection * copyTaskTracks = AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newFrameworkTracks.Data());
    }
  
    // EMCal corrections
    if (newCorrectionFramework == kTRUE)
    {
      AliEmcalCorrectionTask * correctionTask = AddTaskEmcalCorrectionTask();
      correctionTask->SelectCollisionCandidates(kPhysSel);
      correctionTask->SetRunPeriod(sRunPeriod.Data());
      if (bIsRun2) {
        correctionTask->SetUseNewCentralityEstimation(kTRUE);
      }
      // Local configuration
      correctionTask->SetUserConfigurationFilename("userTestConfiguration.yaml");
      // Grid configuration
      //correctionTask->SetUserConfigurationFilename("alien:///alice/cern.ch/user/r/rehlersi/AliEmcalCorrectionConfigurationTestComparison.yaml");
      correctionTask->Initialize();
    }

    if (bDoTender) {
        const char *cPass        = 0;
        Bool_t   bDistBC         = kFALSE; //switch for recalculation cluster position from bad channel
        Bool_t   bRecalibClus    = kFALSE;
        Bool_t   bRecalcClusPos  = kFALSE;
        Bool_t   bNonLinearCorr  = kFALSE;
        Bool_t   bRemExoticCell  = kFALSE;
        Bool_t   bRemExoticClus  = kFALSE;
        Bool_t   bFidRegion      = kFALSE;
        Bool_t   bCalibEnergy    = kTRUE;
        Bool_t   bCalibTime      = kTRUE;
        Bool_t   bRemBC          = kTRUE;
        UInt_t   iNonLinFunct    = AliEMCALRecoUtils::kNoCorrection;
        Bool_t   bReclusterize   = kFALSE;
        Float_t  fSeedThresh     = 0.1;      // 100 MeV
        Float_t  fCellThresh     = 0.05;     // 50 MeV
        UInt_t   iClusterizer    = AliEMCALRecParam::kClusterizerv2;
        Bool_t   bTrackMatch     = kFALSE;
        Bool_t   bUpdateCellOnly = kTRUE;
        Float_t  fEMCtimeMin     = -50e-6;
        Float_t  fEMCtimeMax     =  50e-6;
        Float_t  fEMCtimeCut     =  1e6;
        if (sRunPeriod == "lhc11h") {
            fEMCtimeMin = -50e-9;
            fEMCtimeMax = 100e-9;
        }
        fEMCtimeMin     = -1;
        fEMCtimeMax     =  1;
        fEMCtimeCut     =  1;

        AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
                bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
                fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut, cPass);
        pTenderTask->SelectCollisionCandidates(kPhysSel);

        AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", iClusterizer,
                0.05, 0.1, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut,
                kFALSE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
        pClusterizerTask->SelectCollisionCandidates(kPhysSel);
        bRemExoticClus  = kTRUE;
        iNonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrectedv3;

        AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(iNonLinFunct, bRemExoticClus, "usedefault", "", 0., kTRUE);
        pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
        pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
        pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
    }

    if (bDoHadCorr) {
        // Cluster-track matcher task
        AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher("usedefault", "usedefault", 0.1, kFALSE, kTRUE, kTRUE, kTRUE);
        pMatcherTask->SelectCollisionCandidates(kPhysSel);
        pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
        pMatcherTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
        pMatcherTask->GetClusterContainer(0)->SetClusECut(0.);
        pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);

        if (iDataType == kEsd) {
          pMatcherTask->SetDoPropagation(kTRUE);
        }

        // Hadronic correction task
        AliHadCorrTask *pHadCorrTask = AddTaskHadCorr("usedefault", "usedefault", "",
                kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
        pHadCorrTask->SelectCollisionCandidates(kPhysSel);
        pHadCorrTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
        pHadCorrTask->GetClusterContainer(0)->SetClusECut(0);
        pHadCorrTask->GetClusterContainer(0)->SetClusPtCut(0.);
    }

    // Jet finding
    const AliJetContainer::EJetAlgo_t jetAlgorithm = AliJetContainer::antikt_algorithm;
    const Double_t jetRadius = 0.2;
    const AliJetContainer::EJetType_t jetType = AliJetContainer::kFullJet;
    const Double_t minTrackPt = 0.15;
    const Double_t minClusterPt = 0.30;
    Double_t kGhostArea = 0.01;
    if (iBeamType != AliAnalysisTaskEmcal::kpp) kGhostArea = 0.005;
    const AliJetContainer::ERecoScheme_t recoScheme = AliJetContainer::pt_scheme;
    const char * label = "Jet";
    const Double_t minJetPt = 0;
    const Bool_t lockTask = kTRUE;
    const Bool_t fillGhosts = kFALSE;

    if (newCorrectionFramework)
    {
      AliEmcalJetTask *pFullJet02TaskNew = AddTaskEmcalJet(newFrameworkTracks.Data(), newFrameworkClusters.Data(),
              jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, label, minJetPt, lockTask, fillGhosts);
      pFullJet02TaskNew->SelectCollisionCandidates(kPhysSel);
      pFullJet02TaskNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

      if (newFrameworkTracks != "tracks" && newFrameworkTracks != "Tracks")
      {
        // Remove wrong particle container
        pFullJet02TaskNew->RemoveParticleContainer(0);
        // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
        AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
        newTracks->SetParticlePtCut(minTrackPt);
        pFullJet02TaskNew->AdoptTrackContainer(newTracks);
      }
    }
    if (oldCorrectionFramework)
    {
      AliEmcalJetTask *pFullJet02Task = AddTaskEmcalJet("usedefault", "usedefault",
              jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, label, minJetPt, lockTask, fillGhosts);
      pFullJet02Task->SelectCollisionCandidates(kPhysSel);
      pFullJet02Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

    //////////////////////////////////////////
    // Add you task here!
    //////////////////////////////////////////
    if (newCorrectionFramework) {
      // Your task using newFrameworkTracks, newFrameworkClusters, and newFrameworkCells
    }
    if (oldCorrectionFramework) {
      // Your task using the normal tracks and clusters
    }

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
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCopyCollection.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
}

void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode)
{
    Int_t maxFilesPerWorker = 4;
    Int_t workerTTL = 7200;
    // 000168464
    const char* runNumbers = "168464";
    //const char* runNumbers = "180720";
    // /alice/data/2011/LHC11h_2/000167693/ESDs/pass2/AOD145
    const char* gridDir = "/alice/data/2011/LHC11h_2";
    const char* pattern = "ESDs/pass2/AOD145/*/AliAOD.root";
    //const char* pattern = "pass2/AOD/*/AliAOD.root";
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
