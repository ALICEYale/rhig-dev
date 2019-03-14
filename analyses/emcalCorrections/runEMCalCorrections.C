// runEMCalCorrections.C

class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisGrid;
class AliAnalysisManager;
class AliAnalysisAlien;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliTaskCDBconnect;

class AliClusterContainer;
class AliParticleContainer;
class AliJetContainer;

// Can either do the forward declaration, and then specify throughout
//namespace PWG { namespace EMCAL { class AliEmcalCopyCollection; } }
//namespace PWG { namespace EMCAL { class AliEmcalCorrectionTask; } }
// Or import the classes into the global namespace
//using PWG::EMCAL::AliEmcalCopyCollection;
//using PWG::EMCAL::AliEmcalCorrectionTask;
class AliEmcalCopyCollection;
class AliEmcalCorrectionTask;
class AliEmcalCorrectionComponent;
class AliEmcalJetTask;
class AliAnalysisTaskEmcalJetSample;

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
        /*const Bool_t  bDoChargedJets = kTRUE,
        const Bool_t  bDoFullJets    = kTRUE,
        const Bool_t  bDoMCJets      = kFALSE,*/
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
    const bool testConfigureForLegoTrain = false;
    const Bool_t newCorrectionFramework = kTRUE;
    Bool_t oldCorrectionFramework = kTRUE;
    const Bool_t createReferenceFile = kFALSE;
    if (createReferenceFile) {
      oldCorrectionFramework = kFALSE;
    }
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
    if (!strcmp(cDataType, "ESD")) {
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

    //const AliAnalysisTaskEmcalJetSpectraQA::EHistoType_t kHistoType = AliAnalysisTaskEmcalJetSpectraQA::kTHnSparse;

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

    #ifndef __CLING__
    LoadMacros();
    #endif

    // Analysis manager
    AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

    if (iDataType == kAod) {
        AliAODInputHandler* pAODHandler = AliAnalysisTaskEmcal::AddAODHandler();
    }
    else {  
        AliESDInputHandler* pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
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

    // Copy collections
    // Use newFrameworkTracks as a proxy for if we need to copy the collections
    if (newFrameworkTracks.Contains("New"))
    {
      // Cells
      AliEmcalContainerUtils::InputObject_t inputObject = AliEmcalContainerUtils::kCaloCells;
      bool IsEsd = (iDataType == kEsd);
      TString inputObjectBranchName = AliEmcalContainerUtils::DetermineUseDefaultName(inputObject, IsEsd);
      AliEmcalCopyCollection * copyTaskCells = AliEmcalCopyCollection::AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newFrameworkCells.Data());

      // Clusters
      // We don't need to copy clusters since we are reclusterizing and the new collection will be made automiatcally
      // Tracks
      inputObject = AliEmcalContainerUtils::kTrack;
      inputObjectBranchName = AliEmcalContainerUtils::DetermineUseDefaultName(inputObject, IsEsd);
      AliEmcalCopyCollection * copyTaskTracks = AliEmcalCopyCollection::AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newFrameworkTracks.Data());
    }
  
    // EMCal corrections
    if (newCorrectionFramework == kTRUE)
    {
      // Uncomment for debugging!
      // Mainly used for debugging the main task
      // Must be set before Initialization so that we can see what is occuring
      //AliLog::SetClassDebugLevel("AliYAMLConfiguration", AliLog::kDebug+2);
      //AliLog::SetClassDebugLevel("AliEmcalCorrectionTask", AliLog::kDebug+2);
      //pMgr->AddClassDebug("AliEmcalCorrectionTask", AliLog::kDebug+2);
      // Mainly used for debugging loading values from yaml
      //AliLog::SetClassDebugLevel("AliEmcalCorrectionComponent", AliLog::kDebug+0);
      //pMgr->AddClassDebug("AliEmcalCorrectionComponent", AliLog::kDebug+0);

      AliEmcalCorrectionTask * correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
      correctionTask->SelectCollisionCandidates(kPhysSel);
      if (bIsRun2) {
        correctionTask->SetUseNewCentralityEstimation(kTRUE);
        correctionTask->SetNCentBins(5);
      }
      if (testConfigureForLegoTrain) {
        std::cout << "\tRun Macro: Using configure EMCal Correction Task for lego train. The following print out should have a dummy task.\n";
        correctionTask = AliEmcalCorrectionTask::ConfigureEmcalCorrectionTaskOnLEGOTrain("mySuffix");
        TObjArray * tasks = pMgr->GetTasks();
        std::cout << "Tasks in Analysis Manager:\n";
        for (auto task : *tasks) {
          std::cout <<  "\t-" << task->GetName() << "\n";
        }
      }
      // Use newFrameworkTracks as a proxy for if we need to apply the config
      if (newFrameworkTracks.Contains("New") || createReferenceFile)
      {
        // Local configuration
        correctionTask->SetUserConfigurationFilename("AliEmcalCorrectionConfigurationTestComparison.yaml");
        // Grid configuration
        //correctionTask->SetUserConfigurationFilename("alien:///alice/cern.ch/user/r/rehlersi/AliEmcalCorrectionConfigurationTestComparison.yaml");
      }
      correctionTask->Initialize(testConfigureForLegoTrain);
      if (testConfigureForLegoTrain) {
        std::cout << "\tRun Macro: The dummy task should now be removed!\n";
        TObjArray * tasks = pMgr->GetTasks();
        std::cout << "Tasks in Analysis Manager:\n";
        for (auto task : *tasks) {
          std::cout <<  "\t-" << task->GetName() << "\n";
        }
      }

      // Test retrieving a particular correction component
      AliEmcalCorrectionComponent * badChannel = correctionTask->GetCorrectionComponent("AliEmcalCorrectionCellBadChannel");
      if (badChannel != nullptr) {
        std::cout << "Got back bad channel component!\nComponent name: " << badChannel->GetName() << "\n";
        AliEMCALRecoUtils * recoUtils = badChannel->GetRecoUtils();
        if (recoUtils != nullptr) {
          std::cout << "Got back reco utils!\n";
        }
      }
      else {
        std::cout << "Bad channel component was not retrieved successfully!\n";
        std::exit(0);
      }
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

        #ifdef __CLING__
        std::stringstream emcalTender;
        emcalTender << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/AddTaskEMCALTender.C(";
        emcalTender << bDistBC << ", ";
        emcalTender << bRecalibClus << ", ";
        emcalTender << bRecalcClusPos << ", ";
        emcalTender << bNonLinearCorr << ", ";
        emcalTender << bRemExoticCell << ", ";
        emcalTender << bRemExoticClus << ", ";
        emcalTender << bFidRegion << ", ";
        emcalTender << bCalibEnergy << ", ";
        emcalTender << bCalibTime << ", ";
        emcalTender << bRemBC << ", ";
        emcalTender << iNonLinFunct << ", ";
        emcalTender << bReclusterize << ", ";
        emcalTender << fSeedThresh << ", ";
        emcalTender << fCellThresh << ", ";
        emcalTender << iClusterizer << ", ";
        emcalTender << bTrackMatch << ", ";
        emcalTender << bUpdateCellOnly << ", ";
        emcalTender << fEMCtimeMin << ", ";
        emcalTender << fEMCtimeMax << ", ";
        emcalTender << fEMCtimeCut << ", ";
        emcalTender << "0);"; // Passing cPass when it is 0 will act like dereferencing a null string!
        std::cout << "Calling tender task with " << emcalTender.str().c_str() << std::endl;
        AliAnalysisTaskSE * pTenderTask = reinterpret_cast<AliAnalysisTaskSE *>(gROOT->ProcessLine(emcalTender.str().c_str()));
        #else
        AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
                bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
                fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut, cPass);
        #endif
        pTenderTask->SelectCollisionCandidates(kPhysSel);

        #ifdef __CLING__
        std::stringstream clusterizerTask;
        clusterizerTask << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/AddTaskClusterizerFast.C(";
        clusterizerTask << "\"ClusterizerFast\", ";
        clusterizerTask << "\"\", ";
        clusterizerTask << "\"\", ";
        clusterizerTask << iClusterizer << ", ";
        clusterizerTask << 0.05 << ", ";
        clusterizerTask << 0.1 << ", ";
        clusterizerTask << fEMCtimeMin << ", ";
        clusterizerTask << fEMCtimeMax << ", ";
        clusterizerTask << fEMCtimeCut << ", ";
        clusterizerTask << kFALSE << ", ";
        clusterizerTask << kFALSE << ", ";
        clusterizerTask << static_cast<unsigned int>(AliAnalysisTaskEMCALClusterizeFast::kFEEData) << ");";
        std::cout << "Calling clusterizer task with " << clusterizerTask.str().c_str() << std::endl;
        AliAnalysisTaskEMCALClusterizeFast * pClusterizerTask = reinterpret_cast<AliAnalysisTaskEMCALClusterizeFast *>(gROOT->ProcessLine(clusterizerTask.str().c_str()));
        #else
        AliAnalysisTaskEMCALClusterizeFast * pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", iClusterizer,
                0.05, 0.1, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut,
                kFALSE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
        #endif
        pClusterizerTask->SelectCollisionCandidates(kPhysSel);
        bRemExoticClus  = kTRUE;
        iNonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrectedv3;

        #ifdef __CLING__
        std::stringstream clusterMaker;
        clusterMaker << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C(";
        clusterMaker << iNonLinFunct << ", ";
        clusterMaker << bRemExoticClus << ",";
        clusterMaker << "\"usedefault\", ";
        clusterMaker << "\"\", ";
        clusterMaker << 0. << ", ";
        clusterMaker << kTRUE << ");";
        std::cout << "Calling cluster maker task with " << clusterMaker.str().c_str() << std::endl;
        AliEmcalClusterMaker * pClusterMakerTask = reinterpret_cast<AliEmcalClusterMaker *>(gROOT->ProcessLine(clusterMaker.str().c_str()));
        #else
        AliEmcalClusterMaker * pClusterMakerTask = AddTaskEmcalClusterMaker(iNonLinFunct, bRemExoticClus, "usedefault", "", 0., kTRUE);
        #endif
        pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
        pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
        pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
        if (bIsRun2) {
          pClusterMakerTask->SetUseNewCentralityEstimation(kTRUE);
          pClusterMakerTask->SetNCentBins(5);
        }
    }

    if (bDoHadCorr) {
        // Cluster-track matcher task
        #ifdef __CLING__
        std::stringstream clusterTrackMatcher;
        clusterTrackMatcher << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C(";
        clusterTrackMatcher << "\"usedefault\", ";
        clusterTrackMatcher << "\"usedefault\", ";
        clusterTrackMatcher << 0.1 << ", ";
        clusterTrackMatcher << kFALSE << ", ";
        clusterTrackMatcher << kTRUE << ", ";
        clusterTrackMatcher << kTRUE << ", ";
        clusterTrackMatcher << kTRUE << ");";
        std::cout << "Calling cluster track matcher task with " << clusterTrackMatcher.str().c_str() << std::endl;
        AliEmcalClusTrackMatcherTask * pMatcherTask = reinterpret_cast<AliEmcalClusTrackMatcherTask *>(gROOT->ProcessLine(clusterTrackMatcher.str().c_str()));
        #else
        AliEmcalClusTrackMatcherTask * pMatcherTask = AddTaskEmcalClusTrackMatcher("usedefault", "usedefault", 0.1, kFALSE, kTRUE, kTRUE, kTRUE);
        #endif
        pMatcherTask->SelectCollisionCandidates(kPhysSel);
        pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
        pMatcherTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
        pMatcherTask->GetClusterContainer(0)->SetClusECut(0.);
        pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);

        if (bIsRun2) {
          pMatcherTask->SetUseNewCentralityEstimation(kTRUE);
          pMatcherTask->SetNCentBins(5);
        }

        if (iDataType == kEsd) {
          pMatcherTask->SetDoPropagation(kTRUE);
        }

        // Hadronic correction task
        #ifdef __CLING__
        std::stringstream hadCorr;
        hadCorr << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/AddTaskHadCorr.C(";
        hadCorr << "\"usedefault\", ";
        hadCorr << "\"usedefault\", ";
        hadCorr << "\"\", ";
        hadCorr << kHadCorrF << ", ";
        hadCorr << 0.15 << ", ";
        hadCorr << 0.030 << ", ";
        hadCorr << 0.015 << ", ";
        hadCorr << 0 << ", ";
        hadCorr << kTRUE << ", ";
        hadCorr << kTRUE << ");";
        std::cout << "Calling hadronic correction task with " << hadCorr.str().c_str() << std::endl;
        AliHadCorrTask * pHadCorrTask = reinterpret_cast<AliHadCorrTask *>(gROOT->ProcessLine(hadCorr.str().c_str()));
        #else
        AliHadCorrTask * pHadCorrTask = AddTaskHadCorr("usedefault", "usedefault", "",
                kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
        #endif
        pHadCorrTask->SelectCollisionCandidates(kPhysSel);
        pHadCorrTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.15);
        pHadCorrTask->GetClusterContainer(0)->SetClusECut(0);
        pHadCorrTask->GetClusterContainer(0)->SetClusPtCut(0.);
        if (bIsRun2) {
          pHadCorrTask->SetUseNewCentralityEstimation(kTRUE);
          pHadCorrTask->SetNCentBins(5);
        }
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
      AliEmcalJetTask *pFullJet02TaskNew = AliEmcalJetTask::AddTaskEmcalJet(newFrameworkTracks.Data(), newFrameworkClusters.Data(),
              jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, label, minJetPt, lockTask, fillGhosts);
      pFullJet02TaskNew->SelectCollisionCandidates(kPhysSel);
      pFullJet02TaskNew->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      if (bIsRun2) {
        pFullJet02TaskNew->SetUseNewCentralityEstimation(kTRUE);
        pFullJet02TaskNew->SetNCentBins(5);
      }

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
      AliEmcalJetTask *pFullJet02Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault",
              jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, label, minJetPt, lockTask, fillGhosts);
      pFullJet02Task->SelectCollisionCandidates(kPhysSel);
      pFullJet02Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      pFullJet02Task->SetNCentBins(5);
      if (bIsRun2) {
        pFullJet02Task->SetUseNewCentralityEstimation(kTRUE);
        pFullJet02Task->SetNCentBins(5);
      }
    }

    //////////////////////////////////////////
    // Add you task here!
    //////////////////////////////////////////
    AliEmcalJet::JetAcceptanceType acceptanceType = AliEmcalJet::kEMCALfid;
    TString jetTag = "Jet";
    if (newCorrectionFramework) {
      // Your task using newFrameworkTracks and newFrameworkClusters
      // Sample task
      AliAnalysisTaskEmcalJetSample *sampleTask= 0;
      sampleTask = AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample(newFrameworkTracks.Data(), newFrameworkClusters.Data(), newFrameworkCells.Data());
      if (newFrameworkTracks != "tracks" && newFrameworkTracks != "Tracks")
      {
        // Remove wrong particle container
        sampleTask->RemoveParticleContainer(0);
        // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
        AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
        sampleTask->AdoptTrackContainer(newTracks);
      }
      sampleTask->GetClusterContainer(0)->SetClusECut(0.);
      sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
      sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
      sampleTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
      sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      sampleTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
      sampleTask->SetHistoBins(600, 0, 300);
      sampleTask->SelectCollisionCandidates(kPhysSel);
      if (bIsRun2) {
        sampleTask->SetUseNewCentralityEstimation(kTRUE);
        sampleTask->SetNCentBins(5);
      }

      AliJetContainer* jetCont02 = sampleTask->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, jetTag.Data());
    }
    if (oldCorrectionFramework) {
      // Your task using the normal tracks and clusters
      // Sample task
      AliAnalysisTaskEmcalJetSample *sampleTask= 0;
      sampleTask = AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample("usedefault", "usedefault", "usedefault");
      sampleTask->GetClusterContainer(0)->SetClusECut(0.);
      sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
      sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
      sampleTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
      sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      sampleTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
      sampleTask->SetHistoBins(600, 0, 300);
      sampleTask->SelectCollisionCandidates(kPhysSel);
      if (bIsRun2) {
        sampleTask->SetUseNewCentralityEstimation(kTRUE);
        sampleTask->SetNCentBins(5);
      }

      AliJetContainer* jetCont02 = sampleTask->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, acceptanceType, jetTag.Data());
    }

    // QA task
    /*AliAnalysisTaskEmcalJetQA *pQATask = 0;
    if (bDoFullJets) {
        pQATask = AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault");
        pQATask->GetClusterContainer(0)->SetClusECut(0.);
        pQATask->GetClusterContainer(0)->SetClusPtCut(0.);
        pQATask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
        pQATask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
        pQATask->SetSeparateEMCalDCal(bSeparateEMCalDCal);
    }
    else {
        pQATask = AddTaskEmcalJetQA("usedefault", "", "");
    }
    pQATask->GetParticleContainer(0)->SetParticlePtCut(0.15);
    pQATask->SetHistoBins(300, 0, 150);
    pQATask->SelectCollisionCandidates(kPhysSel);

    // Background
    TString sRhoChName;
    TString sRhoFuName;
    if (iBeamType != AliAnalysisTaskEmcal::kpp) {
        sRhoChName = "Rho";
        sRhoFuName = "Rho_Scaled";

        AliEmcalJetTask *pKtChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, kJetRadius, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
        pKtChJetTask->SelectCollisionCandidates(kPhysSel);

        AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChName, kJetRadius);
        pRhoTask->SetExcludeLeadJets(2);
        pRhoTask->SelectCollisionCandidates(kPhysSel);

        if (bDoFullJets) {
            TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
            TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
            pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
        }
    }

    // Charged jet analysis
    if (bDoChargedJets) {
        AliEmcalJetTask *pChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, kJetRadius, AliJetContainer::kChargedJet, 0.15, 0, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
        pChJetTask->SelectCollisionCandidates(kPhysSel);
    }

    // Full jet analysis
    if (bDoFullJets) {
        AliEmcalJetTask *pFuJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault", AliJetContainer::antikt_algorithm, kJetRadius, AliJetContainer::kFullJet, 0.15, 0.30, kGhostArea, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
        pFuJetTask->SelectCollisionCandidates(kPhysSel);
    }

    AliAnalysisTaskEmcalJetSpectraQA *pSpectraTask = AddTaskEmcalJetSpectraQA("usedefault", "usedefault");

    if (bDoFullJets) {
        AliJetContainer* jetCont = pSpectraTask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliEmcalJet::kEMCALfid, "Jet");
        if (iBeamType != AliAnalysisTaskEmcal::kpp) {
            jetCont->SetRhoName(sRhoFuName);
            jetCont->SetPercAreaCut(0.6);
        }
    }

    if (bDoChargedJets) {
        AliJetContainer* jetCont = pSpectraTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliEmcalJet::kTPCfid, "Jet");
        if (iBeamType != AliAnalysisTaskEmcal::kpp) {
            jetCont->SetRhoName(sRhoChName);
            jetCont->SetPercAreaCut(0.6);
        }
    }

    if (bDoMCJets) {
        AliMCParticleContainer* MCpartCont = pSpectraTask->AddMCParticleContainer("mcparticles");

        if (bDoChargedJets) {
            pSpectraTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliEmcalJet::kTPCfid, MCpartCont, 0);
        }

        if (bDoFullJets) {
            pSpectraTask->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, kJetRadius, AliEmcalJet::kTPCfid, MCpartCont, 0);
        }
    }

    pSpectraTask->SetNLeadingJets(1);
    pSpectraTask->SelectCollisionCandidates(kPhysSel);
    pSpectraTask->SetHistoType(kHistoType);*/

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
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
    /*gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSpectraQA.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");*/
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
