/// \file runJetShapes.C
/// \brief Run the test shapes tasks
///
/// \ingroup EMCALJETFW
/// Run macro for the Jet Shapes tasks
///
/// To select which framework to run, you need to select a few options.
/// To use the new framework, set
///   newFramework = true;
///   oldFramework = false;
///   tag = "new";
/// To use the old framework, set
///   newFramework = false;
///   oldFramework = true;
///   tag = "old";
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
/// \date Mar 10, 2017

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

class AliAnalysisTaskEmcalEmbeddingHelper;
class AliEmcalCopyCollection;
class AliEmcalCorrectionTask;
class AliEmcalJetTask;

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runJetShapes(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC11h",                                // set the run period
    const char   *cLocalFiles    = "aodFiles.txt",                          // set the local list file
    const UInt_t  iNumEvents     = 1000,                                    // number of events to be analyzed.
    const UInt_t  kPhysSel       = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral, //AliVEvent::kAny,                         // physics selection
    const char   *cTaskName      = "JetShapesAnalysis",                     // sets name of analysis manager
    // 0 = only prepare the analysis manager but do not start the analysis
    // 1 = prepare the analysis manager and start the analysis
    // 2 = launch a grid analysis
    Int_t         iStartAnalysis = 1,
    const UInt_t  iNumFiles      = 5,                                     // number of files analyzed locally
    const char   *cGridMode      = "test"
)
{
  // Select which embedding framework to use
  const bool newFramework = true;
  const bool oldFramework = false;
  const bool fullJets = false;
  const bool constituentSubtraction = true;
  const TString tag = "new";

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

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

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

  // Define relevant variables
  const bool IsEsd = (iDataType == kEsd);
  // Note what tag we are using
  TString truthParticlesName = "";
  TString emcalCellsName = "emcalCells";
  TString clustersName = "";
  TString clustersNameCombined = "";
  TString tracksName = "";

  // Define for convenience
  const Double_t kClusPtCut = 0.30;
  const Double_t kTrackPtCut = 0.15;
  // Used to set the new framework to the same mass hypothesis as the PicoTracks (which are hard coded)
  const Double_t massHypothesis = 0.13957;

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

    // Setup embedding task
    AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
    embeddingHelper->SetPtHardBin(4);
    //embeddingHelper->SetFilePattern("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/");
    embeddingHelper->SetFileListFilename("aodFilesEmbed.txt");
    embeddingHelper->SetTriggerMask(AliVEvent::kAny); // Equivalent to 0
    embeddingHelper->SelectCollisionCandidates(kPhysSel);
    embeddingHelper->SetRandomFileAccess(kFALSE);
    embeddingHelper->Initialize();

    // Setup for LEGO train
    /*
    __R_ADDTASK__->SetPtHardBin(kPtHardBin);
    __R_ADDTASK__->SetAnchorRun(169838);
    // Revise to the new function once it is pushed!
    __R_ADDTASK__->SetFilePattern("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/ AliAOD.root");
    //__R_ADDTASK__->SetFilePattern("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/");
    __R_ADDTASK__->SetTriggerMask(AliVEvent::kAny); // Equivalent to 0
    __R_ADDTASK__->SelectCollisionCandidates(kPhysSel);
    __R_ADDTASK__->SetRandomFileAccess(kTRUE);
    __R_ADDTASK__->Initialize();
    */

    // For comparison
    // ZVertex cut should be +/- 10 cm. This is the default

    // TEMP
    //AliLog::SetClassDebugLevel("AliAnalysisTaskEmcalEmbeddingHelper", AliLog::kDebug+0);
  }

  if (oldFramework) {
    // Setup relevant values
    truthParticlesName = "MCSelectedParticles";
    tracksName = "PicoTracks";

    // AODTrackMaker
    // TODO: ROOT6
    // Create AOD tracks. This is using real data
    AliEmcalAodTrackFilterTask *aodTrackMaker = AddTaskEmcalAodTrackFilter("AODFilterTracks", "tracks", sRunPeriod.Data());
    aodTrackMaker->SelectCollisionCandidates(kPhysSel);
    if (fullJets == false) {
      // Disables track prop when no EMCal so the comparison matches up!
      aodTrackMaker->SetAttemptProp(kFALSE);
      aodTrackMaker->SetAttemptPropMatch(kFALSE);
    }

    // PicoTrackMaker
    // TODO: ROOT6
    // Create pico tracks according to the AOD track filter. This uses PbPb data.
    AliEmcalPicoTrackMaker * picoTracks = AddTaskEmcalPicoTrackMaker(tracksName, "AODFilterTracks");
    picoTracks->SelectCollisionCandidates(kPhysSel);

    // EmbeddingFromPythia
    // TODO: ROOT6
    // Extracts information from PYTHIA and injects the tracks into the pico tracks
    // NOTE: Geometry warnings are generated here! The EMCal gemoetry is hard coded in the base class. However, it's not worth my time to fix
    TString mcTracksName = "MCSelectedParticles";
    //TString sPYTHIAPath = "alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root";
    TString sPYTHIAPath = "/Users/re239/code/alice/data/LHC12a15e_fix/%d/%d/AOD149/%04d/root_archive.zip#AliAOD.root";
    Int_t kNpTHardBins = 11;
    //Double_t kPtHardBinsScaling[11] = {0.000000E+00, 5.206101E-05, 5.859497E-06, 4.444755E-07, 4.344664E-08, 5.154750E-09, 6.956634E-10, 1.149828E-10, 2.520137E-11, 6.222240E-12, 2.255832E-12};
    Double_t kPtHardBinsScaling[11] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
    Double_t kMinPythiaJetPt = 0;
    // Convert to AODTask to allow for local files!
    AliJetEmbeddingFromPYTHIATask * pythiaEmbedding = AddTaskJetEmbeddingFromPYTHIA(tracksName, "",
    //AliJetEmbeddingFromAODTask * pythiaEmbedding = AddTaskJetEmbeddingFromAOD(tracksName, "",
                    emcalCellsName,                         // EMCal cells name
                    truthParticlesName,                     // Name of the MC at particle level
                    sPYTHIAPath,                            // PYTHIATask only - Path to PYTHIA files. Seach string
                    kNpTHardBins, kPtHardBinsScaling,       // PYTHIATask only - Pt hard bin settings
                    //"aodFilesEmbed.txt",                    // AODTask only - Path to file list
                    "aodTree", "tracks", "", "emcalCells", "mcparticles",   // AOD standard branch information. The combined tracks and cells will be injected into these branches
                    "lhc12a15e",                            // Run period for MC particles to embed
                    kFALSE,                                 // False includes the ITS
                    -1, -1,                                 // Centrality settings
                    0,                                      // Trigger Mask
                    kMinPythiaJetPt,                        // PYTHIAtask only
                    kFALSE, kTRUE,                          // Copy array and make QA
                    //1234567890,                             // AODTask only - nFiles (Just taking the default so we can set the name argument next)
                    "",                                     // PYTHIATask only - Generates a file table (perhaps use instead of PYTHIApath?)
                    TString::Format("JetEmbeddingFromAODTask_%s", tag.Data())); // Set the file name
    // NOTE: Causes crash due to it running before the tasks above, and therefore creating branches that it shouldn't!
    //pythiaEmbedding->SelectCollisionCandidates(AliVEvent::kAny);
    pythiaEmbedding->SelectCollisionCandidates(kPhysSel);
    // ======== Changes needed for agreement =========
    // NOTE: Disabled random access for a consitent comparison between the two tags!
    //pythiaEmbedding->SetRandomAccess(kTRUE);
    // NOTE: Disable track efficiency for a consistent comparsion between the two tags!
    //pythiaEmbedding->SetTrackEfficiency(0.98);
    // Below much be uncommented for a successful comparsion.
    // NOTE: random access doesn't seem to mean anything in the PYTHIA embedding task. This may be a bug!
    // However, it matters for the AOD task functions that are called by the PYTHIA task
    // Must be explicitly disabled because the task enables it by default
    pythiaEmbedding->SetRandomAccess(kFALSE);
    // Ensure that files are embedded in order. This makes limiting the total number of files obsolete!
    pythiaEmbedding->SetDebugEmbedding(kTRUE);
    // ======== End changes needed for agreement =========
    // Setting to 1 ensures that the same file is always embedded
    //pythiaEmbedding->SetTotalFiles(1);
    pythiaEmbedding->SetTotalFiles(140);
    pythiaEmbedding->SetAODMC(kTRUE);
    pythiaEmbedding->SetIncludeNoITS(kFALSE);
    pythiaEmbedding->SetAODfilterBits(256,512);
    pythiaEmbedding->SetMarkMC(99999);
    // Not compatiable with AODTask
    //pythiaEmbedding->SetMinEntriesPerPtHardBin(999999);

    // TEMP
    //AliLog::SetClassDebugLevel("AliJetEmbeddingFromPYTHIATask", AliLog::kDebug-1);
  }

  AliJetContainer::EJetType_t jetType = AliJetContainer::kChargedJet;
  const AliJetContainer::ERecoScheme_t recoScheme = AliJetContainer::E_scheme;
  AliEmcalJet::JetAcceptanceType acceptanceType = AliEmcalJet::kTPCfid;
  if (fullJets) {
    // TODO: Full implement this
    AliFatalGeneral("EMCal support not fully implemented in this run macro!");
    jetType = AliJetContainer::kFullJet;
    acceptanceType = AliEmcalJet::kEMCALfid;
  }

  // Print status
  Printf("Using the follow input object names:");
  Printf("cells: %s", emcalCellsName.Data());
  Printf("clusters: %s", clustersName.Data());
  Printf("clustersCombined: %s", clustersNameCombined.Data());
  Printf("tracks: %s", tracksName.Data());

  // JetFinderQGAKTCharged_R02_EschemePythiaEmbPart
  // Clusters will always be empty here because this jet finder is for the truth only!
  AliEmcalJetTask *pAKtPythiaEmb = AliEmcalJetTask::AddTaskEmcalJet(truthParticlesName.Data(), "", AliJetContainer::antikt_algorithm, 0.2, jetType, 0, kClusPtCut, kGhostArea, recoScheme, TString::Format("%s_%s","JetMCOnlyPartLevel", tag.Data()), 0, kFALSE, kFALSE);
  pAKtPythiaEmb->SelectCollisionCandidates(AliVEvent::kCentral);
  pAKtPythiaEmb->SetVzRange(-10,10);
  pAKtPythiaEmb->SetNeedEmcalGeom(kFALSE);
  pAKtPythiaEmb->SetZvertexDiffValue(0.1);
  if (newFramework) {
    // Retrieve the MC particles from the embedded event
    pAKtPythiaEmb->GetParticleContainer(0)->SetIsEmbedding(kTRUE);
    // Don't set the mass hypothesis here!! The mass is already known because these are truth level!
  }

  // JetFinderKtTpcQG_ESchemePythiaEmbForR2
  // This includes all tracks
  AliEmcalJetTask *pKtData = AliEmcalJetTask::AddTaskEmcalJet(tracksName.Data(), clustersNameCombined.Data(), AliJetContainer::kt_algorithm, 0.2, AliJetContainer::kChargedJet, kTrackPtCut, 0, kGhostArea, recoScheme, TString::Format("%s_%s", "Jet", tag.Data()),0,kFALSE,kFALSE);
  pKtData->SelectCollisionCandidates(AliVEvent::kCentral);
  pKtData->SetVzRange(-10,10);
  pKtData->SetNeedEmcalGeom(kFALSE);
  pKtData->SetZvertexDiffValue(0.1);
  if (newFramework) {
    // Add the embedded detector tracks
    AliTrackContainer * tracksDetLevel = new AliTrackContainer(tracksName.Data());
    tracksDetLevel->SetIsEmbedding(kTRUE);
    tracksDetLevel->SetParticlePtCut(kTrackPtCut);
    tracksDetLevel->SetMassHypothesis(massHypothesis);
    pKtData->AdoptTrackContainer(tracksDetLevel);

    // Set existing track container mass hypothesis
    pKtData->GetParticleContainer(0)->SetMassHypothesis(massHypothesis);

    if (fullJets) {
      // Nothing explicitly needed, since the cells will be combined together already
    }
  }
  if (oldFramework) {
    // Nothing explicitly needed for the old framework
    if (fullJets) {
      // Nothing explicitly needed for the old framework
      // TODO: Check on this
    }
  }

  // JetFinderQGAKTCharged_R02_EschemePythiaEmbTrue
  // TODO: To be understood, but I think this only the detector level truth
  AliEmcalJetTask *pAKtDataEmb = AliEmcalJetTask::AddTaskEmcalJet(tracksName.Data(), clustersName.Data(), AliJetContainer::antikt_algorithm, 0.2, jetType, kTrackPtCut, kClusPtCut, kGhostArea, recoScheme, TString::Format("%s_%s","JetMCOnly",tag.Data()), 0, kFALSE, kFALSE);
  pAKtDataEmb->SelectCollisionCandidates(AliVEvent::kCentral);
  pAKtDataEmb->SetVzRange(-10,10);
  pAKtDataEmb->SetNeedEmcalGeom(kFALSE);
  pAKtDataEmb->SetZvertexDiffValue(0.1);
  if (newFramework) {
    pAKtDataEmb->GetParticleContainer(0)->SetIsEmbedding(kTRUE);
    pAKtDataEmb->GetParticleContainer(0)->SetMassHypothesis(massHypothesis);

    if (fullJets) {
      pAKtDataEmb->GetClusterContainer(0)->SetIsEmbedding(kTRUE);
    }
  }
  if (oldFramework) {
    pAKtDataEmb->GetParticleContainer(0)->SetMCLabelRange(99999,99999999);
    if (fullJets) {
      // TODO: Is this correct??
      pAKtDataEmb->GetClusterContainer(0)->SetMCLabelRange(99999,99999999);
    }
  }

  // Background
  // AliAnalysisTaskQGRhoTpcExLJ_ESchemePythiaEmbForR2
  // TODO: ROOT 6
  // NOTE: Uses charged jets and TPCfid regardless of whether using charged or full jets!
  TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
  AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew(tracksName.Data(), clustersNameCombined.Data(), "Rho", 0.2, AliJetContainer::kTPCfid, AliJetContainer::kChargedJet, kTRUE, recoScheme, tag.Data());
  //AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("PicoTracks",    "", "Rho", 0.2, AliJetContainer::kTPCfid, AliJetContainer::kChargedJet, kTRUE, AliJetContainer::E_scheme, tag.Data());
  pRhoTask->SetScaleFunction(sfunc);
  pRhoTask->SetExcludeLeadJets(2);
  pRhoTask->SelectCollisionCandidates(AliVEvent::kCentral);
  // NOTE: These were commented in the wagon
  //pRhoTask->GetParticleContainer(0)->SetMCTrackBitMap(TObject::kBitMask);
  //pRhoTask->GetClusterContainer(0)->SetMCClusterBitMap(TObject::kBitMask);
  //pRhoTask->SetHistoBins(250,0,250);
  pRhoTask->SetNeedEmcalGeom(kFALSE);
  pRhoTask->SetZvertexDiffValue(0.1);
  pRhoTask->SetVzRange(-10,10);
  AliJetContainer * jetCont = pRhoTask->GetJetContainer(0);
  // Fix jet container name, since the Rho AddTask doesn't accept a tag for the jet name
  TString jetContName = jetCont->GetName();
  jetContName = jetContName.Insert(jetContName.Index("Jet_") + 4, TString::Format("%s_", tag.Data()));
  jetCont->SetName(jetContName.Data());
  jetCont->SetArrayName(jetContName.Data());
  // Configure the jet conatiner
  jetCont->SetJetRadius(0.2);
  jetCont->SetJetAcceptanceType(AliJetContainer::kTPCfid);
  jetCont->SetMaxTrackPt(100);
  // NOTE: The tracks shouldn't matter much here, but including for completeness
  if (newFramework) {
    // Add the embedded detector tracks
    AliTrackContainer * tracksDetLevel = new AliTrackContainer(tracksName.Data());
    tracksDetLevel->SetIsEmbedding(kTRUE);
    tracksDetLevel->SetParticlePtCut(kTrackPtCut);
    tracksDetLevel->SetMassHypothesis(massHypothesis);
    pRhoTask->AdoptTrackContainer(tracksDetLevel);

    // Set existing track container mass hypothesis
    pRhoTask->GetParticleContainer(0)->SetMassHypothesis(massHypothesis);

    if (fullJets) {
      // Nothing explicitly needed, since the cells will be combined together already
    }
  }

  // Rho Mass
  // RhoMassQGTPCPythiaEmbForR2
  // TODO: ROOT 6
  // Defined for convenience
  TF1* srhomfunc = new TF1("srhomfunc","[0]*x*x+[1]*x+[2]",-1,100);
  srhomfunc->SetParameter(2, 1.68354);
  srhomfunc->SetParameter(1, -2.86991e-03);
  srhomfunc->SetParameter(0, -1.49981e-05);
  // Task
  // PicoTracks should only use a particle container
  AliParticleContainer * partCont = newFramework ? new AliTrackContainer(tracksName.Data()) : new AliParticleContainer(tracksName.Data());
  AliClusterContainer * clusCont = 0;
  if (fullJets) {
    clusCont = new AliClusterContainer(clustersNameCombined.Data());
  }
  // GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag)
  TString rhoMassJetName = AliJetContainer::GenerateJetName(AliJetContainer::kChargedJet, AliJetContainer::kt_algorithm, recoScheme, 0.2, partCont, clusCont, TString::Format("Jet_%s", tag.Data()));
  AliAnalysisTaskRhoMass * pRhoMass = AddTaskRhoMass(rhoMassJetName.Data(), tracksName.Data(), clustersNameCombined.Data(), "Rhomass", 0.2, "TPC", 0.01, 0, 0, 2, kTRUE, TString::Format("%s_%s", "RhoMass", tag.Data()));
  //AliAnalysisTaskRhoMass * pRhoMass = AddTaskRhoMass(TString::Format("Jet_%s_KTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()), "PicoTracks", "", "Rhomass", 0.2, "TPC", 0.01, 0, 0, 2, kTRUE, TString::Format("%s_%s", "RhoMass", tag.Data()));
  pRhoMass->SelectCollisionCandidates(AliVEvent::kCentral);
  pRhoMass->SetHistoBins(250,0,250);
  pRhoMass->SetScaleFunction(srhomfunc);
  pRhoMass->SetNeedEmcalGeom(kFALSE);
  jetCont = pRhoMass->GetJetContainer(0);
  jetCont->SetJetRadius(0.2);
  jetCont->SetJetAcceptanceType(AliJetContainer::kTPCfid);
  jetCont->SetMaxTrackPt(100);
  pRhoMass->SetNeedEmcalGeom(kFALSE);
  pRhoMass->SetZvertexDiffValue(0.1);
  pRhoMass->SetVzRange(-10,10);
  if (newFramework) {
    // Remove first particle container and add it as a track container to ensure that track selection is applied
    // The AddTask does not does this itself.
    // TODO: Consider updating the AddTask to do this...
    pRhoMass->RemoveParticleContainer(0);
    AliTrackContainer * tracks = new AliTrackContainer(tracksName.Data());
    tracks->SetParticlePtCut(kTrackPtCut);
    tracks->SetMassHypothesis(massHypothesis);
    pRhoMass->AdoptTrackContainer(tracks);

    // Add the embedded detector tracks
    AliTrackContainer * tracksDetLevel = new AliTrackContainer(tracksName.Data());
    tracksDetLevel->SetIsEmbedding(kTRUE);
    tracksDetLevel->SetParticlePtCut(kTrackPtCut);
    tracksDetLevel->SetMassHypothesis(massHypothesis);
    pRhoMass->AdoptTrackContainer(tracksDetLevel);

    if (fullJets) {
      // Nothing explicitly needed, since the cells will be combined together already
    }
  }

  // Final jet finder
  // JetFinderQGAKTCharged_R02_EschemePythiaEmb
  // This includes all tracks
  AliEmcalJetTask * pFinalJets = AliEmcalJetTask::AddTaskEmcalJet(tracksName.Data(), clustersNameCombined.Data(), AliJetContainer::antikt_algorithm, 0.2, jetType, kTrackPtCut, kClusPtCut, kGhostArea, recoScheme, TString::Format("%s_%s", "Jet", tag.Data()), 0, kFALSE, kFALSE);
  pFinalJets->SelectCollisionCandidates(AliVEvent::kCentral);
  pFinalJets->SetVzRange(-10,10);
  pFinalJets->SetNeedEmcalGeom(kFALSE);
  pFinalJets->SetZvertexDiffValue(0.1);
  if (newFramework) {
    // Add the embedded particles
    AliTrackContainer * tracksDetLevel = new AliTrackContainer(tracksName.Data());
    tracksDetLevel->SetIsEmbedding(kTRUE);
    tracksDetLevel->SetParticlePtCut(kTrackPtCut);
    tracksDetLevel->SetMassHypothesis(massHypothesis);
    pFinalJets->AdoptTrackContainer(tracksDetLevel);

    // Set existing track container mass hypothesis
    pFinalJets->GetParticleContainer(0)->SetMassHypothesis(massHypothesis);

    if (fullJets) {
      // Nothing explicitly needed, since the cells will be combined together already
    }
  }

  AliEmcalJetUtilityConstSubtractor* constUtil = pFinalJets->AddUtility(new AliEmcalJetUtilityConstSubtractor("ConstSubtractor"));

  AliEmcalJetUtilityGenSubtractor* genSub = pFinalJets->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
  genSub->SetGenericSubtractionJetMass(kTRUE);
  genSub->SetGenericSubtractionExtraJetShapes(kTRUE);
  genSub->SetUseExternalBkg(kTRUE);
  genSub->SetRhoName("Rho");
  genSub->SetRhomName("Rhomass");
  constUtil->SetJetsSubName(Form("%sConstSub", pFinalJets->GetName()));
  //constUtil->SetConstituentSubtraction(kTRUE);
  constUtil->SetParticlesSubName("TracksSub");
  constUtil->SetUseExternalBkg(kTRUE);
  constUtil->SetRhoName("Rho");
  constUtil->SetRhomName("Rhomass");

  // Jet Taggers
  // JetTaggerR02PythiaEmb
  // TODO: ROOT 6
  //
  //TString::Format("Jet_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data())
  //TString::Format("JetMCOnly_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data())
  //
  // GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag)
  TString taggerEmbJetName1 = AliJetContainer::GenerateJetName(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, recoScheme, 0.2, partCont, clusCont, TString::Format("Jet_%s", tag.Data()));
  TString taggerEmbJetName2 = AliJetContainer::GenerateJetName(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, recoScheme, 0.2, partCont, clusCont, TString::Format("JetMCOnly_%s", tag.Data()));
  AliAnalysisTaskEmcalJetTagger * taggerEmb = AddTaskEmcalJetTagger(taggerEmbJetName1.Data(), taggerEmbJetName2, 0.2, "", "", tracksName.Data(), clustersNameCombined.Data(), "TPC", "V0M", kPhysSel, "", "");
  //AliAnalysisTaskEmcalJetTagger * taggerEmb = AddTaskEmcalJetTagger(TString::Format("Jet_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()),TString::Format("JetMCOnly_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()),0.2,"","","PicoTracks","","TPC","V0M",kPhysSel,"","");
  taggerEmb->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest);
  taggerEmb->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
  taggerEmb->SelectCollisionCandidates(AliVEvent::kCentral);
  taggerEmb->SetTypeAcceptance(0);
  cont = taggerEmb->GetJetContainer(0);
  AliJetContainer *cont2 = taggerEmb->GetJetContainer(1);
  cont->SetMaxTrackPt(100);
  cont2->SetMaxTrackPt(100);

  // JetTaggerR02PythiaEmbDetPart
  // TODO: ROOT 6
  //
  //TString::Format("JetMCOnly_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data())
  //TString::Format("JetMCOnlyPartLevel_%s_AKTChargedR020_MCSelectedParticles_pT0000_E_scheme", tag.Data())
  //
  // NOTE: taggerEmbDetPartJetName1 == taggerEmbJetName2
  TString taggerEmbDetPartJetName1 = AliJetContainer::GenerateJetName(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, recoScheme, 0.2, partCont, clusCont, TString::Format("JetMCOnly_%s", tag.Data()));
  // Need the mc particles particle container here
  AliMCParticleContainer * mcPartCont = new AliMCParticleContainer(truthParticlesName.Data());
  mcPartCont->SetParticlePtCut(0);
  if (newFramework) {
    mcPartCont->SetIsEmbedding(kTRUE);
    // Don't set the mass hypothesis here!! The mass is already known because these are truth level!
    //mcPartCont->SetMassHypothesis(massHypothesis);
  }
  TString taggerEmbDetPartJetName2 = AliJetContainer::GenerateJetName(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, recoScheme, 0.2, mcPartCont, clusCont, TString::Format("JetMCOnlyPartLevel_%s", tag.Data()));
  AliAnalysisTaskEmcalJetTagger * taggerEmbDetPart = AddTaskEmcalJetTagger(taggerEmbDetPartJetName1.Data(), taggerEmbDetPartJetName2.Data(), 0.2, "", "", tracksName.Data(), clustersNameCombined.Data(), "TPC", "V0M", kPhysSel, "", "");
  //AliAnalysisTaskEmcalJetTagger * taggerEmbDetPart = AddTaskEmcalJetTagger(TString::Format("JetMCOnly_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()),TString::Format("JetMCOnlyPartLevel_%s_AKTChargedR020_MCSelectedParticles_pT0000_E_scheme", tag.Data()),0.2,"","","PicoTracks","","TPC","V0M",kPhysSel,"","");
  taggerEmbDetPart->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest);
  taggerEmbDetPart->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
  taggerEmbDetPart->SelectCollisionCandidates(AliVEvent::kCentral);
  taggerEmbDetPart->SetTypeAcceptance(0);
  taggerEmbDetPart->SetIsPythia(kTRUE);
  cont = taggerEmbDetPart->GetJetContainer(0);
  cont2 = taggerEmbDetPart->GetJetContainer(1);
  cont->SetMaxTrackPt(100);
  cont2->SetMaxTrackPt(100);

  // QG Tagging Task
  // RhoMassQGTPCPythiaEmbForR2
  // TODO: ROOT 6
  //
  //TString::Format("Jet_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data())
  //TString::Format("Jet_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data())
  //TString::Format("JetMCOnly_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data())
  //TString::Format("JetMCOnlyPartLevel_%s_AKTChargedR020_MCSelectedParticles_pT0000_E_scheme", tag.Data())
  //
  TString tracksSub = tracksName.Data();
  AliAnalysisTaskEmcalQGTagging::JetShapeSub subType = AliAnalysisTaskEmcalQGTagging::kDerivSub;
  TString qgTaggingJetName1 = taggerEmbJetName1.Data();
  if (constituentSubtraction) {
    tracksSub = "TracksSub";
    subType = AliAnalysisTaskEmcalQGTagging::kConstSub;
    qgTaggingJetName1 += "ConstSub";
  }
  AliAnalysisTaskEmcalQGTagging * qgTagging = AddTaskEmcalQGTagging(qgTaggingJetName1.Data(), taggerEmbJetName1.Data(), taggerEmbJetName2.Data(), taggerEmbDetPartJetName2.Data(), 0.2, "Rho", tracksSub.Data(), tracksName.Data(), truthParticlesName.Data(), clustersNameCombined.Data(), truthParticlesName.Data(), "TPC", "V0M", kPhysSel, "", "", "RawConstSub", AliAnalysisTaskEmcalQGTagging::kDetEmbPartPythia, subType, AliAnalysisTaskEmcalQGTagging::kInclusive);
  //AliAnalysisTaskEmcalQGTagging * qgTagging = AddTaskEmcalQGTagging(TString::Format("Jet_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()),TString::Format("Jet_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()),TString::Format("JetMCOnly_%s_AKTChargedR020_PicoTracks_pT0150_E_scheme", tag.Data()),TString::Format("JetMCOnlyPartLevel_%s_AKTChargedR020_MCSelectedParticles_pT0000_E_scheme", tag.Data()),0.2,"Rho","PicoTracks","PicoTracks","MCSelectedParticles","","MCSelectedParticles","TPC","V0M",kPhysSel,"","","RawConstSub",AliAnalysisTaskEmcalQGTagging::kDetEmbPartPythia,AliAnalysisTaskEmcalQGTagging::kDerivSub , AliAnalysisTaskEmcalQGTagging::kInclusive);
  cont = qgTagging->GetJetContainer(0);
  cont->SetRhoName("Rho");
  cont->SetRhoMassName("Rhomass");
  cont->SetJetRadius(0.2);
  cont->SetJetAcceptanceType(AliJetContainer::kTPCfid);
  cont->SetMaxTrackPt(100);

  cont2 = qgTagging->GetJetContainer(1);
  cont2->SetJetRadius(0.2);
  cont2->SetJetAcceptanceType(AliJetContainer::kTPCfid);
  cont2->SetMaxTrackPt(100);
  if (constituentSubtraction) {
    cont2->SetRhoName("Rho");
  }

  qgTagging->SetMinFractionShared(0.5);
  qgTagging->SetJetPtThreshold(20);
  qgTagging->SelectCollisionCandidates(AliVEvent::kCentral);
  qgTagging->SetVzRange(-10,10);
  qgTagging->SetNeedEmcalGeom(kFALSE);
  qgTagging->SetZvertexDiffValue(0.1);

  AliJetContainer *cont3 = qgTagging->GetJetContainer(2);
  cont3->SetJetRadius(0.2);
  cont3->SetJetAcceptanceType(AliJetContainer::kTPCfid);
  cont3->SetMaxTrackPt(100);
  if (constituentSubtraction) {
    cont3->SetRhoName("Rho");
  }

  if (newFramework) {
    // Determine the indices by inspecting the task. This is _not_ robust
    // Don't set the mass hypothesis here!! The mass is already known because both of these containers are truth level! (at least in this configuration)
    // Explict cast to MCParticleContainer so this breaks if the configuration changes
    mcPartCont = qgTagging->GetMCParticleContainer(2);
    mcPartCont->SetIsEmbedding(kTRUE);
    mcPartCont = qgTagging->GetMCParticleContainer(3);
    mcPartCont->SetIsEmbedding(kTRUE);
  }

  ////////////////////////
  //
  ////////////////////////

  AliAnalysisTaskEmcalJetSample * sampleTaskNew = AddTaskEmcalJetSample(tracksName.Data(), "", "", tag.Data());
  /*sampleTaskNew->GetClusterContainer(0)->SetClusECut(0.3);
  sampleTaskNew->GetParticleContainer(0)->SetParticlePtCut(0.15);*/
  if (newFramework) {
    // Add the embedded particles
    AliTrackContainer * tracksDetLevel = new AliTrackContainer(tracksName.Data());
    tracksDetLevel->SetIsEmbedding(kTRUE);
    sampleTaskNew->AdoptTrackContainer(tracksDetLevel);

    if (fullJets) {
      // Nothing explicitly needed, since the cells will be combined together already
    }
  }
  if (oldFramework) {
    sampleTaskNew->GetParticleContainer(0)->SetClassName("AliPicoTrack");
  }
  sampleTaskNew->SetHistoBins(600, 0, 300);
  sampleTaskNew->SelectCollisionCandidates(AliVEvent::kCentral);

  AliJetContainer* jetContAKT = sampleTaskNew->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, 0.2, acceptanceType, TString::Format("Jet_%s", tag.Data()));

  AliJetContainer* jetContKT = sampleTaskNew->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::kt_algorithm, AliJetContainer::E_scheme, 0.2, acceptanceType, TString::Format("Jet_%s", tag.Data()));

  // TEMP
  //AliLog::SetClassDebugLevel("AliFJWrapper", AliLog::kDebug-2);
  //AliLog::SetClassDebugLevel("AliEmcalJetTask", AliLog::kDebug-1);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskRhoMass", AliLog::kDebug+0);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskRhoMassBase", AliLog::kDebug+0);
  //AliLog::SetClassDebugLevel("AliEmcalContainer", AliLog::kDebug+7);
  //AliLog::SetClassDebugLevel("AliParticleContainer", AliLog::kDebug+7);
  //AliLog::SetClassDebugLevel("AliJetContainer", AliLog::kDebug+7);
  // END TEMP

  // Jet finding
  /*const AliJetContainer::EJetAlgo_t jetAlgorithm = AliJetContainer::antikt_algorithm;
  const Double_t jetRadius = 0.2;
  const AliJetContainer::EJetType_t jetType = AliJetContainer::kFullJet;
  const Double_t minTrackPt = 3.0;
  const Double_t minClusterPt = 3.0;
  const AliJetContainer::ERecoScheme_t recoScheme = AliJetContainer::pt_scheme;
  const char * label = "Jet";
  const Double_t minJetPt = 0; // TEMP: Consider at 1!
  const Bool_t lockTask = kTRUE;
  const Bool_t fillGhosts = kFALSE;

  if (newCorrectionFramework)
  {
    AliEmcalJetTask *pFullJet02TaskNew = AliEmcalJetTask::AddTaskEmcalJet(newFrameworkTracks.Data(), newFrameworkClusters.Data(),
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
    AliEmcalJetTask *pFullJet02Task = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault",
            jetAlgorithm, jetRadius, jetType, minTrackPt, minClusterPt, kGhostArea, recoScheme, label, minJetPt, lockTask, fillGhosts);
    pFullJet02Task->SelectCollisionCandidates(kPhysSel);
    pFullJet02Task->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  //////////////////////////////////////////
  // Sample task for QA information
  //////////////////////////////////////////
  #ifdef __CLING__
  std::stringstream sampleTask;
  sampleTask << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C(";
  sampleTask << newFrameworkTracks.Data() << ", ";
  sampleTask << newFrameworkClusters.Data() << ", ";
  sampleTask << newFrameworkCells.Data() << ");";
  std::cout << "Calling sample taks with " << sampleTask.str().c_str() << std::endl;
  AliAnalysisTaskEmcalJetSample * sampleTaskNew = reinterpret_cast<AliTaskCDBconnect *>(gROOT->ProcessLine(sampleTask.str().c_str()));
  #else
  AliAnalysisTaskEmcalJetSample * sampleTaskNew = AddTaskEmcalJetSample(newFrameworkTracks.Data(), newFrameworkClusters.Data(), newFrameworkCells.Data());
  #endif
  if (newFrameworkTracks != "tracks" && newFrameworkTracks != "Tracks")
  {
    // Remove wrong particle container
    sampleTaskNew->RemoveParticleContainer(0);
    // Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
    AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
    sampleTaskNew->AdoptTrackContainer(newTracks);
  }
  sampleTaskNew->GetClusterContainer(0)->SetClusECut(3.);
  sampleTaskNew->GetParticleContainer(0)->SetParticlePtCut(3);
  sampleTaskNew->SetHistoBins(600, 0, 300);
  sampleTaskNew->SelectCollisionCandidates(kPhysSel);

  AliJetContainer* jetCont02 = sampleTaskNew->AddJetContainer(jetType, jetAlgorithm, recoScheme, jetRadius, AliEmcalJet::kEMCALfid);*/

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
      esdChain << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/CreateAODChain.C(";
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
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromPYTHIA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromAOD.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoMass.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalQGTagging.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSample.C");
}

void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode)
{
  Int_t maxFilesPerWorker = 4;
  Int_t workerTTL = 7200;
  // LHC11a
  // /alice/data/2011/LHC11a/000146860/ESDs/pass4_with_SDD/AOD113/0002
  // Run list from EMC_pp
  const char* runNumbers = "146860 146859 146858 146856 146824 146817 146807 146806 146805 146804 146803 146802 146801 146748 146747 146746"
  const char* pattern = "ESDs/pass4_with_SDD/AOD113/*/AliAOD.root";
  const char* gridDir = "/alice/data/2011/LHC11a";
  // LHC11h
  // /alice/data/2011/LHC11h_2/000167693/ESDs/pass2/AOD145
  //const char* runNumbers = "167693";
  //const char* pattern = "ESDs/pass2/AOD145/*/AliAOD.root";
  //const char* gridDir = "/alice/data/2011/LHC11h_2";
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
  plugin->SetAliPhysicsVersion("vAN-20170103-1");

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
