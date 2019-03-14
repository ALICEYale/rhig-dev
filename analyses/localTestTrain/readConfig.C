// Sample code on how to process an MLTrainDefinition.cfg file by hand.
// With it, you can read in the configuration, and then execute AddTask
// macros and their extended configurations.
//
// Can be used to check that __R_ADDTASK__ has been properly defined in the
// macro configuration.
//
// WARNING: This code is very simple, and therfore not very flexible. You
// may need to modify it slightly for your purposes.
//
// author: Raymond Ehlers, Yale University <raymond.ehlers@yale.edu>

#include <TROOT.h>
#include <TMacro.h>
#include <TObjArray.h>

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskCfg.h>
#include <AliAODInputHandler.h>

#include <AliAnalysisTaskEmcal.h>
#else
class AliAnalysisManager;
class AliAnalysisTaskCfg;
class AliAODInputHandler;

class AliAnalysisTaskEmcal;
#endif

void readConfig()
{
  // Read in the global variables so the configuration can be executed
  gROOT->Macro("globalvariables.C");
  // Read the MLTrainDefinition.cfg file
  TObjArray * tempArr = AliAnalysisTaskCfg::ExtractModulesFrom("MLTrainDefinition.cfg");
  // Get the configuration in the first entry. It follows the order of how wagons are defined in MLTrainDefinition
  AliAnalysisTaskCfg * cfg = dynamic_cast<AliAnalysisTaskCfg *>(tempArr->At(0));
  // Create the analysis manager and input handler to allow for the execution of the macro
  AliAnalysisManager * mgr = new AliAnalysisManager("tempMgr");
  AliAODInputHandler * pAODHandler = AliAnalysisTaskEmcal::AddAODHandler();
  // Execute the AddTask macro
  cfg->ExecuteMacro();
  // Execute the configuration of the wagon (where you would use __R_ADDTASK__)
  //cfg->ExecuteConfigMacro();
  cfg->GetConfigMacro()->Print();
}
