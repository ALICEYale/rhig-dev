Running on the GRID via a run macro
===================================

This guide describes a number of steps when using a run macro to submit a grid job

Setup
=====

- Ensure that you have run ``alien-token-init``

Settings
========

- Set iStartAnalysis = 2 and remember to set cGridMode = "full" (although it is probably best to run a test)
- Set the pattern and grid dir for the path to the files on the grid!
- Set whether running on MC in CreateAlienHandler()
- Set the AliPhysics version
- Set how to merge everything. If you want run-by-run, then comment out (?):

  ```
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);
  ```

