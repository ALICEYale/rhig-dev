# IMPORTANT

__THIS GUIDE HAS BEEN SUPERCEDED BY THE OFFICIAL DOCUMENTATION__. See the [main documentation](http://alidoc.cern.ch/AliPhysics/master/_r_e_a_d_m_eemc_corrections.html) and [information on switching](http://alidoc.cern.ch/AliPhysics/master/_r_e_a_d_m_eemc_corrections_change.html)!

Note that the run macro and example configuration are still fine to use!

# Testing the new EMCal correction framework

Thanks for helping out! These instructions are on how to setup and run the new corrections framework. To test it, you should perform a local test using your analysis task. We recommend you do this using the run macro we provide you, which runs the new and old frameworks side-by-side. If you prefer, you alternatively can use your own run macro (it will need a few more modifications outlined below). The entire process shouldn't take too long - likely less than an hour. If you have any questions at any point, please don't hesitate to email us (raymond.ehlers@yale.edu, james.mulligan@yale.edu).

## Prerequisites

You need to be using the "new" EMCal framework developed by Salvatore. See [here](http://alidoc.cern.ch/AliPhysics/master/_r_e_a_d_m_echangefw.html) for instructions to update, if you need. 

## Get the code

To get the main correction code, pull the most recent AliPhysics master and compile. As noted above, we have prepared a run macro and sample configuration for you. The run macro is available [here](https://gitlab.cern.ch/ALICEYale/alice-yale-dev/raw/master/analyses/emcalCorrections/codeForTesters/runEMCalCorrections.C) (you may need to right click -> Save Link As..), and the sample configuration is available [here](https://gitlab.cern.ch/ALICEYale/alice-yale-dev/raw/master/analyses/emcalCorrections/codeForTesters/userTestConfiguration.yaml) (you may need to right click -> Save Link As..).

If you are curious about the code changes (you can use it just fine without looking!), you can see the code at ``$ALICE_PHYSICS/PWG/EMCAL/EMCALtasks``. The steering class is ``AliEmcalCorrectionTask``, and the individual corrections inherit from ``AliEmcalCorrectionComponent``, and are labeled ``AliEmcalCorrectionXXXX``. The default configuration file is at ``$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalConfiguration.yaml``.

## Configure the run macro

We need to configure the run macro to include the correction task and perform some configurations. If you use our sample run macro, some of these steps will be done for you, and will be noted as such. 

### Add the CorrectionTask (skip if you use our run macro)

If you are using your own run macro, we need to enable the new corrections -- add the following lines **before** the EMCal tender task!

```c++
AliEmcalCorrectionTask * correctionTask = AddTaskEmcalCorrectionTask();
correctionTask->SelectCollisionCandidates(kPhysSel);
// Set the run period, same as the track container
// If you derived from the file "runEMCalJetSampleTask.C", then it is likely stored under "sRunPeriod.Data()"
correctionTask->SetRunPeriod("LHC11h");
// Create a copy of cells, clusters, and tracks to compare against the current correction framework
correctionTask->SetCreateNewObjectBranches(true);
// Set the user configuration file
correctionTask->SetUserConfigurationFilename("userTestConfiguration.yaml");
// Initialize the configuration files in the correction task
correctionTask->InitializeConfiguration();
```

Don't forget to also load the macro with:

```c++
gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
```

### Configure your task to run twice

Configure your analysis task to run twice - once with the current corrections and once with the new corrections. To do so, configure your task once with different cluster and track (if needed) containers. Note that if you use jets, you need to run the jet finder twice - once for the old framework, and once for the new framework. For example, if you run on AODs, and set your cluster and track containers in run macro using "usedefault", then it would look something like:

```c++
AliAnalysisTaskMyAnalysis * task =
new AliAnalysisTaskMyAnalysis("usedefault", // clusters
                              "usedefault", // tracks
                              WhateverOptionsYouNormallyUse);
task->OptionOne();
// etc...

AliAnalysisTaskMyAnalysis * taskNew =
new AliAnalysisTaskMyAnalysis("caloClustersNew", // clusters
                              "tracksNew", // tracks
                              WhateverOptionsYouNormallyUse);
taskNew->OptionOne();
// etc...
```

**Note**: If you don't use the "usedefault" pattern, then you'll have to set the names manually. You should just set the tracks and clusters names as usual.

#### Special note on the track containers

As noted above, if you use jets, you will need to configure two jet finders. Due to using nonstandard track and cluster names to compare the two frameworks, the jet finder will not be configured properly if you just pass the names. Instead, you will need to remove the particle container and add a new track container based on the "tracksNew" branch. In code, it will look like,

```c++
// Remove wrong particle container
pFullJet02TaskNew->RemoveParticleContainer(0);
// Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
// Set the cuts for the new track container
newTracks->SetParticlePtCut(minTrackPt);
// Add it to the jet finder
pFullJet02TaskNew->AdoptTrackContainer(newTracks);
```

If you don't want to worry about the details, we suggest using the run macro defined for testing the new framework. It will take care of configuring the jet finders properly.

If your analysis uses a track container, and you use the "usedefault" pattern in your analysis task AddTask, repeat the above lines (adapted for your task) for your task's track container.

## Configure the new correction framework

In the new framework, the corrections will be configured using a YAML file -- see the default at ``$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalConfiguration.yaml``. YAML is a simple, human readable format often used for configuration. If you are interested or if it is unclear, read more about it [here](https://en.wikipedia.org/wiki/YAML). In order to test the new framework against the current framework, we need you to configure the YAML user file (**not the default file!**) to use the same set of parameter values that you currently set in your AddTasks.

First, ensure that the file `userTestConfiguration.yaml` (which was selected above) is in the same directory as you run macro. This will be the basis of your user file.

You can map the current corrections to the new ones with the following table. Note that each new correction is preceded by the name `AliEmcalCorrection`:

| Current Correction Name | New Correction Name       |
| ----------------------- | ------------------------- |
| Tender                  | CellBadChannel            |
|                         | CellEnergy                |
|                         | CellBadChannel            |
| Clusterizer             | Clusterizer               |
| ClusterMaker            | ClusterExotics            |
|                         | ClusterNonLinearity       |
| ClusTrackMatcher        | ClusterTrackMatcher       |
| HadCorr                 | ClusterHadronicCorrection |

Using the table above, we need to make sure that the current and the new corrections are configured the same. Compare the settings in your run macro with those in ``$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalConfiguration.yaml``, known as the default file. In particular, check that each parameter value in the YAML file matches the value for that variable in your run macro (and not the other way around -- it is okay if there is no matching field for every field in the AddTask). The variable names should often be the same. If values are different in the YAML file, then modify the values in **``userTestConfiguration.yaml``**, not in ``AliEmcalConfiguration.yaml``! To make the change, write it in yourself or copy the relevant structure of the default file to your user file and modify there. For instance, in the default file, the ``w0`` parameter in the Clusterizer is to ``4.5``. If I wanted to change that in the user file, I would add:

```
Clusterizer:
    w0: 3 # Changed this value
```

Any setting that you have in the user file will override the default file!

If you have some parameters with the same value and you want to change them in unison (say, if you wanted to change the track name, or a min pt cut), you can set a parameter in the ``sharedParameters`` section, and then set the value in the component to that parameter. For example, if I wanted to set the min pt in the ``ClusterExotics`` and ``ClusterNonLinearity``, I would have:

```
sharedParameters:
    ptMin: 0.15
ClusterExotics:
    clusterPtMin: "sharedParameters:ptMin"
ClusterNonLinearity:
    clusterPtMin: "sharedParameters:ptMin"
```

This is just provided for your convenience!

If there are certain corrections you don't use in your analysis (e.g. ClusterTrackMatcher), and you are using our suggested run macro, then you should remove the unnecessary task from the run macro, and disable it in ``userTestConfiguration.yaml``. 

## Run your run macro

You are all set! Now run the run macro as normal. Once it is finished, we have a python script to compare the output histograms automatically. It is available <a href="https://gitlab.cern.ch/ALICEYale/alice-yale-dev/raw/master/analyses/utilities/compareHistos.py" download>here</a> (you may need to right click -> Save Link As..). To run it, you need to pass the path of the ``AnalysisResults.root`` file, as well as the name of the output list from your task that was generated with the new corrections. For example,

```bash
python compareHistos.py -f ../exampleDir/AnalysisResults.root -n MyAnalysisOutputListFromNewCorrections
```

If you don't know the name of your output list, just pass an invalid name to the ``-n`` argument and it will show you the available options:

```bash
python compareHistos.py -f ../exampleDir/AnalysisResults.root -n invalidName
```

In exceptional cases, you will need to pass the name corresponding to the name of the output list from your task that was run with the current corrections. If needed, the script will inform you. For example,

```bash
python compareHistos.py -f ../exampleDir/AnalysisResults.root -n MyAnalysisOutputListFromNewCorrections -o MyAnalysisOutputListFromOldCorrections
```

If the outputs do not match initially, please double check your settings! If they still do not match, then please let us know!

Help for the script is available with `python compareHistos.py --help`. (If you are using aliBuild, you'll need to set your ``$PYTHONPATH`` variable. In bash, you can do this by setting ``export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH``. Other shells may vary.)

Alternatively, you can just send us your ``AnalysisResults.root`` along with your run macro and user configuration. We can analyze them ourselves and report the results back.
