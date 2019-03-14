# Test train instructions

Make sure to initialize an alien token with `alien-token-init`!

For ROOT6, see below.

Put all of these files in the same directory.

1. Download `runTest.sh` and `generate.C` from [here](https://twiki.cern.ch/twiki/bin/view/ALICE/AnalysisTrains#Local_train_test). Be certain to get the most recent versions.

2. Find the train that you are interested in and open the link to `testing output dir`. 

3. Download the necessary files:
    - In the `config` dir:
         - `MLTrainDefinition.cfg` - Contains the wagons used for the train.
         - `env.sh` - Contains the environmental variables necessary for execution.
         - `globalvariables.C` - Contains the global variables defined for that train.
         - `handlers.C` - Sets up the file type (ie ESD, AOD, MC) handler.

4. All other necessary files should be created by `generate.C`! If it says that it fails to find something like `lego_train.sh`, it is because `generate.C` failed. Check the `generation.log` file for details!

5. To use data local files instead of downloading them again from AliEn, the filename needs to be of a particular form which depends somewhat on the dataset. An example for LHC11h is included here (starting with ``__alice__``). Note that you'll need to adapt it for your local system.

## ROOT6

Markus provided a brief explanation and a set of files needed for using ROOT6 with the test train. See his
email in the chain "root6 test train". The `generate.C` and `handlers.C` files needed are sent there and are
also stored in this directory.

## Checking `MLTrainDefinition.cfg` by Hand

The configuration is read by ``AliAnalysisTaskCfg::ExtractModulesFrom("MLTrainDefinition.cfg")``. It returns a ``TObjArray`` of ``AliAnalysisTaskCfg``, which each one corresponding to a particular entry in the MLTrainDefinition. You can look at this interactively by using the approach shown in ``readConfig.C``.
