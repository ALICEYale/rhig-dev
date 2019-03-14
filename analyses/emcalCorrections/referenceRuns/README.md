# Comparing Against a Reference

To run the comparison script against a reference, a YAML file is needed to define the mapping. In this mode, the reference file takes the place of the "old" corrections, while the file to compare stays as the "new" file. For convenience, a script was created to run the preparation and comparison automatically. It can be run with:

```
$ ./compareToReference.sh (new filename) (reference filename)
```

where `(new filename)` defaults to `AnalysisResults.root` in the `emcalCorrections` directory and `(reference filename)` defaults to `AnalysisResults.root` in the `referenceRuns` directory if the variables are not passed.

You are also welcome to run the preparation and comparison scripts by hand. Note that the option `--overrideCorreationsMap` is most likely necessary, as the YAML file will probably remap the standard map from the new correction tasks to their corresponding reference tasks. For example, it will change the map for `AliEmcalCorrectionClusterHadronicCorrection` from `AliHadCorrTask` to `AliEmcalCorrectionClusterHadronicCorrection_reference`.

# Prepare a Reference File

It is best to prepare a reference with only the results from the new correction framework. One it is created using `runEMCalCorrections.C`, then complete the following steps:

1. Copy to the `referenceRuns` directory. Be certain to tag the file with the date in the filename!
2. Symlink `AnalysisResults.root` to the new reference file.

If you are going to commit it into the get repository, a few additional steps are recommended.

1. Run `prepareReferenceForComparison.py` at least once so the reference will be prepared before being committed.
2. Delete the previous reference.
3. Commit the new reference file (may need to use `-f` with `git add` since `.root` files are usually ignored).
