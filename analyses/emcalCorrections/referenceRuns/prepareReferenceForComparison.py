#!/usr/bin/env python

# Takes a reference file and appends "_reference" to each list name so that
# it doesn't conflict with newly created file which will be compared against this
# reference. In terms of the comparison script, this reference file wiil take the place
# of the "old" corrections file, while the new file that will be compared against the
# reference will continue its place as the "new" file.
#
# author: Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University

import argparse
import os
from pachyderm import histogram
import ROOT
import ruamel.yaml

def rewriteWithDifferentName(element, suffix, elementNamesMap):
    """ Rewrite an element from an open file with a different name.
    It is the user responsibility to delete the element, if so desired.

    Because this is ROOT, it is fine that the file is implicitly open.
    """
    newElementName = "{}_{}".format(element.GetName(), suffix)
    # Same the name in the proper formatt to be used in the YAML map
    elementNamesMap[element.GetName()] = [newElementName]

    # Need to set the name (which will be seen when printing) to be the same as
    # the name that is written
    element.SetName(newElementName)
    element.Write(newElementName, ROOT.TObject.kSingleKey)

    return elementNamesMap

def renameElements(filename, suffix, debug):
    print("Renaming elemnts in filename: \"{}\"".format(filename))
    elementNamesMap = {}
    with histogram.RootOpen(filename, "UPDATE") as f:
        if debug:
            print("f.ls()  pre:")
            f.ls()
        # Cannot just iterate over GetListOfKeys because the hash list is updated when an element is added to a file...
        # Instead, we copy the keys so the iterable is not updated (this should really be treated as a ROOT bug...)
        keys = f.GetListOfKeys()
        simpleListOfKeys = []
        for key in keys:
            simpleListOfKeys.append(key.GetName())

        # Loop over the available keys. If it is the correction task, then unpack the component lists
        for key in simpleListOfKeys:
            if key.endswith("_" + suffix):
                print("Skipping the processing of element {} since it has already been processed!".format(key))
                # Same the name in the proper formatt to be used in the YAML map
                elementNamesMap[key.replace("_" + suffix, "")] = [key]
                # Don't apply the suffix twice
                continue

            element = f.Get(key)
            print("Processing element: {}".format(element.GetName()))
            # Remove the existing element from the file
            f.Delete(element.GetName() + ";*")
            # Rewrite the name
            if "AliEmcalCorrectionTask" in element.GetName():
                for component in element:
                    elementNamesMap = rewriteWithDifferentName(component, suffix, elementNamesMap)
            else:
                elementNamesMap = rewriteWithDifferentName(element, suffix, elementNamesMap)

        if debug:
            print("keys: {}".format(keys.GetEntries()))
            print("f.ls() post:")
            f.ls()

    # Save the map for user with the comparison script
    # Write the reference map to the same directory as the reference file
    yamlReferenceMapLocation = os.path.join(os.path.dirname(filename), "referenceMap.yaml")
    print("Writing yaml reference map to \"{}\"".format(yamlReferenceMapLocation))
    y = ruamel.yaml.YAML()
    with open(yamlReferenceMapLocation, "w") as f:
        y.dump(elementNamesMap, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare root file for use as a reference with the comparison script.")
    parser.add_argument("-f", "--filename", action="store",
                        metavar="filename",
                        default="AnalysisResults.root",
                        help="Path of AnalysisResults.root file.")
    parser.add_argument("-s", "--suffix", action="store",
                        metavar="suffix",
                        default="reference",
                        help="Suffix to append to the elements in the reference file.")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="Print additional debug messages")

    # Parse the arguments
    args = parser.parse_args()

    renameElements(filename = args.filename, suffix = args.suffix, debug = args.debug)
