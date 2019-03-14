#!/usr/bin/env python

#
# Script to compare histograms in at least one .root file for exact or
# statistical equality. The comparison can be between EMCal corrections,
# as well as a user defined task. The user must set the input file.
# If desired, the user can also set the name of the output from an
# additional user defined task to compare. This name should correspond
# to the "new" task. In some cases, the user also needs to specify the
# name of the "old" task (the user will be notified in this case).
#
# NOTE: This script only checks for histograms in the task names that
#       are specified. Additional histograms may exist in "old"
#       components which will never be checked if they do not
#       correspond to "new" histograms
#
# Authors:
#   Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
#   James Mulligan <james.mulligan@yale.edu>, Yale University

# Python 3 compatibility
from __future__ import print_function
# Handle xrange (without using the future package to ensure that users don't need to install packages)
# See: https://stackoverflow.com/a/15015199
try:
    xrange
except NameError:
    xrange = range

# General
import os
import sys
# Identify the system
import socket
# Handle args
import argparse
# ROOT
import ROOT
# Convenience
import pprint

# Force root to ignore command line options...
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Set ROOT into batch mode (as if we ran "root -b")
ROOT.gROOT.SetBatch()

def extractHistsFromList(listOfObjects, outputDict, appendToName = "", correctionListName = "AliEmcalCorrectionTask_histos", hasNotRecursed = True):
    """ Extract histograms from TLists while still saving some structure to allow for comparison.

    Each entry in the dict contains a dict of histograms. The name used to store the histogram is
    a combination of the histogram name and the path to that histogram. A delimiter is used between
    each step in the path. This name is referred to as the "extended name". For a few examples,
    consider the following structure in a root file:

    AnalysisTask/
        fNEvents
        fHist1
        fList1/
            fHist2

    In the main dictionary, there will be an entry for "AnalysisTask", which contains a dict for that
    task. This dict corresponding to "AnalysisTask" will contain:

    >>> outputDict["AnalysisTask"]
    {
        "fNEvents" : fNEvents,
        "fHist1" : fHist1,
        "fList1:fHist2" : fHist2
    }

    (With additional spacing added to the print out for clarity).

    NOTE: ":" is used as the delimiter between one level and the next

    Args:
        listOfObjects (list): List containing all of the top level objects in the root file. Note that
            these must be the actual objects read in from the file - the name of the key is not enough!
        outputDict (dict): Dict where all of the extracted histograms should be stored.
        appendToName (str): Name that should be appended to the histogram name to created the extended
            name. This is primarily used when the function is called recursively. Default: Empty string
        correctionListName (str): Name of list that should be recursed over, but skipped as a dict name
            and as an entry in appendToName. This is used for the main correction task list, since it
            only contains the output lists of other tasks and therefore creates unnecessary structure.
            This could be trivially generalized to support more names if necessary.
            Default: "AliEmcalCorrectionTask_histos"
        hasNotRecursed (bool): True if we have not yet recursed. Generally used to denote that a new dict
            has not been created, and nothing should be appended to the name yet. Default: True

    Returns:
        All of the histograms are stored in outputDict
    """

    # Append to earlier name
    if appendToName != "":
        appendToName = "{0}:".format(appendToName)

    for obj in listOfObjects:
        #print("Processing object {0} of type {1}".format(obj.GetName(), obj.ClassName()))

        # Handle TH1 derived
        if obj.InheritsFrom(ROOT.TH1.Class()):
            #print("Adding object {0} with appendToName \"{1}\" of type {2} as a hist".format(obj.GetName(), appendToName, obj.ClassName()))
            addHistToDict(obj, obj.GetName(), appendToName, outputDict)

        # Handle THnSparse.
        # Same as TH1, but we first need to project the axes
        elif obj.InheritsFrom(ROOT.THnSparse.Class()):
            #print("Adding object {0} of type {1} as a THnSparse".format(obj.GetName(), obj.ClassName()))
            # Loop over the axes of the THnSparse, Project each one, and add them to the dict
            for i in xrange(0, obj.GetNdimensions()):
                histName = "{0}_Axis{1}".format(obj.GetName(), i)
                #print("Projecting axis {0} as {1} from THnSparse".format(i, histName))
                projection = obj.Projection(i)
                projection.SetName(histName)
                addHistToDict(projection, histName, appendToName, outputDict)

        # Handle TList -> Recursively
        elif obj.InheritsFrom(ROOT.TList.Class()):
            # We copy values here because we are in a for loop over all objects. If we changed the
            # assignment of values then they are changed for the rest of the for loop, which is not
            # the desired behavior. Instead, we want the value to be different only once it recurses!
            hasNotRecursedValueToPass = hasNotRecursed
            appendToNameToPass = appendToName
            dictToPass = outputDict

            # If the case of the first recursion, we want to create a new dict to store the hists
            # The exception is if the object is marked as one to skip. In that case, it is ignored.
            # This exception is used primarily for the main new correction list, as we don't want
            # to create a new list since it just adds an unnecessary layer.
            if not obj.GetName() == correctionListName:
                if hasNotRecursed == True:
                    # Things may break from here, so it is good to notify the user. However, this situation
                    # should be extremely rare.
                    if appendToName != "":
                        print("\n\nWARNING: appendToName not empty as is usually expected! Value: {0}".format(appendToName))

                    # We will now recurse with a new dict to hold the contents of the TList.
                    # Any additional TLists stored inside of this TList will be accounted for via appendToName.
                    hasNotRecursedValueToPass = False
                    if obj.GetName() in outputDict.keys():
                        print("\nERROR: Obj name {0} is already in the output dict (which contains keys: {1}. Check that the input files do not have any tasks that are named the same!)".format(obj.GetName(), outputDict.keys()))
                        sys.exit(3)
                    outputDict[obj.GetName()] = dict()
                    dictToPass = outputDict[obj.GetName()]
                else:
                    # Add the object name as we recurse past it to keep track of the structure
                    appendToNameToPass += obj.GetName()
            else:
                # We want to skip over the main new correction first list entry since it adds an unnecessary layer
                pass

            # Call recursively to handle TLists
            #print("Recursing with object {0} of type {1}".format(obj.GetName(), obj.ClassName()))
            extractHistsFromList(obj, dictToPass, appendToName = appendToNameToPass, correctionListName = correctionListName, hasNotRecursed = hasNotRecursedValueToPass)

        # TDirectoryFile - Recursive
        # Handle similarly to TLists
        elif obj.InheritsFrom(ROOT.TDirectoryFile.Class()):
            # TDirectoryFile isn't iterable, so we need to iterate the keys ourselves
            #print("Obj name: {0}".format(obj.GetName()))
            for key in obj.GetListOfKeys():
                #print("Obj key: {0}".format(key.GetName()))
                # We copy values here because we are in a for loop over all objects. If we changed the
                # assignment of values then they are changed for the rest of the for loop, which is not
                # the desired behavior. Instead, we want the value to be different only once it recurses!
                dictToPass = outputDict
                appendToNameToPass = appendToName
                hasNotRecursedValueToPass = hasNotRecursed

                # Same approach as is taken for the TList (however is it not worth refactoring yet - there
                # are a very large number of variables to manage for refactoring just a few lines of code).
                if hasNotRecursed == True:
                    # Things may break from here, so it is good to notify the user. However, this situation
                    # should be extremely rare.
                    if appendToName != "":
                        print("\n\nWARNING: appendToName not empty as is usually expected! Value: {0}".format(appendToName))

                    # We will now recurse with a new dict to hold the contents of the TList. Any TLists below
                    # here will be accounted for via appendToName.
                    hasNotRecursedValueToPass = False
                    if obj.GetName() in outputDict.keys():
                        print("\nERROR: Obj name {0} is already in the output dict (which contains keys: {1}. Check that the input files do not have any tasks that are named the same!)".format(obj.GetName(), outputDict.keys()))
                        sys.exit(3)
                    outputDict[key.GetName()] = dict()
                    dictToPass = outputDict[key.GetName()]
                else:
                    appendToNameToPass += key.GetName()

                # Call recursively to handle TDirectoryFile keys
                extractHistsFromList(key.ReadObj(), dictToPass, appendToName = appendToNameToPass, correctionListName = correctionListName, hasNotRecursed = hasNotRecursedValueToPass)

        # TTree
        # Project each branch as a 1D histogram
        # NOTE: Does not yet handle nested branches!
        elif obj.InheritsFrom(ROOT.TTree.Class()):
            # Setup output
            outputDict[obj.GetName()] = dict()
            dictToPass = outputDict[obj.GetName()]

            # Loop over branches
            for branch in obj.GetListOfLeaves():
                # Need to append here rather than in the addHistToDict function to ensure that the histograms
                # won't have the same names and therefore be overwritten!
                #histName = "{0}_{1}".format(obj.GetName(), branch.GetName())
                histName = "{0}".format(branch.GetName())
                #print("Projecting branch {0} as {1} from TTree".format(branch.GetName(), histName))
                obj.Project(histName, branch.GetName())
                # Get the resulting histogram
                hist = ROOT.gDirectory.Get(histName)
                #print("name: {0}, entries: {1}, hist: {2}".format(hist.GetName(), hist.GetEntries(), hist))
                addHistToDict(hist, histName, appendToName, dictToPass)

        # Unrecognized -> Throw warning and ignore
        else:
            print("Object {0} of type {1} is not of the type TH1, THnSparse, TList, or TDirectoryFile and cannot be handled. Please look into this more closely!".format(obj.GetName(), obj.ClassName()))

def addHistToDict(hist, histName, appendToName, outputDict):
    """ Simple function to add a hist to the input dict.

    It ensures that formatting with appendToName is correct and then stores it
    in the output dict. It works with any histogram type, including TH1 and THnSparse.
    """
    # Make sure that it stays around after closing the file
    hist.SetDirectory(0)
    outputName = "{0}{1}".format(appendToName, histName)
    #print("outputName: {0}".format(outputName))
    outputDict[outputName] = hist

def runComparison(inputFiles, newAnalysisList, oldAnalysisList, correctionListName, newToOldCorrections, normalize, saveHistsToRootFile, printPDF, doStatisticalComparison, thresh):
    """ Steering of the comparison.

    It handles all input, retrieving of histograms, and then doing the initial comparison of the histograms.
    """
    # Setup input objects
    # Holds the input hists
    inputHists = dict()
    # Contains the names of all of the new corrections that are available in the file(s)
    # NOTE: Only histograms that are in the components named in this list will actually be checked!
    newCorrectionsList = list()

    # Loop over all available files.
    # Checks are performed to ensure that objects are not overwritten
    for inputFile in inputFiles:
        # Open the file
        f = ROOT.TFile(inputFile)

        # Get a list of the names all of the new corrections that are present in the input file
        # First check for the correction components in the input
        correctionInput = f.Get(correctionListName)
        # If loading multiple files, it is possible that the new correction output is not in every file
        # Normally, this would be expected to return as an empty list and not be iterated over, but ROOT doesn't make this easy...
        if correctionInput:
            for obj in correctionInput:
                if obj.GetName() not in newCorrectionsList:
                    newCorrectionsList.append(obj.GetName())
                else:
                    print("\nERROR: Attempted to add new correction component {0}, but it is already in the new corrections list. Check that the input files do not have multiple sets of the new correction components! Exiting!".format(obj.GetName()))
                    sys.exit(2)

        # Now check for the rest of the possible tasks
        keys = f.GetListOfKeys()
        for key in keys:
            if key.GetName() in newToOldCorrections.keys():
                if key.GetName() not in newCorrectionsList:
                    newCorrectionsList.append(key.GetName())
                else:
                    print("\nERROR: Attempted to add task {0}, but it is already in the list to check. Check that the input files do not have multiple sets of the task! Exiting!".format(key.GetName()))
                    sys.exit(2)

        # Extract objects from keys to extract and organize hists
        print("\nLoading histograms from file {0}...".format(inputFile))
        inputList = list()
        for key in keys:
            inputList.append(key.ReadObj())

        # Extract histograms from the file and organize into a dict
        extractHistsFromList(inputList, inputHists, correctionListName = correctionListName)

        f.Close()

    # Add user analysis task input list name
    if newAnalysisList != "" and not newAnalysisList in newCorrectionsList:
        newCorrectionsList.append(newAnalysisList)

    # Handle the old analysis name if the new analysis name is meaningful (ie not an empty string)
    if newAnalysisList != "":
        # Guess old analysis name if none was passed
        if oldAnalysisList == "":
            # Often tasks name just differ by changing the tracks and clusters name.
            # The first letter was left off because it can vary for ESDs and AODs, but it is unlikely to match
            # anything else.
            oldAnalysisList = newAnalysisList.replace("racksNew", "racks").replace("aloClustersNew", "aloClusters").replace("ellsNew", "ells")
            print("Guessing old anaylsis output list name: \"{0}\"".format(oldAnalysisList))

        # Check if old and new names are the same (for instance, if the replace didn't work)
        # Without this check, it will appear to work properly, but will simply compare with itself,
        # which should trivially agree
        if oldAnalysisList == newAnalysisList:
            print("\nCould not successfully guess old analysis name! Please pass a name for the analysis with the old framework")
            print("Available names include:")
            print("\t{0}".format("\n\t".join(inputHists.keys())))
            sys.exit(1)

        # Attempt to match new and old analysis
        if oldAnalysisList in inputHists:
            # Store the match for use later
            newToOldCorrections[newAnalysisList] = [oldAnalysisList]
            print("\nMatched new analysis name \"{0}\" to old analysis name \"{1}\"".format(newAnalysisList, oldAnalysisList))
        else:
            print("\nCould not find old analysis list name \"{0}\" in output. Please pass a name for the analysis with the old framework if the guess was wrong or confirm that it is the right value.".format(oldAnalysisList))
            print("Available names include:")
            print("\t{0}".format("\n\t".join(inputHists.keys())))
            sys.exit(1)
    else:
        print("No new analysis list name passed. Only analyzing the new vs old corrections!")

    # Print status if necessary. These are used for debugging
    #print("InputHists: ")
    #pprint.pprint(inputHists)
    #print("newCorrectionsList: ")
    #pprint.pprint(newCorrectionsList)
    #print("newToOldCorrections: ")
    #pprint.pprint(newToOldCorrections)

    # Ability to save out the histograms if desired
    if saveHistsToRootFile:
        fOut = ROOT.TFile("compareHistos.root", "RECREATE")
        for lists in inputHists:
            for k, hist in inputHists[lists].iteritems():
                hist.Write()

        fOut.Close()

    # Print the status for the user
    print("\nProcessing with:")
    print("Tasks using the new correction framework: ")
    print("\t{0}".format("\n\t".join(newCorrectionsList)))
    print("Tasks using the old correction framework: ")
    # Prints all old corrections that are going to be matched based on which new corrections are available
    # See: https://stackoverflow.com/a/952952
    print("\t{0}".format("\n\t".join(set(
        [oldCorrection for newCorrection in newCorrectionsList for oldCorrection in newToOldCorrections.get(newCorrection, [])]
        ))
    ))
    print("")
    if doStatisticalComparison:
        print("This may take a minute...")
        print("")

    # Loop through component histograms, and find matching histogram in old correction framework
    nNewCorrectionHists = 0
    nMatched = 0
    nUnmatched = 0
    nCorrect = 0
    nEmpty = 0
    nCouldNotFindOldComponent = 0
    # Start with each new component
    for newComponent in newCorrectionsList:
        # Iterate over new hists in a new component. The histograms are initially compared by name under which
        # they are stored in the dict. This "extended name" (defined more precisely in extractHistsFromList())
        # includes both the histogram name, as well as the path to where the histogram was found.
        #
        # If the dict corresponding to the new component is not found, then it is skipped 
        for newHistExtendedName in inputHists.get(newComponent, dict()):
            #print("Searching for {2} ({0}) in {1}".format(inputHists[newComponent][newHistExtendedName], newComponent, newHistExtendedName))

            # The old hist may have a slightly different name, so we attempt to handle it here.
            # Possibilities include:
            #  - Removing "New" from cells, clusters, and tracks names.
            #  - Replacing "new" with "old"
            #  - Replacing "New: with "Old"
            #  - Replacing "tracks" with "PicoTracks"
            newHistExtendedNameComparisonName = newHistExtendedName.replace("racksNew", "racks").replace("aloClustersNew", "aloClusters").replace("ellsNew", "ells")
            newHistExtendedNamePossibilities = [newHistExtendedNameComparisonName,
                                               newHistExtendedNameComparisonName.replace("new", "old"),
                                               newHistExtendedNameComparisonName.replace("New", "Old"),
                                               newHistExtendedNameComparisonName.replace("new", "old").replace("tracks", "PicoTracks"),
                                               newHistExtendedNameComparisonName.replace("tracks", "PicoTracks")]

            # Setup for looping over the new hist
            nNewCorrectionHists += 1
            matchFound = False
            couldNotFindOldComponent = False

            # Get each old components histograms associated with the new component by the new to old corrections map
            # In principle there could be multiple components, although in practice there is usually just one
            for oldComponent in newToOldCorrections[newComponent]:
                # Attempt to get the old correction hists.
                # First attempt ESD names instead of AOD names if the AOD name doesn't exist
                if oldComponent not in inputHists:
                    #print("Attempting ESD name instead")
                    oldComponent = oldComponent.replace("caloClusters", "CaloClusters").replace("tracks", "Tracks").replace("emcalCells", "EMCALCells")

                # Only continue if we can find an old component. Note to the user if we cannot
                if oldComponent in inputHists:
                    # Loop over possible new hist extended names
                    for newHistExtendedNameComparisonName in newHistExtendedNamePossibilities:
                        # newHistExtendedNameComparisonName is the name that should be used for comparing against
                        # the old histograms, since it includes additional guesses to match the old hist extended
                        # name.
                        # newHistExtendedName is the actual name of the new histogram in the dict. 
                        #print("Trying name: {0}:{1}".format(oldComponent, newHistExtendedNameComparisonName))
                        if newHistExtendedNameComparisonName in inputHists[oldComponent]:
                            #print("Found match")
                            # Get hists
                            newHist = inputHists[newComponent][newHistExtendedName]
                            # newHistExtendedNameComparisonName must be the name of the object we are looking for
                            # so we can use it directly to retrieve the old hist
                            oldHist = inputHists[oldComponent][newHistExtendedNameComparisonName]

                            # Normalize both histograms for the comparison if desired
                            if normalize:
                                newHistScaleFactor = newHist.Integral()
                                oldHistScaleFactor = oldHist.Integral()
                                if newHistScaleFactor > 0:
                                    newHist.Scale(1/newHistScaleFactor)
                                if oldHistScaleFactor > 0:
                                    oldHist.Scale(1/oldHistScaleFactor)

                            # An additional sanity check that we have the right hist
                            if newHist.GetName() != oldHist.GetName():
                                print("WARNING: Hist names do not match!!\nNew: {0}, {1}\nOld: {2}, {3}".format(newHistExtendedName, newHist.GetName(), newHistExtendedNameComparisonName, oldHist.GetName()))

                            # Note that we found the match
                            matchFound = True
                            nMatched = nMatched + 1

                            # Compare the two histograms together
                            nCorrect = nCorrect + compare(newHist, oldHist, printPDF = printPDF, newComponentName = newComponent, newHistExtendedNameComparisonName = newHistExtendedNameComparisonName, doStatisticalComparison = doStatisticalComparison, thresh = thresh)

                            # Note empty histograms
                            if newHist.GetEntries() < 1e-4 and oldHist.GetEntries() < 1e-4:
                                nEmpty = nEmpty + 1
                            break
                else:
                    # Note that we could not find the old component
                    print("WARNING: For new hist {0}:{1}, could not find old component {2} in inputHists! Skipping!".format(newComponent, newHistExtendedName, newToOldCorrections[newComponent]))
                    couldNotFindOldComponent = True
                    nCouldNotFindOldComponent = nCouldNotFindOldComponent+ 1

            # Note unmatched hists
            if matchFound is False and couldNotFindOldComponent is False:
                nUnmatched = nUnmatched + 1
                print("WARNING: For new hist {0}:{1}, no matching histogram was found! Skipping!".format(newComponent, newHistExtendedName))

    # Print the final summary
    print("\nNumber of histograms in new correction components: {0}".format(nNewCorrectionHists))
    print("Number of matched histograms: {0}".format(nMatched))
    print("Number of unmatched histograms due to not finding old component: {0}".format(nCouldNotFindOldComponent))
    print("Number of unmatched histograms: {0}".format(nUnmatched))
    print("Number of correct histograms: {0}".format(nCorrect))
    print("Number of incorrect histograms: {0}".format(nMatched - nCorrect))
    print("Number of empty histograms: {0}".format(nEmpty))

def compare(h1, h2, printPDF, newComponentName, newHistExtendedNameComparisonName, doStatisticalComparison, thresh):
    """ Compare two histograms to see if they are the same.

    It can do two types of comparisons:
    - Subtract one from the other and compare the difference.
    - Statistical bin-by-bin comparison. See the official documentation and the code for more information.

    Includes options to print the histograms that disagree.
    """
    h1Entries = h1.GetEntries()
    h2Entries = h2.GetEntries()
    # Use this to compare so we can plot both later if desired
    comparisonHist = h1.Clone(h1.GetName() + "_Subtracted")

    # Work around for a ROOT6 bug, which does not work properly with the Add(h2, -1) function when projected from a TTree
    # Below works properly for both ROOT5 and ROOT6. It replaces the contents of comparisonHist, so the result is the same
    # as was done previously - it just uses a new method.
    subtractionResult = comparisonHist.Add(h1, h2, 1, -1)

    if not subtractionResult:
        print("WARNING: Subtraction was not successful! Comparison of histogram {0} is not valid!".format(h1.GetName()))

    match = 1 # flag for whether histograms agree (1) or not (0)
    epsilon = 1e-5
    divByZero = False # flag for whether fractional difference is infinite for some bin

    # Do bin-by-bin statistical comparison, if requested
    if doStatisticalComparison:
        histFractionalDiff = h1.Clone(h1.GetName() + "_FractionalDifference")
        histFractionalDiff.Add(h1, -1) 
        
        # TH1 will always apply
        nBins = h1.GetXaxis().GetNbins() + 2
        # Handle TH2
        if h1.InheritsFrom(ROOT.TH2.Class()):
            nBins *= (h1.GetYaxis().GetNbins() + 2)
        # Handle TH3 (since TH3 doesn't inherit from TH2, but does from TH1)
        if h1.InheritsFrom(ROOT.TH3.Class()):
            nBins *= ((h1.GetYaxis().GetNbins() + 2) * (h1.GetZaxis().GetNbins() + 2))

        # Can only check with GetNcells() on ROOT 6
        if ROOT.gROOT.GetVersionInt() > 60000:
            if nBins != h1.GetNcells():
                print("WARNING! Bin number mismatch!! nBins: {0}, nCells: {1}".format(nBins, h1.GetNcells()))

        for bin in xrange(0, nBins):
            if h1.IsBinUnderflow(bin) or h1.IsBinOverflow(bin): continue

            subtractedVal = comparisonHist.GetBinContent(bin)
            if abs(subtractedVal) < epsilon: continue
            
            originalVal = h1.GetBinContent(bin)
            #print ("bin: " + str(bin) + " / sub: " + str(subtractedVal) + " / orig: " + str(originalVal))
            
            if abs(originalVal) < epsilon: # avoid division by zero
                fractionalDiff = 0
                print("WARNING: bin %d fractional difference contains division by 0!" % bin)
                divByZero = True
            else:
                fractionalDiff = subtractedVal/originalVal
                histFractionalDiff.SetBinContent(bin, fractionalDiff)
                #print("fractional diff: " + str(fractionalDiff))
            
        # print out that histograms disagree if fractional difference exceeds threshold
        maxFracDiff = histFractionalDiff.GetMaximum()
        minFracDiff = histFractionalDiff.GetMinimum()
        if abs(minFracDiff) > maxFracDiff: maxFracDiff = abs(minFracDiff)
        if abs(maxFracDiff) > epsilon:
            print("maximum disagreement of " + h1.GetName() + ": " + str(maxFracDiff))

        if maxFracDiff > thresh:
            print("Histograms {0} in component {1} mismatch!".format(h1.GetName(), newComponentName))
            print("Old entries: {0}, New entries: {1}".format(h2Entries, h1Entries))
            print("max fractional disagreement = {0}".format(maxFracDiff))
            if "area" in h1.GetName().lower():
                print("\tWARNING: \"area\" is in the mismatched histogram name. If this is related to jet area, note that jet area is not reliable to calculate from one run to the next. Consequently, this difference may be ignorable.")
            print("Print the histogram with the \"-p\" option to investigate more closely.")
            match = 0

    # Otherwise, do exact comparison
    else:
        maxValue = comparisonHist.GetBinContent(comparisonHist.GetMaximumBin())
        minValue = comparisonHist.GetBinContent(comparisonHist.GetMinimumBin())
        if abs(maxValue) > epsilon or abs(minValue) > 0 or h1Entries != h2Entries:
            print("Histograms {0} in component {1} mismatch!".format(h1.GetName(), newComponentName))
            print("Old entries: {0}, New entries: {1}".format(h2Entries, h1Entries))
            print("subtracted max = {0}".format(maxValue))
            print("subtracted min = {0}".format(minValue))
            if "area" in h1.GetName().lower():
                print("\tWARNING: \"area\" is in the mismatched histogram name. If this is related to jet area, note that jet area is not reliable to calculate from one run to the next. Consequently, this difference may be ignorable.")
            print("Print the histogram with the \"-p\" option to investigate more closely.")
            match = 0

    # Print histograms, if desired
    if  match == 0 and printPDF:
        ROOT.gStyle.SetOptStat(1)
        # Create canvas
        canvas = ROOT.TCanvas("canvas", "canvas")
        canvas.Divide(2)
        
        # Set options
        # Should be the same for both
        drawOptions = ""
        if h1.InheritsFrom(ROOT.TH2.Class()):
            drawOptions = "colz"
            rangeAxisH2 = h2.GetZaxis()
        elif h1.InheritsFrom(ROOT.TH1.Class()):
            rangeAxisH2 = h2.GetYaxis()
            
        # Get maximum and minimum from h1
        minVal = h1.GetMinimum()
        maxVal = h1.GetMaximum()
        # Set the limits on h2 so that they are the same
        rangeAxisH2.SetRangeUser(minVal, maxVal)

        # Draw
        canvas.cd(1)
        h1.Draw(drawOptions)
        canvas.cd(2)
        h2.Draw(drawOptions)
        
        # Label old and new corrections
        textOld = ROOT.TLatex()
        textOld.SetNDC()
        canvas.cd(1)
        textOld.DrawLatex(0.30,0.92,"New corrections")
        textOld = ROOT.TLatex()
        textOld.SetNDC()
        canvas.cd(2)
        textOld.DrawLatex(0.30,0.92,"Old corrections")

        # Save hists
        outputFilename = os.path.join("output", newComponentName + "__" + newHistExtendedNameComparisonName.replace(":","_") + ".pdf")
        canvas.SaveAs(outputFilename)

        # Draw the subtracted histograms and the fractional difference histograms
        ROOT.gStyle.SetOptStat(0)
        canvasDiff = ROOT.TCanvas("canvasSubtr", "canvasSubtr")
        canvasDiff.Divide(2)
        
        minValSubtr = comparisonHist.GetMinimum()
        maxValSubtr = comparisonHist.GetMaximum()
        rangeAxisH2.SetRangeUser(minValSubtr, maxValSubtr)
        canvasDiff.cd(1)
        comparisonHist.Draw(drawOptions)

        # Label subtracted panel
        textSubtr = ROOT.TLatex()
        textSubtr.SetNDC()
        canvasDiff.cd(1)
        textSubtr.DrawLatex(0.27,0.92,"Subtracted Difference")
        
        if doStatisticalComparison:
          minValFrac = histFractionalDiff.GetMinimum()
          maxValFrac = histFractionalDiff.GetMaximum()
          rangeAxisH2.SetRangeUser(minValFrac, maxValFrac)
          
          if h1.InheritsFrom(ROOT.TH2.Class()):
            rangeAxisFrac = histFractionalDiff.GetZaxis()
          elif h1.InheritsFrom(ROOT.TH1.Class()):
            rangeAxisFrac = histFractionalDiff.GetYaxis()
          rangeAxisFrac.SetTitle("Fractional Difference")
          
          canvasDiff.cd(2)
          histFractionalDiff.Draw(drawOptions)

          # Label fractional difference panel
          textFracDiff = ROOT.TLatex()
          textFracDiff.SetNDC()
          canvasDiff.cd(2)
          textFracDiff.DrawLatex(0.27,0.92,"Fractional Difference")
          if divByZero:
            textDivZero = ROOT.TLatex()
            textDivZero.SetNDC()
            canvasDiff.cd(2)
            textDivZero.DrawLatex(0.39,0.40,"Warning: Div by zero!")
            textDivZero.DrawLatex(0.39,0.36,"See stdout.")

        outputFilenameDiff = os.path.join("output", newComponentName + "__" + h1.GetName() + "_Diff" + ".pdf")
        canvasDiff.SaveAs(outputFilenameDiff)
            
    return match

if __name__ == '__main__':
    # Define arguments
    parser = argparse.ArgumentParser(description="Compare histograms to test the new EMCal corrections framework.")
    parser.add_argument("-f", "--inputFiles", action="store",
                        nargs="*", metavar="filename",
                        default=["AnalysisResults.root"],
                        help="Path of AnalysisResults.root file. Can be more than one file, although the task names cannot be overlap.")
    parser.add_argument("-y", "--yamlNewToOldDictFile", action="store",
                        type=str, metavar="yamlFile",
                        default = None,
                        help="Name of the YAML file which contains additions to the new to old tasks dictionary.")
    parser.add_argument("--overrideCorrectionsMap", action="store_true",
                        help="Allow the YAML file to override the correction framework map.")
    parser.add_argument("-n", "--newAnalysisList", action="store",
                        type=str, metavar="newAnalysisList",
                        default="",
                        help="Name of the list generated by your task run using the new EMCal corrections framework.")
    parser.add_argument("-o", "--oldAnalysisList", action="store",
                        type=str, metavar="oldAnalysisList",
                        default="",
                        help="Name of the list generated by your task run using the old EMCal corrections framework.")
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize all histograms to unity for comparison.")
    parser.add_argument("-r","--saveHistsToRootFile", action="store_true",
                        help="Save out the histograms used in the comparison to a ROOT file.")
    parser.add_argument("-p", "--printPDF", action="store_true",
                        help="Print a pdf of the difference between any mismatched histograms.")
    parser.add_argument("-s", "--doStatisticalComparison", action="store_true", 
                        help="Compare for statistical agreement (rather than exact agreement).")
    parser.add_argument("-t", "--thresh", action="store",
                        type=float, metavar="thresh", default=0.003,
                        help="Set threshold for statistical comparison -- fractional difference.")

    # Parse the arguments
    args = parser.parse_args()

    # User configurable options:
    # AliEmcalCorrections list name
    correctionListName = "AliEmcalCorrectionTask_histos"
    #correctionListName = "AliEmcalCorrectionTask_embed_histos"

    # Set tracks, clusters, cells names
    # (They usually don't have to change. Standard AOD and ESD branches are handled automatically)
    cellsName = "emcalCells"
    clustersName = "caloClusters"
    tracksName = "tracks"
    #clustersName = "EmcCaloClusters"
    #tracksName = "AODFilterTracks"

    # Map from the new corrections to the old corrections to know where to look for a matching histogram
    # Some useful, but nonstandard options are commented out. If they are to be used, be certain to commment
    # out the correspond standard options!
    newToOldCorrections = dict()
    newToOldCorrections["AliEmcalCorrectionCellEnergy"] = ["emcal_tender_plots"]
    newToOldCorrections["AliEmcalCorrectionCellBadChannel"] = ["emcal_tender_plots"]
    newToOldCorrections["AliEmcalCorrectionCellTimeCalib"] = ["emcal_tender_plots"]
    newToOldCorrections["AliEmcalCorrectionClusterizer"] = []
    newToOldCorrections["AliEmcalCorrectionClusterExotics"] = ["EmcalClusterMaker_{clusters}".format(clusters = clustersName)]
    newToOldCorrections["AliEmcalCorrectionClusterNonLinearity"] = ["EmcalClusterMaker_{clusters}".format(clusters = clustersName)]
    #newToOldCorrections["AliEmcalCorrectionClusterExotics"] = ["EmcalClusterMaker_{clusters}".format(clusters = "caloClusters_EmcCaloClusters")]
    #newToOldCorrections["AliEmcalCorrectionClusterNonLinearity"] = ["EmcalClusterMaker_{clusters}".format(clusters = "caloClusters_EmcCaloClusters")]
    newToOldCorrections["AliEmcalCorrectionClusterTrackMatcher"] = ["ClusTrackMatcher_{tracks}_{clusters}_histos".format(clusters = clustersName, tracks = tracksName)]
    newToOldCorrections["AliEmcalCorrectionClusterHadronicCorrection"] = ["HadCorr_{tracks}_{clusters}".format(clusters = clustersName, tracks = tracksName)]
    #newToOldCorrections["AliEmcalCorrectionClusterTrackMatcher_embed"] = ["ClusTrackMatcher_{tracks}_{clusters}_histos".format(clusters = clustersName, tracks = tracksName)]
    #newToOldCorrections["AliEmcalCorrectionClusterHadronicCorrection_embed"] = ["HadCorr_{tracks}_{clusters}".format(clusters = clustersName, tracks = tracksName)]

    # Add newToOldCorrections entires from a YAML file for convenience
    if args.yamlNewToOldDictFile:
        # Import YAML here so that it is not a full dependency.
        try:
            import yaml
        except:
            print("Need yaml installed!")
            sys.exit(1)
        with open(args.yamlNewToOldDictFile, "r") as f:
            yamlTaskCorrespondence = yaml.load(f)

        for key, val in yamlTaskCorrespondence.iteritems():
            if not key in newToOldCorrections.keys() or args.overrideCorrectionsMap:
                if args.overrideCorrectionsMap and key in newToOldCorrections.keys():
                    print("WARNING: Overriding key {0} from {1} to {2}".format(key, newToOldCorrections[key], val))
                newToOldCorrections[key] = val
            else:
                print("WARNING: Key {0} is already defined in newToOldCorrections! Exiting!".format(key))
                sys.exit(1)

    # Parse input file list
    inputFiles = args.inputFiles

    # Convenience for rehlers
    if "ray-MBP" in socket.gethostname() or "ray-desktop" in socket.gethostname():
        inputFiles = [ "../emcalCorrections/{0}".format(f) for f in inputFiles ]

    # Print settings
    print("Comparing histograms with:")
    print("inputFile(s): \n\t- {0}".format("\n\t- ".join(inputFiles)))
    print("newAnalysisList: \"{0}\"".format(args.newAnalysisList))
    if args.oldAnalysisList != "":
        print("oldAnalysisList: \"{0}\"".format(args.oldAnalysisList))
    if args.yamlNewToOldDictFile:
        print("yamlNewToOldCorrectionsFile: \"{0}\"".format(args.yamlNewToOldDictFile))

    # Check whether the input files exist
    for f in inputFiles:
        if not os.path.exists(f):
            print("File \"{0}\" does not exist! Exiting!".format(f))
            sys.exit(1)

    # Create necessary output directory
    if args.printPDF and not os.path.exists("output"):
        os.makedirs("output")

    # Run the actual comparison
    runComparison(inputFiles = inputFiles, newAnalysisList = args.newAnalysisList, oldAnalysisList = args.oldAnalysisList, correctionListName = correctionListName, newToOldCorrections = newToOldCorrections, normalize = args.normalize, saveHistsToRootFile = args.saveHistsToRootFile, printPDF = args.printPDF, doStatisticalComparison = args.doStatisticalComparison, thresh = args.thresh)
