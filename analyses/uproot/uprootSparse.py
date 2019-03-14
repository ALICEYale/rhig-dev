#!/usr/bin/env python

# Tests for reading a sparse using uproot.
#
# author: Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
# date: 18 July 2018

import numpy as np
import collections
import ROOT

import uproot

def createSparseWithRoot(filename):
    """ Create a THnSparseF using ROOT. We fill some axis we are interested in (2, 4, 5) and leave others empty. """
    # namedtuple is just for convenience
    sparseAxis = collections.namedtuple("sparseAxis", ["nBins", "min", "max"])
    ignoredAxis   = sparseAxis(nBins =  1, min =   0.0, max =  1.0)
    selectedAxis1 = sparseAxis(nBins = 10, min =   0.0, max = 20.0)
    selectedAxis2 = sparseAxis(nBins = 20, min = -10.0, max = 10.0)
    selectedAxis3 = sparseAxis(nBins = 30, min =   0.0, max = 30.0)
    # We want to select axes 2, 4, 5
    axes = [ignoredAxis, ignoredAxis, selectedAxis1, ignoredAxis, selectedAxis2, selectedAxis3, ignoredAxis]

    # Create the actual sparse
    # dtype is required here for the sparse to be created successfully.
    bins = np.array([el.nBins for el in axes], dtype=np.int32)
    mins = np.array([el.min for el in axes])
    maxes = np.array([el.max for el in axes])
    print("bins: {}, mins: {}, maxs: {}".format(bins, mins, maxes))
    sparse = ROOT.THnSparseF("testSparse", "testSparse", len(axes),
                bins, mins, maxes)

    # Fill a few values in the axes of interest
    arr1 = np.array([0, 0, 10, 0,  0, 10, 0], dtype = np.float64)
    val = sparse.Fill(arr1)
    print("val: {}, Nbins: {}, arr: {}".format(val, sparse.GetNbins(), arr1))
    sparse.PrintBin(val, "")

    arr2 = np.array([0, 1,  5, 0,  -5, 20, 0], dtype = np.float64)
    val = sparse.Fill(arr2)
    print("val: {}, Nbins: {}, arr: {}".format(val, sparse.GetNbins(), arr2))
    sparse.PrintBin(val, "")

    arr3 = np.array([0, 0,  5, 0,  2, 20, 0], dtype = np.float64)
    val = sparse.Fill(arr3)
    print("val: {}, Nbins: {}, arr: {}".format(val, sparse.GetNbins(), arr3))
    sparse.PrintBin(val, "")

    arr4 = np.array([0, 0,  5, 0, -5, 20, 0], dtype = np.float64)
    val = sparse.Fill(arr4)
    print("val: {}, Nbins: {}, arr: {}".format(val, sparse.GetNbins(), arr4))
    sparse.PrintBin(val, "")

    print("Writing test sparse to filename {}.".format(filename))
    f = ROOT.TFile(filename, "RECREATE")
    sparse.Write()
    f.Close()

def readWithUproot(filename):
    """ Read the sparse with uproot, printing the available attributes. """
    print("Test reading sparse with uproot.")
    f = uproot.open(filename)
    sparse = f["testSparse"]
    # We expect to see the members.
    print("THnSparse: {!r}".format(sparse))
    print("THnSparse available attributes: {}".format(dir(sparse)))

    # Check for which members are available and whether the values are correct.
    # These currently fail
    # For example, the TObjArray of axes.
    assert "fAxes" in dir(sparse)
    # The dimensionality should match the number of entries in the axes array.
    assert sparse.fAxes.GetEntries() == sparse.fNdimensions

def readWithRoot(filename):
    """ Reading back the sparse with ROOT to check the validity of the assertions. """
    f = ROOT.TFile(filename, "READ")
    sparse = f.Get("testSparse")

    # We expect to see the members.
    print("THnSparse: {!r}".format(sparse))
    print("THnSparse available attributes: {}".format(dir(sparse)))

    # Check for which members are available and whether the values are correct.
    # These currently fail
    # For example, the TObjArray of axes.
    #assert "fAxes" in dir(sparse)
    # The dimensionality should match the number of entries in the axes array.
    assert sparse.GetListOfAxes().GetEntries() == sparse.GetNdimensions()

if __name__ == "__main__":
    filename = "testSparseThreeEntry.root"
    createSparseWithRoot(filename = filename)
    readWithRoot(filename = filename)
    readWithUproot(filename = filename)

