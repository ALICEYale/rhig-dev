#!/usr/bin/env python

import ROOT

def extractTree():
    tags = [ "new", "old" ]
    filename = "../emcalEmbedding/AnalysisResults.{0}.constSub.pythia.root";
    treeName = "JetQGTaggings_Jet_{0}_AKTChargedR020_PicoTracks_pT0150_E_schemeConstSub_TCRawConstSubTree_ConstSub_Incl"

    for tag in tags:
        fIn = ROOT.TFile(filename.format(tag), "READ")
        print("Tree name: {0}".format(treeName.format(tag)))
        tree = fIn.Get(treeName.format(tag))
        print("Tree entries: {0}".format(tree.GetEntries()))

        fOut = ROOT.TFile("{0}.root".format(tag), "CREATE")
        fOut.cd()
        tree.Write()

        # Close files
        fOut.ls()
        fOut.Close()
        fIn.Close()

if __name__ == "__main__":
    extractTree()

