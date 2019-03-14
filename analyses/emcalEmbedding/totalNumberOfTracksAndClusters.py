
from __future__ import print_function
import ROOT

def getNumberFromHist(sampleTask, histName):
    hist = sampleTask.FindObject(histName)

    val = 0;
    #print("{0}".format(histName))
    for i in xrange(1, hist.GetXaxis().GetNbins()):
        #print("binNumber: {0}, binCenter: {1}, tempVal: {2}".format(i, hist.GetBinCenter(i), hist.GetBinContent(i)*hist.GetBinCenter(i)))
        val += hist.GetBinContent(i)*hist.GetBinCenter(i)

    return val

def getNumber():
    fIn = ROOT.TFile("AnalysisResults.root", "READ")

    sampleTask = fIn.Get("AliAnalysisTaskEmcalJetSample_tracks_caloClustersCombined_emcalCellsCombined_hybridJets_histos")

    nClusters = getNumberFromHist(sampleTask, "fHistSumNClusters")
    nTracks = getNumberFromHist(sampleTask, "fHistSumNTracks")

    print("nClusters: {0}, nTracks: {1}, n*m: {2}".format(nClusters, nTracks, nClusters*nTracks))

if __name__ == "__main__":
    getNumber()
