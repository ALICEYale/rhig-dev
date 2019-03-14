#!/usr/bin/env python

import rootpy
import rootpy.io

def checkSpectra():
    fSeparate = rootpy.io.File.Open("AnalysisResults.separate.root", "READ")
    fCombined = rootpy.io.File.Open("AnalysisResults.combined.root", "READ")

    fSepTracksList = fSeparate.AliAnalysisTaskEmcalJetSample_tracks_caloClustersCombined_emcalCellsCombined_hybridJets_histos.FindObject("tracks")
    fSepTracksEmbList = fSeparate.AliAnalysisTaskEmcalJetSample_tracks_caloClustersCombined_emcalCellsCombined_hybridJets_histos.FindObject("tracks_Emb")
    fCombinedList = fCombined.AliAnalysisTaskEmcalJetSample_tracks_caloClustersCombined_emcalCellsCombined_hybridJets_histos.FindObject("tracks")

    print("Starting comparison...")

    if len(fSepTracksList) != len(fSepTracksEmbList) or len(fSepTracksList) != len(fCombinedList):
        print("ERROR: List length mismatch: sepTracks: {0}, sepTracksEmb: {1}, combined: {2}".format(len(fSepTracksList), len(fSepTracksEmbList), len(fCombinedList)))

    count = 0
    for i, (track, trackEmb, combined) in enumerate(zip(fSepTracksList, fSepTracksEmbList, fCombinedList)):
        if track.GetName() != combined.GetName():
            print("ERROR: Not comparing the same hists! Have sep {0} and combined {1}".format(track.GetName(), combined.GetName()))
            count += 1
            continue

        sepHist = track.Clone("{0}_{1}".format(track.GetName(), "combined"))
        #print("Pre-Add entries: {0} + {1} = {2}".format(track.GetEntries(), trackEmb.GetEntries(), sepHist.GetEntries()))
        sepHist.Add(trackEmb)
        #print("PostAdd entries: {0}".format(sepHist.GetEntries()))
        
        if sepHist.GetEntries() != combined.GetEntries():
            print("ERROR: Number of entries mismatch between hists {0} (nEntries:{1}) and {2} (nEntries: {3})".format(sepHist.GetName(), sepHist.GetEntries(), combined.GetName(), combined.GetEntries()))
            count += 1
            continue

        sepHist.Add(combined, -1)

        if sepHist.GetEntries() != 0:
            print("ERROR: Mismatch for hist {0} with entries {1}".format(sepHist.GetName(), sepHist.GetEntries()))
            count += 1
            continue

    print("Completed comparison. Number of mismatches: {0}/{1}".format(count, len(fSepTracksList)))

if __name__ == "__main__":
    checkSpectra()
