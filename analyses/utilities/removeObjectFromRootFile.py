#!/usr/bin/env python

import ROOT
import os

def removeElement(elementName, inNewCorrectionTask = False):
    # See: https://stackoverflow.com/a/3964691
    for root, dirs, files in os.walk("data"):
        for f in files:
            if f.endswith(".root"):
                print("Filename: {0}".format(os.path.join(root,f)))
                fIn = ROOT.TFile(os.path.join(root, f), "UPDATE")
                if inNewCorrectionTask:
                    correctionTaskListName = "AliEmcalCorrectionTask_histos"
                    correctionTask = fIn.Get(correctionTaskListName)
                    print("Correction task GetName(): {0}".format(correctionTask.GetName()))
                    print("Correction task: {0}".format(", ".join([el.GetName() for el in correctionTask])))
                    print("Deleting \"{0}\" from AliEmcalCorrectionTask list".format(elementName))
                    correctionTask.Remove(correctionTask.FindObject(elementName))
                    print("Correction task GetName(): {0}, Class(): {1}".format(correctionTask.GetName(), correctionTask.Class()))
                    print("Correction task: {0}".format(", ".join([el.GetName() for el in correctionTask])))

                    # Remove old list from file and add new list
                    fIn.Delete(correctionTaskListName + ";1")
                    #ROOT.gDirectory.Delete(correctionTaskListName)
                    # kSingleKey ensures that it is written as a list rather than each individual hist
                    correctionTask.Write(correctionTask.GetName(), ROOT.TObject.kSingleKey)
                else:
                    fIn.ls()
                    print("Deleting \"{0}\"".format(elementName))
                    fIn.Delete(elementName)
                    ROOT.gDirectory.Delete(elementName)
                    fIn.ls()
                fIn.Close()
                #element = fIn.Get(elementName)
                #print(os.path.join(root, f))

if __name__ == "__main__":
    #removeElement("HadCorr_tracks_caloClusters;*")
    #removeElement("AliEmcalCorrectionClusterHadronicCorrection", inNewCorrectionTask = True)

    # Sample tasks
    # Old
    removeElement("AliAnalysisTaskEmcalJetSample_tracks_caloClusters_emcalCells_histos;*")
    # New
    removeElement("AliAnalysisTaskEmcalJetSample_tracksNew_caloClustersNew_emcalCellsNew_histos;*")
