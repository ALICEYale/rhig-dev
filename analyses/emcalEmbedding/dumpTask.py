#!/usr/bin/env python
import ROOT

fIn = ROOT.TFile("train.root")

# Get the analysis manager
jetShapes = fIn.Get("JetShapesAnalysis")

# Retrieve new task
jetTask = jetShapes.GetTask("Jet_new_KTChargedR020_tracks_pT0150_E_scheme")
if not jetTask:
    jetTask = jetShapes.GetTask("Jet_old_KTChargedR020_PicoTracks_pT0150_E_scheme")

# Dump the task
jetTask.Dump()
