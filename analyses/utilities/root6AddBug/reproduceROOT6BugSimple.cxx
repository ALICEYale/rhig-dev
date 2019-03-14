#include <vector>
#include <string>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1F.h>

TTree * generateTree(std::string treeName);

void reproduceROOT6BugSimple()
{
  std::vector <TH1F*> hists;
  TFile * fIn = nullptr;
  TTree * tree = nullptr;

  // Setup tree and project
  std::string treeName = "tree";
  tree = generateTree(treeName);

  // Project
  tree->Project("weightPythia", "weightPythia");

  // Get the histogram
  TH1F * hist = dynamic_cast<TH1F *>(gDirectory->Get("weightPythia"));
  if (!hist) {
    Printf("Invalid pointer for hist");
    return;
  }
  Printf("hist: %p", hist);
  hists.push_back(hist);
  hist->SetDirectory(0);

  // Create the second histogram
  hists.push_back(dynamic_cast<TH1F *>(hists.at(0)->Clone("tempCloneTest")));
  hists.at(1)->SetDirectory(0);

  // Ensure that we have both histograms
  for (auto hist : hists)
  {
    std::cout << "histName: " << hist->GetName() << " entries: " << hist->GetEntries() << " address: " << hist << "\n";
  }

  // Create a clone to use for subtraction
  TH1F * clone = dynamic_cast<TH1F*>(hists.at(0)->Clone("testClone"));
  clone->Add(hists.at(1), -1);

  // Should be 0!
  std::cout << "Clone entries (should be 0): " << clone->GetEntries() << "\n";
  clone->Draw();

  // Try ading the two
  TH1F * tempHist = dynamic_cast<TH1F*>(hists.at(0)->Clone("testClone2"));
  tempHist->Add(hists.at(0), hists.at(1), 1, -1);
  std::cout << "Adding the two after cloning hist entries (should be 0): " << tempHist->GetEntries() << "\n";
  //tempHist->Draw();
}

TTree * generateTree(std::string treeName)
{
  TTree * tree = new TTree(treeName.c_str(), treeName.c_str());

  Double_t weightPythia = 1;
  TBranch * branch = tree->Branch("weightPythia", &weightPythia);

  for (int i = 0; i < 22; i++)
  {
    // Just use 1
    weightPythia = 1;
    tree->Fill();
  }

  // Avoid crashes
  tree->ResetBranchAddress(branch);

  return tree;
}
