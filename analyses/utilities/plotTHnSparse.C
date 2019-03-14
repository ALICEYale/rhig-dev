// plotTHnSparse.C
// Small macro to plot THnSparse and its error propagation
//
// There is a bug in ROOT (has since been fixed ROOT-4453, but is not in ALICE's ROOT),
// in which Sumw2() does not work when called on THnSparse after filling.
//
// Author: James Mulligan <james.mulligan@yale.edu>

#include "THnSparse.h"
#include "TFile.h"

//______________________________________________________________________________
void plotTHnSparse()
{
  
  TFile* f = new TFile("myFile.root", "RECREATE");

  // Create a THnSparse and fill some content
  Int_t nbins[20] = {0};
  Double_t min[20] = {0.};
  Double_t max[20] = {0.};
  for(Int_t dim = 0; dim < 3; dim++) {
    nbins[dim] = 50;
    min[dim] = 0;
    max[dim] = 100;
  }

  THnSparse* thn = new THnSparseF("thn","thn",3,nbins,min,max);
  //thn->Sumw2();
  Double_t* content = new Double_t[3];
  for (Int_t i = 0; i < 100; i++) {
    content[0] = 10;
    content[1] = 20;
    content[2] = 30;
    thn->Fill(content);
  }
  
  printf("\nScale THnSparse by 10...");
  
  // Project to a TH1 and print the bin content and error for the bin we filled
  TH1* h1 = thn->Projection(0);
  h1->SetName("Proj Before Scaling");
  printf("\nbin content: %d", h1->GetBinContent(6));
  printf("\nbin error: %d", h1->GetBinError(6));
  
  // Scale the THnSparse
  thn->Sumw2();
  thn->Scale(10);
  
  // Project the scaled THnSparse to a TH1 and again print the bin content and error
  TH1* h2 = thn->Projection(0);
  h2->SetName("Proj After Scaling");
  printf("\nscaled bin content: %d", h2->GetBinContent(6));
  printf("\nscaled bin error: %d", h2->GetBinError(6));
  
  // Check that the same procedure for TH1 works...
  printf("\n\nDo the same procedure for a TH1...");
  TH1* h1Scaled = h1->Clone();
  h1Scaled->SetName("TH1 After Scaling");
  printf("\nbin content: %d", h1Scaled->GetBinContent(6));
  printf("\nbin error: %d", h1Scaled->GetBinError(6));
  h1Scaled->Sumw2();
  h1Scaled->Scale(10);
  printf("\nscaled bin content: %d", h1Scaled->GetBinContent(6));
  printf("\nscaled bin error: %d\n", h1Scaled->GetBinError(6));
  
  h1->Write();
  h2->Write();
  thn->Write();
  h1Scaled->Write();

}
