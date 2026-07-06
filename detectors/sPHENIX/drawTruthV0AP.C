// =============================================================================
// drawTruthV0AP.C
// -----------------------------------------------------------------------------
// Draw generator-level AP plots from truthV0Tree written by TpcTruthV0DecayTree.
//
// Usage:
//   root -l -b -q 'drawTruthV0AP.C("truth_v0_decays.root")'
//   root -l -b -q 'drawTruthV0AP.C("truth_v0_decays.root","truth_kshort","parent_pid == 310")'
// =============================================================================

#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>

void drawTruthV0AP(
    const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/truth_v0_decays_100k_zero_vtx.root",
    const char* outprefix = "AP_truth_v0",
    const char* cuts = "parent_pid == 310 || abs(parent_pid) == 3122")
{
  const char* style = "/sphenix/user/dcxchenxi/develope/macros/TrackingProduction/sPhenixStyle.C";
  if (!gSystem->AccessPathName(style))
  {
    gROOT->LoadMacro(style);
    gROOT->ProcessLine("SetsPhenixStyle();");
  }

  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "[drawTruthV0AP] cannot open " << infile << std::endl;
    return;
  }

  TTree* t = nullptr;
  f->GetObject("truthV0Tree", t);
  if (!t)
  {
    std::cerr << "[drawTruthV0AP] truthV0Tree not found in " << infile << std::endl;
    f->Close();
    return;
  }

  auto* hAP = new TH2F("hAP_truthV0",
                       "Generator-level V0 AP;#alpha;q_{T}^{+}  (GeV/c)",
                       240, -1.5, 1.5,
                       240, 0.0, 0.6);

  auto* c1 = new TCanvas("c_truth_v0_ap", "Generator-level V0 AP", 1100, 800);
  c1->cd();
  t->Draw("qT:alpha>>hAP_truthV0", TCut(cuts), "colz");
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();

  const TString pdf = TString(outprefix) + ".pdf";
  const TString root = TString(outprefix) + ".root";
  c1->SaveAs(pdf);

  TFile fout(root, "RECREATE");
  hAP->Write();
  fout.Close();

  std::cout << "[drawTruthV0AP] truthV0Tree entries: " << t->GetEntries() << std::endl;
  std::cout << "[drawTruthV0AP] histogram entries after cuts: " << hAP->GetEntries() << std::endl;
  std::cout << "[drawTruthV0AP] wrote " << pdf << " and " << root << std::endl;

  f->Close();
}
