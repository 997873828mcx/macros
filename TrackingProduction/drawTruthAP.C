// =============================================================================
// drawTruthAP.C
// -----------------------------------------------------------------------------
// Draw an Armenteros-Podolanski plot from the pairTree written by
// make_truth_ap.py.
//
// Usage:
//   root -l -q 'drawTruthAP.C("ap_truthpoints.root")'
//   root -l -q 'drawTruthAP.C("ap_truthpoints.root","ap.pdf","charge1 != charge2 && qT < 0.3")'
// =============================================================================

#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <string>

void drawTruthAP(const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/ap_truthpoints_10evt_truepv.root",
                 const char* outfile = "AP_truthpoints.pdf",
                 const char* cuts = "charge1 != charge2")
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
    std::cerr << "[drawTruthAP] cannot open " << infile << std::endl;
    return;
  }

  TTree* t = nullptr;
  f->GetObject("pairTree", t);
  if (!t)
  {
    std::cerr << "[drawTruthAP] pairTree not found in " << infile << std::endl;
    f->Close();
    return;
  }

  auto* c1 = new TCanvas("c_truth_ap", "Truth-point Armenteros-Podolanski", 1100, 800);
  c1->cd();

  auto* hAP = new TH2F("h_ArmPod_TruthPoints",
                       "Armenteros-Podolanski from ideal TPC truth points;#alpha;q_{T}^{+}  (GeV/c)",
                       240, -1.5, 1.5,
                       240, 0.0, 0.6);

  t->Draw("qT:alpha>>h_ArmPod_TruthPoints", TCut(cuts), "colz");
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  c1->SaveAs(outfile);

  std::string root_out = outfile;
  const auto dot = root_out.find_last_of('.');
  if (dot != std::string::npos)
  {
    root_out = root_out.substr(0, dot);
  }
  root_out += ".root";

  TFile fout(root_out.c_str(), "RECREATE");
  hAP->Write();
  fout.Close();

  std::cout << "[drawTruthAP] pairTree entries: " << t->GetEntries() << std::endl;
  std::cout << "[drawTruthAP] histogram entries after cuts: " << hAP->GetEntries() << std::endl;
  std::cout << "[drawTruthAP] wrote " << outfile << " and " << root_out << std::endl;

  f->Close();
}
