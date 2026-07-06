// =============================================================================
// drawRecoTruthAP.C
// -----------------------------------------------------------------------------
// Compare reconstructed AP variables with truth AP variables stored in pairTree.
//
// Usage:
//   root -l -b -q 'drawRecoTruthAP.C("v0_candidates.root")'
//   root -l -b -q 'drawRecoTruthAP.C("v0_candidates.root","ap_kshort","parent_id1==parent_id2 && parent_pid==310")'
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

void drawRecoTruthAP(
    const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_v0_candidates_100k_zero_vtx_merged.root",
    const char* outprefix = "AP_reco_truth",
    const char* cuts = "charge1 != charge2 && parent_id1 == parent_id2 && parent_id1 != 0 && (abs(parent_pid) == 310 || abs(parent_pid) == 3122) && truth_alpha == truth_alpha && truth_qT == truth_qT")
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
    std::cerr << "[drawRecoTruthAP] cannot open " << infile << std::endl;
    return;
  }

  TTree* t = nullptr;
  f->GetObject("pairTree", t);
  if (!t)
  {
    std::cerr << "[drawRecoTruthAP] pairTree not found in " << infile << std::endl;
    f->Close();
    return;
  }

  TCut allCuts(cuts);

  auto* hReco = new TH2F("hAP_reco",
                         "Reconstructed AP;#alpha_{reco};q_{T,reco}^{+}  (GeV/c)",
                         240, -1.5, 1.5,
                         240, 0.0, 0.6);
  auto* hTruth = new TH2F("hAP_truth",
                          "Truth AP for same candidates;#alpha_{truth};q_{T,truth}^{+}  (GeV/c)",
                          240, -1.5, 1.5,
                          240, 0.0, 0.6);
  auto* hDelta = new TH2F("hAP_delta",
                          "Reco - truth;#Delta#alpha;#Delta q_{T}  (GeV/c)",
                          240, -0.1, 0.1,
                          240, -0.05, 0.05);

  auto* c = new TCanvas("c_reco_truth_ap", "Reco vs Truth AP", 1800, 650);
  c->Divide(3, 1);

  c->cd(1);
  t->Draw("qT:alpha>>hAP_reco", allCuts, "colz");
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();

  c->cd(2);
  t->Draw("truth_qT:truth_alpha>>hAP_truth", allCuts, "colz");
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();

  c->cd(3);
  t->Draw("delta_qT:delta_alpha>>hAP_delta", allCuts, "colz");
  gPad->SetRightMargin(0.15);

  const TString pdf = TString(outprefix) + ".pdf";
  const TString root = TString(outprefix) + ".root";
  c->SaveAs(pdf);

  TFile fout(root, "RECREATE");
  hReco->Write();
  hTruth->Write();
  hDelta->Write();
  fout.Close();

  std::cout << "[drawRecoTruthAP] input entries: " << t->GetEntries() << std::endl;
  std::cout << "[drawRecoTruthAP] reco entries after cuts: " << hReco->GetEntries() << std::endl;
  std::cout << "[drawRecoTruthAP] truth entries after cuts: " << hTruth->GetEntries() << std::endl;
  std::cout << "[drawRecoTruthAP] wrote " << pdf << " and " << root << std::endl;

  f->Close();
}
