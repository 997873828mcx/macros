// Scatter plot: deiondx_kev_per_cm vs ionized_electrons_dedx_kev_per_cm
// Tree: driftStepQA in merged_driftstepqa.root

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"

void plot_dedx_scatter(
    const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/iondedx_gain1400_sig004/completed/merged_driftstepqa.root",
    const char* outfile = "dedx_scatter.pdf")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile *f = TFile::Open(infile);
  if (!f || f->IsZombie()) {
    printf("ERROR: cannot open %s\n", infile);
    return;
  }

  TTree *t = (TTree*)f->Get("driftStepQA");
  if (!t) {
    printf("ERROR: tree 'driftStepQA' not found\n");
    f->Close();
    return;
  }

  // 2D histogram: x = ionized_electrons_dedx_kev_per_cm, y = deiondx_kev_per_cm
  // Axis ranges cover ~99th percentile of each variable
  TH2F *h = new TH2F("h_dedx_scatter",
                      ";Ionized electrons dE/dx (keV/cm);dE_{ion}/dx (keV/cm)",
                      500, 0, 20,
                      500, 0, 900);

  t->Draw("deiondx_kev_per_cm:ionized_electrons_dedx_kev_per_cm>>h_dedx_scatter",
          "deiondx_kev_per_cm>0 && ionized_electrons_dedx_kev_per_cm>0",
          "goff");

  TCanvas *c = new TCanvas("c", "dE/dx scatter", 900, 700);
  c->SetLogz();
  c->SetRightMargin(0.13);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleOffset(1.1);

  h->Draw("COLZ");

  c->SaveAs(outfile);
  printf("Saved: %s\n", outfile);

  f->Close();
}
