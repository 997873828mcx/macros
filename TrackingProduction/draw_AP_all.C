#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

void draw_AP_all(const char* infile = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/myKShortReco/inReconstruction/ap.root",
                 const char* histname = "h_AP_all",
                 const char* outfile = "ArmenterosPodolanski_all.pdf")
{
  // Open the file and retrieve the 2D-histogram
  TFile f(infile, "READ");
  if (!f.IsOpen()) {
    fprintf(stderr, "ERROR: could not open %s\n", infile);
    return;
  }
  TH2D* h = dynamic_cast<TH2D*>(f.Get(histname));
  if (!h) {
    fprintf(stderr, "ERROR: no TH2D called %s in %s\n", histname, infile);
    return;
  }

  // Optional: adjust style
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.2, "Z");

  // Create canvas
  TCanvas* c = new TCanvas("c","Armenteros-Podolanski (all)", 800, 600);
  c->SetRightMargin(0.15);      // give room for color palette
  c->SetLogz();                 // turn on log-scale for Z axis

  // Draw
  h->GetXaxis()->SetTitle("#alpha = (p^{+}_{L} - p^{-}_{L}) / (p^{+}_{L} + p^{-}_{L})");
  h->GetYaxis()->SetTitle("p_{T}^{+} [GeV/c]");
  h->Draw("COLZ");

  // Save
  c->SaveAs(outfile);

  // clean up
  delete c;
  f.Close();
}