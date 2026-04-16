#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>

void DrawRecoPtTrackID(
    const char* infile = "output/pionplus_pt10/completed/merged_reco_pt.root",
    const unsigned int selected_track_id = 0,
    const int nbins = 120,
    const double pt_min = 0.0,
    const double pt_max = 20.0,
    const double fill_pt_min = 7.0,
    const double fill_pt_max = 14.0)
{
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "Failed to open input file: " << infile << std::endl;
    return;
  }

  TTree* reco_tracks = dynamic_cast<TTree*>(fin->Get("reco_tracks"));
  if (!reco_tracks)
  {
    std::cerr << "Could not find TTree 'reco_tracks' in " << infile << std::endl;
    fin->Close();
    return;
  }

  TH1D* h_pt = new TH1D(
      "h_pt_track_id",
      Form("Reco p_{T} distribution for track_id = %u, %.2f < p_{T} < %.2f;Reco p_{T} [GeV/c];Counts",
           selected_track_id, fill_pt_min, fill_pt_max),
      nbins, pt_min, pt_max);
  h_pt->SetDirectory(nullptr);
  h_pt->SetLineWidth(2);

  Float_t pt = 0.0;
  UInt_t track_id = 0;

  reco_tracks->SetBranchStatus("*", 0);
  reco_tracks->SetBranchStatus("pt", 1);
  reco_tracks->SetBranchStatus("track_id", 1);
  reco_tracks->SetBranchAddress("pt", &pt);
  reco_tracks->SetBranchAddress("track_id", &track_id);

  Long64_t nfilled = 0;
  const Long64_t nentries = reco_tracks->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    reco_tracks->GetEntry(i);
    if (track_id != selected_track_id)
    {
      continue;
    }
    if (!std::isfinite(pt))
    {
      continue;
    }
    if (pt < fill_pt_min || pt > fill_pt_max)
    {
      continue;
    }
    h_pt->Fill(pt);
    ++nfilled;
  }

  TCanvas* c1 = new TCanvas("c_reco_pt_track_id", "Reco pT", 900, 700);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);
  h_pt->Draw("hist");

  TF1* f_gaus = nullptr;
  if (nfilled > 0)
  {
    f_gaus = new TF1(
        "f_gaus_track_id",
        "gaus",
        fill_pt_min,
        fill_pt_max);
    f_gaus->SetLineColor(2);
    f_gaus->SetLineWidth(2);
    f_gaus->SetParameters(h_pt->GetMaximum(), h_pt->GetMean(), h_pt->GetRMS());
    h_pt->Fit(f_gaus, "RQ0");
    f_gaus->Draw("same");
  }
  c1->Update();

  TString outbase = gSystem->BaseName(infile);
  outbase.ReplaceAll(".root", "");
  outbase += Form("_trackid%u_pt", selected_track_id);

  const TString outdir = gSystem->DirName(infile);
  TString outprefix = outbase;
  if (!outdir.IsNull() && outdir != ".")
  {
    outprefix = outdir + "/" + outbase;
  }

  c1->SaveAs(outprefix + ".png");
  c1->SaveAs(outprefix + ".pdf");

  std::cout << "Filled " << nfilled << " entries for track_id=" << selected_track_id
            << std::endl;
  if (f_gaus)
  {
    std::cout << "Gaussian fit mean = " << f_gaus->GetParameter(1)
              << " GeV/c, sigma = " << f_gaus->GetParameter(2) << " GeV/c"
              << std::endl;
  }
  std::cout << "Saved plots to " << outprefix << ".png and " << outprefix << ".pdf"
            << std::endl;

  fin->Close();
}
