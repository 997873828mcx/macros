#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>

namespace
{
  bool fill_reco_hist(TTree* tree,
                      TH1D* hist,
                      const unsigned int selected_track_id,
                      const double fill_pt_min,
                      const double fill_pt_max,
                      Long64_t& nfilled)
  {
    if (!tree || !hist) return false;

    Float_t pt = 0.0;
    UInt_t track_id = 0;

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("pt", 1);
    tree->SetBranchStatus("track_id", 1);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("track_id", &track_id);

    nfilled = 0;
    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
      tree->GetEntry(i);
      if (track_id != selected_track_id) continue;
      if (!std::isfinite(pt)) continue;
      if (pt < fill_pt_min || pt > fill_pt_max) continue;
      hist->Fill(pt);
      ++nfilled;
    }

    return true;
  }

  bool fill_circle_hist(TTree* tree,
                        TH1D* hist,
                        const std::string& pt_branch_name,
                        const double fill_pt_min,
                        const double fill_pt_max,
                        Long64_t& nfilled)
  {
    if (!tree || !hist) return false;
    if (!tree->GetBranch(pt_branch_name.c_str())) return false;

    Double_t pt = 0.0;

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(pt_branch_name.c_str(), 1);
    tree->SetBranchAddress(pt_branch_name.c_str(), &pt);

    nfilled = 0;
    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
      tree->GetEntry(i);
      if (!std::isfinite(pt)) continue;
      if (pt < fill_pt_min || pt > fill_pt_max) continue;
      hist->Fill(pt);
      ++nfilled;
    }

    return true;
  }

  TF1* fit_gaussian(TH1D* hist,
                    const char* name,
                    const double fit_min,
                    const double fit_max,
                    const int line_color)
  {
    if (!hist || hist->GetEntries() <= 0) return nullptr;

    auto* f = new TF1(name, "gaus", fit_min, fit_max);
    f->SetLineColor(line_color);
    f->SetLineWidth(2);
    f->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());

    const int status = hist->Fit(f, "RQ0");
    if (status != 0 || !std::isfinite(f->GetParameter(1)) || !std::isfinite(f->GetParameter(2)))
    {
      if (hist->GetListOfFunctions()) hist->GetListOfFunctions()->Remove(f);
      delete f;
      return nullptr;
    }

    return f;
  }
}  // namespace

void CompareRecoVsCircleFitPt(
    const char* reco_file = "output/pionplus_pt10/completed/merged_reco_pt.root",
    const char* circlefit_file = "tpc_dnl_pt_resolution.root",
    const char* circlefit_branch = "pt_after",
    const unsigned int selected_track_id = 0,
    const int nbins = 120,
    const double pt_min = 0.0,
    const double pt_max = 20.0,
    const double fill_pt_min = 7.0,
    const double fill_pt_max = 14.0,
    const bool normalize = true)
{
  TFile* f_reco = TFile::Open(reco_file, "READ");
  if (!f_reco || f_reco->IsZombie())
  {
    std::cerr << "Failed to open reco file: " << reco_file << std::endl;
    return;
  }

  TFile* f_circle = TFile::Open(circlefit_file, "READ");
  if (!f_circle || f_circle->IsZombie())
  {
    std::cerr << "Failed to open circle-fit file: " << circlefit_file << std::endl;
    f_reco->Close();
    return;
  }

  TTree* reco_tracks = dynamic_cast<TTree*>(f_reco->Get("reco_tracks"));
  if (!reco_tracks)
  {
    std::cerr << "Could not find TTree 'reco_tracks' in " << reco_file << std::endl;
    f_circle->Close();
    f_reco->Close();
    return;
  }

  TTree* ptfits = dynamic_cast<TTree*>(f_circle->Get("ptfits"));
  if (!ptfits)
  {
    std::cerr << "Could not find TTree 'ptfits' in " << circlefit_file << std::endl;
    f_circle->Close();
    f_reco->Close();
    return;
  }

  auto* h_reco = new TH1D(
      "h_reco_track0_compare",
      Form("Reco vs circle-fit p_{T};p_{T} [GeV/c];%s",
           normalize ? "Normalized counts" : "Counts"),
      nbins, pt_min, pt_max);
  auto* h_circle = new TH1D(
      "h_circlefit_compare",
      Form("Reco vs circle-fit p_{T};p_{T} [GeV/c];%s",
           normalize ? "Normalized counts" : "Counts"),
      nbins, pt_min, pt_max);
  h_reco->SetDirectory(nullptr);
  h_circle->SetDirectory(nullptr);

  Long64_t nreco = 0;
  Long64_t ncircle = 0;
  if (!fill_reco_hist(reco_tracks, h_reco, selected_track_id, fill_pt_min, fill_pt_max, nreco))
  {
    std::cerr << "Failed to fill reco histogram" << std::endl;
    f_circle->Close();
    f_reco->Close();
    return;
  }

  if (!fill_circle_hist(ptfits, h_circle, circlefit_branch, fill_pt_min, fill_pt_max, ncircle))
  {
    std::cerr << "Failed to fill circle-fit histogram from branch '" << circlefit_branch << "'" << std::endl;
    f_circle->Close();
    f_reco->Close();
    return;
  }

  if (normalize)
  {
    const double reco_int = h_reco->Integral("width");
    const double circle_int = h_circle->Integral("width");
    if (reco_int > 0.0) h_reco->Scale(1.0 / reco_int);
    if (circle_int > 0.0) h_circle->Scale(1.0 / circle_int);
  }

  h_reco->SetLineColor(kBlue + 1);
  h_reco->SetLineWidth(2);
  h_reco->SetStats(0);

  h_circle->SetLineColor(kRed + 1);
  h_circle->SetLineWidth(2);
  h_circle->SetStats(0);

  auto* f_reco_gaus = fit_gaussian(h_reco, "f_reco_gaus_compare", fill_pt_min, fill_pt_max, kBlue + 1);
  auto* f_circle_gaus = fit_gaussian(h_circle, "f_circle_gaus_compare", fill_pt_min, fill_pt_max, kRed + 1);

  const double ymax = 1.15 * std::max(h_reco->GetMaximum(), h_circle->GetMaximum());
  h_reco->SetMaximum(ymax > 0.0 ? ymax : 1.0);

  auto* c1 = new TCanvas("c_reco_vs_circlefit_pt", "Reco vs circle fit pT", 950, 750);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);
  h_reco->Draw("hist");
  h_circle->Draw("hist same");
  if (f_reco_gaus) f_reco_gaus->Draw("same");
  if (f_circle_gaus) f_circle_gaus->Draw("same");

  auto* leg = new TLegend(0.52, 0.68, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_reco, Form("sPHENIX reco track_id=%u", selected_track_id), "l");
  leg->AddEntry(h_circle, Form("Circle fit (%s)", circlefit_branch), "l");
  leg->Draw();

  c1->Update();

  TString outbase = TString::Format("compare_reco_trackid%u_vs_circle_%s",
                                    selected_track_id,
                                    circlefit_branch);
  const TString outdir = gSystem->DirName(circlefit_file);
  TString outprefix = outbase;
  if (!outdir.IsNull() && outdir != ".")
  {
    outprefix = outdir + "/" + outbase;
  }

  c1->SaveAs(outprefix + ".png");
  c1->SaveAs(outprefix + ".pdf");

  std::cout << "Reco entries filled: " << nreco << std::endl;
  std::cout << "Circle-fit entries filled: " << ncircle << std::endl;
  if (f_reco_gaus)
  {
    std::cout << "Reco gaus mean = " << f_reco_gaus->GetParameter(1)
              << " GeV/c, sigma = " << std::abs(f_reco_gaus->GetParameter(2)) << " GeV/c"
              << std::endl;
  }
  if (f_circle_gaus)
  {
    std::cout << "Circle-fit gaus mean = " << f_circle_gaus->GetParameter(1)
              << " GeV/c, sigma = " << std::abs(f_circle_gaus->GetParameter(2)) << " GeV/c"
              << std::endl;
  }
  std::cout << "Saved comparison plots to " << outprefix << ".png and " << outprefix << ".pdf"
            << std::endl;

  f_circle->Close();
  f_reco->Close();
}
