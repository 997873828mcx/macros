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

namespace
{
  TF1* fit_gaussian_reco_vs_tpc_circle(TH1D* hist,
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

void CompareRecoVsTpcCircleTrackID(
    const char* infile = "output/pionplus_pt10/completed/merged_reco_pt.root",
    const int selected_track_id = -1,
    const int selected_truth_pid = 0,
    const int nbins = 120,
    const double pt_min = 0.0,
    const double pt_max = 20.0,
    const double fill_pt_min = 7.0,
    const double fill_pt_max = 14.0,
    const bool normalize = false)
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

  if (!reco_tracks->GetBranch("pt_tpc_circle"))
  {
    std::cerr << "Tree 'reco_tracks' does not contain branch 'pt_tpc_circle'. "
              << "You need to rerun with the updated G4_User.C first." << std::endl;
    fin->Close();
    return;
  }
  const bool has_truth_pid = reco_tracks->GetBranch("best_truth_pid");
  const bool has_truth_match_ok = reco_tracks->GetBranch("best_truth_match_ok");
  const bool has_tpc_circle_ok = reco_tracks->GetBranch("tpc_circle_ok");
  const bool use_truth_selection = (selected_truth_pid != 0);

  if (use_truth_selection && (!has_truth_pid || !has_truth_match_ok))
  {
    std::cerr << "Tree 'reco_tracks' does not contain truth-match branches "
              << "('best_truth_pid' and 'best_truth_match_ok'). "
              << "Either rerun with the updated G4_User.C or call "
              << "CompareRecoVsTpcCircleTrackID(..., selected_truth_pid = 0, ...)." << std::endl;
    fin->Close();
    return;
  }

  TString selection_label = "all track_id values";
  if (selected_track_id >= 0)
  {
    selection_label = Form("track_id = %d", selected_track_id);
  }
  if (use_truth_selection)
  {
    selection_label += Form("%struth pid = %d",
                            selection_label.IsNull() ? "" : ", ",
                            selected_truth_pid);
  }

  auto* h_reco = new TH1D(
      "h_reco_pt_compare",
      Form("Reco vs TPC circle-fit p_{T} for %s;p_{T} [GeV/c];%s",
           selection_label.Data(),
           normalize ? "Normalized counts" : "Counts"),
      nbins, pt_min, pt_max);
  auto* h_circle = new TH1D(
      "h_tpc_circle_pt_compare",
      Form("Reco vs TPC circle-fit p_{T} for %s;p_{T} [GeV/c];%s",
           selection_label.Data(),
           normalize ? "Normalized counts" : "Counts"),
      nbins, pt_min, pt_max);
  h_reco->SetDirectory(nullptr);
  h_circle->SetDirectory(nullptr);

  Float_t pt = 0.0;
  Float_t pt_tpc_circle = 0.0;
  UInt_t track_id = 0;
  Int_t best_truth_pid = 0;
  Int_t best_truth_match_ok = 0;
  Int_t tpc_circle_ok = 0;

  reco_tracks->SetBranchStatus("*", 0);
  reco_tracks->SetBranchStatus("pt", 1);
  reco_tracks->SetBranchStatus("pt_tpc_circle", 1);
  reco_tracks->SetBranchStatus("track_id", 1);
  reco_tracks->SetBranchAddress("pt", &pt);
  reco_tracks->SetBranchAddress("pt_tpc_circle", &pt_tpc_circle);
  reco_tracks->SetBranchAddress("track_id", &track_id);
  if (has_truth_pid)
  {
    reco_tracks->SetBranchStatus("best_truth_pid", 1);
    reco_tracks->SetBranchAddress("best_truth_pid", &best_truth_pid);
  }
  if (has_truth_match_ok)
  {
    reco_tracks->SetBranchStatus("best_truth_match_ok", 1);
    reco_tracks->SetBranchAddress("best_truth_match_ok", &best_truth_match_ok);
  }
  if (has_tpc_circle_ok)
  {
    reco_tracks->SetBranchStatus("tpc_circle_ok", 1);
    reco_tracks->SetBranchAddress("tpc_circle_ok", &tpc_circle_ok);
  }

  Long64_t nreco = 0;
  Long64_t ncircle = 0;
  const Long64_t nentries = reco_tracks->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    reco_tracks->GetEntry(i);
    if (selected_track_id >= 0 && static_cast<int>(track_id) != selected_track_id)
    {
      continue;
    }

    if (use_truth_selection && (!best_truth_match_ok || best_truth_pid != selected_truth_pid))
    {
      continue;
    }

    if (std::isfinite(pt) && pt >= fill_pt_min && pt <= fill_pt_max)
    {
      h_reco->Fill(pt);
      ++nreco;
    }

    if ((!has_tpc_circle_ok || tpc_circle_ok) && std::isfinite(pt_tpc_circle) &&
        pt_tpc_circle >= fill_pt_min && pt_tpc_circle <= fill_pt_max)
    {
      h_circle->Fill(pt_tpc_circle);
      ++ncircle;
    }
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

  auto* f_reco = fit_gaussian_reco_vs_tpc_circle(
      h_reco, "f_reco_pt_compare", fill_pt_min, fill_pt_max, kBlue + 1);
  auto* f_circle = fit_gaussian_reco_vs_tpc_circle(
      h_circle, "f_tpc_circle_pt_compare", fill_pt_min, fill_pt_max, kRed + 1);

  const double ymax = 1.15 * std::max(h_reco->GetMaximum(), h_circle->GetMaximum());
  h_reco->SetMaximum(ymax > 0.0 ? ymax : 1.0);

  auto* c1 = new TCanvas("c_reco_vs_tpc_circle_pt", "Reco vs TPC circle-fit pT", 950, 750);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);
  h_reco->Draw("hist");
  h_circle->Draw("hist same");
  if (f_reco) f_reco->Draw("same");
  if (f_circle) f_circle->Draw("same");

  auto* leg = new TLegend(0.50, 0.68, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_reco, "sPHENIX TrackingReco p_{T}", "l");
  leg->AddEntry(h_circle, "TPC circle-fit p_{T}", "l");
  leg->Draw();

  c1->Update();

  TString outbase = gSystem->BaseName(infile);
  outbase.ReplaceAll(".root", "");
  if (use_truth_selection)
  {
    outbase += Form("_truthpid%d", selected_truth_pid);
  }
  if (selected_track_id >= 0)
  {
    outbase += Form("_trackid%d", selected_track_id);
  }
  else
  {
    outbase += "_alltrackids";
  }
  outbase += "_reco_vs_tpc_circle_pt";

  const TString outdir = gSystem->DirName(infile);
  TString outprefix = outbase;
  if (!outdir.IsNull() && outdir != ".")
  {
    outprefix = outdir + "/" + outbase;
  }

  c1->SaveAs(outprefix + ".png");
  c1->SaveAs(outprefix + ".pdf");

  std::cout << "Selection: " << selection_label.Data() << std::endl;
  std::cout << "Reco entries filled: " << nreco << std::endl;
  std::cout << "TPC circle-fit entries filled: " << ncircle << std::endl;
  if (f_reco)
  {
    std::cout << "Reco gaus mean = " << f_reco->GetParameter(1)
              << " GeV/c, sigma = " << std::abs(f_reco->GetParameter(2)) << " GeV/c"
              << std::endl;
  }
  if (f_circle)
  {
    std::cout << "TPC circle-fit gaus mean = " << f_circle->GetParameter(1)
              << " GeV/c, sigma = " << std::abs(f_circle->GetParameter(2)) << " GeV/c"
              << std::endl;
  }
  std::cout << "Saved comparison plots to " << outprefix << ".png and " << outprefix << ".pdf"
            << std::endl;

  fin->Close();
}
