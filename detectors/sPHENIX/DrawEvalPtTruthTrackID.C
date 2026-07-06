#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>

void DrawEvalPtTruthTrackID(
    const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pionplus_pt10/completed/merged_eval.root",
    const char* tree_name = "ntp_track",
    const int selected_truth_track_id = -1,
    const int selected_truth_pid = 0,
    const bool require_primary = true,
    const bool require_singlematch = false,
    const int nbins = 120,
    const double pt_min = 0.0,
    const double pt_max = 20.0,
    const double fill_pt_min = 7.0,
    const double fill_pt_max = 14.0,
    const int rel_nbins = 120,
    const double rel_min = -0.5,
    const double rel_max = 0.5,
    const double rel_fit_min = -0.3,
    const double rel_fit_max = 0.3)
{
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "Failed to open input file: " << infile << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(fin->Get(tree_name));
  if (!tree)
  {
    std::cerr << "Could not find tree '" << tree_name << "' in " << infile << std::endl;
    fin->Close();
    return;
  }

  if (!tree->GetBranch("pt") || !tree->GetBranch("gtrackID") || !tree->GetBranch("gpt"))
  {
    std::cerr << "Tree '" << tree_name
              << "' is missing required branches 'pt', 'gtrackID', or 'gpt'" << std::endl;
    fin->Close();
    return;
  }

  const bool use_truth_track_cut = (selected_truth_track_id >= 0);
  const bool use_truth_pid_cut = (selected_truth_pid != 0);
  const bool has_gflavor = tree->GetBranch("gflavor");
  const bool has_gprimary = tree->GetBranch("gprimary");
  const bool has_singlematch = tree->GetBranch("singlematch");

  if (use_truth_pid_cut && !has_gflavor)
  {
    std::cerr << "Tree '" << tree_name << "' is missing branch 'gflavor'" << std::endl;
    fin->Close();
    return;
  }
  if (require_primary && !has_gprimary)
  {
    std::cerr << "Tree '" << tree_name << "' is missing branch 'gprimary'" << std::endl;
    fin->Close();
    return;
  }
  if (require_singlematch && !has_singlematch)
  {
    std::cerr << "Tree '" << tree_name << "' is missing branch 'singlematch'" << std::endl;
    fin->Close();
    return;
  }

  TH1D* h_pt = new TH1D(
      "h_eval_pt_truth_track",
      Form("%s reco p_{T}%s%s;Reco p_{T} [GeV/c];Counts",
           tree_name,
           use_truth_track_cut ? Form(" for gtrackID=%d", selected_truth_track_id) : " for all truth IDs",
           use_truth_pid_cut ? Form(", pid=%d", selected_truth_pid) : ", all pids"),
      nbins, pt_min, pt_max);
  h_pt->SetDirectory(nullptr);
  h_pt->SetLineWidth(2);
  h_pt->SetLineColor(kBlue + 1);
  h_pt->SetStats(0);

  TH1D* h_rel = new TH1D(
      "h_eval_rel_pt_truth_track",
      Form("%s (p_{T}/g_{T} - 1)%s%s;(p_{T}/g_{T} - 1);Counts",
           tree_name,
           use_truth_track_cut ? Form(" for gtrackID=%d", selected_truth_track_id) : " for all truth IDs",
           use_truth_pid_cut ? Form(", pid=%d", selected_truth_pid) : ", all pids"),
      rel_nbins, rel_min, rel_max);
  h_rel->SetDirectory(nullptr);
  h_rel->SetLineWidth(2);
  h_rel->SetLineColor(kBlue + 1);
  h_rel->SetStats(0);

  Float_t pt = 0.0f;
  Float_t gpt = 0.0f;
  Float_t gtrackID = 0.0f;
  Float_t gflavor = 0.0f;
  Float_t gprimary = 0.0f;
  Float_t singlematch = 0.0f;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("gpt", 1);
  tree->SetBranchStatus("gtrackID", 1);
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("gpt", &gpt);
  tree->SetBranchAddress("gtrackID", &gtrackID);

  if (use_truth_pid_cut)
  {
    tree->SetBranchStatus("gflavor", 1);
    tree->SetBranchAddress("gflavor", &gflavor);
  }
  if (require_primary)
  {
    tree->SetBranchStatus("gprimary", 1);
    tree->SetBranchAddress("gprimary", &gprimary);
  }
  if (require_singlematch)
  {
    tree->SetBranchStatus("singlematch", 1);
    tree->SetBranchAddress("singlematch", &singlematch);
  }

  Long64_t nfilled = 0;
  Long64_t nrelfilled = 0;
  const Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    tree->GetEntry(i);
    if (!std::isfinite(pt))
    {
      continue;
    }
    if (use_truth_track_cut && (!std::isfinite(gtrackID) ||
        static_cast<int>(std::lround(gtrackID)) != selected_truth_track_id))
    {
      continue;
    }
    if (use_truth_pid_cut && (!std::isfinite(gflavor) ||
        static_cast<int>(std::lround(gflavor)) != selected_truth_pid))
    {
      continue;
    }
    if (require_primary && (!std::isfinite(gprimary) || gprimary < 0.5f))
    {
      continue;
    }
    if (require_singlematch && (!std::isfinite(singlematch) || singlematch < 0.5f))
    {
      continue;
    }
    if (pt < fill_pt_min || pt > fill_pt_max)
    {
      continue;
    }

    h_pt->Fill(pt);
    ++nfilled;

    if (!std::isfinite(gpt) || std::abs(gpt) < 1e-9f)
    {
      continue;
    }

    const double rel = static_cast<double>(pt) / static_cast<double>(gpt) - 1.0;
    if (!std::isfinite(rel))
    {
      continue;
    }
    if (rel < rel_min || rel > rel_max)
    {
      continue;
    }

    h_rel->Fill(rel);
    ++nrelfilled;
  }

  TCanvas* c1 = new TCanvas("c_eval_pt_truth_track", "Evaluator pT", 900, 700);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);
  h_pt->Draw("hist");

  TF1* f_gaus = nullptr;
  if (nfilled > 0)
  {
    f_gaus = new TF1("f_eval_pt_truth_track", "gaus", fill_pt_min, fill_pt_max);
    f_gaus->SetLineColor(kRed + 1);
    f_gaus->SetLineWidth(2);
    f_gaus->SetParameters(h_pt->GetMaximum(), h_pt->GetMean(), h_pt->GetRMS());
    h_pt->Fit(f_gaus, "RQ0");
    f_gaus->Draw("same");
  }

  auto* leg1 = new TLegend(0.54, 0.76, 0.89, 0.89);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_pt, tree_name, "l");
  if (f_gaus) leg1->AddEntry(f_gaus, "Gaussian fit", "l");
  leg1->Draw();
  c1->Update();

  TCanvas* c2 = new TCanvas("c_eval_rel_pt_truth_track", "Evaluator relative pT", 900, 700);
  c2->SetMargin(0.12, 0.04, 0.12, 0.06);
  h_rel->Draw("hist");

  TF1* f_rel = nullptr;
  if (nrelfilled > 0)
  {
    f_rel = new TF1("f_eval_rel_pt_truth_track", "gaus", rel_fit_min, rel_fit_max);
    f_rel->SetLineColor(kRed + 1);
    f_rel->SetLineWidth(2);
    f_rel->SetParameters(h_rel->GetMaximum(), h_rel->GetMean(), h_rel->GetRMS());
    h_rel->Fit(f_rel, "RQ0");
    f_rel->Draw("same");
  }

  auto* leg2 = new TLegend(0.54, 0.76, 0.89, 0.89);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry(h_rel, tree_name, "l");
  if (f_rel) leg2->AddEntry(f_rel, "Gaussian fit", "l");
  leg2->Draw();
  c2->Update();

  TString outbase = gSystem->BaseName(infile);
  outbase.ReplaceAll(".root", "");
  outbase += Form("_%s", tree_name);
  if (use_truth_track_cut) outbase += Form("_gtrack%d", selected_truth_track_id);
  else outbase += "_allgtracks";
  if (use_truth_pid_cut) outbase += Form("_pid%d", selected_truth_pid);
  else outbase += "_allpid";
  if (require_primary) outbase += "_primary";
  if (require_singlematch) outbase += "_singlematch";
  outbase += "_pt";

  const TString outdir = gSystem->DirName(infile);
  TString outprefix = outbase;
  if (!outdir.IsNull() && outdir != ".")
  {
    outprefix = outdir + "/" + outbase;
  }

  c1->SaveAs(outprefix + ".png");
  c1->SaveAs(outprefix + ".pdf");
  c2->SaveAs(outprefix + "_relpt.png");
  c2->SaveAs(outprefix + "_relpt.pdf");

  std::cout << "Filled " << nfilled << " entries from tree '" << tree_name << "'" << std::endl;
  if (nfilled > 0)
  {
    std::cout << "Histogram mean in [" << fill_pt_min << ", " << fill_pt_max
              << "] GeV/c = " << h_pt->GetMean()
              << " GeV/c, stddev = " << h_pt->GetStdDev() << " GeV/c" << std::endl;
  }
  if (f_gaus)
  {
    std::cout << "Gaussian fit mean = " << f_gaus->GetParameter(1)
              << " GeV/c, sigma = " << f_gaus->GetParameter(2) << " GeV/c" << std::endl;
  }
  std::cout << "Filled " << nrelfilled << " entries for (pt/gpt - 1)" << std::endl;
  if (nrelfilled > 0)
  {
    std::cout << "Relative histogram mean = " << h_rel->GetMean()
              << ", stddev = " << h_rel->GetStdDev() << std::endl;
  }
  if (f_rel)
  {
    std::cout << "Relative Gaussian mean = " << f_rel->GetParameter(1)
              << ", sigma = " << f_rel->GetParameter(2) << std::endl;
  }
  std::cout << "Saved plots to " << outprefix << ".png/.pdf and "
            << outprefix << "_relpt.png/.pdf" << std::endl;

  fin->Close();
}
