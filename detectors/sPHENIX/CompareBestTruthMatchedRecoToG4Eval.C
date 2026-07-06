#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace
{
  struct RecoCandidate
  {
    bool valid = false;
    float pt = std::numeric_limits<float>::quiet_NaN();
    float truth_weight = -1.0f;
    unsigned int reco_track_id = 0;
  };

  TF1* fit_gaussian_best_truth_match(TH1D* hist,
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

  bool fill_best_truth_reco_hist(TChain* reco_chain,
                                 TChain* event_chain,
                                 TH1D* hist,
                                 const int selected_truth_track_id,
                                 const int selected_truth_pid,
                                 const bool require_truth_pid,
                                 const double fill_pt_min,
                                 const double fill_pt_max,
                                 Long64_t& nmatching_candidates,
                                 Long64_t& nevent_best_tracks,
                                 Long64_t& nfilled)
  {
    if (!reco_chain || !event_chain || !hist) return false;

    if (!reco_chain->GetBranch("pt") ||
        !reco_chain->GetBranch("best_truth_match_ok") ||
        !reco_chain->GetBranch("best_truth_track_id") ||
        !reco_chain->GetBranch("best_truth_weight") ||
        !event_chain->GetBranch("nreco"))
    {
      return false;
    }

    Float_t pt = std::numeric_limits<float>::quiet_NaN();
    Int_t best_truth_match_ok = 0;
    Int_t best_truth_track_id = std::numeric_limits<int>::min();
    Int_t best_truth_pid = std::numeric_limits<int>::min();
    Float_t best_truth_weight = std::numeric_limits<float>::quiet_NaN();
    UInt_t reco_track_id = 0;
    Int_t nreco_in_event = 0;

    reco_chain->SetBranchStatus("*", 0);
    reco_chain->SetBranchStatus("pt", 1);
    reco_chain->SetBranchStatus("best_truth_match_ok", 1);
    reco_chain->SetBranchStatus("best_truth_track_id", 1);
    reco_chain->SetBranchStatus("best_truth_weight", 1);
    reco_chain->SetBranchStatus("track_id", 1);
    if (require_truth_pid)
    {
      if (!reco_chain->GetBranch("best_truth_pid")) return false;
      reco_chain->SetBranchStatus("best_truth_pid", 1);
      reco_chain->SetBranchAddress("best_truth_pid", &best_truth_pid);
    }

    reco_chain->SetBranchAddress("pt", &pt);
    reco_chain->SetBranchAddress("best_truth_match_ok", &best_truth_match_ok);
    reco_chain->SetBranchAddress("best_truth_track_id", &best_truth_track_id);
    reco_chain->SetBranchAddress("best_truth_weight", &best_truth_weight);
    reco_chain->SetBranchAddress("track_id", &reco_track_id);

    event_chain->SetBranchStatus("*", 0);
    event_chain->SetBranchStatus("nreco", 1);
    event_chain->SetBranchAddress("nreco", &nreco_in_event);

    nmatching_candidates = 0;
    nevent_best_tracks = 0;
    nfilled = 0;

    const Long64_t nreco_entries_total = reco_chain->GetEntries();
    const Long64_t nevent_entries_total = event_chain->GetEntries();
    Long64_t reco_entry = 0;

    for (Long64_t ievent = 0; ievent < nevent_entries_total; ++ievent)
    {
      event_chain->GetEntry(ievent);
      if (nreco_in_event < 0)
      {
        return false;
      }

      RecoCandidate best;
      for (Int_t i = 0; i < nreco_in_event; ++i)
      {
        if (reco_entry >= nreco_entries_total)
        {
          return false;
        }

        reco_chain->GetEntry(reco_entry++);

        if (!best_truth_match_ok) continue;
        if (best_truth_track_id != selected_truth_track_id) continue;
        if (require_truth_pid && best_truth_pid != selected_truth_pid) continue;
        if (!std::isfinite(pt)) continue;
        if (!std::isfinite(best_truth_weight)) continue;

        ++nmatching_candidates;

        if (!best.valid || best_truth_weight > best.truth_weight)
        {
          best.valid = true;
          best.pt = pt;
          best.truth_weight = best_truth_weight;
          best.reco_track_id = reco_track_id;
        }
      }

      if (!best.valid) continue;
      ++nevent_best_tracks;
      if (!std::isfinite(best.pt)) continue;
      if (best.pt < fill_pt_min || best.pt > fill_pt_max) continue;
      hist->Fill(best.pt);
      ++nfilled;
    }

    return reco_entry == nreco_entries_total;
  }

  bool fill_g4eval_hist(TChain* chain,
                        TH1D* hist,
                        const int selected_truth_track_id,
                        const int selected_truth_pid,
                        const bool require_truth_pid,
                        const bool require_primary,
                        const double fill_pt_min,
                        const double fill_pt_max,
                        Long64_t& nfilled)
  {
    if (!chain || !hist) return false;
    if (!chain->GetBranch("gtrackID") || !chain->GetBranch("pt")) return false;

    Float_t gtrackID = std::numeric_limits<float>::quiet_NaN();
    Float_t gflavor = std::numeric_limits<float>::quiet_NaN();
    Float_t gprimary = std::numeric_limits<float>::quiet_NaN();
    Float_t pt = std::numeric_limits<float>::quiet_NaN();

    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("gtrackID", 1);
    chain->SetBranchStatus("pt", 1);
    chain->SetBranchAddress("gtrackID", &gtrackID);
    chain->SetBranchAddress("pt", &pt);

    if (require_truth_pid)
    {
      if (!chain->GetBranch("gflavor")) return false;
      chain->SetBranchStatus("gflavor", 1);
      chain->SetBranchAddress("gflavor", &gflavor);
    }

    if (require_primary)
    {
      if (!chain->GetBranch("gprimary")) return false;
      chain->SetBranchStatus("gprimary", 1);
      chain->SetBranchAddress("gprimary", &gprimary);
    }

    nfilled = 0;
    const Long64_t nentries = chain->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
      chain->GetEntry(i);
      if (!std::isfinite(gtrackID) || static_cast<int>(std::lround(gtrackID)) != selected_truth_track_id) continue;
      if (require_truth_pid && (!std::isfinite(gflavor) || static_cast<int>(std::lround(gflavor)) != selected_truth_pid)) continue;
      if (require_primary && (!std::isfinite(gprimary) || gprimary < 0.5f)) continue;
      if (!std::isfinite(pt)) continue;
      if (pt < fill_pt_min || pt > fill_pt_max) continue;
      hist->Fill(pt);
      ++nfilled;
    }

    return true;
  }
}  // namespace

void CompareBestTruthMatchedRecoToG4Eval(
    const char* reco_input = "output/pionplus_pt10/completed/general_*_reco_pt.root",
    const char* g4eval_input = "",
    const int selected_truth_track_id = 0,
    const int selected_truth_pid = 211,
    const int nbins = 120,
    const double pt_min = 0.0,
    const double pt_max = 20.0,
    const double fill_pt_min = 7.0,
    const double fill_pt_max = 14.0,
    const bool normalize = false,
    const bool require_truth_pid = true,
    const bool require_primary = true)
{
  auto* reco_chain = new TChain("reco_tracks");
  const int nrecofiles = reco_chain->Add(reco_input);
  if (nrecofiles <= 0)
  {
    std::cerr << "Failed to add reco input: " << reco_input << std::endl;
    delete reco_chain;
    return;
  }

  auto* event_chain = new TChain("event_info");
  const int neventfiles = event_chain->Add(reco_input);
  if (neventfiles <= 0)
  {
    std::cerr << "Failed to add event_info from reco input: " << reco_input << std::endl;
    delete event_chain;
    delete reco_chain;
    return;
  }

  if (nrecofiles == 1 || neventfiles == 1)
  {
    std::cout << "Reco input resolved to one file. Merged files are supported here by using "
              << "the event_info tree to reconstruct event boundaries." << std::endl;
  }

  auto* h_reco = new TH1D(
      "h_best_truth_reco_pt",
      Form("Best truth-matched reco p_{T};p_{T} [GeV/c];%s",
           normalize ? "Normalized counts" : "Counts"),
      nbins, pt_min, pt_max);
  h_reco->SetDirectory(nullptr);
  h_reco->SetLineColor(kBlue + 1);
  h_reco->SetLineWidth(2);
  h_reco->SetStats(0);

  Long64_t nmatching_candidates = 0;
  Long64_t nevent_best_tracks = 0;
  Long64_t nrecofilled = 0;
  if (!fill_best_truth_reco_hist(
          reco_chain,
          event_chain,
          h_reco,
          selected_truth_track_id,
          selected_truth_pid,
          require_truth_pid,
          fill_pt_min,
          fill_pt_max,
          nmatching_candidates,
          nevent_best_tracks,
          nrecofilled))
  {
    std::cerr << "Failed to fill reco histogram from reco_tracks. "
              << "Make sure the truth-match branches are present in your G4_User output." << std::endl;
    delete event_chain;
    delete h_reco;
    delete reco_chain;
    return;
  }

  TH1D* h_eval = nullptr;
  Long64_t nevalfilled = 0;
  if (g4eval_input && std::string(g4eval_input).size() > 0)
  {
    auto* eval_chain = new TChain("ntp_gtrack");
    const int nevalfiles = eval_chain->Add(g4eval_input);
    if (nevalfiles <= 0)
    {
      std::cerr << "Failed to add g4svtx evaluator input: " << g4eval_input << std::endl;
      delete eval_chain;
      delete event_chain;
      delete h_reco;
      delete reco_chain;
      return;
    }

    h_eval = new TH1D(
        "h_g4eval_best_track_pt",
        Form("g4svtx_eval ntp_gtrack p_{T};p_{T} [GeV/c];%s",
             normalize ? "Normalized counts" : "Counts"),
        nbins, pt_min, pt_max);
    h_eval->SetDirectory(nullptr);
    h_eval->SetLineColor(kRed + 1);
    h_eval->SetLineWidth(2);
    h_eval->SetStats(0);

    if (!fill_g4eval_hist(
            eval_chain,
            h_eval,
            selected_truth_track_id,
            selected_truth_pid,
            require_truth_pid,
            require_primary,
            fill_pt_min,
            fill_pt_max,
            nevalfilled))
    {
      std::cerr << "Failed to fill evaluator histogram from ntp_gtrack." << std::endl;
      delete eval_chain;
      delete h_eval;
      delete event_chain;
      delete h_reco;
      delete reco_chain;
      return;
    }

    delete eval_chain;
  }

  if (normalize)
  {
    const double reco_int = h_reco->Integral("width");
    if (reco_int > 0.0) h_reco->Scale(1.0 / reco_int);

    if (h_eval)
    {
      const double eval_int = h_eval->Integral("width");
      if (eval_int > 0.0) h_eval->Scale(1.0 / eval_int);
    }
  }

  auto* f_reco = fit_gaussian_best_truth_match(
      h_reco, "f_best_truth_reco", fill_pt_min, fill_pt_max, kBlue + 1);
  auto* f_eval = fit_gaussian_best_truth_match(
      h_eval, "f_g4eval_best_track", fill_pt_min, fill_pt_max, kRed + 1);

  double ymax = h_reco->GetMaximum();
  if (h_eval) ymax = std::max(ymax, h_eval->GetMaximum());
  h_reco->SetMaximum(ymax > 0.0 ? 1.15 * ymax : 1.0);

  auto* c1 = new TCanvas("c_best_truth_match_compare", "Best truth-matched reco pT", 950, 750);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);
  h_reco->SetTitle(
      Form("Truth track %d%s%s",
           selected_truth_track_id,
           require_truth_pid ? Form(", pid = %d", selected_truth_pid) : "",
           require_primary ? ", primary only in g4eval" : ""));
  h_reco->Draw("hist");
  if (h_eval) h_eval->Draw("hist same");
  if (f_reco) f_reco->Draw("same");
  if (f_eval) f_eval->Draw("same");

  auto* leg = new TLegend(0.47, 0.68, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_reco, "reco_tracks best truth-matched reco track", "l");
  if (h_eval) leg->AddEntry(h_eval, "g4svtx_eval ntp_gtrack", "l");
  leg->Draw();

  c1->Update();

  TString outdir = gSystem->DirName(reco_input);
  if (outdir.IsNull() || outdir == ".")
  {
    outdir = gSystem->WorkingDirectory();
  }

  TString outbase = Form("truthtrack%d_best_truth_match_pt", selected_truth_track_id);
  if (require_truth_pid)
  {
    outbase += Form("_pid%d", selected_truth_pid);
  }
  if (h_eval)
  {
    outbase += "_vs_g4svtx_eval";
  }

  const TString outprefix = outdir + "/" + outbase;
  c1->SaveAs(outprefix + ".png");
  c1->SaveAs(outprefix + ".pdf");

  std::cout << "Reco matching candidates found: " << nmatching_candidates << std::endl;
  std::cout << "Reco best-track candidates kept: " << nevent_best_tracks << std::endl;
  std::cout << "Reco histogram entries filled: " << nrecofilled << std::endl;
  if (h_eval)
  {
    std::cout << "g4svtx_eval histogram entries filled: " << nevalfilled << std::endl;
  }
  if (f_reco)
  {
    std::cout << "Reco Gaussian mean = " << f_reco->GetParameter(1)
              << " GeV/c, sigma = " << f_reco->GetParameter(2) << " GeV/c" << std::endl;
  }
  if (f_eval)
  {
    std::cout << "g4svtx_eval Gaussian mean = " << f_eval->GetParameter(1)
              << " GeV/c, sigma = " << f_eval->GetParameter(2) << " GeV/c" << std::endl;
  }
  std::cout << "Saved plots to " << outprefix << ".png and " << outprefix << ".pdf" << std::endl;

  delete h_eval;
  delete event_chain;
  delete h_reco;
  delete reco_chain;
}
