#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>

#include <iostream>
#include <string>

void draw_track_pt(
    const char *infile = "output/pp_minbias_sim_helix_auto_trackqa_new15_merged.root",
    const char *outbase = "output/track_pt_all",
    const char *tree_name = "trackTree",
    const double pt_min = 0.0,
    const double pt_max = 10.0,
    const int nbins = 250,
    const char *extra_cut = "",
    const bool draw_truth = false)
{
  TFile *file = TFile::Open(infile, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "[draw_track_pt] cannot open " << infile << std::endl;
    return;
  }

  TTree *tree = nullptr;
  file->GetObject(tree_name, tree);
  if (!tree)
  {
    std::cerr << "[draw_track_pt] tree " << tree_name << " not found in " << infile << std::endl;
    file->Close();
    return;
  }

  gStyle->SetOptStat(0);

  const std::string extra = extra_cut ? extra_cut : "";
  auto cut = [&extra](const std::string &base)
  {
    if (extra.empty())
    {
      return base;
    }
    return "(" + base + ") && (" + extra + ")";
  };

  TH1F *h_pt_all = new TH1F(
      "h_pt_all",
      ";p_{T} [GeV/c];Tracks",
      nbins, pt_min, pt_max);
  TH1F *h_pt_pos = new TH1F(
      "h_pt_pos",
      ";p_{T} [GeV/c];Tracks",
      nbins, pt_min, pt_max);
  TH1F *h_pt_neg = new TH1F(
      "h_pt_neg",
      ";p_{T} [GeV/c];Tracks",
      nbins, pt_min, pt_max);
  TH1F *h_pt_zero = new TH1F(
      "h_pt_zero",
      ";p_{T} [GeV/c];Tracks",
      nbins, pt_min, pt_max);

  tree->Draw("pt>>h_pt_all", cut("pt > 0").c_str(), "goff");
  tree->Draw("pt>>h_pt_pos", cut("pt > 0 && charge > 0").c_str(), "goff");
  tree->Draw("pt>>h_pt_neg", cut("pt > 0 && charge < 0").c_str(), "goff");
  tree->Draw("pt>>h_pt_zero", cut("pt > 0 && charge == 0").c_str(), "goff");

  h_pt_all->SetLineColor(kBlack);
  h_pt_all->SetLineWidth(2);
  h_pt_pos->SetLineColor(kRed + 1);
  h_pt_pos->SetLineWidth(2);
  h_pt_neg->SetLineColor(kBlue + 1);
  h_pt_neg->SetLineWidth(2);
  h_pt_zero->SetLineColor(kGreen + 2);
  h_pt_zero->SetLineWidth(2);

  TCanvas *c_pt = new TCanvas("c_track_pt", "track_pt", 1000, 800);
  c_pt->SetLogy();
  h_pt_all->Draw("hist");
  h_pt_pos->Draw("hist same");
  h_pt_neg->Draw("hist same");
  if (h_pt_zero->GetEntries() > 0)
  {
    h_pt_zero->Draw("hist same");
  }

  TLegend *leg = new TLegend(0.58, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_pt_all, "all tracks", "l");
  leg->AddEntry(h_pt_pos, "q > 0", "l");
  leg->AddEntry(h_pt_neg, "q < 0", "l");
  if (h_pt_zero->GetEntries() > 0)
  {
    leg->AddEntry(h_pt_zero, "q = 0", "l");
  }
  leg->Draw();
  c_pt->SaveAs((std::string(outbase) + ".pdf").c_str());
  c_pt->SaveAs((std::string(outbase) + ".png").c_str());

  const bool has_truth =
      draw_truth &&
      tree->GetBranch("truth_px") &&
      tree->GetBranch("truth_py");

  TH1F *h_truth_pt_all = nullptr;
  TH2F *h_reco_vs_truth_pt = nullptr;
  TH1F *h_delta_pt = nullptr;
  TCanvas *c_truth = nullptr;
  TCanvas *c_reco_vs_truth = nullptr;
  TCanvas *c_delta = nullptr;

  if (has_truth)
  {
    h_truth_pt_all = new TH1F(
        "h_truth_pt_all",
        ";true p_{T} [GeV/c];Tracks",
        nbins, pt_min, pt_max);
    h_reco_vs_truth_pt = new TH2F(
        "h_reco_vs_truth_pt",
        ";true p_{T} [GeV/c];reco/fitted p_{T} [GeV/c]",
        nbins, pt_min, pt_max,
        nbins, pt_min, pt_max);
    h_delta_pt = new TH1F(
        "h_delta_pt",
        ";p_{T}^{reco} - p_{T}^{true} [GeV/c];Tracks",
        240, -2.0, 2.0);

    const std::string truth_pt = "sqrt(truth_px*truth_px + truth_py*truth_py)";
    tree->Draw((truth_pt + ">>h_truth_pt_all").c_str(),
               cut("pt > 0 && truth_px == truth_px && truth_py == truth_py").c_str(),
               "goff");
    tree->Draw(("pt:" + truth_pt + ">>h_reco_vs_truth_pt").c_str(),
               cut("pt > 0 && truth_px == truth_px && truth_py == truth_py").c_str(),
               "goff");
    tree->Draw(("pt - " + truth_pt + ">>h_delta_pt").c_str(),
               cut("pt > 0 && truth_px == truth_px && truth_py == truth_py").c_str(),
               "goff");

    h_truth_pt_all->SetLineColor(kMagenta + 2);
    h_truth_pt_all->SetLineWidth(2);

    c_truth = new TCanvas("c_truth_track_pt", "truth_track_pt", 1000, 800);
    c_truth->SetLogy();
    h_pt_all->Draw("hist");
    h_truth_pt_all->Draw("hist same");
    TLegend *truth_leg = new TLegend(0.58, 0.74, 0.88, 0.88);
    truth_leg->SetBorderSize(0);
    truth_leg->SetFillStyle(0);
    truth_leg->AddEntry(h_pt_all, "reco/fitted p_{T}", "l");
    truth_leg->AddEntry(h_truth_pt_all, "truth p_{T}", "l");
    truth_leg->Draw();
    c_truth->SaveAs((std::string(outbase) + "_truth_overlay.pdf").c_str());
    c_truth->SaveAs((std::string(outbase) + "_truth_overlay.png").c_str());

    c_reco_vs_truth = new TCanvas("c_reco_vs_truth_pt", "reco_vs_truth_pt", 1000, 800);
    c_reco_vs_truth->SetRightMargin(0.15);
    c_reco_vs_truth->SetLogz();
    h_reco_vs_truth_pt->Draw("colz");
    c_reco_vs_truth->SaveAs((std::string(outbase) + "_reco_vs_truth.pdf").c_str());
    c_reco_vs_truth->SaveAs((std::string(outbase) + "_reco_vs_truth.png").c_str());

    c_delta = new TCanvas("c_delta_pt", "delta_pt", 1000, 800);
    c_delta->SetLogy();
    h_delta_pt->Draw("hist");
    c_delta->SaveAs((std::string(outbase) + "_delta_pt.pdf").c_str());
    c_delta->SaveAs((std::string(outbase) + "_delta_pt.png").c_str());
  }

  const std::string root_output = std::string(outbase) + ".root";
  TFile output_file(root_output.c_str(), "RECREATE");
  h_pt_all->Write();
  h_pt_pos->Write();
  h_pt_neg->Write();
  h_pt_zero->Write();
  c_pt->Write("c_track_pt");
  if (has_truth)
  {
    h_truth_pt_all->Write();
    h_reco_vs_truth_pt->Write();
    h_delta_pt->Write();
    c_truth->Write("c_truth_track_pt");
    c_reco_vs_truth->Write("c_reco_vs_truth_pt");
    c_delta->Write("c_delta_pt");
  }
  output_file.Close();

  std::cout << "[draw_track_pt] wrote " << outbase << ".pdf and .png" << std::endl;
  if (has_truth)
  {
    std::cout << "[draw_track_pt] wrote truth overlay/reco-vs-truth/delta plots" << std::endl;
  }
  std::cout << "[draw_track_pt] wrote " << root_output << std::endl;

  file->Close();
}
