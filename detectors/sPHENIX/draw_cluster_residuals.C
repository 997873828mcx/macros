#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <tuple>

void draw_cluster_residuals(
    const char *infile = "output/tpc_pattern_v0_79513_keff1_kalman_samesign_v1_merged.root",
    const char *outbase = "cluster_residuals",
    const int select_run = -999999,
    const int select_evt = -999999,
    const int select_track_id = -999999,
    const int select_side = -1,
    const int select_charge = 0,
    const double min_track_pt = 2.0,
    const int min_npoints = 30,
    const double max_fit_chi2_ndf = -1.0,
    const double z_min = -15.0,
    const double z_max = 0.0,
    const double residual_z_min = -2.0,
    const double residual_z_max = 2.0,
    const double r_min = 20.0,
    const double r_max = 75.0,
    const double residual_rphi_min = -0.8,
    const double residual_rphi_max = 0.8,
    const bool split_charge = false)
{
  TFile *file = TFile::Open(infile, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "[draw_cluster_residuals] cannot open " << infile << std::endl;
    return;
  }

  TTree *cluster_tree = nullptr;
  file->GetObject("clusterResidualTree", cluster_tree);
  if (!cluster_tree)
  {
    std::cerr << "[draw_cluster_residuals] clusterResidualTree not found in "
              << infile << std::endl;
    file->Close();
    return;
  }

  using TrackKey = std::tuple<int, int, int>;
  std::set<TrackKey> accepted_tracks;
  const bool use_track_tree_selection = min_track_pt >= 0.0;

  if (use_track_tree_selection)
  {
    TTree *track_tree = nullptr;
    file->GetObject("trackTree", track_tree);
    if (!track_tree)
    {
      std::cerr << "[draw_cluster_residuals] requested min_track_pt but trackTree is missing"
                << std::endl;
      file->Close();
      return;
    }

    int run = 0;
    int evt = 0;
    int track_id = 0;
    int npoints = 0;
    float pt = 0.0F;
    track_tree->SetBranchAddress("run", &run);
    track_tree->SetBranchAddress("evt", &evt);
    track_tree->SetBranchAddress("track_id", &track_id);
    track_tree->SetBranchAddress("npoints", &npoints);
    track_tree->SetBranchAddress("pt", &pt);

    const Long64_t ntracks = track_tree->GetEntries();
    for (Long64_t i = 0; i < ntracks; ++i)
    {
      track_tree->GetEntry(i);
      if (select_run != -999999 && run != select_run)
      {
        continue;
      }
      if (select_evt != -999999 && evt != select_evt)
      {
        continue;
      }
      if (select_track_id != -999999 && track_id != select_track_id)
      {
        continue;
      }
      if (pt < min_track_pt)
      {
        continue;
      }
      if (min_npoints > 0 && npoints < min_npoints)
      {
        continue;
      }
      accepted_tracks.emplace(run, evt, track_id);
    }

    std::cout << "[draw_cluster_residuals] accepted " << accepted_tracks.size()
              << " tracks from trackTree" << std::endl;
  }

  int run = 0;
  int evt = 0;
  int track_id = 0;
  int side = 0;
  int charge = 0;
  int npoints = 0;
  float cluster_z = 0.0F;
  float cluster_r = 0.0F;
  float residual_z = 0.0F;
  float residual_rphi = 0.0F;
  float fit_chi2_ndf = 0.0F;

  cluster_tree->SetBranchAddress("run", &run);
  cluster_tree->SetBranchAddress("evt", &evt);
  cluster_tree->SetBranchAddress("track_id", &track_id);
  cluster_tree->SetBranchAddress("side", &side);
  cluster_tree->SetBranchAddress("charge", &charge);
  cluster_tree->SetBranchAddress("npoints", &npoints);
  cluster_tree->SetBranchAddress("cluster_z", &cluster_z);
  cluster_tree->SetBranchAddress("cluster_r", &cluster_r);
  cluster_tree->SetBranchAddress("residual_z", &residual_z);
  cluster_tree->SetBranchAddress("residual_rphi", &residual_rphi);
  cluster_tree->SetBranchAddress("fit_chi2_ndf", &fit_chi2_ndf);

  TH2F *h_residual_z = new TH2F(
      "h_residual_z_vs_cluster_z",
      ";cluster z [cm];z residual = z_{cluster} - z_{fit} [cm]",
      240, z_min, z_max,
      240, residual_z_min, residual_z_max);

  TH2F *h_residual_rphi = new TH2F(
      "h_residual_rphi_vs_cluster_r",
      ";cluster r [cm];r#phi residual [cm]",
      240, r_min, r_max,
      240, residual_rphi_min, residual_rphi_max);

  TH2F *h_residual_z_pos = nullptr;
  TH2F *h_residual_z_neg = nullptr;
  TH2F *h_residual_rphi_pos = nullptr;
  TH2F *h_residual_rphi_neg = nullptr;
  if (split_charge)
  {
    h_residual_z_pos = new TH2F(
        "h_residual_z_vs_cluster_z_pos",
        "positive charge;cluster z [cm];z residual = z_{cluster} - z_{fit} [cm]",
        240, z_min, z_max,
        240, residual_z_min, residual_z_max);
    h_residual_z_neg = new TH2F(
        "h_residual_z_vs_cluster_z_neg",
        "negative charge;cluster z [cm];z residual = z_{cluster} - z_{fit} [cm]",
        240, z_min, z_max,
        240, residual_z_min, residual_z_max);
    h_residual_rphi_pos = new TH2F(
        "h_residual_rphi_vs_cluster_r_pos",
        "positive charge;cluster r [cm];r#phi residual [cm]",
        240, r_min, r_max,
        240, residual_rphi_min, residual_rphi_max);
    h_residual_rphi_neg = new TH2F(
        "h_residual_rphi_vs_cluster_r_neg",
        "negative charge;cluster r [cm];r#phi residual [cm]",
        240, r_min, r_max,
        240, residual_rphi_min, residual_rphi_max);
  }

  Long64_t filled = 0;
  Long64_t filled_positive = 0;
  Long64_t filled_negative = 0;
  const Long64_t nclusters = cluster_tree->GetEntries();
  for (Long64_t i = 0; i < nclusters; ++i)
  {
    cluster_tree->GetEntry(i);

    if (select_run != -999999 && run != select_run)
    {
      continue;
    }
    if (select_evt != -999999 && evt != select_evt)
    {
      continue;
    }
    if (select_track_id != -999999 && track_id != select_track_id)
    {
      continue;
    }
    if (select_side != 0 && side != select_side)
    {
      continue;
    }
    if (select_charge != 0 && charge != select_charge)
    {
      continue;
    }
    if (min_npoints > 0 && npoints < min_npoints)
    {
      continue;
    }
    if (max_fit_chi2_ndf >= 0.0 &&
        (!std::isfinite(fit_chi2_ndf) || fit_chi2_ndf > max_fit_chi2_ndf))
    {
      continue;
    }
    if (use_track_tree_selection &&
        accepted_tracks.find(TrackKey(run, evt, track_id)) == accepted_tracks.end())
    {
      continue;
    }

    h_residual_z->Fill(cluster_z, residual_z);
    h_residual_rphi->Fill(cluster_r, residual_rphi);
    if (split_charge && charge > 0)
    {
      h_residual_z_pos->Fill(cluster_z, residual_z);
      h_residual_rphi_pos->Fill(cluster_r, residual_rphi);
      ++filled_positive;
    }
    else if (split_charge && charge < 0)
    {
      h_residual_z_neg->Fill(cluster_z, residual_z);
      h_residual_rphi_neg->Fill(cluster_r, residual_rphi);
      ++filled_negative;
    }
    ++filled;
  }

  gStyle->SetOptStat(0);

  TCanvas *canvas = new TCanvas("c_cluster_residuals", "cluster residuals", 1600, 750);
  canvas->Divide(2, 1);
  canvas->cd(1);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  h_residual_z->Draw("colz");
  canvas->cd(2);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  h_residual_rphi->Draw("colz");

  const std::string base = outbase ? outbase : "cluster_residuals";
  canvas->SaveAs((base + ".pdf").c_str());
  canvas->SaveAs((base + ".png").c_str());

  TCanvas *charge_canvas = nullptr;
  if (split_charge)
  {
    charge_canvas = new TCanvas("c_cluster_residuals_by_charge", "cluster residuals by charge", 1600, 1400);
    charge_canvas->Divide(2, 2);
    charge_canvas->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    h_residual_z_pos->Draw("colz");
    charge_canvas->cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    h_residual_z_neg->Draw("colz");
    charge_canvas->cd(3);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    h_residual_rphi_pos->Draw("colz");
    charge_canvas->cd(4);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    h_residual_rphi_neg->Draw("colz");
    charge_canvas->SaveAs((base + "_by_charge.pdf").c_str());
    charge_canvas->SaveAs((base + "_by_charge.png").c_str());
  }

  TFile output_file((base + ".root").c_str(), "RECREATE");
  h_residual_z->Write();
  h_residual_rphi->Write();
  if (split_charge)
  {
    h_residual_z_pos->Write();
    h_residual_z_neg->Write();
    h_residual_rphi_pos->Write();
    h_residual_rphi_neg->Write();
  }
  canvas->Write();
  if (charge_canvas)
  {
    charge_canvas->Write();
  }
  output_file.Close();

  std::cout << "[draw_cluster_residuals] filled " << filled
            << " clusters from " << infile << std::endl;
  if (split_charge)
  {
    std::cout << "[draw_cluster_residuals] positive clusters=" << filled_positive
              << ", negative clusters=" << filled_negative << std::endl;
  }
  std::cout << "[draw_cluster_residuals] wrote " << base << ".pdf, "
            << base << ".png, and " << base << ".root" << std::endl;
  if (split_charge)
  {
    std::cout << "[draw_cluster_residuals] wrote " << base
              << "_by_charge.pdf and " << base << "_by_charge.png" << std::endl;
  }

  file->Close();
}
