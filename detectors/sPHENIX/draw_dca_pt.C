#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TTree.h>

#include <iostream>
#include <string>

void draw_dca_pt(
    const char *infile = "output/tpc_pattern_v0_kalman_79513_good_auto_order_partial177_merged.root",
    const char *outbase = "dca_xy_vs_pt_kalman",
    const char *tree_name = "trackTree",
    const double bfield_t = 1.4,
    const double rdca_y_min = -30.0,
    const double rdca_y_max = 30.0,
    const double qoverpt_x_min = -20.0,
    const double qoverpt_x_max = 20.0,
    const char *extra_cut = "",
    const char *rdca_mode = "auto")
{
  TFile *file = TFile::Open(infile, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "[draw_dca_pt] cannot open " << infile << std::endl;
    return;
  }

  TTree *tree = nullptr;
  file->GetObject(tree_name, tree);
  if (!tree)
  {
    std::cerr << "[draw_dca_pt] tree " << tree_name << " not found in " << infile << std::endl;
    file->Close();
    return;
  }

  gStyle->SetOptStat(0);

  const std::string extra_cut_string = extra_cut ? extra_cut : "";
  auto with_extra_cut = [&extra_cut_string](const std::string &base_cut)
  {
    if (extra_cut_string.empty())
    {
      return base_cut;
    }
    return "(" + base_cut + ") && (" + extra_cut_string + ")";
  };

  TCanvas *canvas = new TCanvas("canvas", "dca_xy_vs_pt", 1000, 800);
  TH2F *hist = new TH2F(
      "h_dca_xy_vs_pt",
      ";p_{T} [GeV/c];DCA_{xy} [cm]",
      240, 0.0, 6.0,
      300, 0.0, 200.0);

  tree->Draw("dca_xy:pt>>h_dca_xy_vs_pt", with_extra_cut("pt > 0 && dca_xy >= 0").c_str(), "colz");
  canvas->SetRightMargin(0.15);
  canvas->SetLogz();
  canvas->SaveAs((std::string(outbase) + ".pdf").c_str());
  canvas->SaveAs((std::string(outbase) + ".png").c_str());

  const std::string qpt_base = std::string(outbase) + "_qpt";
  TCanvas *canvas_qpt = new TCanvas("canvas_qpt", "dca_xy_vs_qpt", 1000, 800);
  TH2F *hist_qpt = new TH2F(
      "h_dca_xy_vs_qpt",
      ";q #times p_{T} [GeV/c];DCA_{xy} [cm]",
      240, -6.0, 6.0,
      300, 0.0, 200.0);

  tree->Draw("dca_xy:(charge*pt)>>h_dca_xy_vs_qpt",
             with_extra_cut("pt > 0 && charge != 0 && dca_xy >= 0").c_str(), "colz");
  canvas_qpt->SetRightMargin(0.15);
  canvas_qpt->SetLogz();
  canvas_qpt->SaveAs((qpt_base + ".pdf").c_str());
  canvas_qpt->SaveAs((qpt_base + ".png").c_str());

  const std::string rdca_base = std::string(outbase) + "_signed_rDCA";
  TCanvas *canvas_rdca = new TCanvas("canvas_rdca", "rdca_vs_pt", 1000, 800);
  TH2F *hist_rdca = new TH2F(
      "h_rdca_vs_pt",
      ";p_{T} [GeV/c];rDCA = |C_{xy}| - R [cm]",
      240, 0.0, 6.0,
      300, rdca_y_min, rdca_y_max);

  const double kappa = 0.003 * bfield_t;
  const std::string rdca_mode_string = rdca_mode ? rdca_mode : "auto";
  const bool has_helix_branches =
      tree->GetBranch("helix_cx") &&
      tree->GetBranch("helix_cy") &&
      tree->GetBranch("helix_radius");
  const bool has_kalman_branches =
      tree->GetBranch("kalman_cx") &&
      tree->GetBranch("kalman_cy") &&
      tree->GetBranch("kalman_radius");
  const std::string state_rdca_expr =
      "sqrt((x + py/(" + std::to_string(kappa) + "*charge))*(x + py/(" + std::to_string(kappa) + "*charge))"
      " + (y - px/(" + std::to_string(kappa) + "*charge))*(y - px/(" + std::to_string(kappa) + "*charge)))"
      " - pt/" + std::to_string(kappa);
  std::string rdca_expr = state_rdca_expr;
  if (rdca_mode_string == "kalman")
  {
    if (!has_kalman_branches)
    {
      std::cerr << "[draw_dca_pt] requested rdca_mode=kalman but Kalman circle branches are missing" << std::endl;
      file->Close();
      return;
    }
    rdca_expr = "sqrt(kalman_cx*kalman_cx + kalman_cy*kalman_cy) - kalman_radius";
  }
  else if (rdca_mode_string == "helix" || (rdca_mode_string == "auto" && has_helix_branches))
  {
    rdca_expr = "sqrt(helix_cx*helix_cx + helix_cy*helix_cy) - helix_radius";
  }

  tree->Draw((rdca_expr + ":pt>>h_rdca_vs_pt").c_str(),
             with_extra_cut("pt > 0 && charge != 0").c_str(), "colz");
  canvas_rdca->SetRightMargin(0.15);
  canvas_rdca->SetLogz();
  canvas_rdca->SaveAs((rdca_base + ".pdf").c_str());
  canvas_rdca->SaveAs((rdca_base + ".png").c_str());

  const std::string rdca_qpt_base = std::string(outbase) + "_signed_rDCA_qpt";
  TCanvas *canvas_rdca_qpt = new TCanvas("canvas_rdca_qpt", "rdca_vs_qpt", 1000, 800);
  TH2F *hist_rdca_qpt = new TH2F(
      "h_rdca_vs_qpt",
      ";q #times p_{T} [GeV/c];rDCA = |C_{xy}| - R [cm]",
      240, -6.0, 6.0,
      300, rdca_y_min, rdca_y_max);

  tree->Draw((rdca_expr + ":(charge*pt)>>h_rdca_vs_qpt").c_str(),
             with_extra_cut("pt > 0 && charge != 0").c_str(), "colz");
  canvas_rdca_qpt->SetRightMargin(0.15);
  canvas_rdca_qpt->SetLogz();
  canvas_rdca_qpt->SaveAs((rdca_qpt_base + ".pdf").c_str());
  canvas_rdca_qpt->SaveAs((rdca_qpt_base + ".png").c_str());

  const std::string rdca_qoverpt_base = std::string(outbase) + "_signed_rDCA_qoverpt";
  TCanvas *canvas_rdca_qoverpt = new TCanvas("canvas_rdca_qoverpt", "rdca_vs_qoverpt", 1000, 800);
  TH2F *hist_rdca_qoverpt = new TH2F(
      "h_rdca_vs_qoverpt",
      ";q / p_{T} [(GeV/c)^{-1}];rDCA = |C_{xy}| - R [cm]",
      240, qoverpt_x_min, qoverpt_x_max,
      300, rdca_y_min, rdca_y_max);

  tree->Draw((rdca_expr + ":(charge/pt)>>h_rdca_vs_qoverpt").c_str(),
             with_extra_cut("pt > 0 && charge != 0").c_str(), "colz");
  canvas_rdca_qoverpt->SetRightMargin(0.15);
  canvas_rdca_qoverpt->SetLogz();
  canvas_rdca_qoverpt->SaveAs((rdca_qoverpt_base + ".pdf").c_str());
  canvas_rdca_qoverpt->SaveAs((rdca_qoverpt_base + ".png").c_str());

  const std::string root_output = std::string(outbase) + ".root";
  TFile output_file(root_output.c_str(), "RECREATE");
  hist->Write();
  hist_qpt->Write();
  hist_rdca->Write();
  hist_rdca_qpt->Write();
  hist_rdca_qoverpt->Write();
  canvas->Write("c_dca_xy_vs_pt");
  canvas_qpt->Write("c_dca_xy_vs_qpt");
  canvas_rdca->Write("c_rdca_vs_pt");
  canvas_rdca_qpt->Write("c_rdca_vs_qpt");
  canvas_rdca_qoverpt->Write("c_rdca_vs_qoverpt");
  output_file.Close();

  std::cout << "[draw_dca_pt] wrote " << outbase << ".pdf and " << outbase << ".png" << std::endl;
  std::cout << "[draw_dca_pt] wrote " << qpt_base << ".pdf and " << qpt_base << ".png" << std::endl;
  std::cout << "[draw_dca_pt] wrote " << rdca_base << ".pdf and " << rdca_base << ".png" << std::endl;
  std::cout << "[draw_dca_pt] wrote " << rdca_qpt_base << ".pdf and " << rdca_qpt_base << ".png" << std::endl;
  std::cout << "[draw_dca_pt] wrote " << rdca_qoverpt_base << ".pdf and " << rdca_qoverpt_base << ".png" << std::endl;
  std::cout << "[draw_dca_pt] wrote " << root_output << std::endl;
  file->Close();
}
