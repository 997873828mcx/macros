#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>
#include <string>

void draw_kshort_mass_pt_ratio(
    const char *input_pattern =
        "output/tpc_pattern_v0_79513_new_8_380_no_mv_rot_kalman_analytic_selected_vtx0p158_0p285_v1_all1000_shards/part_*.root",
    const char *outbase = "output/kshort_mass_daughter_pt_79513_new_8_380",
    const char *extra_cut = "",
    const char *topology = "all",
    const double bfield_z_t = 1.4,
    const char *selection = "auto")
{
  TChain pair_tree("pairTree");
  const int nfiles = pair_tree.Add(input_pattern);
  if (nfiles <= 0)
  {
    std::cerr << "[draw_kshort_mass_pt_ratio] no files matched "
              << input_pattern << std::endl;
    return;
  }

  const Long64_t nentries = pair_tree.GetEntries();
  if (nentries <= 0)
  {
    std::cerr << "[draw_kshort_mass_pt_ratio] pairTree is empty" << std::endl;
    return;
  }
  if (pair_tree.LoadTree(0) < 0)
  {
    std::cerr << "[draw_kshort_mass_pt_ratio] cannot load the first pairTree entry"
              << std::endl;
    return;
  }

  const char *required_branches[] = {
      "charge1", "charge2", "px1", "py1", "px2", "py2", "pairDCA", "cosThetaReco"};
  for (const char *branch : required_branches)
  {
    if (!pair_tree.GetBranch(branch))
    {
      std::cerr << "[draw_kshort_mass_pt_ratio] required branch missing: "
                << branch << std::endl;
      return;
    }
  }

  const bool has_stored_mass = pair_tree.GetBranch("mass_Kshort");
  if (!has_stored_mass &&
      (!pair_tree.GetBranch("pz1") || !pair_tree.GetBranch("pz2")))
  {
    std::cerr << "[draw_kshort_mass_pt_ratio] mass_Kshort is absent and pz1/pz2 "
                 "are unavailable for calculating it"
              << std::endl;
    return;
  }

  pair_tree.SetAlias("pt1", "sqrt(px1*px1+py1*py1)");
  pair_tree.SetAlias("pt2", "sqrt(px2*px2+py2*py2)");
  pair_tree.SetAlias("ptplus", "(charge1>0)*pt1+(charge2>0)*pt2");
  pair_tree.SetAlias("ptminus", "(charge1<0)*pt1+(charge2<0)*pt2");
  pair_tree.SetAlias("ptsum", "ptplus+ptminus");
  pair_tree.SetAlias("ptratio", "ptplus/ptminus");
  pair_tree.SetAlias("ptasymmetry", "(ptplus-ptminus)/ptsum");
  if (has_stored_mass)
  {
    pair_tree.SetAlias("mass_kshort_plot", "mass_Kshort");
  }
  else
  {
    pair_tree.SetAlias(
        "pion_energy1",
        "sqrt(px1*px1+py1*py1+pz1*pz1+0.13957039*0.13957039)");
    pair_tree.SetAlias(
        "pion_energy2",
        "sqrt(px2*px2+py2*py2+pz2*pz2+0.13957039*0.13957039)");
    pair_tree.SetAlias(
        "mass_kshort_plot",
        "sqrt((pion_energy1+pion_energy2)*(pion_energy1+pion_energy2)"
        "-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2))");
  }
  pair_tree.SetAlias("pxplus", "(charge1>0)*px1+(charge2>0)*px2");
  pair_tree.SetAlias("pyplus", "(charge1>0)*py1+(charge2>0)*py2");
  pair_tree.SetAlias("pxminus", "(charge1<0)*px1+(charge2<0)*px2");
  pair_tree.SetAlias("pyminus", "(charge1<0)*py1+(charge2<0)*py2");
  const std::string sailor_metric =
      std::to_string(bfield_z_t) +
      "*(pxplus*pyminus-pyplus*pxminus)";
  pair_tree.SetAlias("sailor_metric", sailor_metric.c_str());

  TCut cuts =
      "charge1*charge2 < 0"
      " && pt1 > 0.2"
      " && pt2 > 0.2"
      " && abs(pairDCA) < 1.0"
      " && cosThetaReco > 0.88"
      " && ptsum > 0";

  const bool has_modern_selection =
      pair_tree.GetBranch("pca_x") && pair_tree.GetBranch("pca_y") &&
      pair_tree.GetBranch("pca_z") && pair_tree.GetBranch("pca1_z") &&
      pair_tree.GetBranch("pca2_z") && pair_tree.GetBranch("npoints1") &&
      pair_tree.GetBranch("npoints2");
  const std::string schema_mode = has_modern_selection ? "modern" : "legacy";
  std::string selection_mode = selection ? selection : "auto";
  if (selection_mode == "auto")
  {
    selection_mode = has_modern_selection ? "modern" : "common";
  }

  if (selection_mode == "modern")
  {
    if (!has_modern_selection)
    {
      std::cerr << "[draw_kshort_mass_pt_ratio] modern selection requested, but "
                   "the required PCA/npoints branches are absent"
                << std::endl;
      return;
    }
    cuts = cuts && TCut(
                       "abs(pca_z) < 15.0"
                       " && sqrt(pca_x*pca_x+pca_y*pca_y) > 2.0"
                       " && npoints1 > 30"
                       " && npoints2 > 30"
                       " && abs(pca1_z-pca2_z) < 0.5");
  }
  else if (selection_mode == "common")
  {
    if (!pair_tree.GetBranch("cross1") || !pair_tree.GetBranch("cross2") ||
        !pair_tree.GetBranch("Lproj"))
    {
      std::cerr << "[draw_kshort_mass_pt_ratio] common selection requires "
                   "cross1, cross2, and Lproj"
                << std::endl;
      return;
    }
    cuts = cuts && TCut("cross1 == cross2 && Lproj > 2.0");
  }
  else
  {
    std::cerr << "[draw_kshort_mass_pt_ratio] selection must be auto, modern, or common"
              << std::endl;
    return;
  }

  const std::string topology_mode = topology ? topology : "all";
  if (topology_mode == "sailor")
  {
    if (bfield_z_t == 0.0)
    {
      std::cerr << "[draw_kshort_mass_pt_ratio] sailor topology requires nonzero Bz"
                << std::endl;
      return;
    }
    cuts = cuts && TCut("sailor_metric > 0");
  }
  else if (topology_mode == "cowboy")
  {
    if (bfield_z_t == 0.0)
    {
      std::cerr << "[draw_kshort_mass_pt_ratio] cowboy topology requires nonzero Bz"
                << std::endl;
      return;
    }
    cuts = cuts && TCut("sailor_metric < 0");
  }
  else if (topology_mode != "all")
  {
    std::cerr << "[draw_kshort_mass_pt_ratio] topology must be all, sailor, or cowboy"
              << std::endl;
    return;
  }

  if (extra_cut && std::string(extra_cut).size() > 0)
  {
    cuts = cuts && TCut(extra_cut);
  }

  gROOT->cd();
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);

  auto *asymmetry_hist = new TH2D(
      "h_kshort_mass_vs_pt_asymmetry",
      ";(p_{T}^{+}-p_{T}^{-})/(p_{T}^{+}+p_{T}^{-});M_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
      120, -1.0, 1.0,
      120, 0.35, 0.65);

  auto *ratio_hist = new TH2D(
      "h_kshort_mass_vs_pt_ratio",
      ";p_{T}^{+}/p_{T}^{-};M_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
      120, 0.0, 6.0,
      120, 0.35, 0.65);

  auto *asymmetry_canvas = new TCanvas(
      "c_kshort_mass_vs_pt_asymmetry",
      "K-short mass versus daughter pT asymmetry",
      1100, 850);
  asymmetry_canvas->SetRightMargin(0.15);

  const Long64_t selected_asymmetry = pair_tree.Draw(
      "mass_kshort_plot:ptasymmetry>>h_kshort_mass_vs_pt_asymmetry",
      cuts,
      "colz");
  asymmetry_canvas->SetLogz();

  auto *ratio_canvas = new TCanvas(
      "c_kshort_mass_vs_pt_ratio",
      "K-short mass versus daughter pT ratio",
      1100, 850);
  ratio_canvas->SetRightMargin(0.15);

  const Long64_t selected_ratio = pair_tree.Draw(
      "mass_kshort_plot:ptratio>>h_kshort_mass_vs_pt_ratio",
      cuts,
      "colz");
  ratio_canvas->SetLogz();

  const std::string output_base = outbase;
  asymmetry_canvas->SaveAs((output_base + "_pt_asymmetry.pdf").c_str());
  asymmetry_canvas->SaveAs((output_base + "_pt_asymmetry.png").c_str());
  ratio_canvas->SaveAs((output_base + "_ptplus_over_ptminus.pdf").c_str());
  ratio_canvas->SaveAs((output_base + "_ptplus_over_ptminus.png").c_str());

  TFile output_file((output_base + ".root").c_str(), "RECREATE");
  asymmetry_hist->Write();
  ratio_hist->Write();
  asymmetry_canvas->Write();
  ratio_canvas->Write();
  output_file.Close();

  std::cout << "[draw_kshort_mass_pt_ratio] files: " << nfiles << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] pairTree entries: " << nentries << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] schema: " << schema_mode
            << (has_stored_mass ? " (stored mass)" : " (calculated pion-pair mass)")
            << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] selection: " << selection_mode << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] selected pairs (asymmetry): "
            << selected_asymmetry << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] selected pairs (ratio): "
            << selected_ratio << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] topology: " << topology_mode
            << " (Bz = " << bfield_z_t << " T)" << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] cuts: " << cuts.GetTitle() << std::endl;
  std::cout << "[draw_kshort_mass_pt_ratio] wrote " << output_base
            << "_pt_asymmetry.pdf/.png, " << output_base
            << "_ptplus_over_ptminus.pdf/.png, and " << output_base
            << ".root" << std::endl;
}
