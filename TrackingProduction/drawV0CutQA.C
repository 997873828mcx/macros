#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TStyle.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  struct HistDef
  {
    TString name;
    TString title;
    TString expr;
    int nbins{100};
    double xmin{0.0};
    double xmax{1.0};
    bool logy{false};
  };

  TH1F* draw_1d(TTree* tree, const HistDef& def, const TCut& cut)
  {
    auto* hist = new TH1F(def.name, def.title, def.nbins, def.xmin, def.xmax);
    hist->SetLineWidth(2);
    tree->Draw(Form("%s>>%s", def.expr.Data(), def.name.Data()), cut, "hist");
    return hist;
  }
}  // namespace

void drawV0CutQA(
    const char* infile = "output/tpc_pattern_v0_79513_keff1_kalman_merged.root",
    const char* outtag = "V0CutQA_keff1_kalman",
    const char* baseCut = "charge1 != charge2",
    const char* treeName = "pairTree",
    const double vertex_x = 0.0,
    const double vertex_y = 0.0,
    const double vertex_z = 0.0)
{
  std::unique_ptr<TFile> file(TFile::Open(infile, "READ"));
  if (!file || file->IsZombie())
  {
    std::cerr << "[drawV0CutQA] cannot open " << infile << std::endl;
    return;
  }

  TTree* tree = nullptr;
  file->GetObject(treeName, tree);
  if (!tree)
  {
    std::cerr << "[drawV0CutQA] tree " << treeName << " not found in " << infile << std::endl;
    return;
  }

  gStyle->SetOptStat(1110);

  const TString cut_string = baseCut ? baseCut : "";
  const TCut cut(cut_string);
  const TString dx = Form("(pca_x-(%.9g))", vertex_x);
  const TString dy = Form("(pca_y-(%.9g))", vertex_y);
  const TString dz = Form("(pca_z-(%.9g))", vertex_z);
  const TString lxy = Form("sqrt(%s*%s + %s*%s)", dx.Data(), dx.Data(), dy.Data(), dy.Data());
  const TString lxyz = Form("sqrt(%s*%s + %s*%s + %s*%s)",
                            dx.Data(), dx.Data(), dy.Data(), dy.Data(), dz.Data(), dz.Data());
  const TString v0pt = "sqrt(v0_px*v0_px + v0_py*v0_py)";
  const TString cosxy = Form("((%s)*v0_px + (%s)*v0_py)/(%s*%s)",
                             dx.Data(), dy.Data(), lxy.Data(), v0pt.Data());

  std::vector<HistDef> hists = {
      {"h_pt1", ";daughter 1 p_{T} [GeV/c];pairs", "sqrt(px1*px1 + py1*py1)", 120, 0.0, 3.0, true},
      {"h_pt2", ";daughter 2 p_{T} [GeV/c];pairs", "sqrt(px2*px2 + py2*py2)", 120, 0.0, 3.0, true},
      {"h_pt_min", ";min daughter p_{T} [GeV/c];pairs", "min(sqrt(px1*px1 + py1*py1), sqrt(px2*px2 + py2*py2))", 120, 0.0, 3.0, true},
      {"h_pairDCA", ";pair DCA [cm];pairs", "abs(pairDCA)", 120, 0.0, 10.0, true},
      {"h_pca_rxy", ";PCA R_{xy} [cm];pairs", lxy, 160, 0.0, 120.0, true},
      {"h_pca_z", ";PCA z [cm];pairs", "pca_z", 160, -120.0, 120.0, true},
      {"h_abs_pca_z", ";|PCA z| [cm];pairs", "abs(pca_z)", 120, 0.0, 120.0, true},
      {"h_dira3d", ";DIRA = cos#theta_{3D};pairs", "cosThetaReco", 120, -1.0, 1.0, false},
      {"h_dira3d_zoom", ";DIRA = cos#theta_{3D};pairs", "cosThetaReco", 120, 0.8, 1.0, false},
      {"h_lproj3d", ";3D decay length L [cm];pairs", "Lproj", 160, 0.0, 160.0, true},
      {"h_lprojxy", ";XY decay length L_{xy} [cm];pairs", lxy, 160, 0.0, 120.0, true},
      {"h_dira_xy", ";XY pointing cosine;pairs", cosxy, 120, -1.0, 1.0, false},
      {"h_dira_xy_zoom", ";XY pointing cosine;pairs", cosxy, 120, 0.8, 1.0, false},
      {"h_mass_kshort", ";M_{#pi#pi} [GeV/c^{2}];pairs", "mass_Kshort", 160, 0.3, 0.8, true},
      {"h_mass_lambda", ";M_{p#pi^{-}} [GeV/c^{2}];pairs", "mass_Lambda", 160, 1.0, 1.3, true},
      {"h_mass_antilambda", ";M_{#bar{p}#pi^{+}} [GeV/c^{2}];pairs", "mass_AntiLambda", 160, 1.0, 1.3, true},
  };

  const TString pdf = Form("%s.pdf", outtag);
  const TString root = Form("%s.root", outtag);
  std::unique_ptr<TFile> output(TFile::Open(root, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "[drawV0CutQA] cannot create " << root << std::endl;
    return;
  }

  TCanvas canvas("c_v0_cutqa", "V0 cut QA", 1000, 800);
  canvas.Print(pdf + "[");

  std::vector<TH1F*> written_hists;
  written_hists.reserve(hists.size());
  for (const auto& def : hists)
  {
    canvas.Clear();
    canvas.SetLogy(def.logy);
    auto* hist = draw_1d(tree, def, cut);
    hist->Draw("hist");
    canvas.Print(pdf);
    output->cd();
    hist->Write();
    written_hists.push_back(hist);
  }

  canvas.Clear();
  canvas.SetLogy(false);
  canvas.SetLogz(true);
  auto* h_pca_rz = new TH2F("h_pca_rxy_vs_z",
                            ";PCA z [cm];PCA R_{xy} [cm]",
                            160, -120.0, 120.0,
                            160, 0.0, 120.0);
  tree->Draw(Form("%s:pca_z>>h_pca_rxy_vs_z", lxy.Data()), cut, "colz");
  canvas.Print(pdf);
  output->cd();
  h_pca_rz->Write();

  canvas.Clear();
  canvas.SetLogz(true);
  auto* h_dira_pairdca = new TH2F("h_dira_vs_pairDCA",
                                  ";pair DCA [cm];DIRA = cos#theta_{3D}",
                                  120, 0.0, 10.0,
                                  120, -1.0, 1.0);
  tree->Draw("cosThetaReco:abs(pairDCA)>>h_dira_vs_pairDCA", cut, "colz");
  canvas.Print(pdf);
  output->cd();
  h_dira_pairdca->Write();

  canvas.Print(pdf + "]");
  output->Close();

  std::cout << "[drawV0CutQA] entries passing base cut: " << tree->GetEntries(cut_string) << std::endl;
  std::cout << "[drawV0CutQA] wrote " << pdf << std::endl;
  std::cout << "[drawV0CutQA] wrote " << root << std::endl;
}
