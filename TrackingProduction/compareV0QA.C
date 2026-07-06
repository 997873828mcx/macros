// =============================================================================
// compareV0QA.C
// -----------------------------------------------------------------------------
// Compare two V0/AP pairTree ROOT files, for example helix vs Kalman.
//
// Usage:
//   root -l -b -q 'compareV0QA.C("helix.root", "kalman.root", "helix", "kalman", "v0qa_compare")'
//
// Optional final argument is a ROOT TTree selection string applied to both files:
//   root -l -b -q 'compareV0QA.C("h.root","k.root","helix","kalman","qa","pairDCA<0.15 && cosThetaReco>0.94")'
//
// Outputs:
//   <outPrefix>.pdf
//   <outPrefix>.root
//   <outPrefix>_summary.txt
// =============================================================================

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  struct Sample
  {
    std::unique_ptr<TFile> file;
    TTree* tree = nullptr;
    TString label;
    int color = kBlack;
  };

  struct HistSpec
  {
    TString key;
    TString expression;
    TString finiteExpression;
    TString title;
    int bins = 100;
    double xmin = 0.0;
    double xmax = 1.0;
    bool logy = false;
  };

  TString sanitize(TString text)
  {
    text.ReplaceAll("/", "_");
    text.ReplaceAll(" ", "_");
    text.ReplaceAll("-", "_");
    text.ReplaceAll("+", "p");
    text.ReplaceAll("(", "");
    text.ReplaceAll(")", "");
    text.ReplaceAll("[", "");
    text.ReplaceAll("]", "");
    text.ReplaceAll(":", "_");
    return text;
  }

  bool has_branch(TTree* tree, const TString& branch)
  {
    return tree && tree->GetBranch(branch.Data());
  }

  bool has_needed_branch(TTree* tree, const HistSpec& spec)
  {
    // This intentionally checks the simple finite branch, not the full draw
    // expression, so expressions like abs(delta_alpha) work.
    return has_branch(tree, spec.finiteExpression);
  }

  TString make_selection(const TString& finiteExpression, const TString& userCut)
  {
    TString selection = Form("(%s==%s)", finiteExpression.Data(), finiteExpression.Data());
    if (!userCut.IsNull())
    {
      selection += "&&(";
      selection += userCut;
      selection += ")";
    }
    return selection;
  }

  Sample open_sample(const char* filename, const char* label, const int color)
  {
    Sample sample;
    sample.file.reset(TFile::Open(filename, "READ"));
    sample.label = label;
    sample.color = color;
    if (!sample.file || sample.file->IsZombie())
    {
      std::cerr << "[compareV0QA] cannot open " << filename << std::endl;
      return sample;
    }
    sample.file->GetObject("pairTree", sample.tree);
    if (!sample.tree)
    {
      std::cerr << "[compareV0QA] pairTree not found in " << filename << std::endl;
    }
    return sample;
  }

  TH1D* draw_hist(const Sample& sample, const HistSpec& spec, const TString& userCut, const bool normalize)
  {
    const TString name = "h_" + sanitize(sample.label) + "_" + sanitize(spec.key);
    if (gDirectory && gDirectory->Get(name))
    {
      delete gDirectory->Get(name);
    }
    auto* hist = new TH1D(name, spec.title, spec.bins, spec.xmin, spec.xmax);
    hist->Sumw2();
    hist->SetLineColor(sample.color);
    hist->SetMarkerColor(sample.color);
    hist->SetLineWidth(2);
    hist->SetStats(false);

    if (!sample.tree || !has_needed_branch(sample.tree, spec))
    {
      hist->SetDirectory(nullptr);
      return hist;
    }

    const TString draw = spec.expression + ">>" + name;
    sample.tree->Draw(draw, make_selection(spec.finiteExpression, userCut), "goff");
    hist->SetDirectory(nullptr);
    if (normalize && hist->Integral() > 0.0)
    {
      hist->Scale(1.0 / hist->Integral());
      hist->GetYaxis()->SetTitle("Normalized candidates");
    }
    else
    {
      hist->GetYaxis()->SetTitle("Candidates");
    }
    return hist;
  }

  void write_summary(std::ofstream& out, const TString& variable, const Sample& sample, TH1D* hist)
  {
    double probs[5] = {0.16, 0.50, 0.68, 0.84, 0.95};
    double q[5] = {0, 0, 0, 0, 0};
    if (hist && hist->GetEntries() > 0)
    {
      hist->GetQuantiles(5, q, probs);
    }
    const double half68 = 0.5 * (q[3] - q[0]);
    out << std::left << std::setw(22) << variable
        << std::setw(16) << sample.label
        << std::right << std::setw(12) << (hist ? hist->GetEntries() : 0)
        << std::setw(14) << (hist ? hist->GetMean() : 0.0)
        << std::setw(14) << (hist ? hist->GetRMS() : 0.0)
        << std::setw(14) << q[1]
        << std::setw(14) << half68
        << std::setw(14) << q[2]
        << std::setw(14) << q[4]
        << "\n";
  }

  void draw_overlay_page(TCanvas& canvas,
                         const HistSpec& spec,
                         const Sample& sampleA,
                         const Sample& sampleB,
                         TH1D* histA,
                         TH1D* histB,
                         const TString& pdf)
  {
    canvas.Clear();
    canvas.SetLogy(spec.logy);
    histA->SetTitle(spec.title);
    histA->GetXaxis()->SetTitle(spec.title.Contains(";") ? "" : spec.key);

    const double ymax = std::max(histA->GetMaximum(), histB->GetMaximum());
    histA->SetMaximum(ymax > 0.0 ? ymax * (spec.logy ? 10.0 : 1.25) : 1.0);
    if (spec.logy)
    {
      histA->SetMinimum(1.0e-6);
    }
    histA->Draw("hist");
    histB->Draw("hist same");

    TLegend legend(0.58, 0.72, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(histA, Form("%s, N=%.0f", sampleA.label.Data(), histA->GetEntries()), "l");
    legend.AddEntry(histB, Form("%s, N=%.0f", sampleB.label.Data(), histB->GetEntries()), "l");
    legend.Draw();

    canvas.Print(pdf);
  }

  TH2D* draw_ap_hist(const Sample& sample, const TString& userCut)
  {
    const TString name = "h_ap_" + sanitize(sample.label);
    if (gDirectory && gDirectory->Get(name))
    {
      delete gDirectory->Get(name);
    }
    auto* hist = new TH2D(name, Form("%s;#alpha;q_{T} (GeV/c)", sample.label.Data()),
                          200, -1.0, 1.0, 200, 0.0, 0.25);
    hist->SetStats(false);
    if (sample.tree && has_branch(sample.tree, "alpha") && has_branch(sample.tree, "qT"))
    {
      sample.tree->Draw("qT:alpha>>" + name, make_selection("alpha", userCut) + "&&(qT==qT)", "goff");
    }
    hist->SetDirectory(nullptr);
    return hist;
  }

  void draw_ap_page(TCanvas& canvas,
                    const Sample& sampleA,
                    const Sample& sampleB,
                    TH2D* apA,
                    TH2D* apB,
                    const TString& pdf)
  {
    canvas.Clear();
    canvas.Divide(2, 1);
    canvas.cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    apA->Draw("colz");
    canvas.cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    apB->Draw("colz");
    canvas.Print(pdf);
    (void) sampleA;
    (void) sampleB;
  }
}  // namespace

void compareV0QA(const char* fileA,
                 const char* fileB,
                 const char* labelA = "helix",
                 const char* labelB = "kalman",
                 const char* outPrefix = "v0qa_compare",
                 const char* cut = "")
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  auto sampleA = open_sample(fileA, labelA, kBlue + 1);
  auto sampleB = open_sample(fileB, labelB, kRed + 1);
  if (!sampleA.tree || !sampleB.tree)
  {
    return;
  }

  const TString userCut = cut;
  const TString pdf = TString(outPrefix) + ".pdf";
  const TString rootName = TString(outPrefix) + ".root";
  const TString summaryName = TString(outPrefix) + "_summary.txt";

  std::vector<HistSpec> specs = {
      {"pca_to_true_3d", "pca_to_true_3d", "pca_to_true_3d", "PCA to true decay vertex;|#Delta r_{3D}| (cm);", 120, 0.0, 5.0, true},
      {"pca_to_true_xy", "pca_to_true_xy", "pca_to_true_xy", "PCA to true decay vertex, xy;|#Delta r_{xy}| (cm);", 120, 0.0, 5.0, true},
      {"pca_to_true_z", "pca_to_true_z", "pca_to_true_z", "PCA to true decay vertex, z;|#Delta z| (cm);", 120, 0.0, 5.0, true},
      {"pairDCA", "pairDCA", "pairDCA", "Daughter pair DCA;pair DCA (cm);", 120, 0.0, 2.0, true},
      {"abs_delta_alpha", "abs(delta_alpha)", "delta_alpha", "AP residual;|#Delta#alpha|;", 120, 0.0, 0.5, true},
      {"abs_delta_qT", "abs(delta_qT)", "delta_qT", "AP residual;|#Delta q_{T}| (GeV/c);", 120, 0.0, 0.10, true},
      {"cos_mom1_truth", "cos_mom1_truth", "cos_mom1_truth", "Momentum direction agreement, daughter 1;cos(#vec{p}_{reco},#vec{p}_{true});", 120, 0.90, 1.001, true},
      {"cos_mom2_truth", "cos_mom2_truth", "cos_mom2_truth", "Momentum direction agreement, daughter 2;cos(#vec{p}_{reco},#vec{p}_{true});", 120, 0.90, 1.001, true},
      {"mass_Kshort", "mass_Kshort", "mass_Kshort", "K^{0}_{S} hypothesis mass;m_{#pi#pi} (GeV/c^{2});", 140, 0.35, 0.65, false},
      {"mass_Lambda", "mass_Lambda", "mass_Lambda", "#Lambda hypothesis mass;m_{p#pi} (GeV/c^{2});", 140, 1.05, 1.25, false},
      {"mass_AntiLambda", "mass_AntiLambda", "mass_AntiLambda", "#bar{#Lambda} hypothesis mass;m_{#bar{p}#pi} (GeV/c^{2});", 140, 1.05, 1.25, false},
  };

  TCanvas canvas("c_v0qa", "V0 QA comparison", 1100, 800);
  canvas.Print(pdf + "[");

  TFile output(rootName, "RECREATE");
  std::ofstream summary(summaryName.Data());
  summary << "# cut: " << (userCut.IsNull() ? "<none>" : userCut.Data()) << "\n";
  summary << std::left << std::setw(22) << "# variable"
          << std::setw(16) << "sample"
          << std::right << std::setw(12) << "entries"
          << std::setw(14) << "mean"
          << std::setw(14) << "rms"
          << std::setw(14) << "median"
          << std::setw(14) << "half68"
          << std::setw(14) << "q68"
          << std::setw(14) << "q95"
          << "\n";

  std::vector<TH1*> writtenHists;
  for (const auto& spec : specs)
  {
    if (!has_needed_branch(sampleA.tree, spec) || !has_needed_branch(sampleB.tree, spec))
    {
      std::cout << "[compareV0QA] skip missing branch/expression: " << spec.key << std::endl;
      continue;
    }

    TH1D* histA = draw_hist(sampleA, spec, userCut, true);
    TH1D* histB = draw_hist(sampleB, spec, userCut, true);
    draw_overlay_page(canvas, spec, sampleA, sampleB, histA, histB, pdf);
    write_summary(summary, spec.key, sampleA, histA);
    write_summary(summary, spec.key, sampleB, histB);
    writtenHists.push_back(histA);
    writtenHists.push_back(histB);
  }

  TH2D* apA = draw_ap_hist(sampleA, userCut);
  TH2D* apB = draw_ap_hist(sampleB, userCut);
  draw_ap_page(canvas, sampleA, sampleB, apA, apB, pdf);
  writtenHists.push_back(apA);
  writtenHists.push_back(apB);

  canvas.Print(pdf + "]");

  output.cd();
  for (auto* hist : writtenHists)
  {
    if (hist)
    {
      hist->Write();
    }
  }
  output.Close();
  summary.close();

  std::cout << "[compareV0QA] wrote " << pdf << std::endl;
  std::cout << "[compareV0QA] wrote " << rootName << std::endl;
  std::cout << "[compareV0QA] wrote " << summaryName << std::endl;
}
