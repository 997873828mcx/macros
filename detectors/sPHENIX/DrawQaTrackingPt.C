#include <TAxis.h>
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>
#include <memory>

namespace
{
  void draw_stats(TH1* hist, double x = 0.58, double y = 0.88)
  {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(x, y, Form("Entries: %.0f", hist->GetEntries()));
    latex.DrawLatex(x, y - 0.05, Form("Mean: %.4g", hist->GetMean()));
    latex.DrawLatex(x, y - 0.10, Form("RMS: %.4g", hist->GetRMS()));
  }

  TString build_output_stem(const TString& infile, const TString& histname)
  {
    TString stem = gSystem->ConcatFileName(gSystem->DirName(infile), gSystem->BaseName(infile));
    if (stem.EndsWith(".root"))
    {
      stem.Resize(stem.Length() - 5);
    }
    TString clean_hist = histname;
    clean_hist.ReplaceAll("/", "_");
    clean_hist.ReplaceAll(";", "_");
    clean_hist.ReplaceAll(" ", "_");
    return stem + "_" + clean_hist;
  }

  void draw_and_save(TH1* hist, TH1* hist_cuts, const TString& title, const TString& png, const TString& pdf, bool logy)
  {
    TCanvas canvas("c_qa_tracking_pt", "QA Tracking pT", 1100, 800);
    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.05);
    canvas.SetBottomMargin(0.12);
    canvas.SetTopMargin(0.08);
    canvas.SetLogy(logy);

    hist->SetTitle(title);
    hist->SetLineColor(kBlue + 1);
    hist->SetMarkerColor(kBlue + 1);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.1);
    hist->SetLineWidth(2);
    hist->GetXaxis()->SetTitle("Reco p_{T} [GeV/c]");
    hist->GetYaxis()->SetTitle("Counts");
    hist->Draw("E1");

    if (hist_cuts)
    {
      hist_cuts->SetLineColor(kRed + 1);
      hist_cuts->SetMarkerColor(kRed + 1);
      hist_cuts->SetMarkerStyle(24);
      hist_cuts->SetMarkerSize(1.0);
      hist_cuts->SetLineWidth(2);
      hist_cuts->Draw("E1 SAME");
    }

    draw_stats(hist);
    canvas.SaveAs(png);
    canvas.SaveAs(pdf);
  }
}

void DrawQaTrackingPt(
    const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pionplus_pt10/completed/general_00004_qa.root",
    const char* histname = "h_QAG4SimulationTracking_nReco_pTReco",
    const bool overlay_cuts = true,
    const double zoom_xmin = 0.0,
    const double zoom_xmax = 20.0,
    const bool logy = false)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  std::unique_ptr<TFile> infile_ptr(TFile::Open(infile, "READ"));
  if (!infile_ptr || infile_ptr->IsZombie())
  {
    ::Error("DrawQaTrackingPt", "Failed to open input file %s", infile);
    return;
  }

  auto* hist_in = dynamic_cast<TH1*>(infile_ptr->Get(histname));
  if (!hist_in)
  {
    ::Error("DrawQaTrackingPt", "Histogram %s not found in %s", histname, infile);
    infile_ptr->ls();
    return;
  }

  std::unique_ptr<TH1> hist(dynamic_cast<TH1*>(hist_in->Clone(Form("%s_clone", hist_in->GetName()))));
  hist->SetDirectory(nullptr);

  std::unique_ptr<TH1> hist_cuts;
  if (overlay_cuts && TString(histname) == "h_QAG4SimulationTracking_nReco_pTReco")
  {
    auto* hist_cuts_in = dynamic_cast<TH1*>(infile_ptr->Get("h_QAG4SimulationTracking_nReco_pTReco_cuts"));
    if (hist_cuts_in)
    {
      hist_cuts.reset(dynamic_cast<TH1*>(hist_cuts_in->Clone(Form("%s_clone", hist_cuts_in->GetName()))));
      hist_cuts->SetDirectory(nullptr);
    }
  }

  const TString stem = build_output_stem(infile, histname);
  draw_and_save(
      hist.get(),
      hist_cuts.get(),
      Form("%s from %s", hist->GetName(), gSystem->BaseName(infile)),
      stem + ".png",
      stem + ".pdf",
      logy);

  if (zoom_xmax > zoom_xmin)
  {
    const double hist_xmin = hist->GetXaxis()->GetXmin();
    const double hist_xmax = hist->GetXaxis()->GetXmax();
    const double draw_xmin = std::max(zoom_xmin, hist_xmin);
    const double draw_xmax = std::min(zoom_xmax, hist_xmax);

    std::unique_ptr<TH1> hist_zoom(dynamic_cast<TH1*>(hist->Clone(Form("%s_zoom", hist->GetName()))));
    hist_zoom->SetDirectory(nullptr);
    hist_zoom->GetXaxis()->SetRangeUser(draw_xmin, draw_xmax);

    std::unique_ptr<TH1> hist_cuts_zoom;
    if (hist_cuts)
    {
      hist_cuts_zoom.reset(dynamic_cast<TH1*>(hist_cuts->Clone(Form("%s_zoom", hist_cuts->GetName()))));
      hist_cuts_zoom->SetDirectory(nullptr);
      hist_cuts_zoom->GetXaxis()->SetRangeUser(draw_xmin, draw_xmax);
    }

    draw_and_save(
        hist_zoom.get(),
        hist_cuts_zoom.get(),
        Form("%s from %s (zoom %.0f-%.0f GeV/c)", hist->GetName(), gSystem->BaseName(infile), draw_xmin, draw_xmax),
        stem + Form("_zoom_%g_%g.png", zoom_xmin, zoom_xmax),
        stem + Form("_zoom_%g_%g.pdf", zoom_xmin, zoom_xmax),
        logy);
  }

  std::cout << "Saved plots to " << stem << ".png/.pdf";
  if (zoom_xmax > zoom_xmin)
  {
    std::cout << " and " << stem << Form("_zoom_%g_%g.png/.pdf", zoom_xmin, zoom_xmax);
  }
  std::cout << std::endl;
}
