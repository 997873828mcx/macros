#ifndef DRAWPHASE1GENERATORQA_C
#define DRAWPHASE1GENERATORQA_C

#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>

void drawPhase1GeneratorQA(
    const char *eventFile = "output/phase1_truth_events_100k_zero_vtx.root",
    const char *v0File = "output/phase1_truth_v0_decays_100k_zero_vtx.root",
    const char *outtag = "Phase1_generatorQA")
{
  TFile *fevent = TFile::Open(eventFile, "READ");
  if (!fevent || fevent->IsZombie())
  {
    std::cerr << "[drawPhase1GeneratorQA] cannot open " << eventFile << std::endl;
    return;
  }

  TFile *fv0 = TFile::Open(v0File, "READ");
  if (!fv0 || fv0->IsZombie())
  {
    std::cerr << "[drawPhase1GeneratorQA] cannot open " << v0File << std::endl;
    return;
  }

  TTree *eventTree = nullptr;
  TTree *v0Tree = nullptr;
  fevent->GetObject("truthEventTree", eventTree);
  fv0->GetObject("truthV0Tree", v0Tree);
  if (!eventTree)
  {
    std::cerr << "[drawPhase1GeneratorQA] truthEventTree not found in " << eventFile << std::endl;
    return;
  }
  if (!v0Tree)
  {
    std::cerr << "[drawPhase1GeneratorQA] truthV0Tree not found in " << v0File << std::endl;
    return;
  }

  const TString pdf = Form("%s.pdf", outtag);
  const TString rootName = Form("%s.root", outtag);
  TFile fout(rootName, "RECREATE");

  auto *c = new TCanvas("c_phase1", "Phase 1 generator QA", 1100, 850);
  c->Print(pdf + "[");

  auto draw1 = [&](TTree *tree, const char *expr, TH1 *hist, const char *cut = "", const char *opt = "hist") {
    fout.cd();
    hist->SetDirectory(&fout);
    tree->Draw(Form("%s>>%s", expr, hist->GetName()), cut, opt);
    hist->Write();
    c->Clear();
    hist->Draw(opt);
    c->Print(pdf);
  };

  auto draw2 = [&](TTree *tree, const char *expr, TH2 *hist, const char *cut = "", const char *opt = "colz") {
    fout.cd();
    hist->SetDirectory(&fout);
    tree->Draw(Form("%s>>%s", expr, hist->GetName()), cut, opt);
    hist->Write();
    c->Clear();
    hist->Draw(opt);
    c->Print(pdf);
  };

  draw1(eventTree, "n_charged_primary_eta_pt",
        new TH1F("h_nch_eta_pt", ";N_{ch}^{primary} (|#eta|<1.1, p_{T}>0.2 GeV/c);Events", 120, -0.5, 119.5));
  draw1(eventTree, "primary_z",
        new TH1F("h_primary_z", ";Primary vertex z (cm);Events", 120, -6.0, 6.0));
  draw1(eventTree, "n_kshort_charged_decay",
        new TH1F("h_n_kshort_event", ";Generated K_{S}^{0}#rightarrow#pi^{+}#pi^{-} per event;Events", 12, -0.5, 11.5));
  draw1(eventTree, "n_lambda_charged_decay",
        new TH1F("h_n_lambda_event", ";Generated #Lambda#rightarrowp#pi^{-} per event;Events", 12, -0.5, 11.5));
  draw1(eventTree, "n_antilambda_charged_decay",
        new TH1F("h_n_antilambda_event", ";Generated #bar{#Lambda}#rightarrow#bar{p}#pi^{+} per event;Events", 12, -0.5, 11.5));

  draw2(eventTree, "n_kshort_charged_decay:n_charged_primary_eta_pt",
        new TH2F("h_kshort_vs_nch", ";N_{ch}^{primary} (|#eta|<1.1, p_{T}>0.2);K_{S}^{0} per event", 80, -0.5, 79.5, 10, -0.5, 9.5));
  draw2(eventTree, "n_lambda_charged_decay:n_charged_primary_eta_pt",
        new TH2F("h_lambda_vs_nch", ";N_{ch}^{primary} (|#eta|<1.1, p_{T}>0.2);#Lambda per event", 80, -0.5, 79.5, 10, -0.5, 9.5));
  draw2(eventTree, "n_antilambda_charged_decay:n_charged_primary_eta_pt",
        new TH2F("h_antilambda_vs_nch", ";N_{ch}^{primary} (|#eta|<1.1, p_{T}>0.2);#bar{#Lambda} per event", 80, -0.5, 79.5, 10, -0.5, 9.5));

  const TCut kshort = "parent_pid == 310";
  const TCut lambda = "parent_pid == 3122";
  const TCut antilambda = "parent_pid == -3122";

  auto overlay3 = [&](const char *name, const char *title, const char *expr, int bins, double lo, double hi, const char *drawOpt = "hist") {
    fout.cd();
    auto *hK = new TH1F(Form("%s_kshort", name), title, bins, lo, hi);
    auto *hL = new TH1F(Form("%s_lambda", name), title, bins, lo, hi);
    auto *hA = new TH1F(Form("%s_antilambda", name), title, bins, lo, hi);
    hK->SetDirectory(&fout);
    hL->SetDirectory(&fout);
    hA->SetDirectory(&fout);
    v0Tree->Draw(Form("%s>>%s", expr, hK->GetName()), kshort, "goff");
    v0Tree->Draw(Form("%s>>%s", expr, hL->GetName()), lambda, "goff");
    v0Tree->Draw(Form("%s>>%s", expr, hA->GetName()), antilambda, "goff");
    hK->SetLineColor(kBlue + 1);
    hL->SetLineColor(kRed + 1);
    hA->SetLineColor(kGreen + 2);
    hK->SetLineWidth(2);
    hL->SetLineWidth(2);
    hA->SetLineWidth(2);
    hK->SetStats(false);
    hL->SetStats(false);
    hA->SetStats(false);
    hK->Write();
    hL->Write();
    hA->Write();

    c->Clear();
    const double ymax = std::max({hK->GetMaximum(), hL->GetMaximum(), hA->GetMaximum()});
    hK->SetMinimum(0.0);
    hK->SetMaximum(ymax > 0.0 ? 1.2 * ymax : 1.0);
    hK->Draw(drawOpt);
    hL->Draw(Form("%s same", drawOpt));
    hA->Draw(Form("%s same", drawOpt));
    auto *leg = new TLegend(0.62, 0.70, 0.88, 0.88);
    leg->AddEntry(hK, "K_{S}^{0}", "l");
    leg->AddEntry(hL, "#Lambda", "l");
    leg->AddEntry(hA, "#bar{#Lambda}", "l");
    leg->Draw();
    c->Print(pdf);
  };

  overlay3("h_v0_pt", ";Generated V0 p_{T} (GeV/c);Counts", "parent_pt", 120, 0.0, 8.0);
  overlay3("h_v0_eta", ";Generated V0 #eta;Counts", "parent_eta", 120, -4.0, 4.0);
  overlay3("h_v0_y", ";Generated V0 y;Counts", "0.5*log((parent_e+parent_pz)/(parent_e-parent_pz))", 120, -4.0, 4.0);
  overlay3("h_v0_lxy", ";Generated decay radius L_{xy} (cm);Counts", "Lxy", 120, 0.0, 120.0);
  overlay3("h_v0_ct", ";Generated proper decay length mL/p (cm);Counts",
           "((parent_pid==310)*0.497611+(abs(parent_pid)==3122)*1.115683)*Lxyz/sqrt(parent_px*parent_px+parent_py*parent_py+parent_pz*parent_pz)", 120, 0.0, 30.0);

  draw2(v0Tree, "qT:alpha",
        new TH2F("h_truth_ap_all", ";#alpha;q_{T}^{+} (GeV/c)", 200, -1.0, 1.0, 200, 0.0, 0.25));
  draw2(v0Tree, "qT:alpha",
        new TH2F("h_truth_ap_kshort", "K_{S}^{0};#alpha;q_{T}^{+} (GeV/c)", 200, -1.0, 1.0, 200, 0.0, 0.25),
        "parent_pid == 310");
  draw2(v0Tree, "qT:alpha",
        new TH2F("h_truth_ap_lambda", "#Lambda + #bar{#Lambda};#alpha;q_{T}^{+} (GeV/c)", 200, -1.0, 1.0, 200, 0.0, 0.25),
        "abs(parent_pid) == 3122");

  c->Print(pdf + "]");
  fout.Close();

  std::cout << "[drawPhase1GeneratorQA] wrote " << pdf << " and " << rootName << std::endl;
}

#endif
