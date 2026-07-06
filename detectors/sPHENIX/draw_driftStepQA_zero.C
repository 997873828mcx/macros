// Draw quick diagnostics for entries where ionized_electrons_dedx_kev_per_cm is zero.
//
// Example:
//   root -l -b -q 'draw_driftStepQA_zero.C("output/iondedx_test/completed/general_00000_DriftStepQA.root")'
//   root -l -b -q 'draw_driftStepQA_zero.C("driftStepQA_merged.root","plots_zero")'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <cstdio>
#include <memory>

namespace
{
  void save_canvas(TCanvas& canvas, const TString& outdir, const TString& name)
  {
    canvas.SaveAs(outdir + "/" + name + ".png");
    canvas.SaveAs(outdir + "/" + name + ".pdf");
  }
}

void draw_driftStepQA_zero(
    const char* input_file = "output/iondedx_test/completed/general_00000_DriftStepQA.root",
    const char* output_dir = "driftStepQA_zero_plots")
{
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie())
  {
    std::printf("draw_driftStepQA_zero: could not open %s\n", input_file);
    return;
  }

  auto* tree = dynamic_cast<TTree*>(file->Get("driftStepQA"));
  if (!tree)
  {
    std::printf("draw_driftStepQA_zero: no driftStepQA tree in %s\n", input_file);
    return;
  }

  const TString outdir(output_dir);
  gSystem->mkdir(outdir, true);

  const char* zero_cut = "ionized_electrons_dedx_kev_per_cm==0";
  const char* nonzero_cut = "ionized_electrons_dedx_kev_per_cm>0";

  std::printf("Entries: %lld\n", tree->GetEntries());
  std::printf("Zero ionized dE/dx entries: %lld\n", tree->GetEntries(zero_cut));
  std::printf("Nonzero ionized dE/dx entries: %lld\n", tree->GetEntries(nonzero_cut));

  TCanvas c1("c1", "ionized electron dE/dx", 900, 700);
  c1.SetLogy();
  tree->Draw("ionized_electrons_dedx_kev_per_cm>>hIonDedx(200,0,20)");
  if (auto* h = dynamic_cast<TH1D*>(gDirectory->Get("hIonDedx")))
  {
    h->SetTitle("Ionized-electron dE/dx;ionized-electron dE/dx (keV/cm);steps");
    h->SetLineWidth(2);
  }
  save_canvas(c1, outdir, "ionized_electron_dedx_all");

  TCanvas c2("c2", "zero ionized dE/dx diagnostics", 1400, 900);
  c2.Divide(2, 2);

  c2.cd(1);
  gPad->SetLogy();
  tree->Draw("step_length_cm>>hStepZero(120,0,1.2)", zero_cut);
  if (auto* h = dynamic_cast<TH1D*>(gDirectory->Get("hStepZero")))
  {
    h->SetTitle("Zero ionized dE/dx: step length;step length (cm);steps");
    h->SetLineWidth(2);
  }

  c2.cd(2);
  gPad->SetLogy();
  tree->Draw("eion_kev>>hEionZero(120,0,20)", zero_cut);
  if (auto* h = dynamic_cast<TH1D*>(gDirectory->Get("hEionZero")))
  {
    h->SetTitle("Zero ionized dE/dx: Geant4 ionized energy;eion (keV);steps");
    h->SetLineWidth(2);
  }

  c2.cd(3);
  gPad->SetLogy();
  tree->Draw("n_primary_clusters>>hPrimaryZero(80,-0.5,79.5)", zero_cut);
  if (auto* h = dynamic_cast<TH1D*>(gDirectory->Get("hPrimaryZero")))
  {
    h->SetTitle("Zero ionized dE/dx: primary clusters;N primary clusters;steps");
    h->SetLineWidth(2);
  }

  c2.cd(4);
  gPad->SetLogy();
  tree->Draw("n_ionization_electrons>>hElectronZero(120,-0.5,119.5)", zero_cut);
  if (auto* h = dynamic_cast<TH1D*>(gDirectory->Get("hElectronZero")))
  {
    h->SetTitle("Zero ionized dE/dx: ionization electrons;N ionization electrons;steps");
    h->SetLineWidth(2);
  }

  save_canvas(c2, outdir, "zero_ionized_dedx_diagnostics");

  TCanvas c3("c3", "Geant4 vs ionized-electron dE/dx", 900, 700);
  c3.SetLogz();
  tree->Draw("ionized_electrons_dedx_kev_per_cm:deiondx_kev_per_cm>>hDedxCompare(160,0,20,160,0,20)", "", "colz");
  if (auto* h = dynamic_cast<TH2D*>(gDirectory->Get("hDedxCompare")))
  {
    h->SetTitle("Ionized-electron dE/dx vs Geant4 ionized dE/dx;Geant4 eion/step (keV/cm);N_{e} #times W / step (keV/cm)");
  }
  save_canvas(c3, outdir, "ionized_vs_geant4_ionized_dedx");
}
