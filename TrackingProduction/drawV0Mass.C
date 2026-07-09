#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <TStyle.h>

#include <iostream>
#include <memory>
#include <string>

void drawV0Mass(
    const char* infile = "output/tpc_pattern_v0_79513_keff1_kalman_merged.root",
    const char* outtag = "KshortMass",
    const char* massBranch = "mass_Kshort",
    const char* cutString = "charge1 != charge2",
    const int nbins = 200,
    const double xmin = 0.3,
    const double xmax = 0.8,
    const char* treeName = "pairTree")
{
  std::unique_ptr<TFile> file(TFile::Open(infile, "READ"));
  if (!file || file->IsZombie())
  {
    std::cerr << "[drawV0Mass] cannot open " << infile << std::endl;
    return;
  }

  TTree* tree = nullptr;
  file->GetObject(treeName, tree);
  if (!tree)
  {
    std::cerr << "[drawV0Mass] tree " << treeName << " not found in " << infile << std::endl;
    return;
  }

  gStyle->SetOptStat(1110);

  const TString histName = Form("h_%s", massBranch);
  auto* hist = new TH1F(histName, Form(";%s [GeV/c^{2}];candidates", massBranch),
                        nbins, xmin, xmax);
  hist->SetLineWidth(2);

  TCanvas canvas("c_mass", "V0 invariant mass", 1000, 800);
  const TCut cuts(cutString ? cutString : "");
  tree->Draw(Form("%s>>%s", massBranch, histName.Data()), cuts, "hist");
  canvas.SetLogy();

  const TString pdfName = Form("%s.pdf", outtag);
  const TString pngName = Form("%s.png", outtag);
  const TString rootName = Form("%s.root", outtag);
  canvas.SaveAs(pdfName);
  canvas.SaveAs(pngName);

  TFile output(rootName, "RECREATE");
  hist->Write();
  canvas.Write("c_mass");
  output.Close();

  std::cout << "[drawV0Mass] entries passing cuts: " << tree->GetEntries(cuts) << std::endl;
  std::cout << "[drawV0Mass] wrote " << pdfName << std::endl;
  std::cout << "[drawV0Mass] wrote " << pngName << std::endl;
  std::cout << "[drawV0Mass] wrote " << rootName << std::endl;
}
