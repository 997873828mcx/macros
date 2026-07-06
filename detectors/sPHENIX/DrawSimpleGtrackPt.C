#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>

#include <iostream>

void DrawSimpleGtrackPt(
    const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pionplus_pt10_1000evt_split_emcoff_hcaloff_1/completed/merged_eval_2.root",
    const char* tree_name = "ntp_gtrack",
    const char* branch_name = "pt",
    const int nbins = 120,
    const double xmin = 0.0,
    const double xmax = 20.0)
{
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "Failed to open input file: " << infile << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(fin->Get(tree_name));
  if (!tree)
  {
    std::cerr << "Could not find tree '" << tree_name << "' in " << infile << std::endl;
    fin->Close();
    return;
  }

  if (!tree->GetBranch(branch_name))
  {
    std::cerr << "Could not find branch '" << branch_name << "' in tree '" << tree_name << "'" << std::endl;
    fin->Close();
    return;
  }

  gROOT->cd();
  if (gROOT->FindObject("h_simple_gtrack_pt"))
  {
    delete gROOT->FindObject("h_simple_gtrack_pt");
  }

  tree->Draw(Form("%s>>h_simple_gtrack_pt(%d,%g,%g)", branch_name, nbins, xmin, xmax), "", "goff");
  TH1* h = dynamic_cast<TH1*>(gROOT->FindObject("h_simple_gtrack_pt"));
  if (!h)
  {
    std::cerr << "Failed to create histogram from branch '" << branch_name << "'" << std::endl;
    fin->Close();
    return;
  }
  h->SetTitle(Form("%s:%s;%s;Counts", tree_name, branch_name, branch_name));

  TCanvas* c = new TCanvas("c_simple_gtrack_pt", "Simple gtrack pt", 900, 700);
  c->SetMargin(0.12, 0.04, 0.12, 0.06);
  h->SetLineWidth(2);
  h->Draw("hist");

  TString outbase = gSystem->BaseName(infile);
  outbase.ReplaceAll(".root", "");
  outbase += Form("_%s_%s_simple", tree_name, branch_name);

  const TString outdir = gSystem->DirName(infile);
  TString outprefix = outbase;
  if (!outdir.IsNull() && outdir != ".")
  {
    outprefix = outdir + "/" + outbase;
  }

  c->SaveAs(outprefix + ".png");
  c->SaveAs(outprefix + ".pdf");

  std::cout << "Filled entries: " << h->GetEntries() << std::endl;
  std::cout << "Histogram mean: " << h->GetMean() << std::endl;
  std::cout << "Histogram stddev: " << h->GetStdDev() << std::endl;
  std::cout << "Saved plots to " << outprefix << ".png/.pdf" << std::endl;

  fin->Close();
}
