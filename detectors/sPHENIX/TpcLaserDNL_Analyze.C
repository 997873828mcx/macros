// Simple analysis macro for TpcLaserDNL output
// Usage:
//   root -l -b -q 'develope/macros/analysis/TpcLaserDNL_Analyze.C("laser_dnl.root","laser_dnl_summary.root", -1, -1, 72)'

#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TColor.h>
#include <TString.h>
#include <vector>
#include <memory>
#include <cstdio>

void TpcLaserDNL_Analyze(const char* infile = "laser_sim_dnl_total.root",
                         const char* outfile = "laser_dnl_summary_total_more_bins.root",
                         int layer_filter = -1,   // -1 for all
                         int side_filter = -1,    // -1 for both; 0 (south), 1 (north)
                         int nphi_bins = 100,      // requested ~30 phi bins
                         double phi_min = 1.54,   // requested phi range
                         double phi_max = 1.58)   // requested phi range
{
  // open input
  std::unique_ptr<TFile> fin(TFile::Open(infile, "READ"));
  if (!fin || fin->IsZombie()) { printf("Cannot open %s\n", infile); return; }
  TTree* tt = dynamic_cast<TTree*>(fin->Get("dnl"));
  if (!tt) { printf("Tree 'dnl' not found in %s\n", infile); return; }

  // branches from TpcLaserDNL
  int event=0, trkid=0, side=0, nused=0;
  unsigned int layer=0; // match tree type (UInt_t)
  double r=0, phi_true=0, phi_reco=0, dphi=0, dRphi=0, adcsum=0;
  double xtrue=0, ytrue=0, ztrue=0, xreco=0, yreco=0, zreco=0;

  tt->SetBranchAddress("event",   &event);
  tt->SetBranchAddress("trkid",   &trkid);
  tt->SetBranchAddress("layer",   &layer);
  tt->SetBranchAddress("side",    &side);
  tt->SetBranchAddress("r",       &r);
  tt->SetBranchAddress("phi_true",&phi_true);
  tt->SetBranchAddress("phi_reco",&phi_reco);
  tt->SetBranchAddress("dphi",    &dphi);
  tt->SetBranchAddress("dRphi",   &dRphi);
  tt->SetBranchAddress("nused",   &nused);
  tt->SetBranchAddress("adcsum",  &adcsum);
  tt->SetBranchAddress("xtrue",   &xtrue);
  tt->SetBranchAddress("ytrue",   &ytrue);
  tt->SetBranchAddress("ztrue",   &ztrue);
  tt->SetBranchAddress("xreco",   &xreco);
  tt->SetBranchAddress("yreco",   &yreco);
  tt->SetBranchAddress("zreco",   &zreco);

  // output
  std::unique_ptr<TFile> fout(TFile::Open(outfile, "RECREATE"));
  if (!fout || fout->IsZombie()) { printf("Cannot create %s\n", outfile); return; }

  // global 2D histogram: dRphi vs phi_true (counts heatmap, layer-agnostic)
  // keep as a wide-range overview
  TH2D* h2_dRphi_vs_phi = new TH2D("h2_dRphi_vs_phi",
                                   "dRphi vs phi_{true} (all layers);phi_{true} [rad];dRphi [cm]",
                                   180, -M_PI, M_PI,
                                   400, -1.0, 1.0);

  // global 2D histogram: dRphi vs phi_reco for comparison against reconstruction
  TH2D* h2_dRphi_vs_phiReco = new TH2D("h2_dRphi_vs_phiReco",
                                       "dRphi vs phi_{reco} (all layers);phi_{reco} [rad];dRphi [cm]",
                                       180, -M_PI, M_PI,
                                       400, -1.0, 1.0);

  // per-layer directories
  const int maxLayers = 64; // safety
  std::vector<TProfile*> prof_phi_all(maxLayers, nullptr);
  std::vector<TProfile*> prof_phi_side0(maxLayers, nullptr);
  std::vector<TProfile*> prof_phi_side1(maxLayers, nullptr);
  std::vector<TH1D*> h_dRphi_all(maxLayers, nullptr);
  // requested: per-layer 2D histograms of dRphi vs phi_true, in a narrow phi window
  std::vector<TH2D*> h2_dRphi_vs_phi_layer(maxLayers, nullptr);
  std::vector<TH2D*> h2_dRphi_vs_phiReco_layer(maxLayers, nullptr);
  std::vector<TGraph*> gr_phiReco_vs_phiTrue(maxLayers, nullptr);

  auto mkProfName = [](const char* base, int lyr, int side)->TString{
    if (side < 0) return TString::Format("%s_layer%02d", base, lyr);
    else return TString::Format("%s_layer%02d_side%d", base, lyr, side);
  };

  // create lazily
  auto getProf = [&](int lyr, int sideSel) -> TProfile* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    TProfile** slot = nullptr;
    if (sideSel < 0) slot = &prof_phi_all[lyr];
    else if (sideSel == 0) slot = &prof_phi_side0[lyr];
    else slot = &prof_phi_side1[lyr];
    if (!*slot)
    {
      const char* base = "prof_dRphi_vs_phi";
      auto name = mkProfName(base, lyr, sideSel);
      *slot = new TProfile(name, name, nphi_bins, -M_PI, M_PI);
      (*slot)->GetXaxis()->SetTitle("phi_true [rad]");
      (*slot)->GetYaxis()->SetTitle("<dRphi> [cm]");
    }
    return *slot;
  };

  auto getHist = [&](int lyr) -> TH1D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h_dRphi_all[lyr])
    {
      auto name = TString::Format("h_dRphi_layer%02d", lyr);
      h_dRphi_all[lyr] = new TH1D(name, name, 200, -1.0, 1.0);
      h_dRphi_all[lyr]->GetXaxis()->SetTitle("dRphi [cm]");
    }
    return h_dRphi_all[lyr];
  };

  auto getH2Layer = [&](int lyr) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h2_dRphi_vs_phi_layer[lyr])
    {
      auto name = TString::Format("h2_dRphi_vs_phi_layer%02d", lyr);
      auto title = TString::Format("dRphi vs phi_{true};phi_{true} [rad];dRphi [cm]");
      // X: phi_true in requested range [phi_min, phi_max] with nphi_bins
      // Y: dRphi with a symmetric range around 0
      h2_dRphi_vs_phi_layer[lyr] = new TH2D(name, title, nphi_bins, phi_min, phi_max, 400, -0.03, 0.03);
    }
    return h2_dRphi_vs_phi_layer[lyr];
  };

  auto getH2PhiRecoLayer = [&](int lyr) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h2_dRphi_vs_phiReco_layer[lyr])
    {
      auto name = TString::Format("h2_dRphi_vs_phiReco_layer%02d", lyr);
      auto title = TString::Format("dRphi vs phi_{reco};phi_{reco} [rad];dRphi [cm]");
      h2_dRphi_vs_phiReco_layer[lyr] = new TH2D(name, title, nphi_bins, phi_min, phi_max, 400, -0.03, 0.03);
    }
    return h2_dRphi_vs_phiReco_layer[lyr];
  };

  auto getGraph = [&](int lyr) -> TGraph* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!gr_phiReco_vs_phiTrue[lyr])
    {
      auto name = TString::Format("gr_phiReco_vs_phiTrue_layer%02d", lyr);
      auto title = TString::Format("phi_{reco} vs phi_{true} (layer %d);phi_{true} [rad];phi_{reco} [rad]", lyr);
      auto* gr = new TGraph();
      gr->SetName(name);
      gr->SetTitle(title);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.4);
      gr->SetMarkerColor(kBlue+1);
      gr_phiReco_vs_phiTrue[lyr] = gr;
    }
    return gr_phiReco_vs_phiTrue[lyr];
  };

  const Long64_t nentries = tt->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    tt->GetEntry(i);
    if (layer_filter >= 0 && layer != layer_filter) continue;
    if (side_filter >= 0 && side  != side_filter) continue;

    // profiles vs phi_true
    if (auto p = getProf(layer, -1)) p->Fill(phi_true, dRphi);
    if (side == 0) { if (auto p = getProf(layer,0)) p->Fill(phi_true, dRphi); }
    else           { if (auto p = getProf(layer,1)) p->Fill(phi_true, dRphi); }

    // overall distribution
    if (auto h = getHist(layer)) h->Fill(dRphi);

    // per-layer 2D in the specified phi window
    if (auto h2 = getH2Layer(layer)) h2->Fill(phi_true, dRphi);
    if (auto h2 = getH2PhiRecoLayer(layer)) h2->Fill(phi_reco, dRphi);

    // per-layer graph of phi_reco vs phi_true
    if (auto gr = getGraph(layer)) gr->SetPoint(gr->GetN(), phi_true, phi_reco);

    // global 2D counts map (ignores layer dimension)
    h2_dRphi_vs_phi->Fill(phi_true, dRphi);
    h2_dRphi_vs_phiReco->Fill(phi_reco, dRphi);
  }

  // write out
  fout->cd();
  if (h2_dRphi_vs_phi) h2_dRphi_vs_phi->Write(); // top-level
  if (h2_dRphi_vs_phiReco) h2_dRphi_vs_phiReco->Write();
  TDirectory* d_all = fout->mkdir("by_layer"); d_all->cd();
  for (auto* h : h_dRphi_all) if (h) h->Write();
  for (auto* p : prof_phi_all) if (p) p->Write();
  TDirectory* d_s0 = fout->mkdir("by_layer_side0"); d_s0->cd();
  for (auto* p : prof_phi_side0) if (p) p->Write();
  TDirectory* d_s1 = fout->mkdir("by_layer_side1"); d_s1->cd();
  for (auto* p : prof_phi_side1) if (p) p->Write();

  // write per-layer 2D histograms in their own directory
  TDirectory* d_2d = fout->mkdir("by_layer_2d"); d_2d->cd();
  for (auto* h2 : h2_dRphi_vs_phi_layer) if (h2) h2->Write();
  for (auto* h2 : h2_dRphi_vs_phiReco_layer) if (h2) h2->Write();

  // write per-layer phi_reco vs phi_true graphs
  fout->cd();
  TDirectory* d_gr = fout->mkdir("by_layer_phi_graph"); d_gr->cd();
  for (auto* gr : gr_phiReco_vs_phiTrue) if (gr) gr->Write();

  // also store canvases with log Z so visualization defaults are preserved
  fout->cd();
  TDirectory* d_cv = fout->mkdir("canvases"); d_cv->cd();

  {
    auto c = new TCanvas("c_h2_dRphi_vs_phi", "c_h2_dRphi_vs_phi", 900, 700);
    c->SetLogz(1);
    h2_dRphi_vs_phi->Draw("COLZ");
    c->Write();
    c->SaveAs(TString::Format("%s.pdf", c->GetName()));
  }
  {
    auto c = new TCanvas("c_h2_dRphi_vs_phiReco", "c_h2_dRphi_vs_phiReco", 900, 700);
    c->SetLogz(1);
    h2_dRphi_vs_phiReco->Draw("COLZ");
    c->Write();
    c->SaveAs(TString::Format("%s.pdf", c->GetName()));
  }
  for (int lyr = 0; lyr < (int)h2_dRphi_vs_phi_layer.size(); ++lyr)
  {
    auto* h2 = h2_dRphi_vs_phi_layer[lyr];
    if (!h2) continue;
    auto cname = TString::Format("c_h2_dRphi_vs_phi_layer%02d", lyr);
    auto c = new TCanvas(cname, cname, 900, 700);
    c->SetLogz(1);
    h2->Draw("COLZ");
    c->Write();
    c->SaveAs(TString::Format("%s.pdf", c->GetName()));
  }
  for (int lyr = 0; lyr < (int)h2_dRphi_vs_phiReco_layer.size(); ++lyr)
  {
    auto* h2 = h2_dRphi_vs_phiReco_layer[lyr];
    if (!h2) continue;
    auto cname = TString::Format("c_h2_dRphi_vs_phiReco_layer%02d", lyr);
    auto c = new TCanvas(cname, cname, 900, 700);
    c->SetLogz(1);
    h2->Draw("COLZ");
    c->Write();
    c->SaveAs(TString::Format("%s.pdf", c->GetName()));
  }
  for (int lyr = 0; lyr < (int)gr_phiReco_vs_phiTrue.size(); ++lyr)
  {
    auto* gr = gr_phiReco_vs_phiTrue[lyr];
    if (!gr) continue;
    auto cname = TString::Format("c_gr_phiReco_vs_phiTrue_layer%02d", lyr);
    auto c = new TCanvas(cname, cname, 800, 600);
    c->SetGrid();
    gr->Draw("AP");
    c->Write();
    c->SaveAs(TString::Format("%s.pdf", c->GetName()));
  }

  fout->Write();
  //fout->Close();
}
