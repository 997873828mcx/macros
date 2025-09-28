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
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <TLine.h>

void TpcLaserDNL_Analyze(const char* infile = "laser_sim_dnl_total.root",
                         const char* outfile = "laser_dnl_summary_total_more_bins.root",
                         int layer_filter = -1,   // -1 for all
                         int side_filter = -1,    // -1 for both; 0 (south), 1 (north)
                         int nphi_bins = 1000,      // requested ~30 phi bins
                         double phi_min = 1.548,   // requested phi range
                         double phi_max = 1.566)   // requested phi range
{
  // open input
  std::unique_ptr<TFile> fin(TFile::Open(infile, "READ"));
  if (!fin || fin->IsZombie()) { printf("Cannot open %s\n", infile); return; }
  TTree* tt = dynamic_cast<TTree*>(fin->Get("dnl"));
  if (!tt) { printf("Tree 'dnl' not found in %s\n", infile); return; }

  // branches from TpcLaserDNL
  int event=0, trkid=0, side=0, nused=0;
  int npad_used = -1;
  unsigned int layer=0; // match tree type (UInt_t)
  double r=0, phi_true=0, phi_reco=0, dphi=0, dRphi=0, adcsum=0;
  double xtrue=0, ytrue=0, ztrue=0, xreco=0, yreco=0, zreco=0;
  double phase = std::numeric_limits<double>::quiet_NaN();
  double phase_reco = std::numeric_limits<double>::quiet_NaN();
  std::vector<double>* pad_phi_centers = nullptr;

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

  const bool has_phase_branch = (tt->GetBranch("phase") != nullptr);
  const bool has_phase_reco_branch = (tt->GetBranch("phase_reco") != nullptr);
  const bool has_npad_branch = (tt->GetBranch("npad_used") != nullptr);
  const bool has_pad_phi_center_branch = (tt->GetBranch("pad_phi_center") != nullptr);
  if (has_phase_branch) tt->SetBranchAddress("phase", &phase);
  if (has_phase_reco_branch) tt->SetBranchAddress("phase_reco", &phase_reco);
  if (has_npad_branch) tt->SetBranchAddress("npad_used", &npad_used);
  if (has_pad_phi_center_branch) tt->SetBranchAddress("pad_phi_center", &pad_phi_centers);

  // output
  std::unique_ptr<TFile> fout(TFile::Open(outfile, "RECREATE"));
  if (!fout || fout->IsZombie()) { printf("Cannot create %s\n", outfile); return; }

  // per-layer directories
  const int maxLayers = 64; // safety
  /* std::vector<TProfile*> prof_phi_all(maxLayers, nullptr);
  std::vector<TProfile*> prof_phi_side0(maxLayers, nullptr);
  std::vector<TProfile*> prof_phi_side1(maxLayers, nullptr); */
  std::vector<TH1D*> h_dRphi_all(maxLayers, nullptr);
  // requested: per-layer 2D histograms of dRphi vs phi_true, in a narrow phi window
  std::vector<TH2D*> h2_dRphi_vs_phi_layer(maxLayers, nullptr);
  std::vector<TH2D*> h2_dRphi_vs_phiReco_layer(maxLayers, nullptr);
  std::vector<TH2D*> h2_dRphi_vs_phase_layer(maxLayers, nullptr);
  std::vector<TH2D*> h2_dRphi_vs_phaseReco_layer(maxLayers, nullptr);
  std::vector<TGraph*> gr_phiReco_vs_phiTrue(maxLayers, nullptr);
  std::vector<TGraph*> gr_phaseReco_vs_phase(maxLayers, nullptr);
  std::vector<TH1D*> h_phase_layer(maxLayers, nullptr);
  std::vector<TH1D*> h_phaseReco_layer(maxLayers, nullptr);
  std::vector<TH2D*> h2_npad_vs_phase_layer(maxLayers, nullptr);
  std::vector<TH2D*> h2_npad_vs_phaseReco_layer(maxLayers, nullptr);

  const int nphase_bins = 120;
  const double phase_min = -0.6;
  const double phase_max = 0.6;
  const double dRphi_min_phase = -0.03;
  const double dRphi_max_phase = 0.03;

  std::vector<std::map<int, TH2D*>> h2_dRphi_vs_phi_layer_byPads(maxLayers);
  std::vector<std::map<int, TH2D*>> h2_dRphi_vs_phiReco_layer_byPads(maxLayers);
  std::vector<std::map<int, TH2D*>> h2_dRphi_vs_phase_layer_byPads(maxLayers);
  std::vector<std::map<int, TH2D*>> h2_dRphi_vs_phaseReco_layer_byPads(maxLayers);
  std::vector<std::map<int, TGraph*>> gr_phiReco_vs_phiTrue_byPads(maxLayers);
  std::vector<std::map<int, TGraph*>> gr_phaseReco_vs_phase_byPads(maxLayers);
  std::vector<std::map<int, TH1D*>> h_phase_layer_byPads(maxLayers);
  std::vector<std::map<int, TH1D*>> h_phaseReco_layer_byPads(maxLayers);
  std::vector<std::map<int, TH1D*>> h_dRphi_layer_byPads(maxLayers);
  std::vector<std::set<double>> pad_phi_centers_by_layer(maxLayers);
  std::vector<double> layer_phiReco_min(maxLayers, std::numeric_limits<double>::infinity());
  std::vector<double> layer_phiReco_max(maxLayers, -std::numeric_limits<double>::infinity());

  std::vector<std::vector<TLine*>> lines_phiTrue_overlays(maxLayers);
  std::vector<std::vector<TLine*>> lines_phiGraph_overlays(maxLayers);

  /* auto mkProfName = [](const char* base, int lyr, int side)->TString{
    if (side < 0) return TString::Format("%s_layer%02d", base, lyr);
    else return TString::Format("%s_layer%02d_side%d", base, lyr, side);
  }; */

  // create lazily
  /* auto getProf = [&](int lyr, int sideSel) -> TProfile* {
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
  }; */

  auto getHist = [&](int lyr) -> TH1D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h_dRphi_all[lyr])
    {
      auto name = TString::Format("h_dRphi_layer%02d", lyr);
      h_dRphi_all[lyr] = new TH1D(name, name, 200, -0.03, 0.03);
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

  auto getPhaseGraph = [&](int lyr) -> TGraph* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!gr_phaseReco_vs_phase[lyr])
    {
      auto name = TString::Format("gr_phaseReco_vs_phase_layer%02d", lyr);
      auto title = TString::Format("phase_{reco} vs phase (layer %d);phase;phase_{reco}", lyr);
      auto* gr = new TGraph();
      gr->SetName(name);
      gr->SetTitle(title);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.4);
      gr->SetMarkerColor(kRed+1);
      gr_phaseReco_vs_phase[lyr] = gr;
    }
    return gr_phaseReco_vs_phase[lyr];
  };

  auto getH2PhaseLayer = [&](int lyr) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h2_dRphi_vs_phase_layer[lyr])
    {
      auto name = TString::Format("h2_dRphi_vs_phase_layer%02d", lyr);
      auto title = TString::Format("dRphi vs phase;phase;dRphi [cm]");
      h2_dRphi_vs_phase_layer[lyr] = new TH2D(name, title, nphase_bins, phase_min, phase_max, 400, dRphi_min_phase, dRphi_max_phase);
    }
    return h2_dRphi_vs_phase_layer[lyr];
  };

  auto getH2PhaseRecoLayer = [&](int lyr) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h2_dRphi_vs_phaseReco_layer[lyr])
    {
      auto name = TString::Format("h2_dRphi_vs_phaseReco_layer%02d", lyr);
      auto title = TString::Format("dRphi vs phase_{reco};phase_{reco};dRphi [cm]");
      h2_dRphi_vs_phaseReco_layer[lyr] = new TH2D(name, title, nphase_bins, phase_min, phase_max, 400, dRphi_min_phase, dRphi_max_phase);
    }
    return h2_dRphi_vs_phaseReco_layer[lyr];
  };

  auto getPhaseHist1D = [&](int lyr) -> TH1D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h_phase_layer[lyr])
    {
      auto name = TString::Format("h_phase_layer%02d", lyr);
      auto title = TString::Format("Layer %d: phase;phase;entries", lyr);
      h_phase_layer[lyr] = new TH1D(name, title, nphase_bins, -1.0, 1.0);
    }
    return h_phase_layer[lyr];
  };

  auto getPhaseRecoHist1D = [&](int lyr) -> TH1D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h_phaseReco_layer[lyr])
    {
      auto name = TString::Format("h_phaseReco_layer%02d", lyr);
      auto title = TString::Format("Layer %d: phase_{reco};phase_{reco};entries", lyr);
      h_phaseReco_layer[lyr] = new TH1D(name, title, nphase_bins, -1.0, 1.0);
    }
    return h_phaseReco_layer[lyr];
  };

  auto getNpadVsPhaseLayer = [&](int lyr) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h2_npad_vs_phase_layer[lyr])
    {
      auto name = TString::Format("h2_npad_vs_phase_layer%02d", lyr);
      auto title = TString::Format("Layer %d: npad vs phase;npad_used;phase", lyr);
      h2_npad_vs_phase_layer[lyr] = new TH2D(name, title, 80, 1, 7, nphase_bins, phase_min, phase_max);
    }
    return h2_npad_vs_phase_layer[lyr];
  };

  auto getNpadVsPhaseRecoLayer = [&](int lyr) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0) return nullptr;
    if (!h2_npad_vs_phaseReco_layer[lyr])
    {
      auto name = TString::Format("h2_npad_vs_phaseReco_layer%02d", lyr);
      auto title = TString::Format("Layer %d: npad vs phase_{reco};npad_used;phase_{reco}", lyr);
      h2_npad_vs_phaseReco_layer[lyr] = new TH2D(name, title, 80, 1, 7, nphase_bins, phase_min, phase_max);
    }
    return h2_npad_vs_phaseReco_layer[lyr];
  };

  auto getPadHist = [&](std::vector<std::map<int, TH2D*>>& store,
                        const char* nameFmt,
                        const char* titleFmt,
                        int lyr,
                        int nPads,
                        int nbinsX,
                        double xmin,
                        double xmax,
                        int nbinsY,
                        double ymin,
                        double ymax) -> TH2D* {
    if (lyr >= maxLayers || lyr < 0 || nPads < 0) return nullptr;
    auto& layerMap = store[lyr];
    auto& hist = layerMap[nPads];
    if (!hist)
    {
      auto name = TString::Format(nameFmt, lyr, nPads);
      auto title = TString::Format(titleFmt, lyr, nPads);
      hist = new TH2D(name, title, nbinsX, xmin, xmax, nbinsY, ymin, ymax);
    }
    return hist;
  };

  auto getPadGraph = [&](std::vector<std::map<int, TGraph*>>& store,
                         const char* nameFmt,
                         const char* titleFmt,
                         int lyr,
                         int nPads,
                         int color) -> TGraph* {
    if (lyr >= maxLayers || lyr < 0 || nPads < 0) return nullptr;
    auto& layerMap = store[lyr];
    auto& gr = layerMap[nPads];
    if (!gr)
    {
      auto name = TString::Format(nameFmt, lyr, nPads);
      auto title = TString::Format(titleFmt, lyr, nPads);
      auto* g = new TGraph();
      g->SetName(name);
      g->SetTitle(title);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.4);
      g->SetMarkerColor(color);
      gr = g;
    }
    return gr;
  };

  auto getPadHist1D = [&](std::vector<std::map<int, TH1D*>>& store,
                          const char* nameFmt,
                          const char* titleFmt,
                          int lyr,
                          int nPads,
                          int nbins,
                          double xmin,
                          double xmax) -> TH1D* {
    if (lyr >= maxLayers || lyr < 0 || nPads < 0) return nullptr;
    auto& layerMap = store[lyr];
    auto& hist = layerMap[nPads];
    if (!hist)
    {
      auto name = TString::Format(nameFmt, lyr, nPads);
      auto title = TString::Format(titleFmt, lyr, nPads);
      hist = new TH1D(name, title, nbins, xmin, xmax);
    }
    return hist;
  };

  const Long64_t nentries = tt->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    npad_used = -1;
    tt->GetEntry(i);
    if (layer_filter >= 0 && layer != layer_filter) continue;
    if (side_filter >= 0 && side  != side_filter) continue;

    // profiles vs phi_true
    /* if (auto p = getProf(layer, -1)) p->Fill(phi_true, dRphi);
    if (side == 0) { if (auto p = getProf(layer,0)) p->Fill(phi_true, dRphi); }
    else           { if (auto p = getProf(layer,1)) p->Fill(phi_true, dRphi); } */

    const int layerIndex = static_cast<int>(layer);

    if (std::isfinite(phi_reco))
    {
      layer_phiReco_min[layerIndex] = std::min(layer_phiReco_min[layerIndex], phi_reco);
      layer_phiReco_max[layerIndex] = std::max(layer_phiReco_max[layerIndex], phi_reco);
    }

    if (has_pad_phi_center_branch && pad_phi_centers)
    {
      auto& phiSet = pad_phi_centers_by_layer[layerIndex];
      for (const double phi_center_value : *pad_phi_centers)
      {
        if (!std::isfinite(phi_center_value)) continue;
        const double phi_key = std::round(phi_center_value * 1e6) / 1e6;
        phiSet.insert(phi_key);
      }
    }

    // overall distribution
    if (auto h = getHist(layer)) h->Fill(dRphi);

    // per-layer 2D in the specified phi window
    if (auto h2 = getH2Layer(layer)) h2->Fill(phi_true, dRphi);
    if (auto h2 = getH2PhiRecoLayer(layer)) h2->Fill(phi_reco, dRphi);

    if (has_phase_branch && std::isfinite(phase))
    {
      if (auto h2 = getH2PhaseLayer(layer)) h2->Fill(phase, dRphi);
      if (auto h1 = getPhaseHist1D(layerIndex)) h1->Fill(phase);
    }
    if (has_phase_reco_branch && std::isfinite(phase_reco))
    {
      if (auto h2 = getH2PhaseRecoLayer(layer)) h2->Fill(phase_reco, dRphi);
      if (auto h1 = getPhaseRecoHist1D(layerIndex)) h1->Fill(phase_reco);
    }

    // per-layer graph of phi_reco vs phi_true
    if (auto gr = getGraph(layer)) gr->SetPoint(gr->GetN(), phi_true, phi_reco);

    if (has_phase_branch && has_phase_reco_branch && std::isfinite(phase) && std::isfinite(phase_reco))
    {
      if (auto gr = getPhaseGraph(layer)) gr->SetPoint(gr->GetN(), phase, phase_reco);
    }

    if (has_npad_branch && npad_used >= 0 && std::isfinite(dRphi))
    {
      const int pad_count = npad_used;

      if (has_phase_branch && std::isfinite(phase))
      {
        if (auto* h2 = getNpadVsPhaseLayer(layerIndex)) h2->Fill(pad_count, phase);
      }
      if (has_phase_reco_branch && std::isfinite(phase_reco))
      {
        if (auto* h2 = getNpadVsPhaseRecoLayer(layerIndex)) h2->Fill(pad_count, phase_reco);
      }

      if (auto* h2 = getPadHist(h2_dRphi_vs_phi_layer_byPads,
                                "h2_dRphi_vs_phi_layer%02d_pads%02d",
                                "Layer %d: dRphi vs phi_{true} (pads=%d);phi_{true} [rad];dRphi [cm]",
                                layerIndex,
                                pad_count,
                                nphi_bins,
                                phi_min,
                                phi_max,
                                400,
                                -0.03,
                                0.03))
      {
        h2->Fill(phi_true, dRphi);
      }

      if (auto* h2 = getPadHist(h2_dRphi_vs_phiReco_layer_byPads,
                                "h2_dRphi_vs_phiReco_layer%02d_pads%02d",
                                "Layer %d: dRphi vs phi_{reco} (pads=%d);phi_{reco} [rad];dRphi [cm]",
                                layerIndex,
                                pad_count,
                                nphi_bins,
                                phi_min,
                                phi_max,
                                400,
                                -0.03,
                                0.03))
      {
        h2->Fill(phi_reco, dRphi);
      }

      if (auto* h1 = getPadHist1D(h_dRphi_layer_byPads,
                                   "h_dRphi_layer%02d_pads%02d",
                                   "Layer %d (pads=%d): dRphi;dRphi [cm];entries",
                                   layerIndex,
                                   pad_count,
                                   200,
                                   -1.0,
                                   1.0))
      {
        h1->Fill(dRphi);
      }

      if (has_phase_branch && std::isfinite(phase))
      {
        if (auto* h2 = getPadHist(h2_dRphi_vs_phase_layer_byPads,
                                  "h2_dRphi_vs_phase_layer%02d_pads%02d",
                                  "Layer %d: dRphi vs phase (pads=%d);phase;dRphi [cm]",
                                  layerIndex,
                                  pad_count,
                                  nphase_bins,
                                  phase_min,
                                  phase_max,
                                  400,
                                  dRphi_min_phase,
                                  dRphi_max_phase))
        {
          h2->Fill(phase, dRphi);
        }
        if (auto* h1 = getPadHist1D(h_phase_layer_byPads,
                                     "h_phase_layer%02d_pads%02d",
                                     "Layer %d (pads=%d): phase;phase;entries",
                                     layerIndex,
                                     pad_count,
                                     nphase_bins,
                                     phase_min,
                                     phase_max))
        {
          h1->Fill(phase);
        }
      }

      if (has_phase_reco_branch && std::isfinite(phase_reco))
      {
        if (auto* h2 = getPadHist(h2_dRphi_vs_phaseReco_layer_byPads,
                                  "h2_dRphi_vs_phaseReco_layer%02d_pads%02d",
                                  "Layer %d: dRphi vs phase_{reco} (pads=%d);phase_{reco};dRphi [cm]",
                                  layerIndex,
                                  pad_count,
                                  nphase_bins,
                                  phase_min,
                                  phase_max,
                                  400,
                                  dRphi_min_phase,
                                  dRphi_max_phase))
        {
          h2->Fill(phase_reco, dRphi);
        }
        if (auto* h1 = getPadHist1D(h_phaseReco_layer_byPads,
                                     "h_phaseReco_layer%02d_pads%02d",
                                     "Layer %d (pads=%d): phase_{reco};phase_{reco};entries",
                                     layerIndex,
                                     pad_count,
                                     nphase_bins,
                                     phase_min,
                                     phase_max))
        {
          h1->Fill(phase_reco);
        }
      }

      if (auto* gr = getPadGraph(gr_phiReco_vs_phiTrue_byPads,
                                  "gr_phiReco_vs_phiTrue_layer%02d_pads%02d",
                                  "Layer %d (pads=%d): phi_{reco} vs phi_{true};phi_{true} [rad];phi_{reco} [rad]",
                                  layerIndex,
                                  pad_count,
                                  kBlue+2))
      {
        gr->SetPoint(gr->GetN(), phi_true, phi_reco);
      }

      if (has_phase_branch && has_phase_reco_branch && std::isfinite(phase) && std::isfinite(phase_reco))
      {
        if (auto* gr = getPadGraph(gr_phaseReco_vs_phase_byPads,
                                    "gr_phaseReco_vs_phase_layer%02d_pads%02d",
                                    "Layer %d (pads=%d): phase_{reco} vs phase;phase;phase_{reco}",
                                    layerIndex,
                                    pad_count,
                                    kRed+2))
        {
          gr->SetPoint(gr->GetN(), phase, phase_reco);
        }
      }
    }
  }

  // build phi-center overlay helpers
  for (int lyr = 0; lyr < maxLayers; ++lyr)
  {
    const auto& phiCenters = pad_phi_centers_by_layer[lyr];
    if (phiCenters.empty()) continue;

    if (auto* h2 = h2_dRphi_vs_phi_layer[lyr])
    {
      const double y1 = h2->GetYaxis()->GetXmin();
      const double y2 = h2->GetYaxis()->GetXmax();
      int idx = 0;
      for (const double phi_center : phiCenters)
      {
        auto* line = new TLine(phi_center, y1, phi_center, y2);
        line->SetLineStyle(2);
        line->SetLineColor(kGray+2);
        line->SetLineWidth(1);
        lines_phiTrue_overlays[lyr].push_back(line);
      }
    }

    if (!lines_phiTrue_overlays[lyr].empty())
    {
      if (gr_phiReco_vs_phiTrue[lyr])
      {
        const double yMinCandidate = layer_phiReco_min[lyr];
        const double yMaxCandidate = layer_phiReco_max[lyr];
        double y1 = std::isfinite(yMinCandidate) ? yMinCandidate : -M_PI;
        double y2 = std::isfinite(yMaxCandidate) ? yMaxCandidate : M_PI;
        if (!(y2 > y1))
        {
          y1 = std::min(y1, -M_PI);
          y2 = std::max(y2, M_PI);
          if (y2 <= y1) { y1 = -M_PI; y2 = M_PI; }
        }

        int idx = 0;
        for (const double phi_center : phiCenters)
        {
          auto* line = new TLine(phi_center, y1, phi_center, y2);
          line->SetLineStyle(2);
          line->SetLineColor(kGray+3);
          line->SetLineWidth(1);
          lines_phiGraph_overlays[lyr].push_back(line);
        }
      }
    }
  }

  // write out
  auto writeTH1Vector = [&](TDirectory* dir, const std::vector<TH1D*>& vec)
  {
    if (!dir) return;
    dir->cd();
    for (auto* h : vec) if (h) h->Write();
  };

  auto writeTH2Vector = [&](TDirectory* dir, const std::vector<TH2D*>& vec)
  {
    if (!dir) return;
    dir->cd();
    for (auto* h : vec) if (h) h->Write();
  };

  auto writeTH2Maps = [&](TDirectory* dir, const std::vector<std::map<int, TH2D*>>& store)
  {
    if (!dir) return;
    dir->cd();
    for (const auto& m : store)
    {
      for (const auto& kv : m) if (kv.second) kv.second->Write();
    }
  };

  auto writeTH1Maps = [&](TDirectory* dir, const std::vector<std::map<int, TH1D*>>& store)
  {
    if (!dir) return;
    dir->cd();
    for (const auto& m : store)
    {
      for (const auto& kv : m) if (kv.second) kv.second->Write();
    }
  };

  auto writeLineCollections = [&](TDirectory* dir, const std::vector<std::vector<TLine*>>& store)
  {
    if (!dir) return;
    for (size_t lyr = 0; lyr < store.size(); ++lyr)
    {
      const auto& lines = store[lyr];
      if (lines.empty()) continue;
      auto* subdir = dir->mkdir(TString::Format("layer%02zu", lyr));
      if (!subdir) continue;
      subdir->cd();
      for (auto* line : lines) if (line) line->Write();
      dir->cd();
    }
  };

  auto drawLines = [](const std::vector<TLine*>& lines)
  {
    for (auto* line : lines) if (line) line->Draw("same");
  };

  auto writeGraphVector = [&](TDirectory* dir, const std::vector<TGraph*>& vec)
  {
    if (!dir) return;
    dir->cd();
    for (auto* g : vec) if (g) g->Write();
  };

  auto writeGraphMaps = [&](TDirectory* dir, const std::vector<std::map<int, TGraph*>>& store)
  {
    if (!dir) return;
    dir->cd();
    for (const auto& m : store)
    {
      for (const auto& kv : m) if (kv.second) kv.second->Write();
    }
  };

  auto hasAnyTH2Map = [&](const std::vector<std::map<int, TH2D*>>& store)
  {
    for (const auto& m : store) if (!m.empty()) return true;
    return false;
  };

  auto hasAnyTH1Map = [&](const std::vector<std::map<int, TH1D*>>& store)
  {
    for (const auto& m : store) if (!m.empty()) return true;
    return false;
  };

  auto hasAnyGraphMap = [&](const std::vector<std::map<int, TGraph*>>& store)
  {
    for (const auto& m : store) if (!m.empty()) return true;
    return false;
  };

  fout->cd();
  TDirectory* d_all = fout->mkdir("by_layer"); d_all->cd();
  for (auto* h : h_dRphi_all) if (h) h->Write();
  //for (auto* p : prof_phi_all) if (p) p->Write();
  TDirectory* d_s0 = fout->mkdir("by_layer_side0"); d_s0->cd();
  //for (auto* p : prof_phi_side0) if (p) p->Write();
  TDirectory* d_s1 = fout->mkdir("by_layer_side1"); d_s1->cd();
  //for (auto* p : prof_phi_side1) if (p) p->Write();

  // write per-layer 2D histograms in their own directory
  fout->cd();
  if (auto* d_2d = fout->mkdir("by_layer_2d"))
  {
    if (auto* d_phi_true = d_2d->mkdir("dRphi_vs_phi_true"))
    {
      writeTH2Vector(d_phi_true, h2_dRphi_vs_phi_layer);
    }
    if (auto* d_phi_reco = d_2d->mkdir("dRphi_vs_phi_reco"))
    {
      writeTH2Vector(d_phi_reco, h2_dRphi_vs_phiReco_layer);
    }
    if (auto* d_phase = d_2d->mkdir("dRphi_vs_phase"))
    {
      writeTH2Vector(d_phase, h2_dRphi_vs_phase_layer);
    }
    if (auto* d_phase_reco = d_2d->mkdir("dRphi_vs_phase_reco"))
    {
      writeTH2Vector(d_phase_reco, h2_dRphi_vs_phaseReco_layer);
    }
    if (auto* d_npad_phase = d_2d->mkdir("npad_vs_phase"))
    {
      writeTH2Vector(d_npad_phase, h2_npad_vs_phase_layer);
    }
    if (auto* d_npad_phaseReco = d_2d->mkdir("npad_vs_phase_reco"))
    {
      writeTH2Vector(d_npad_phaseReco, h2_npad_vs_phaseReco_layer);
    }
  }

  const bool hasPadTH2 = hasAnyTH2Map(h2_dRphi_vs_phi_layer_byPads) ||
                         hasAnyTH2Map(h2_dRphi_vs_phiReco_layer_byPads) ||
                         hasAnyTH2Map(h2_dRphi_vs_phase_layer_byPads) ||
                         hasAnyTH2Map(h2_dRphi_vs_phaseReco_layer_byPads);
  const bool hasPadTH1 = hasAnyTH1Map(h_dRphi_layer_byPads) ||
                         hasAnyTH1Map(h_phase_layer_byPads) ||
                         hasAnyTH1Map(h_phaseReco_layer_byPads);

  if (hasPadTH2 || hasPadTH1)
  {
    if (auto* d_byPads = fout->mkdir("by_layer_by_npads"))
    {
      if (hasAnyTH2Map(h2_dRphi_vs_phi_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("dRphi_vs_phi_true");
        writeTH2Maps(dir, h2_dRphi_vs_phi_layer_byPads);
      }
      if (hasAnyTH2Map(h2_dRphi_vs_phiReco_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("dRphi_vs_phi_reco");
        writeTH2Maps(dir, h2_dRphi_vs_phiReco_layer_byPads);
      }
      if (hasAnyTH2Map(h2_dRphi_vs_phase_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("dRphi_vs_phase");
        writeTH2Maps(dir, h2_dRphi_vs_phase_layer_byPads);
      }
      if (hasAnyTH2Map(h2_dRphi_vs_phaseReco_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("dRphi_vs_phase_reco");
        writeTH2Maps(dir, h2_dRphi_vs_phaseReco_layer_byPads);
      }
      if (hasAnyTH1Map(h_dRphi_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("dRphi_1d");
        writeTH1Maps(dir, h_dRphi_layer_byPads);
      }
      if (hasAnyTH1Map(h_phase_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("phase_1d");
        writeTH1Maps(dir, h_phase_layer_byPads);
      }
      if (hasAnyTH1Map(h_phaseReco_layer_byPads))
      {
        auto* dir = d_byPads->mkdir("phase_reco_1d");
        writeTH1Maps(dir, h_phaseReco_layer_byPads);
      }
    }
  }

  if (auto* d_phase1d = fout->mkdir("phase_1d_by_layer"))
  {
    writeTH1Vector(d_phase1d, h_phase_layer);
  }

  if (auto* d_phaseReco1d = fout->mkdir("phase_reco_1d_by_layer"))
  {
    writeTH1Vector(d_phaseReco1d, h_phaseReco_layer);
  }

  // write graphs
  fout->cd();
  if (auto* d_gr_phi = fout->mkdir("graphs_phi_true_vs_phi_reco"))
  {
    writeGraphVector(d_gr_phi, gr_phiReco_vs_phiTrue);
    if (hasAnyGraphMap(gr_phiReco_vs_phiTrue_byPads))
    {
      auto* dir = d_gr_phi->mkdir("by_npads");
      writeGraphMaps(dir, gr_phiReco_vs_phiTrue_byPads);
    }
  }

  if (auto* d_gr_phase = fout->mkdir("graphs_phase_vs_phase_reco"))
  {
    writeGraphVector(d_gr_phase, gr_phaseReco_vs_phase);
    if (hasAnyGraphMap(gr_phaseReco_vs_phase_byPads))
    {
      auto* dir = d_gr_phase->mkdir("by_npads");
      writeGraphMaps(dir, gr_phaseReco_vs_phase_byPads);
    }
  }

  if (auto* d_lines = fout->mkdir("phi_center_lines"))
  {
    writeLineCollections(d_lines->mkdir("dRphi_vs_phi_true"), lines_phiTrue_overlays);
    writeLineCollections(d_lines->mkdir("graph_phi_true_vs_phi_reco"), lines_phiGraph_overlays);
  }

  // also store canvases with log Z so visualization defaults are preserved
  fout->cd();
  if (auto* d_cv = fout->mkdir("canvases"))
  {
    d_cv->cd();

    for (int lyr = 0; lyr < (int)h2_dRphi_vs_phi_layer.size(); ++lyr)
    {
      auto* h2 = h2_dRphi_vs_phi_layer[lyr];
      if (!h2) continue;
      auto cname = TString::Format("c_h2_dRphi_vs_phi_layer%02d", lyr);
      auto* c = new TCanvas(cname, cname, 900, 700);
      c->SetLogz(1);
      h2->Draw("COLZ");
      if (lyr < (int)lines_phiTrue_overlays.size()) drawLines(lines_phiTrue_overlays[lyr]);
      c->Write();
      c->SaveAs(TString::Format("%s.pdf", c->GetName()));
    }

    for (int lyr = 0; lyr < (int)h2_dRphi_vs_phiReco_layer.size(); ++lyr)
    {
      auto* h2 = h2_dRphi_vs_phiReco_layer[lyr];
      if (!h2) continue;
      auto cname = TString::Format("c_h2_dRphi_vs_phiReco_layer%02d", lyr);
      auto* c = new TCanvas(cname, cname, 900, 700);
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
      auto* c = new TCanvas(cname, cname, 800, 600);
      c->SetGrid();
      gr->Draw("AP");
      if (lyr < (int)lines_phiGraph_overlays.size()) drawLines(lines_phiGraph_overlays[lyr]);
      c->Write();
      c->SaveAs(TString::Format("%s.pdf", c->GetName()));
    }
  }

  fout->Write();
  //fout->Close();
}
