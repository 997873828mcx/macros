#ifndef FUN4ALL_RAWFINALTRACKPT_C
#define FUN4ALL_RAWFINALTRACKPT_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/FinalTrack.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/FinalTrackContainer.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/lib/libInModuleTracks.so)

class RawFinalTrackPtDumper : public SubsysReco
{
 public:
  RawFinalTrackPtDumper(const std::string &name,
                        const std::string &outfile,
                        const std::string &node_name,
                        const double pt_min,
                        const double pt_max,
                        const int nbins)
    : SubsysReco(name)
    , m_outfile(outfile)
    , m_node_name(node_name)
    , m_pt_min(pt_min)
    , m_pt_max(pt_max)
    , m_nbins(nbins)
  {
  }

  int Init(PHCompositeNode *) override
  {
    m_file = TFile::Open(m_outfile.c_str(), "RECREATE");
    if (!m_file || m_file->IsZombie())
    {
      std::cerr << "[RawFinalTrackPtDumper] cannot create " << m_outfile << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_tree = new TTree("rawFinalTrackTree", "Raw FINALTRACKS momentum from DST");
    m_tree->Branch("event", &m_event, "event/I");
    m_tree->Branch("track_index", &m_track_index, "track_index/I");
    m_tree->Branch("track_id", &m_track_id, "track_id/i");
    m_tree->Branch("source_full_track_id", &m_source_full_track_id, "source_full_track_id/i");
    m_tree->Branch("fit_status", &m_fit_status, "fit_status/I");
    m_tree->Branch("nclusters", &m_nclusters, "nclusters/i");
    m_tree->Branch("charge", &m_charge, "charge/F");
    m_tree->Branch("x", &m_x, "x/F");
    m_tree->Branch("y", &m_y, "y/F");
    m_tree->Branch("z", &m_z, "z/F");
    m_tree->Branch("px", &m_px, "px/F");
    m_tree->Branch("py", &m_py, "py/F");
    m_tree->Branch("pz", &m_pz, "pz/F");
    m_tree->Branch("pt", &m_pt, "pt/F");
    m_tree->Branch("p", &m_p, "p/F");
    m_tree->Branch("chi2", &m_chi2, "chi2/F");
    m_tree->Branch("ndf", &m_ndf, "ndf/F");

    m_h_pt_all = new TH1F("h_raw_finaltrack_pt_all",
                          ";raw FINALTRACKS p_{T} [GeV/c];Tracks",
                          m_nbins, m_pt_min, m_pt_max);
    m_h_pt_good = new TH1F("h_raw_finaltrack_pt_fitstatus_nonzero",
                           ";raw FINALTRACKS p_{T} [GeV/c];Tracks",
                           m_nbins, m_pt_min, m_pt_max);
    m_h_pt_pos = new TH1F("h_raw_finaltrack_pt_charge_pos",
                          ";raw FINALTRACKS p_{T} [GeV/c];Tracks",
                          m_nbins, m_pt_min, m_pt_max);
    m_h_pt_neg = new TH1F("h_raw_finaltrack_pt_charge_neg",
                          ";raw FINALTRACKS p_{T} [GeV/c];Tracks",
                          m_nbins, m_pt_min, m_pt_max);

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode *topNode) override
  {
    auto *tracks = findNode::getClass<FinalTrackContainer>(topNode, m_node_name);
    if (!tracks)
    {
      if (!m_warned_missing_node)
      {
        std::cerr << "[RawFinalTrackPtDumper] missing node " << m_node_name << std::endl;
        m_warned_missing_node = true;
      }
      ++m_event;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    for (unsigned int i = 0; i < tracks->size(); ++i)
    {
      const FinalTrack *track = tracks->get_track(i);
      if (!track)
      {
        continue;
      }

      m_track_index = static_cast<int>(i);
      m_track_id = track->get_track_id();
      m_source_full_track_id = track->get_source_full_track_id();
      m_fit_status = track->get_fit_status();
      m_nclusters = track->get_nclusters();
      m_charge = static_cast<float>(track->get_charge());
      m_x = static_cast<float>(track->get_x());
      m_y = static_cast<float>(track->get_y());
      m_z = static_cast<float>(track->get_z());
      m_px = static_cast<float>(track->get_px());
      m_py = static_cast<float>(track->get_py());
      m_pz = static_cast<float>(track->get_pz());
      m_pt = static_cast<float>(std::sqrt(m_px * m_px + m_py * m_py));
      m_p = static_cast<float>(std::sqrt(m_px * m_px + m_py * m_py + m_pz * m_pz));
      m_chi2 = static_cast<float>(track->get_chi2());
      m_ndf = static_cast<float>(track->get_ndf());

      m_tree->Fill();
      m_h_pt_all->Fill(m_pt);
      if (m_fit_status != 0)
      {
        m_h_pt_good->Fill(m_pt);
      }
      if (m_charge > 0)
      {
        m_h_pt_pos->Fill(m_pt);
      }
      else if (m_charge < 0)
      {
        m_h_pt_neg->Fill(m_pt);
      }
    }

    ++m_event;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int End(PHCompositeNode *) override
  {
    if (!m_file)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    gStyle->SetOptStat(0);
    m_h_pt_all->SetLineColor(kBlack);
    m_h_pt_all->SetLineWidth(2);
    m_h_pt_good->SetLineColor(kMagenta + 2);
    m_h_pt_good->SetLineWidth(2);
    m_h_pt_pos->SetLineColor(kRed + 1);
    m_h_pt_pos->SetLineWidth(2);
    m_h_pt_neg->SetLineColor(kBlue + 1);
    m_h_pt_neg->SetLineWidth(2);

    auto *canvas = new TCanvas("c_raw_finaltrack_pt", "raw_finaltrack_pt", 1000, 800);
    canvas->SetLogy();
    m_h_pt_all->Draw("hist");
    m_h_pt_good->Draw("hist same");
    m_h_pt_pos->Draw("hist same");
    m_h_pt_neg->Draw("hist same");

    auto *legend = new TLegend(0.54, 0.66, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(m_h_pt_all, "all FINALTRACKS", "l");
    legend->AddEntry(m_h_pt_good, "fit_status != 0", "l");
    legend->AddEntry(m_h_pt_pos, "charge > 0", "l");
    legend->AddEntry(m_h_pt_neg, "charge < 0", "l");
    legend->Draw();

    const std::string basename = m_outfile.substr(0, m_outfile.rfind(".root"));
    canvas->SaveAs((basename + ".pdf").c_str());
    canvas->SaveAs((basename + ".png").c_str());

    m_file->cd();
    m_tree->Write();
    m_h_pt_all->Write();
    m_h_pt_good->Write();
    m_h_pt_pos->Write();
    m_h_pt_neg->Write();
    canvas->Write();
    m_file->Close();

    std::cout << "[RawFinalTrackPtDumper] wrote " << m_outfile << std::endl;
    std::cout << "[RawFinalTrackPtDumper] wrote " << basename << ".pdf and .png" << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_outfile;
  std::string m_node_name;
  double m_pt_min = 0.0;
  double m_pt_max = 10.0;
  int m_nbins = 250;
  bool m_warned_missing_node = false;

  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;
  TH1F *m_h_pt_all = nullptr;
  TH1F *m_h_pt_good = nullptr;
  TH1F *m_h_pt_pos = nullptr;
  TH1F *m_h_pt_neg = nullptr;

  int m_event = 0;
  int m_track_index = 0;
  unsigned int m_track_id = 0;
  unsigned int m_source_full_track_id = 0;
  int m_fit_status = 0;
  unsigned int m_nclusters = 0;
  float m_charge = 0.0;
  float m_x = 0.0;
  float m_y = 0.0;
  float m_z = 0.0;
  float m_px = 0.0;
  float m_py = 0.0;
  float m_pz = 0.0;
  float m_pt = 0.0;
  float m_p = 0.0;
  float m_chi2 = 0.0;
  float m_ndf = 0.0;
};

int Fun4All_RawFinalTrackPt(
    const int nEvents = 10,
    const std::string &inputDst = "/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513/output_DST/HITS_clusters_seeds_79513_0.root",
    const std::string &outputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/raw_finaltrack_pt.root",
    const std::string &nodeName = "FINALTRACKS",
    const double ptMin = 0.0,
    const double ptMax = 10.0,
    const int nbins = 250)
{
  const std::string pattern_reco_libdir = "/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/lib";
  gSystem->AddDynamicPath(pattern_reco_libdir.c_str());
  const int load_pattern = gSystem->Load((pattern_reco_libdir + "/libInModuleTracks.so").c_str());
  if (load_pattern < 0)
  {
    std::cerr << "Failed to load libInModuleTracks.so" << std::endl;
    gSystem->Exit(1);
    return 1;
  }

  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto *dumper = new RawFinalTrackPtDumper("RawFinalTrackPtDumper", outputFile, nodeName, ptMin, ptMax, nbins);
  se->registerSubsystem(dumper);

  auto *input = new Fun4AllDstInputManager("RawFinalTrackPtInput");
  if (inputDst.size() >= 5 && inputDst.substr(inputDst.size() - 5) == ".list")
  {
    input->AddListFile(inputDst);
  }
  else
  {
    input->fileopen(inputDst);
  }
  se->registerInputManager(input);

  se->run(nEvents);
  se->End();
  se->PrintTimer();
  delete se;

  std::cout << "Finished raw FINALTRACKS pT dump: " << outputFile << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif
