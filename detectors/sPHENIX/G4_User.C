#ifndef MACRO_G4USER_C
#define MACRO_G4USER_C

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <Rtypes.h>  // for R__LOAD_LIBRARY
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libtrackbase_historic.so)

class PHG4Reco;

namespace Enable
{
  bool USER = false;
  int USER_VERBOSITY = 0;
}  // namespace Enable

namespace G4USER
{
  bool WRITE_RECO_PT_TREE = false;
  std::string RECO_PT_TREE_OUTPUT = "reco_pt.root";
  std::string TRACK_MAP_NAME = "SvtxTrackMap";
}  // namespace G4USER

class RecoPtTreeWriter : public SubsysReco
{
 public:
  RecoPtTreeWriter(const std::string& name,
                   const std::string& outputfile,
                   const std::string& track_map_name)
    : SubsysReco(name)
    , m_outputfile(outputfile)
    , m_track_map_name(track_map_name)
  {
  }

  int Init(PHCompositeNode* /*topNode*/) override
  {
    m_outfile = TFile::Open(m_outputfile.c_str(), "RECREATE");
    if (!m_outfile)
    {
      std::cout << Name() << ": failed to open output file " << m_outputfile << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_track_tree = new TTree("reco_tracks", "Slim reconstructed track tree");
    m_track_tree->Branch("event", &m_track_event, "event/I");
    m_track_tree->Branch("track_id", &m_track_id, "track_id/i");
    m_track_tree->Branch("charge", &m_charge, "charge/I");
    m_track_tree->Branch("crossing", &m_crossing, "crossing/I");
    m_track_tree->Branch("nclusters", &m_nclusters, "nclusters/I");
    m_track_tree->Branch("quality", &m_quality, "quality/F");
    m_track_tree->Branch("pt", &m_pt, "pt/F");
    m_track_tree->Branch("eta", &m_eta, "eta/F");
    m_track_tree->Branch("phi", &m_phi, "phi/F");

    m_event_tree = new TTree("event_info", "Per-event reconstructed track counts");
    m_event_tree->Branch("event", &m_event, "event/I");
    m_event_tree->Branch("nreco", &m_event_nreco, "nreco/I");

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    auto* trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name.c_str());
    if (!trackmap && !m_warned_missing_track_map)
    {
      std::cout << Name() << ": missing track map " << m_track_map_name << std::endl;
      m_warned_missing_track_map = true;
    }

    m_event_nreco = 0;
    if (trackmap)
    {
      for (auto it = trackmap->begin(); it != trackmap->end(); ++it)
      {
        auto* track = it->second;
        if (!track)
        {
          continue;
        }

        m_pt = track->get_pt();
        if (!std::isfinite(m_pt))
        {
          continue;
        }

        m_track_event = m_event;
        m_track_id = track->get_id();
        m_charge = track->get_charge();
        m_crossing = track->get_crossing();
        m_nclusters = static_cast<int>(track->size_cluster_keys());
        if (m_nclusters == 0)
        {
          m_nclusters = static_cast<int>(track->size_clusters());
        }
        m_quality = track->get_quality();
        m_eta = track->get_eta();
        m_phi = track->get_phi();

        m_track_tree->Fill();
        ++m_event_nreco;
        ++m_total_tracks;
      }
    }

    m_event_tree->Fill();
    ++m_total_events;
    ++m_event;

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int End(PHCompositeNode* /*topNode*/) override
  {
    if (!m_outfile)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    m_outfile->cd();
    if (m_track_tree)
    {
      m_track_tree->Write();
    }
    if (m_event_tree)
    {
      m_event_tree->Write();
    }
    m_outfile->Close();

    std::cout << Name() << ": wrote " << m_total_tracks
              << " reconstructed tracks from " << m_total_events
              << " events to " << m_outputfile << std::endl;

    delete m_outfile;
    m_outfile = nullptr;
    m_track_tree = nullptr;
    m_event_tree = nullptr;

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_outputfile;
  std::string m_track_map_name;

  TFile* m_outfile{nullptr};
  TTree* m_track_tree{nullptr};
  TTree* m_event_tree{nullptr};

  int m_event{0};
  int m_track_event{0};
  unsigned int m_track_id{0};
  int m_charge{0};
  int m_crossing{0};
  int m_nclusters{0};
  int m_event_nreco{0};

  float m_quality{NAN};
  float m_pt{NAN};
  float m_eta{NAN};
  float m_phi{NAN};

  std::size_t m_total_tracks{0};
  std::size_t m_total_events{0};
  bool m_warned_missing_track_map{false};
};

void UserInit()
{
}

void UserDetector(PHG4Reco* /*g4Reco*/)
{
}

void UserAnalysisInit()
{
  if (!G4USER::WRITE_RECO_PT_TREE)
  {
    return;
  }

  if (G4USER::RECO_PT_TREE_OUTPUT.empty())
  {
    std::cout << "UserAnalysisInit: empty RECO_PT_TREE_OUTPUT, skipping reco pt tree writer" << std::endl;
    return;
  }

  Fun4AllServer* se = Fun4AllServer::instance();
  auto* reco_pt_tree_writer =
      new RecoPtTreeWriter("RecoPtTreeWriter", G4USER::RECO_PT_TREE_OUTPUT, G4USER::TRACK_MAP_NAME);
  reco_pt_tree_writer->Verbosity(Enable::USER_VERBOSITY);
  se->registerSubsystem(reco_pt_tree_writer);
}

#endif
