#ifndef MACRO_FUN4ALLREADTRKRHITSET_C
#define MACRO_FUN4ALLREADTRKRHITSET_C

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffaobjects.so)
R__LOAD_LIBRARY(libtrack_io.so)
R__LOAD_LIBRARY(libtrack.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4tpc.so)
R__LOAD_LIBRARY(libmicromegas.so)

class TrkrHitSetDumper : public SubsysReco
{
 public:
  TrkrHitSetDumper(const std::string& flatOutputFile, const int maxPrint)
    : SubsysReco("TrkrHitSetDumper")
    , m_flatOutputFile(flatOutputFile)
    , m_maxPrint(maxPrint)
  {
  }

  int Init(PHCompositeNode*) override
  {
    if (!m_flatOutputFile.empty())
    {
      m_file = TFile::Open(m_flatOutputFile.c_str(), "RECREATE");
      if (!m_file || m_file->IsZombie())
      {
        std::cout << "TrkrHitSetDumper: could not open output file "
                  << m_flatOutputFile << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      m_tree = new TTree("trkr_hits", "Flat dump of TRKR_HITSET");
      m_tree->Branch("event", &m_event, "event/I");
      m_tree->Branch("trkrid", &m_trkrid, "trkrid/I");
      m_tree->Branch("layer", &m_layer, "layer/I");
      m_tree->Branch("sector", &m_sector, "sector/I");
      m_tree->Branch("side", &m_side, "side/I");
      m_tree->Branch("hitsetkey", &m_hitsetkey, "hitsetkey/i");
      m_tree->Branch("hitkey", &m_hitkey, "hitkey/i");
      m_tree->Branch("pad", &m_pad, "pad/I");
      m_tree->Branch("tbin", &m_tbin, "tbin/I");
      m_tree->Branch("adc", &m_adc, "adc/i");
      m_tree->Branch("energy", &m_energy, "energy/F");
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    auto* hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (!hitsets)
    {
      std::cout << "TrkrHitSetDumper: missing TRKR_HITSET" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    const auto* eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
    m_event = eventHeader ? eventHeader->get_EvtSequence() : m_eventIndex;

    uint64_t eventHitsets = 0;
    uint64_t eventHits = 0;
    uint64_t eventTpcHitsets = 0;
    uint64_t eventTpcHits = 0;

    const auto hitsetRange = hitsets->getHitSets();
    for (auto hitsetIter = hitsetRange.first; hitsetIter != hitsetRange.second; ++hitsetIter)
    {
      ++eventHitsets;
      ++m_totalHitsets;

      m_hitsetkey = hitsetIter->first;
      m_trkrid = TrkrDefs::getTrkrId(m_hitsetkey);
      m_layer = TrkrDefs::getLayer(m_hitsetkey);
      m_sector = -1;
      m_side = -1;

      const bool isTpc = (m_trkrid == TrkrDefs::tpcId);
      if (isTpc)
      {
        ++eventTpcHitsets;
        m_sector = TpcDefs::getSectorId(m_hitsetkey);
        m_side = TpcDefs::getSide(m_hitsetkey);
      }

      const auto hitRange = hitsetIter->second->getHits();
      for (auto hitIter = hitRange.first; hitIter != hitRange.second; ++hitIter)
      {
        ++eventHits;
        ++m_totalHits;
        if (isTpc)
        {
          ++eventTpcHits;
          ++m_totalTpcHits;
        }

        m_hitkey = hitIter->first;
        m_pad = isTpc ? TpcDefs::getPad(m_hitkey) : -1;
        m_tbin = isTpc ? TpcDefs::getTBin(m_hitkey) : -1;

        const auto* hit = hitIter->second;
        m_adc = hit ? hit->getAdc() : 0;
        m_energy = hit ? hit->getEnergy() : 0;

        if (m_printed < m_maxPrint)
        {
          std::cout << "event " << m_event
                    << " trkrid " << m_trkrid
                    << " layer " << m_layer
                    << " sector " << m_sector
                    << " side " << m_side
                    << " hitsetkey " << m_hitsetkey
                    << " hitkey " << m_hitkey
                    << " pad " << m_pad
                    << " tbin " << m_tbin
                    << " adc " << m_adc
                    << " energy " << m_energy
                    << std::endl;
          ++m_printed;
        }

        if (m_tree)
        {
          m_tree->Fill();
        }
      }
    }

    std::cout << "TrkrHitSetDumper event " << m_event
              << ": hitsets=" << eventHitsets
              << ", hits=" << eventHits
              << ", tpc_hitsets=" << eventTpcHitsets
              << ", tpc_hits=" << eventTpcHits
              << std::endl;

    ++m_eventIndex;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int End(PHCompositeNode*) override
  {
    if (m_file)
    {
      m_file->cd();
      if (m_tree)
      {
        m_tree->Write();
      }
      m_file->Close();
      std::cout << "TrkrHitSetDumper wrote " << m_flatOutputFile << std::endl;
    }

    std::cout << "TrkrHitSetDumper totals: hitsets=" << m_totalHitsets
              << ", hits=" << m_totalHits
              << ", tpc_hits=" << m_totalTpcHits
              << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_flatOutputFile;
  int m_maxPrint = 20;
  int m_printed = 0;
  int m_eventIndex = 0;

  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;

  int m_event = 0;
  int m_trkrid = -1;
  int m_layer = -1;
  int m_sector = -1;
  int m_side = -1;
  unsigned int m_hitsetkey = 0;
  unsigned int m_hitkey = 0;
  int m_pad = -1;
  int m_tbin = -1;
  unsigned int m_adc = 0;
  float m_energy = 0;

  uint64_t m_totalHitsets = 0;
  uint64_t m_totalHits = 0;
  uint64_t m_totalTpcHits = 0;
};

int Fun4All_ReadTrkrHitSet(
    const std::string& inputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/DST_trkrhitset.root",
    const int nEvents = 1,
    const std::string& flatOutputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/DST_trkrhitset_flat.root",
    const int maxPrint = 20)
{
  auto* se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto* in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(inputFile);
  se->registerInputManager(in);

  auto* dumper = new TrkrHitSetDumper(flatOutputFile, maxPrint);
  se->registerSubsystem(dumper);

  se->run(nEvents);
  se->End();
  delete se;

  std::cout << "Finished reading " << inputFile << std::endl;
  return 0;
}

#endif
