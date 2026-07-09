#include "PHHitsModule.h"


#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <trackbase/TrackFitUtils.h>

//#include <tpc/TpcDistortionCorrectionContainer.h>
//#include <tpc/TpcGlobalPositionWrapper.h>




#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1RawHit.h>

//#include <trackbase_historic/SvtxTrack.h>
//#include <trackbase_historic/SvtxTrackMap.h>
//#include <trackbase_historic/TrackSeedContainer.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxAlignmentStateMap.h>


#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <cmath>
#include <string>

#include <iostream>
#include <sstream>
#include "Riostream.h"
#include <vector>  

namespace
{

    // square
    template <class T>
    inline constexpr T square(const T& x)
    {
        return x * x;
    }

    // radius
    template <class T>
    inline T r(const T& x, const T& y)
    {
        return std::sqrt(square(x) + square(y));
    }

    //sector edge from in the range -pi to pi
/*    inline constexpr float phi_edge(int isec)
    {
        return M_PI - (2*M_PI/TpcDefs::NSectors)/2 - isec*(2*M_PI/TpcDefs::NSectors);
    }
    //sector phi edge in the range from 0 to 2pi
    inline constexpr float phi_edge_pos(int isec)
    {
        return ( isec < 6 ) ? (M_PI - (2*M_PI/TpcDefs::NSectors)/2 - isec*(2*M_PI/TpcDefs::NSectors)) : (M_PI - (2*M_PI/TpcDefs::NSectors)/2 - isec*(2*M_PI/TpcDefs::NSectors) + 2*M_PI);
    }

    //! get cluster keys from a given track
    std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
    {
        std::vector<TrkrDefs::cluskey> out;
        for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
        {
            if (seed)
            {
                std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
            }
        }
        return out;
    }
*/

}  // namespace

//___________________________________________________________________________________
PHHitsModule::PHHitsModule(const std::string& name)
        : SubsysReco(name)
{
    std::cout << "PHHitsModule::PHHitsModule" <<std::endl;
}

//___________________________________________________________________________________
int PHHitsModule::Init(PHCompositeNode* /*topNode*/)
{
    // configuration printout
    outfile = new TFile(m_outputfile.c_str(), "RECREATE");
    m_tree = new TTree("residualtree", "A tree with TPC hits and clusters");
 
    m_tree->Branch("event", &m_event, "m_event/I");
    m_tree->Branch("gl1bco", &m_bco, "m_bco/l");
    m_tree->Branch("trbco", &m_bcotr, "m_bcotr/l");

    m_tree->Branch("hitsetkey", &m_hitsetkey, "m_hitsetkey/i");
    m_tree->Branch("hitgx", &m_hitgx, "m_hitgx/F");
    m_tree->Branch("hitgy", &m_hitgy, "m_hitgy/F");
    m_tree->Branch("hitgz", &m_hitgz, "m_hitgz/F");
    m_tree->Branch("layer", &m_hitlayer, "m_hitlayer/I");
    m_tree->Branch("sector", &m_sector, "m_sector/I");
    m_tree->Branch("side", &m_side, "m_side/I");
    m_tree->Branch("pad", &m_hitpad, "m_hitpad/I");
    m_tree->Branch("tbin", &m_hittbin, "m_hittbin/I");
    m_tree->Branch("adc", &m_adc, "m_adc/F");
    m_tree->Branch("zdriftlength", &m_zdriftlength, "m_zdriftlength/F");

    m_tree->Branch("lx", &m_scluslx, "m_scluslx/F");
    m_tree->Branch("lz", &m_scluslz, "m_scluslz/F");
    m_tree->Branch("gx", &m_sclusgx, "m_sclusgx/F");
    m_tree->Branch("gy", &m_sclusgy, "m_sclusgy/F");
    m_tree->Branch("gz", &m_sclusgz, "m_sclusgz/F");
    m_tree->Branch("r", &m_sclusgr, "m_sclusgr/F");
    m_tree->Branch("phi", &m_sclusphi, "m_sclusphi/F");
    m_tree->Branch("eta", &m_scluseta, "m_scluseta/F");
    m_tree->Branch("clussize", &m_sclusize, "m_sclusize/I");
    m_tree->Branch("clusoverlap", &m_sclusoverlap, "m_sclusoverlap/F");
    m_tree->Branch("clusedge", &m_sclusedge, "m_sclusedge/F");
    m_tree->Branch("phisize", &m_phisize, "m_phisize/I");
    m_tree->Branch("zsize", &m_zsize, "m_zsize/I");
    m_tree->Branch("erphi", &m_scluselx, "m_scluselx/F");
    m_tree->Branch("ez", &m_scluselz, "m_scluselz/F");
    m_tree->Branch("maxadc", &m_clusmaxadc, "m_clusmaxadc/F");
    m_tree->Branch("sector", &m_clussector, "m_clussector/I");




    std::cout << "PHHitsModule::Init - EVENT_OK " << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHHitsModule::InitRun(PHCompositeNode* topNode)
{
    std::cout << "PHHitsModule::InitRun " << std::endl;
    std::cout << "PHHitsModule::InitRun - EVENT_OK " << std::endl;

    auto tpccellgeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    m_clusterMover.initialize_geometry(tpccellgeo);
    m_clusterMover.set_verbosity(0);
    auto se = Fun4AllServer::instance();
    m_runnumber = se->RunNumber();

    return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHHitsModule::process_event(PHCompositeNode* topNode)
{
    std::cout << "PHHitsModule::process_event - Event " <<m_event<< std::endl;

    auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    auto hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    auto tpcGeom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  auto gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    m_bco = gl1->get_bco();
    auto lbshift = m_bco << 24U;
    m_bcotr = lbshift >> 24U;
  }
  else
  {
    Gl1Packet* gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1PacketInfo)
    {
      m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
      m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
    }
    m_firedTriggers.clear();

    if (gl1PacketInfo)
    {
      m_gl1BunchCrossing = gl1PacketInfo->getBunchNumber();
      uint64_t triggervec = gl1PacketInfo->getTriggerVector();
      m_bco = gl1PacketInfo->getBCO();
      auto lbshift = m_bco << 24U;
      m_bcotr = lbshift >> 24U;
      for (int i = 0; i < 64; i++)
      {
        bool trig_decision = ((triggervec & 0x1U) == 0x1U);
        if (trig_decision)
        {
          m_firedTriggers.push_back(i);
        }
        triggervec = (triggervec >> 1U) & 0xffffffffU;
      }
    }
  }

    m_ntpcclus = 0;

    fillHitTree(hitmap, geometry, tpcGeom);

    fillClusterTree(clustermap, geometry);


  
    m_event++;


    return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHHitsModule::End(PHCompositeNode* /*topNode*/)
{
    std::cout << "PHHitsModule::End - writing residuals to " << m_outputfile << std::endl;
    outfile->cd();
    m_tree->Write();

    // print counters

    return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
void PHHitsModule::fillHitTree(TrkrHitSetContainer* hitmap,
                                 ActsGeometry* geometry,
                                 PHG4TpcCylinderGeomContainer* tpcGeom)
{
    if (Verbosity())
    {
        std::cout << "PHHitsModule::fillHitTree - Hitmap size " << hitmap->size() << std::endl;
    }

   if (!tpcGeom)
   {
       std::cout << PHWHERE << "missing hit map, can't continue with hit tree" << std::endl;
    return;
   }

  TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
  TrkrHitSet* hitset;
  unsigned char det;
  for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
       hitsetiter != all_hitsets.second;
       ++hitsetiter)
  {
    m_hitsetkey = hitsetiter->first;
    //TrkrHitSet* hitset = hitsetiter->second;
    hitset = hitsetiter->second;

    m_hitlayer = TrkrDefs::getLayer(m_hitsetkey);
    det = TrkrDefs::getTrkrId(m_hitsetkey);
    if (det == TrkrDefs::TrkrId::tpcId)
    {
      m_sector = TpcDefs::getSectorId(m_hitsetkey);
      m_side = TpcDefs::getSide(m_hitsetkey);
    }
  

  auto hitrangei = hitset->getHits();

  for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
  {
      auto hitkey = hitr->first;
      auto hit = hitr->second;
      m_adc = hit->getAdc();
    if (det == TrkrDefs::TrkrId::tpcId)
    {
                m_hitpad = TpcDefs::getPad(hitkey);
        m_hittbin = TpcDefs::getTBin(hitkey);

        auto geoLayer = tpcGeom->GetLayerCellGeom(m_hitlayer);
        auto phi = geoLayer->get_phicenter(m_hitpad);
        auto radius = geoLayer->get_radius();
        float AdcClockPeriod = geoLayer->get_zstep();
        auto glob = geometry->getGlobalPositionTpc(m_hitsetkey, hitkey, phi, radius, AdcClockPeriod);
        m_hitgx = glob.x();
        m_hitgy = glob.y();
        m_hitgz = glob.z();
    }
    m_tree->Fill();
  }
}
}

//_______________________________________________________________________________


void PHHitsModule::fillClusterTree(TrkrClusterContainer* clusters,
                                     ActsGeometry* geometry)
{
  if (clusters->size()< m_min_cluster_size)
  {
    return;
  }
  auto det = TrkrDefs::TrkrId::tpcId;
    for (const auto& hitsetkey : clusters->getHitSetKeys(det))
    {
      m_scluslayer = TrkrDefs::getLayer(hitsetkey);
      auto range = clusters->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        auto key = iter->first;
        auto cluster = clusters->findCluster(key);
        Acts::Vector3 glob;
        glob = geometry->getGlobalPosition(key, cluster);
        m_sclusize = cluster->getSize();
        m_sclusoverlap = cluster->getOverlap();
        m_sclusedge = cluster->getEdge();
        m_sclusgx = glob.x();
        m_sclusgy = glob.y();
        m_sclusgz = glob.z();
        m_sclusgr = r(m_sclusgx, m_sclusgy);
        m_sclusphi = atan2(glob.y(), glob.x());
        m_scluseta = acos(glob.z() / std::sqrt(square(glob.x()) + square(glob.y()) + square(glob.z())));
        m_clusadc = cluster->getAdc();
        m_clusmaxadc = cluster->getMaxAdc();
        m_scluslx = cluster->getLocalX();
        m_scluslz = cluster->getLocalY();
        auto para_errors = m_clusErrPara.get_clusterv5_modified_error(cluster, m_sclusgr, key);
        m_phisize = cluster->getPhiSize();
        m_zsize = cluster->getZSize();
        m_scluselx = std::sqrt(para_errors.first);
        m_scluselz = std::sqrt(para_errors.second);

          m_clussector = TpcDefs::getSectorId(key);
          m_side = TpcDefs::getSide(key);

       m_tree->Fill();
      }
    }
}



//____________________________________________________________________________..
void PHHitsModule::clearClusterStateVectors()
{
}





