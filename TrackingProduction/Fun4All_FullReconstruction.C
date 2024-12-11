/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <eventdisplay/TrackerEventDisplay.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <tpcqa/TpcRawHitQA.h>
#include <trackingqa/TpcSeedsQA.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
#include <trackingdiagnostics/KshortReconstruction.h>
#include <stdio.h>

/*#include <float.h>

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#pragma GCC diagnostic pop*/

R__LOAD_LIBRARY(libkfparticle_sphenix.so)

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libEventDisplay.so)
R__LOAD_LIBRARY(libtpcqa.so)

void Fun4All_FullReconstruction(
    const int nEvents = 10,
    const std::string tpcfilename = "DST_STREAMING_EVENT_run2pp_new_2024p002-00053217-00000.root",
    const std::string tpcdir = "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/new_2024p002/run_00053200_00053300/",
    const std::string outfilename = "clusters_seeds",
    const bool convertSeeds = false)
{
  std::string inputtpcRawHitFile = tpcdir + tpcfilename;

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
  std::pair<int, int>
      runseg = Fun4AllUtils::GetRunSegment(tpcfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;

  std::cout<< " run: " << runnumber
	   << " samples: " << TRACKING::reco_tpc_maxtime_sample
	   << " pre: " << TRACKING::reco_tpc_time_presample
	   << " vdrift: " << G4TPC::tpc_drift_velocity_reco
	   << std::endl;

 TRACKING::pp_mode = true;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  G4TRACKING::SC_CALIBMODE = false;

  ACTSGEOM::mvtxMisalignment = 100;
  ACTSGEOM::inttMisalignment = 100.;
  ACTSGEOM::tpotMisalignment = 100.;
  TString outfile = outfilename + "_" + runnumber + "-" + segment + ".root";
  std::string theOutfile = outfile.Data();

    string outDir = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/myKShortReco/";
    string outputFileName = "outputFile_kso_" + to_string(runnumber) + "_" + to_string(segment);

    string outputRecoDir = outDir + "inReconstruction/";
    string makeDirectory = "mkdir -p " + outputRecoDir;
    system(makeDirectory.c_str());
    string outputRecoFile = outputRecoDir + outputFileName + ".root";

  auto se = Fun4AllServer::instance();
  se->Verbosity(2);
  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");//准备载入几何信息

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TpcReadoutInit( runnumber );

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = false;
  //Flag for running the tpc hit unpacker with zero suppression on
  TRACKING::tpc_zero_supp = true;

  //to turn on the default static corrections, enable the two lines below
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections derived from simulation, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  //G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
  //G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
  //G4TPC:average_correction_filename = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/average_minus_static_distortion_inverted_10-new.root";

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto hitsin = new Fun4AllDstInputManager("InputManager");
  hitsin->fileopen(inputtpcRawHitFile);
  // hitsin->AddFile(inputMbd);
  se->registerInputManager(hitsin);

  Mvtx_HitUnpacking();
  Intt_HitUnpacking();
  Tpc_HitUnpacking();
 Micromegas_HitUnpacking();

  Mvtx_Clustering();
  Intt_Clustering();


  Tpc_LaserEventIdentifying();

  auto tpcclusterizer = new TpcClusterizer;//把hit打包成cluster

  tpcclusterizer->Verbosity(0);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  se->registerSubsystem(tpcclusterizer);


  Micromegas_Clustering();

  /*
   * Begin Track Seeding
   */

  /*
   * Silicon Seeding
   */


  /*
  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(0);
  silicon_Seeding->searchInIntt();
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(1.6);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);
  */

  auto silicon_Seeding = new PHActsSiliconSeeding;//种子并合并种子
  silicon_Seeding->Verbosity(0);
  // these get us to about 83% INTT > 1
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(2.0);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);

  auto merger = new PHSiliconSeedMerger;
  merger->Verbosity(0);
  se->registerSubsystem(merger);

  /*
   * Tpc Seeding
   */

    TString outfileSeed = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/" + outputFileName + "_clusters_edgeOn_staticOff_1120_0.root";
    std::string outstring(outfileSeed.Data());

  auto seeder = new PHCASeeding("PHCASeeding");
  double fieldstrength = std::numeric_limits<double>::quiet_NaN();  // set by isConstantField if constant
  bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);
  if (ConstField)
  {
    seeder->useConstBField(true);
    seeder->constBField(fieldstrength);
  }
  else
  {
    seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
    seeder->useConstBField(false);
    seeder->magFieldFile(G4MAGNET::magfield_tracking);  // to get charge sign right
  }
  seeder->Verbosity(0);
  seeder->SetLayerRange(7, 55);
  seeder->SetSearchWindow(2.,0.05); // z-width and phi-width, default in macro at 1.5 and 0.05
  seeder->SetClusAdd_delta_window(3.0,0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
  //seeder->SetNClustersPerSeedRange(4,60); // default is 6, 6
  seeder->SetMinHitsPerCluster(0);
  seeder->SetMinClustersPerTrack(3);
  seeder->useFixedClusterError(true);
  seeder->set_pp_mode(true);
  se->registerSubsystem(seeder);

  // expand stubs in the TPC using simple kalman filter
  auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
  cprop->set_field_dir(G4MAGNET::magfield_rescale);
  if (ConstField)
  {
    cprop->useConstBField(true);
    cprop->setConstBField(fieldstrength);
  }
  else
  {
    cprop->magFieldFile(G4MAGNET::magfield_tracking);
    cprop->set_field_dir(-1 * G4MAGNET::magfield_rescale);
  }
  cprop->useFixedClusterError(true);
  cprop->set_max_window(5.);
  cprop->Verbosity(0);
  cprop->set_pp_mode(true);
  se->registerSubsystem(cprop);

  // Always apply preliminary distortion corrections to TPC clusters before silicon matching
  // and refit the trackseeds. Replace KFProp fits with the new fit parameters in the TPC seeds.
  auto prelim_distcorr = new PrelimDistortionCorrection;
  prelim_distcorr->set_pp_mode(true);
  prelim_distcorr->Verbosity(0);
  se->registerSubsystem(prelim_distcorr);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
    //为什么silicon的种子没有拓展成stubs?
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  silicon_match->set_x_search_window(2.);
  silicon_match->set_y_search_window(2.);
  silicon_match->set_z_search_window(5.);
  silicon_match->set_phi_search_window(0.2);
  silicon_match->set_eta_search_window(0.1);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_rphi_search_window_lyr1(3.);
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);

  /*
   * End Track Seeding
   */

  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("TpcTrackSeedContainer");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    se->registerSubsystem(converter);
  }
  else
  {
    auto deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(0);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(cleaner);

    if (G4TRACKING::SC_CALIBMODE)
    {
      /*
      * in calibration mode, calculate residuals between TPC and fitted tracks,
      * store in dedicated structure for distortion correction
      */
      auto residuals = new PHTpcResiduals;
      const TString tpc_residoutfile = theOutfile + "_PhTpcResiduals.root";
      residuals->setOutputfile(tpc_residoutfile.Data());
      residuals->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

      // matches Tony's analysis
      residuals->setMinPt( 0.2 );

      // reconstructed distortion grid size (phi, r, z)
      residuals->setGridDimensions(36, 48, 80);
      se->registerSubsystem(residuals);
    }

  }

  auto finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(-99999.);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000000000);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);


    /*    //KFParticle setup
      KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX("myKShortReco");
      kfparticle->Verbosity(1);
      kfparticle->setDecayDescriptor("K_S0 -> pi^+ pi^-");

      //Basic node selection and configuration
      kfparticle->magFieldFile("FIELDMAP_TRACKING");
      kfparticle->getAllPVInfo(false);
      kfparticle->allowZeroMassTracks(true);
      kfparticle->useFakePrimaryVertex(true);

      kfparticle->constrainToPrimaryVertex(false);
      kfparticle->setMotherIPchi2(FLT_MAX);                               //不对IP chi2进行任何限制
      kfparticle->setFlightDistancechi2(-1.);
      kfparticle->setMinDIRA(-1.1);
      kfparticle->setDecayLengthRange(0., FLT_MAX);
      kfparticle->setDecayTimeRange(-1 * FLT_MAX, FLT_MAX);

      //Track parameters
      kfparticle->setMinMVTXhits(0);
      kfparticle->setMinTPChits(20);
      kfparticle->setMinimumTrackPT(-1.);
      kfparticle->setMaximumTrackPTchi2(FLT_MAX);
      kfparticle->setMinimumTrackIPchi2(-1.);
      kfparticle->setMinimumTrackIP(-1.);
      kfparticle->setMaximumTrackchi2nDOF(20.);

      //Vertex parameters
      kfparticle->setMaximumVertexchi2nDOF(50);
      kfparticle->setMaximumDaughterDCA(1.);

      //Parent parameters
      kfparticle->setMotherPT(0);
      kfparticle->setMinimumMass(0.300);
      kfparticle->setMaximumMass(0.700);
      kfparticle->setMaximumMotherVertexVolume(0.1);

      kfparticle->setOutputName(outputRecoFile);

      se->registerSubsystem(kfparticle);
      std::cout << "KFParticle output file: " << outputRecoFile << std::endl;*/

  TString residoutfile = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/" + outputFileName + "_resid.root";
  std::string residstring(residoutfile.Data());

  auto resid = new TrackResiduals("TrackResiduals");
  resid->outfileName(residstring);
  resid->alignment(false);

  // adjust track map name
  if(G4TRACKING::SC_CALIBMODE && !G4TRACKING::convert_seeds_to_svtxtracks)
  {
    resid->trackmapName("SvtxSiliconMMTrackMap");
    if( G4TRACKING::SC_USE_MICROMEGAS )
    { resid->set_doMicromegasOnly(true); }
  }

  resid->clusterTree();
  resid->hitTree();
  resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);

  resid->Verbosity(0);
  se->registerSubsystem(resid);

    auto ks0reco = new KshortReconstruction("KshortReconstruction");
  ks0reco->Verbosity(5);

  ks0reco->setPtCut(0.000000001);
    /*ks0reco->setApplyInvariantPtCut(false);
    ks0reco->setApplyQualityCut(false);
    ks0reco->setApplyDCACut(false);
    ks0reco->setApplyPairDCACut(false);*/
    ks0reco->setRequireMVTX(false);
    ks0reco->setTrackQualityCut(10000000000000000);
    ks0reco->setPairDCACut(10000000000000);
    ks0reco->setTrackDCACut(0.0000000000000001);
    ks0reco->set_output_file(outputRecoFile);
    se->registerSubsystem(ks0reco);
    //tree->Draw("pt", "pdg==11");

  //auto ntuplizer = new TrkrNtuplizer("TrkrNtuplizer");
  //se->registerSubsystem(ntuplizer);

  // Fun4AllOutputManager *out = new Fun4AllDstOutputManager("out", "/sphenix/tg/tg01/hf/jdosbo/tracking_development/Run24/Beam/41626/hitsets.root");
  // se->registerOutputManager(out);
  if (Enable::QA)
  {
    se->registerSubsystem(new TpcRawHitQA);
    se->registerSubsystem(new MvtxClusterQA);
    se->registerSubsystem(new InttClusterQA);
    se->registerSubsystem(new TpcClusterQA);
    se->registerSubsystem(new MicromegasClusterQA);
    auto tpcqa = new TpcSeedsQA;
    tpcqa->setTrackMapName("TpcSvtxTrackMap");
    tpcqa->setVertexMapName("TpcSvtxVertexMap");
    tpcqa->setSegment(rc->get_IntFlag("RUNSEGMENT"));
    se->registerSubsystem(tpcqa);

  }
  se->run(nEvents);
  se->End();
  se->PrintTimer();

       ifstream file(outputRecoFile.c_str());
      if (file.good())
      {
          std::cout << "Output file found: " << outputRecoFile << std::endl;
          string moveOutput = "mv " + outputRecoFile + " " + outDir;
          system(moveOutput.c_str());
          std::cout << "Moved output file to: " << outDir << std::endl;
      } else {
          std::cout << "Output file not found: " << outputRecoFile << std::endl;
      }

  if (Enable::QA)
  {
    TString qaname = outputRecoFile + "_qa.root";
    std::string qaOutputFileName(qaname.Data());
    QAHistManagerDef::saveQARootFile(qaOutputFileName);
  }

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}





/*
 * This file is part of KFParticle package
 * Copyright ( C ) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * ( at your option ) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "KFParticle_sPHENIX.h"

#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TEntryList.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <TTree.h>  // for getting the TTree from the file

#include <KFPVertex.h>
#include <KFParticle.h>           // for KFParticle
#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY...
#include <fun4all/SubsysReco.h>   // for SubsysReco

#include <ffamodules/CDBInterface.h>  // for accessing the field map file from the CDB
#include <cctype>                     // for toupper
#include <cmath>                      // for sqrt
#include <cstdlib>                    // for size_t, exit
#include <filesystem>
#include <iostream>  // for operator<<, endl, basi...
#include <map>       // for map
#include <tuple>     // for tie, tuple

class PHCompositeNode;

namespace TMVA
{
    class Reader;
}

int candidateCounter = 0;

/// KFParticle constructor
KFParticle_sPHENIX::KFParticle_sPHENIX()
        : SubsysReco("KFPARTICLE")
        , m_has_intermediates_sPHENIX(false)
        , m_constrain_to_vertex_sPHENIX(false)
        , m_require_mva(false)
        , m_save_dst(false)
        , m_save_output(true)
        , m_outfile_name("outputData.root")
        , m_outfile(nullptr)
{
}

KFParticle_sPHENIX::KFParticle_sPHENIX(const std::string &name)
        : SubsysReco(name)
        , m_has_intermediates_sPHENIX(false)
        , m_constrain_to_vertex_sPHENIX(false)
        , m_require_mva(false)
        , m_save_dst(false)
        , m_save_output(true)
        , m_outfile_name("outputData.root")
        , m_outfile(nullptr)
{
}

int KFParticle_sPHENIX::Init(PHCompositeNode *topNode)
{
    if (m_save_output && Verbosity() >= VERBOSITY_SOME)
    {
        std::cout << "Output nTuple: " << m_outfile_name << std::endl;
    }

    if (m_save_dst)
    {
        createParticleNode(topNode);
    }

    if (m_require_mva)
    {
        TMVA::Reader *reader;
        std::vector<Float_t> MVA_parValues;
        tie(reader, MVA_parValues) = initMVA();
    }

    int returnCode = 0;
    if (!m_decayDescriptor.empty())
    {
        returnCode = parseDecayDescriptor();
    }

    return returnCode;
}

int KFParticle_sPHENIX::InitRun(PHCompositeNode *topNode)
{
    assert(topNode);

    getField();

    return 0;
}

int KFParticle_sPHENIX::process_event(PHCompositeNode *topNode)
{
    std::vector<KFParticle> mother, vertex_kfparticle;
    std::vector<std::vector<KFParticle>> daughters, intermediates;
    int nPVs, multiplicity;

    if (!m_use_fake_pv)
    {
        SvtxVertexMap *check_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name);
        if (check_vertexmap->size() == 0)
        {
            if (Verbosity() >= VERBOSITY_SOME)
            {
                std::cout << "KFParticle: Event skipped as there are no vertices" << std::endl;
            }
            return Fun4AllReturnCodes::ABORTEVENT;
        }
    }

    SvtxTrackMap *check_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name);
    if (check_trackmap->size() == 0)
    {
        if (Verbosity() >= VERBOSITY_SOME)
        {
            std::cout << "KFParticle: Event skipped as there are no tracks" << std::endl;
        }
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    multiplicity = check_trackmap->size();

    createDecay(topNode, mother, vertex_kfparticle, daughters, intermediates, nPVs);

    if (!m_has_intermediates_sPHENIX)
    {
        intermediates = daughters;
    }
    if (!m_constrain_to_vertex_sPHENIX)
    {
        vertex_kfparticle = mother;
    }

    if (mother.size() != 0)
    {
        for (unsigned int i = 0; i < mother.size(); ++i)
        {
            if (m_save_output && candidateCounter == 0)
            {
                m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
                initializeBranches();
            }

            candidateCounter += 1;

            if (m_save_output)
            {
                fillBranch(topNode, mother[i], vertex_kfparticle[i], daughters[i], intermediates[i], nPVs, multiplicity);
            }
            if (m_save_dst)
            {
                fillParticleNode(topNode, mother[i], daughters[i], intermediates[i]);
            }

            if (Verbosity() >= VERBOSITY_SOME)
            {
                printParticles(mother[i], vertex_kfparticle[i], daughters[i], intermediates[i], nPVs, multiplicity);
            }
            if (Verbosity() >= VERBOSITY_MORE)
            {
                if (m_save_dst)
                {
                    printNode(topNode);
                }
            }
        }
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int KFParticle_sPHENIX::End(PHCompositeNode * /*topNode*/)
{
    std::cout << "KFParticle_sPHENIX object " << Name() << " finished. Number of candidates: " << candidateCounter << std::endl;

    if (m_save_output && candidateCounter != 0)
    {
        m_outfile->Write();
        m_outfile->Close();
        delete m_outfile;
    }

    return 0;
}

void KFParticle_sPHENIX::printParticles(const KFParticle &motherParticle,
                                        const KFParticle &chosenVertex,
                                        const std::vector<KFParticle> &daughterParticles,
                                        const std::vector<KFParticle> &intermediateParticles,
                                        const int numPVs, const int numTracks)
{
    std::cout << "\n---------------KFParticle candidate information---------------" << std::endl;

    std::cout << "Mother information:" << std::endl;
    identify(motherParticle);

    if (m_has_intermediates_sPHENIX)
    {
        std::cout << "Intermediate state information:" << std::endl;
        for (const auto &intermediateParticle : intermediateParticles)
        {
            identify(intermediateParticle);
        }
    }

    std::cout << "Final track information:" << std::endl;
    for (const auto &daughterParticle : daughterParticles)
    {
        identify(daughterParticle);
    }

    if (m_constrain_to_vertex_sPHENIX)
    {
        std::cout << "Primary vertex information:" << std::endl;
        std::cout << "(x,y,z) = (" << chosenVertex.GetX() << " +/- " << std::sqrt(chosenVertex.GetCovariance(0, 0)) << ", ";
        std::cout << chosenVertex.GetY() << " +/- " << std::sqrt(chosenVertex.GetCovariance(1, 1)) << ", ";
        std::cout << chosenVertex.GetZ() << " +/- " << std::sqrt(chosenVertex.GetCovariance(2, 2)) << ") cm\n"
                  << std::endl;
    }

    std::cout << "The number of primary vertices is: " << numPVs << std::endl;
    std::cout << "The number of tracks in the event is: " << numTracks << std::endl;

    std::cout << "------------------------------------------------------------\n"
              << std::endl;
}

int KFParticle_sPHENIX::parseDecayDescriptor()
{
    bool ddCanBeParsed = true;

    size_t daughterLocator;

    std::string mother;
    std::string intermediate;
    std::string daughter;

    std::vector<std::pair<std::string, int>> intermediate_list;
    std::vector<std::string> intermediates_name;
    std::vector<int> intermediates_charge;

    std::vector<std::pair<std::string, int>> daughter_list;
    std::vector<std::string> daughters_name;
    std::vector<int> daughters_charge;

    int nTracks = 0;
    std::vector<int> m_nTracksFromIntermediates;

    std::string decayArrow = "->";
    std::string chargeIndicator = "^";
    std::string startIntermediate = "{";
    std::string endIntermediate = "}";

    // These tracks require a + or - after their name for TDatabasePDG
    std::string specialTracks[] = {"e", "mu", "pi", "K"};

    std::string manipulateDecayDescriptor = m_decayDescriptor;

    // Remove all white space before we begin
    size_t pos;
    while ((pos = manipulateDecayDescriptor.find(' ')) != std::string::npos)
    {
        manipulateDecayDescriptor.replace(pos, 1, "");
    }

    // Check for charge conjugate requirement
    std::string checkForCC = manipulateDecayDescriptor.substr(0, 1) + manipulateDecayDescriptor.substr(manipulateDecayDescriptor.size() - 3, 3);
    std::for_each(checkForCC.begin(), checkForCC.end(), [](char &c)
    { c = ::toupper(c); });

    // Remove the CC check if needed
    if (checkForCC == "[]CC")
    {
        manipulateDecayDescriptor = manipulateDecayDescriptor.substr(1, manipulateDecayDescriptor.size() - 4);
        getChargeConjugate();
    }

    // Find the initial particle
    size_t findMotherEndPoint = manipulateDecayDescriptor.find(decayArrow);
    mother = manipulateDecayDescriptor.substr(0, findMotherEndPoint);
    if (!findParticle(mother))
    {
        ddCanBeParsed = false;
    }
    manipulateDecayDescriptor.erase(0, findMotherEndPoint + decayArrow.length());

    // Try and find the intermediates
    while ((pos = manipulateDecayDescriptor.find(startIntermediate)) != std::string::npos)
    {
        size_t findIntermediateStartPoint = manipulateDecayDescriptor.find(startIntermediate, pos);
        size_t findIntermediateEndPoint = manipulateDecayDescriptor.find(endIntermediate, pos);
        std::string intermediateDecay = manipulateDecayDescriptor.substr(pos + 1, findIntermediateEndPoint - (pos + 1));

        intermediate = intermediateDecay.substr(0, intermediateDecay.find(decayArrow));
        if (findParticle(intermediate))
        {
            intermediates_name.emplace_back(intermediate.c_str());
        }
        else
        {
            ddCanBeParsed = false;
        }

        // Now find the daughters associated to this intermediate
        int nDaughters = 0;
        intermediateDecay.erase(0, intermediateDecay.find(decayArrow) + decayArrow.length());
        while ((daughterLocator = intermediateDecay.find(chargeIndicator)) != std::string::npos)
        {
            daughter = intermediateDecay.substr(0, daughterLocator);
            std::string daughterChargeString = intermediateDecay.substr(daughterLocator + 1, 1);
            if (std::find(std::begin(specialTracks), std::end(specialTracks), daughter) != std::end(specialTracks))
            {
                daughter += daughterChargeString;
            }
            if (findParticle(daughter))
            {
                daughters_name.emplace_back(daughter.c_str());

                if (daughterChargeString == "+")
                {
                    daughters_charge.push_back(+1);
                }
                else if (daughterChargeString == "-")
                {
                    daughters_charge.push_back(-1);
                }
                else if (daughterChargeString == "0")
                {
                    daughters_charge.push_back(0);
                }
                else
                {
                    if (Verbosity() >= VERBOSITY_MORE)
                    {
                        std::cout << "The charge of " << daughterChargeString << " was not known" << std::endl;
                    }
                    ddCanBeParsed = false;
                }
            }
            else
            {
                ddCanBeParsed = false;
            }
            intermediateDecay.erase(0, daughterLocator + 2);
            ++nDaughters;
        }
        manipulateDecayDescriptor.erase(findIntermediateStartPoint, findIntermediateEndPoint + 1 - findIntermediateStartPoint);
        m_nTracksFromIntermediates.push_back(nDaughters);
        nTracks += nDaughters;
    }

    // Now find any remaining reconstructable tracks from the mother
    while ((daughterLocator = manipulateDecayDescriptor.find(chargeIndicator)) != std::string::npos)
    {
        daughter = manipulateDecayDescriptor.substr(0, daughterLocator);
        std::string daughterChargeString = manipulateDecayDescriptor.substr(daughterLocator + 1, 1);
        if (std::find(std::begin(specialTracks), std::end(specialTracks), daughter) != std::end(specialTracks))
        {
            daughter += daughterChargeString;
        }
        if (findParticle(daughter))
        {
            daughters_name.emplace_back(daughter.c_str());
            if (daughterChargeString == "+")
            {
                daughters_charge.push_back(+1);
            }
            else if (daughterChargeString == "-")
            {
                daughters_charge.push_back(-1);
            }
            else if (daughterChargeString == "0")
            {
                daughters_charge.push_back(0);
            }
            else
            {
                if (Verbosity() >= VERBOSITY_MORE)
                {
                    std::cout << "The charge of " << daughterChargeString << " was not known" << std::endl;
                }
                ddCanBeParsed = false;
            }
        }
        else
        {
            ddCanBeParsed = false;
        }
        manipulateDecayDescriptor.erase(0, daughterLocator + 2);
        nTracks += 1;
    }

    int trackEnd = 0;
    for (unsigned int i = 0; i < intermediates_name.size(); ++i)
    {
        int trackStart = trackEnd;
        trackEnd = m_nTracksFromIntermediates[i] + trackStart;

        int vtxCharge = 0;

        for (int j = trackStart; j < trackEnd; ++j)
        {
            vtxCharge += daughters_charge[j];
        }

        intermediates_charge.push_back(vtxCharge);

        intermediate_list.emplace_back(intermediates_name[i], intermediates_charge[i]);
    }

    daughter_list.reserve(nTracks);
    for (int i = 0; i < nTracks; ++i)
    {
        daughter_list.emplace_back(daughters_name[i], daughters_charge[i]);
    }

    setMotherName(mother);
    setNumberOfTracks(nTracks);
    setDaughters(daughter_list);

    if (intermediates_name.size() > 0)
    {
        hasIntermediateStates();
        setIntermediateStates(intermediate_list);
        setNumberOfIntermediateStates(intermediates_name.size());
        setNumberTracksFromIntermeditateState(m_nTracksFromIntermediates);
    }

    if (ddCanBeParsed)
    {
        if (Verbosity() >= VERBOSITY_MORE)
        {
            std::cout << "Your decay descriptor can be parsed" << std::endl;
        }
        return 0;
    }
    else
    {
        if (Verbosity() >= VERBOSITY_SOME)
        {
            std::cout << "KFParticle: Your decay descriptor, " << Name() << " cannot be parsed"
                      << "\nExiting!" << std::endl;
        }
        return Fun4AllReturnCodes::ABORTRUN;
    }
}

void KFParticle_sPHENIX::getField()
{
    //This sweeps the sPHENIX magnetic field map from some point radially then grabs the first event that passes the selection
    m_magField = std::filesystem::exists(m_magField) ? m_magField : CDBInterface::instance()->getUrl(m_magField);

    if (Verbosity() > 0)
    {
        std::cout << PHWHERE << ": using fieldmap : " << m_magField << std::endl;
    }

    TFile *fin = new TFile(m_magField.c_str());
    TTree *fieldmap = (TTree *) fin->Get("fieldmap");

    float Bz = 0.;
    unsigned int r = 0.;
    float z = 0.;

    double arc = M_PI/2;
    unsigned int n = 0;

    while (Bz == 0)
    {
        if (n == 4)
        {
            ++r;
        }

        if (r == 3) //Dont go too far out radially
        {
            ++z;
        }

        n = n & 0x3U; //Constrains n from 0 to 3
        r = r & 0x2U;

        double x = r*std::cos(n*arc);
        double y = r*std::sin(n*arc);

        std::string sweep = "x == " + std::to_string(x) + " && y == " + std::to_string(y) + " && z == " + std::to_string(z);

        fieldmap->Draw(">>elist", sweep.c_str(), "entrylist");
        TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
        if (elist->GetEntry(0))
        {
            TLeaf *fieldValue = fieldmap->GetLeaf("bz");
            fieldValue->GetBranch()->GetEntry(elist->GetEntry(0));
            Bz = fieldValue->GetValue();
        }

        ++n;

        if (r == 0) // No point in rescanning (0,0)
        {
            ++r;
            n = 0;
        }
    }
    // The actual unit of KFParticle is in kilo Gauss (kG), which is equivalent to 0.1 T, instead of Tesla (T). The positive value indicates the B field is in the +z direction
    Bz *= 10;  // Factor of 10 to convert the B field unit from kG to T
    KFParticle::SetField((double) Bz);

    fieldmap->Delete();
    fin->Close();
}
