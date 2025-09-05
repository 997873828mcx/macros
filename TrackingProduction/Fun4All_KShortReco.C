/*
 * This macro shows a working example of running TrackSeeding over the cluster DST
 * This has track residuals as default output but has KFParticle set up with a togglable flag
 * with the default set up for K Short reconstruction
 */

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_Clustering.C>
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
#include <globalvertex/GlobalVertexReco.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/TpcSiliconQA.h>

#include <trackingdiagnostics/KshortReconstruction.h>
#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>

// #include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#include <stdio.h>

// R__LOAD_LIBRARY(libkfparticle_sphenix.so)

R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
void Fun4All_KShortReco(
    const int nEvents = 5,
    const std::string trackfilename = "DST_TRKR_TRACKS_run2pp_ana475_2024p017_v001-00053877-00000.root",
    const std::string dir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana475_2024p017_v001/DST_TRKR_TRACKS/run_00053800_00053900/dst/",
    const std::string outfilename = "clusters_seeds")
    // const bool convertSeeds = false,
    // const bool doKFParticle = false

{
  // std::string inputseedRawHitFile = dir + seedfilename;
  std::string inputTrackFile = dir + trackfilename;

  std::pair<int, int>
      runseg = Fun4AllUtils::GetRunSegment(trackfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;

  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  TpcReadoutInit(runnumber);
  // these lines show how to override the drift velocity and time offset values set in TpcReadoutInit
  // G4TPC::tpc_drift_velocity_reco = 0.0073844; // cm/ns
  // TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
  // G4TPC::tpc_tzero_reco = -5*50;  // ns
  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  string outDir = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/myKShortReco/";
  string outputFileName = "outputKShortReco_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
  string outputRecoDir = outDir + "inReconstruction/";
  string outputRecoFile = outputRecoDir + outputFileName;

  string makeMainDirectory = "mkdir -p " + outDir;
  system(makeMainDirectory.c_str());
  string makeDirectory = "mkdir -p " + outputRecoDir;
  system(makeDirectory.c_str());

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  G4TRACKING::SC_CALIBMODE = false;
  TRACKING::pp_mode = true;

  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  auto se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  std::time_t now = std::time(nullptr);
  std::tm tm = *std::localtime(&now);
  char dateStr[9];
  std::strftime(dateStr, sizeof(dateStr), "%m%d", &tm);

  TrackingInit();

  auto tracks = new Fun4AllDstInputManager("TrackInputManager");
  tracks->fileopen(inputTrackFile);
  se->registerInputManager(tracks);

  GlobalVertexReco *gblvertex = new GlobalVertexReco();
  gblvertex->Verbosity(0);
  se->registerSubsystem(gblvertex);

  auto ks0reco = new KshortReconstruction("KshortReconstruction");
  ks0reco->Verbosity(5);

  ks0reco->setPtCut(0.000000001);
  ks0reco->setApplyInvariantPtCut(false);
  ks0reco->setApplyQualityCut(false);
  ks0reco->setApplyDCACut(false);
  ks0reco->setApplyPairDCACut(false);
  ks0reco->setRequireMVTX(false);
  ks0reco->setApplyTrackPtCut(false);
  ks0reco->setApplyChargeCut(false);
  // ks0reco->setTrackQualityCut(100000);
  // ks0reco->setPairDCACut(0.3);
  // ks0reco->setTrackDCACut(0.01);
  ks0reco->set_output_file(outputRecoFile);
  se->registerSubsystem(ks0reco);

  se->run(nEvents);
  se->End();
  se->PrintTimer();
  // CDBInterface::instance()->Print();

  ifstream file(outputRecoFile.c_str());
  if (file.good())
  {
    string moveOutput = "mv " + outputRecoFile + " " + outDir;
    system(moveOutput.c_str());
  }

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
