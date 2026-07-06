#ifndef FUN4ALL_TPCPATTERNTRACKPOINTS_C
#define FUN4ALL_TPCPATTERNTRACKPOINTS_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/TpcPatternTrackPointTree.h>

#include <TSystem.h>

#include <iostream>
#include <string>

int Fun4All_TpcPatternTrackPoints(
    const int nEvents = 10,
    const std::string &inputDst = "/sphenix/user/dcxchenxi/develope/TPC_pattern_reco_Analyse/macro/good_dst_files.list",
    const std::string &outputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_pattern_track_points.root")
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

  const int load_g4tpc = gSystem->Load("libg4tpc.so");
  if (load_g4tpc < 0)
  {
    std::cerr << "Failed to load libg4tpc.so" << std::endl;
    gSystem->Exit(1);
    return 1;
  }

  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto *points = new TpcPatternTrackPointTree("TpcPatternTrackPointTree", outputFile);
  points->set_pattern_cluster_track_node("TPCPOLYCLUSTERTRACKS");
  points->set_pattern_final_track_node("FINALTRACKS");
  points->set_primary_vertex(0.0, 0.0, 0.0);
  points->set_min_points(5);
  points->Verbosity(1);
  se->registerSubsystem(points);

  auto *input = new Fun4AllDstInputManager("TPCPatternRecoInput");
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

  std::cout << "Finished TPC pattern track-point tree: " << outputFile << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif
