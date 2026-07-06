#ifndef FUN4ALL_TPCTRUTHV0_C
#define FUN4ALL_TPCTRUTHV0_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/TpcV0CandidateTree.h>

#include <TSystem.h>

#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4tpc.so)

int Fun4All_TpcTruthV0(
    const int nEvents = 10,
    const std::string &inputDst = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_truthpoints_dst.root",
    const std::string &outputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_truth_v0_candidates.root")
{
  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto *v0 = new TpcV0CandidateTree("TpcV0CandidateTree", outputFile);
  v0->set_truth_point_node("G4HIT_TPC_TRUECLUSTER");
  v0->set_truth_info_node("G4TruthInfo");
  v0->set_use_truth_primary_vertex(true);
  v0->set_fit_helix(true);
  v0->set_fit_first_points(8);
  v0->set_theta_extension(2.0);
  v0->set_coarse_steps(64);
  v0->set_pca_candidates(32);

  // Loose preselection before the expensive helix-helix PCA.
  v0->set_pre_track_pt_min(0.05);
  v0->set_pre_track_dca_xy_min(0.01);
  v0->set_pre_pair_dca_max(5.0);
  v0->set_pre_lproj_min(0.1);
  v0->Verbosity(1);
  se->registerSubsystem(v0);

  auto *input = new Fun4AllDstInputManager("TPCTruthPointInput");
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

  std::cout << "Finished TPC truth V0 candidate tree: " << outputFile << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif
