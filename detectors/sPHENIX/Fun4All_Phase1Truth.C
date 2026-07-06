#ifndef FUN4ALL_PHASE1TRUTH_C
#define FUN4ALL_PHASE1TRUTH_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/TpcTruthEventTree.h>
#include <g4tpc/TpcTruthV0DecayTree.h>

#include <TSystem.h>

#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4tpc.so)

int Fun4All_Phase1Truth(
    const int nEvents = 100000,
    const std::string &inputDst = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pp_minbias_truthdst_100k_zero_vtx/dst_files.list",
    const std::string &eventOutputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/phase1_truth_events_100k_zero_vtx.root",
    const std::string &v0OutputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/phase1_truth_v0_decays_100k_zero_vtx.root")
{
  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto *truthEvent = new TpcTruthEventTree("TpcTruthEventTree", eventOutputFile);
  truthEvent->set_truth_info_node("G4TruthInfo");
  truthEvent->set_use_truth_primary_vertex(true);
  truthEvent->set_charged_eta_max(1.1);
  truthEvent->set_charged_pt_min(0.2);
  truthEvent->set_v0_abs_y_max(1.1);
  truthEvent->set_v0_pt_min(0.0);
  truthEvent->Verbosity(1);
  se->registerSubsystem(truthEvent);

  auto *truthV0 = new TpcTruthV0DecayTree("TpcTruthV0DecayTree", v0OutputFile);
  truthV0->set_truth_info_node("G4TruthInfo");
  truthV0->set_use_truth_primary_vertex(true);
  truthV0->Verbosity(1);
  se->registerSubsystem(truthV0);

  auto *input = new Fun4AllDstInputManager("Phase1TruthInput");
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

  std::cout << "Finished Phase 1 truth trees:" << std::endl;
  std::cout << "  event tree: " << eventOutputFile << std::endl;
  std::cout << "  V0 tree:    " << v0OutputFile << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif
