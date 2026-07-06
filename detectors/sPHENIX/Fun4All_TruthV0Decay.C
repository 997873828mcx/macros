#ifndef FUN4ALL_TRUTHV0DECAY_C
#define FUN4ALL_TRUTHV0DECAY_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/TpcTruthV0DecayTree.h>

#include <TSystem.h>

#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4tpc.so)

int Fun4All_TruthV0Decay(
    const int nEvents = 100000,
    const std::string &inputDst = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pp_minbias_truthdst_100k_zero_vtx/dst_files.list",
    const std::string &outputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/truth_v0_decays_100k_zero_vtx.root")
{
  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto *truthV0 = new TpcTruthV0DecayTree("TpcTruthV0DecayTree", outputFile);
  truthV0->set_truth_info_node("G4TruthInfo");
  truthV0->set_use_truth_primary_vertex(true);
  truthV0->Verbosity(1);
  se->registerSubsystem(truthV0);

  auto *input = new Fun4AllDstInputManager("TruthV0Input");
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

  std::cout << "Finished truth V0 decay tree: " << outputFile << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif
