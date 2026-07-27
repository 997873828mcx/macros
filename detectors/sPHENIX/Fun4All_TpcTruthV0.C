#ifndef FUN4ALL_TPCTRUTHV0_C
#define FUN4ALL_TPCTRUTHV0_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/TpcV0CandidateTree.h>

#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4tpc.so)

int Fun4All_TpcTruthV0(
    const int nEvents = 10,
    const std::string &inputDst = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_truthpoints_dst.root",
    const std::string &outputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_truth_v0_candidates.root",
    const std::string &fitMethod = "helix",
    const std::string &pointOrder = "path",
    const double kalmanSigmaRPhiCm = 0.03,
    const double kalmanSigmaRCm = 0.03,
    const double kalmanSigmaZCm = 0.05,
    const bool writeKalmanInnovationDiagnostics = false,
    const bool kalmanAnalyticUniform = true,
    const double bfieldT = 1.4,
    const int coarseSteps = 64,
    const int pcaCandidates = 32)
{
  if (nEvents <= 0 || outputFile.empty() ||
      kalmanSigmaRPhiCm < 0.0 || kalmanSigmaRCm < 0.0 || kalmanSigmaZCm < 0.0 ||
      !std::isfinite(bfieldT) || bfieldT == 0.0 ||
      coarseSteps < 8 || pcaCandidates <= 0)
  {
    std::cout << "Fun4All_TpcTruthV0: invalid configuration" << std::endl;
    return 1;
  }

  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto *v0 = new TpcV0CandidateTree("TpcV0CandidateTree", outputFile);
  v0->set_truth_point_node("G4HIT_TPC_TRUECLUSTER");
  v0->set_truth_info_node("G4TruthInfo");
  v0->set_use_truth_primary_vertex(true);
  if (!v0->set_track_fit_method(fitMethod) || !v0->set_point_order(pointOrder))
  {
    std::cout << "Fun4All_TpcTruthV0: unsupported fit method or point order: fitMethod="
              << fitMethod << " pointOrder=" << pointOrder << std::endl;
    delete v0;
    delete se;
    return 1;
  }
  v0->set_fit_first_points(8);
  v0->set_bfield(bfieldT);
  v0->use_kalman_field_map(false);
  v0->set_kalman_analytic_uniform_propagation(kalmanAnalyticUniform);
  v0->set_kalman_measurement_sigmas(kalmanSigmaRPhiCm, kalmanSigmaRCm,
                                    kalmanSigmaZCm);
  v0->set_write_kalman_innovation_diagnostics(writeKalmanInnovationDiagnostics);
  v0->set_write_cluster_residual_tree(false);
  v0->set_theta_extension(2.0);
  v0->set_kalman_search(80.0, 5.0);
  v0->set_coarse_steps(coarseSteps);
  v0->set_pca_candidates(pcaCandidates);

  // Loose preselection before the expensive helix-helix PCA.
  v0->set_pre_track_pt_min(0.05);
  v0->set_pre_track_dca_xy_min(0.01);
  v0->set_pre_pair_dca_max(5.0);
  v0->set_pre_lproj_min(0.1);
  v0->Verbosity(1);
  se->registerSubsystem(v0);

  std::cout << "Fun4All_TpcTruthV0 configuration:"
            << " fit_method=" << fitMethod
            << " point_order=" << pointOrder
            << " bfield_t=" << bfieldT
            << " kalman_analytic_uniform=" << kalmanAnalyticUniform
            << " sigma_rphi_cm=" << kalmanSigmaRPhiCm
            << " sigma_r_cm=" << kalmanSigmaRCm
            << " sigma_z_cm=" << kalmanSigmaZCm
            << " innovation_diagnostics=" << writeKalmanInnovationDiagnostics
            << " coarse_steps=" << coarseSteps
            << " pca_candidates=" << pcaCandidates
            << std::endl;

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
