#ifndef FUN4ALL_TPCPATTERNRECOV0_C
#define FUN4ALL_TPCPATTERNRECOV0_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4tpc/TpcV0CandidateTree.h>

#include <TSystem.h>

#include <iostream>
#include <string>

int Fun4All_TpcPatternRecoV0(
    const int nEvents = 10,
    const std::string &inputDst = "/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513/output_DST/HITS_clusters_seeds_79513_0.root",
    const std::string &outputFile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/tpc_pattern_v0_candidates.root",
    const double preTrackPtMin = 0.05,
    const double preTrackDcaXyMin = -1.0,
    const double prePairDcaMax = 10.0,
    const double preLprojMin = -1.0,
    const double preCosThetaMin = -2.0,
    const bool useFinalTrackHelix = false,
    const std::string &pointOrder = "auto",
    const std::string &fitMethod = "helix",
    const double kalmanSigmaRphiCm = 0.03,
    const double kalmanSigmaRCm = 0.03,
    const double kalmanSigmaZCm = 0.05,
    const bool writeSameSignPairs = false,
    const bool writeClusterResidualTree = false)
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

  auto *v0 = new TpcV0CandidateTree("TpcPatternRecoV0CandidateTree", outputFile);
  v0->use_pattern_cluster_tracks(true);
  v0->set_pattern_cluster_track_node("TPCPOLYCLUSTERTRACKS");
  v0->set_pattern_final_track_node("FINALTRACKS");
  v0->set_use_truth_primary_vertex(false);
  v0->set_primary_vertex(0.0, 0.0, 0.0);
  v0->set_fit_helix(true);
  v0->set_use_final_track_helix(useFinalTrackHelix);
  if (!v0->set_point_order(pointOrder))
  {
    gSystem->Exit(1);
    return 1;
  }
  if (!v0->set_track_fit_method(fitMethod))
  {
    gSystem->Exit(1);
    return 1;
  }
  v0->set_kalman_measurement_sigmas(kalmanSigmaRphiCm, kalmanSigmaRCm, kalmanSigmaZCm);
  v0->set_fit_first_points(0);
  v0->set_theta_extension(2.0);
  v0->set_coarse_steps(64);
  v0->set_pca_candidates(32);

  // Loose preselection before the expensive helix-helix PCA.  Negative values
  // disable cuts; for TPC-only data the primary vertex is not known reliably,
  // so the default avoids PV-dependent DCA/Lproj/pointing cuts.
  v0->set_pre_track_pt_min(preTrackPtMin);
  v0->set_pre_track_dca_xy_min(preTrackDcaXyMin);
  v0->set_pre_pair_dca_max(prePairDcaMax);
  v0->set_pre_lproj_min(preLprojMin);
  v0->set_pre_cos_theta_min(preCosThetaMin);
  v0->set_write_same_sign_pairs(writeSameSignPairs);
  v0->set_write_cluster_residual_tree(writeClusterResidualTree);
  v0->Verbosity(1);
  se->registerSubsystem(v0);

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

  std::cout << "Finished TPC pattern-reco V0 candidate tree: " << outputFile << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif
