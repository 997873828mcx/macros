#ifndef MACRO_G4USER_C
#define MACRO_G4USER_C

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4main/PHG4Particle.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <g4detectors/PHG4TpcGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <tpc/TpcClusterMover.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <Rtypes.h>  // for R__LOAD_LIBRARY
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libtrack.so)
R__LOAD_LIBRARY(libtrackbase_historic.so)
R__LOAD_LIBRARY(libtpc.so)

class PHG4Reco;

namespace Enable
{
  bool USER = false;
  int USER_VERBOSITY = 0;
}  // namespace Enable

namespace G4USER
{
  bool WRITE_RECO_PT_TREE = false;
  bool WRITE_TPC_CIRCLEFIT = false;
  bool TPC_CIRCLEFIT_USE_CLUSTER_MOVER = true;
  int TPC_CIRCLEFIT_MIN_CLUSTERS = 8;
  double TPC_CIRCLEFIT_BFIELD_TESLA = 1.4;
  int SELECTED_TRUTH_TRACK_ID = 0;
  int SELECTED_TRUTH_PID = 211;
  bool REQUIRE_SELECTED_TRUTH_PID = true;

  std::string RECO_PT_TREE_OUTPUT = "reco_pt.root";
  std::string TRACK_MAP_NAME = "SvtxTrackMap";
  std::string CLUSTER_CONTAINER_NAME = "TRKR_CLUSTER";
}  // namespace G4USER

namespace G4UserDetail
{
  struct FitPoint
  {
    double x{std::numeric_limits<double>::quiet_NaN()};
    double y{std::numeric_limits<double>::quiet_NaN()};
    unsigned int layer{0};
    unsigned int crossing{0};
  };

  bool solve3x3(double a[3][3], double b[3], double x[3])
  {
    for (int i = 0; i < 3; ++i)
    {
      int pivot = i;
      double pivot_abs = std::abs(a[i][i]);
      for (int row = i + 1; row < 3; ++row)
      {
        const double value_abs = std::abs(a[row][i]);
        if (value_abs > pivot_abs)
        {
          pivot = row;
          pivot_abs = value_abs;
        }
      }

      if (pivot_abs < 1e-16) return false;

      if (pivot != i)
      {
        for (int col = i; col < 3; ++col) std::swap(a[i][col], a[pivot][col]);
        std::swap(b[i], b[pivot]);
      }

      const double diag = a[i][i];
      for (int col = i; col < 3; ++col) a[i][col] /= diag;
      b[i] /= diag;

      for (int row = 0; row < 3; ++row)
      {
        if (row == i) continue;
        const double factor = a[row][i];
        if (std::abs(factor) < 1e-20) continue;
        for (int col = i; col < 3; ++col) a[row][col] -= factor * a[i][col];
        b[row] -= factor * b[i];
      }
    }

    x[0] = b[0];
    x[1] = b[1];
    x[2] = b[2];
    return true;
  }

  bool fit_circle_xy(const std::vector<FitPoint>& input_points,
                     double& xc,
                     double& yc,
                     double& radius_cm,
                     int& npoints_used)
  {
    xc = std::numeric_limits<double>::quiet_NaN();
    yc = std::numeric_limits<double>::quiet_NaN();
    radius_cm = std::numeric_limits<double>::quiet_NaN();
    npoints_used = 0;

    std::vector<FitPoint> points;
    points.reserve(input_points.size());
    std::set<std::pair<unsigned int, unsigned int>> seen_keys;
    for (const auto& point : input_points)
    {
      if (!std::isfinite(point.x) || !std::isfinite(point.y)) continue;
      const auto key = std::make_pair(point.layer, point.crossing);
      if (!seen_keys.insert(key).second) continue;
      points.push_back(point);
    }

    if (points.size() < 3) return false;

    const double dN = static_cast<double>(points.size());

    double mean_x = 0.0;
    double mean_y = 0.0;
    for (const auto& p : points)
    {
      mean_x += p.x;
      mean_y += p.y;
    }
    mean_x /= dN;
    mean_y /= dN;

    double Cuu = 0.0;
    double Cvv = 0.0;
    double Cuv = 0.0;
    for (const auto& p : points)
    {
      const double u = p.x - mean_x;
      const double v = p.y - mean_y;
      Cuu += u * u;
      Cvv += v * v;
      Cuv += u * v;
    }
    double phi = 0.5 * std::atan2(2.0 * Cuv, Cuu - Cvv);

    const double cphi = std::cos(phi);
    const double sphi = std::sin(phi);
    double sum_s = 0.0;
    double sum_s2 = 0.0;
    double sum_s3 = 0.0;
    double sum_s4 = 0.0;
    double sum_n = 0.0;
    double sum_sn = 0.0;
    double sum_s2n = 0.0;

    for (const auto& p : points)
    {
      const double u = p.x - mean_x;
      const double v = p.y - mean_y;
      const double s = cphi * u + sphi * v;
      const double n = -sphi * u + cphi * v;
      sum_s += s;
      sum_s2 += s * s;
      sum_s3 += s * s * s;
      sum_s4 += s * s * s * s;
      sum_n += n;
      sum_sn += s * n;
      sum_s2n += s * s * n;
    }

    double matrix[3][3] = {
        {dN, sum_s, sum_s2},
        {sum_s, sum_s2, sum_s3},
        {sum_s2, sum_s3, sum_s4}};
    double rhs[3] = {sum_n, sum_sn, sum_s2n};
    double coeffs[3] = {0.0, 0.0, 0.0};
    if (!solve3x3(matrix, rhs, coeffs)) return false;

    const double a = coeffs[0];
    const double b = coeffs[1];
    const double c = coeffs[2];
    const double denom = std::pow(1.0 + b * b, 1.5);
    if (!(denom > 0.0) || !std::isfinite(denom)) return false;
    const double kappa = 2.0 * c / denom;
    if (!std::isfinite(kappa) || !(std::abs(kappa) > 1e-14)) return false;

    const auto circle_from_local = [&](const double a_par,
                                       const double b_par,
                                       const double kappa_par,
                                       double& xc_out,
                                       double& yc_out,
                                       double& radius_out) -> bool
    {
      if (!std::isfinite(a_par) || !std::isfinite(b_par) || !std::isfinite(kappa_par)) return false;
      if (!(std::abs(kappa_par) > 1e-18)) return false;

      const double slope_norm_local = std::sqrt(1.0 + b_par * b_par);
      if (!(slope_norm_local > 0.0) || !std::isfinite(slope_norm_local)) return false;

      const double center_s_local = -b_par / (slope_norm_local * kappa_par);
      const double center_n_local = a_par + 1.0 / (slope_norm_local * kappa_par);
      const double radius_local = std::abs(1.0 / kappa_par);
      if (!(radius_local > 0.0) || !std::isfinite(radius_local)) return false;

      xc_out = mean_x + center_s_local * cphi - center_n_local * sphi;
      yc_out = mean_y + center_s_local * sphi + center_n_local * cphi;
      radius_out = radius_local;
      return std::isfinite(xc_out) && std::isfinite(yc_out);
    };

    const auto residuals_from_local = [&](const double a_par,
                                          const double b_par,
                                          const double kappa_par,
                                          std::vector<double>& residuals,
                                          double& chi2_out,
                                          double* xc_out = nullptr,
                                          double* yc_out = nullptr,
                                          double* radius_out = nullptr) -> bool
    {
      double xc_fit = 0.0;
      double yc_fit = 0.0;
      double radius_fit = 0.0;
      if (!circle_from_local(a_par, b_par, kappa_par, xc_fit, yc_fit, radius_fit)) return false;

      residuals.resize(points.size());
      chi2_out = 0.0;
      for (std::size_t ipoint = 0; ipoint < points.size(); ++ipoint)
      {
        const double dx = points[ipoint].x - xc_fit;
        const double dy = points[ipoint].y - yc_fit;
        const double distance = std::sqrt(dx * dx + dy * dy);
        if (!std::isfinite(distance)) return false;
        const double resid = distance - radius_fit;
        residuals[ipoint] = resid;
        chi2_out += resid * resid;
      }

      if (xc_out) *xc_out = xc_fit;
      if (yc_out) *yc_out = yc_fit;
      if (radius_out) *radius_out = radius_fit;
      return std::isfinite(chi2_out);
    };

    double a_cur = a;
    double b_cur = b;
    double kappa_cur = kappa;
    double xc_cur = 0.0;
    double yc_cur = 0.0;
    double r_cur = 0.0;
    std::vector<double> residuals_cur;
    double chi2_cur = 0.0;
    if (!residuals_from_local(a_cur, b_cur, kappa_cur, residuals_cur, chi2_cur, &xc_cur, &yc_cur, &r_cur))
    {
      return false;
    }

    double lambda = 1e-3;
    for (int iter = 0; iter < 30; ++iter)
    {
      const double step_a = std::max(1e-8, 1e-6 * std::max({1.0, std::abs(a_cur), 0.01 * r_cur}));
      const double step_b = std::max(1e-8, 1e-6 * std::max(1.0, std::abs(b_cur)));
      const double step_kappa = std::max(1e-12, 1e-6 * std::max(std::abs(kappa_cur), 1.0 / std::max(r_cur, 1.0)));
      const double steps[3] = {step_a, step_b, step_kappa};

      std::vector<double> jac_cols[3];
      bool jacobian_ok = true;
      for (int ipar = 0; ipar < 3; ++ipar)
      {
        jac_cols[ipar].assign(points.size(), 0.0);

        double a_plus = a_cur;
        double b_plus = b_cur;
        double kappa_plus = kappa_cur;
        double a_minus = a_cur;
        double b_minus = b_cur;
        double kappa_minus = kappa_cur;
        if (ipar == 0)
        {
          a_plus += steps[ipar];
          a_minus -= steps[ipar];
        }
        else if (ipar == 1)
        {
          b_plus += steps[ipar];
          b_minus -= steps[ipar];
        }
        else
        {
          kappa_plus += steps[ipar];
          kappa_minus -= steps[ipar];
        }

        std::vector<double> residuals_plus;
        std::vector<double> residuals_minus;
        double chi2_dummy = 0.0;
        const bool have_plus = residuals_from_local(a_plus, b_plus, kappa_plus, residuals_plus, chi2_dummy);
        const bool have_minus = residuals_from_local(a_minus, b_minus, kappa_minus, residuals_minus, chi2_dummy);
        if (!have_plus && !have_minus)
        {
          jacobian_ok = false;
          break;
        }

        for (std::size_t ipoint = 0; ipoint < points.size(); ++ipoint)
        {
          if (have_plus && have_minus)
          {
            jac_cols[ipar][ipoint] = (residuals_plus[ipoint] - residuals_minus[ipoint]) / (2.0 * steps[ipar]);
          }
          else if (have_plus)
          {
            jac_cols[ipar][ipoint] = (residuals_plus[ipoint] - residuals_cur[ipoint]) / steps[ipar];
          }
          else
          {
            jac_cols[ipar][ipoint] = (residuals_cur[ipoint] - residuals_minus[ipoint]) / steps[ipar];
          }
        }
      }
      if (!jacobian_ok) break;

      double jtj[3][3] = {};
      double jtr[3] = {};
      for (std::size_t ipoint = 0; ipoint < points.size(); ++ipoint)
      {
        for (int i = 0; i < 3; ++i)
        {
          const double ji = jac_cols[i][ipoint];
          jtr[i] += ji * residuals_cur[ipoint];
          for (int j = 0; j < 3; ++j) jtj[i][j] += ji * jac_cols[j][ipoint];
        }
      }

      bool accepted = false;
      double delta[3] = {0.0, 0.0, 0.0};
      for (int trial = 0; trial < 8; ++trial)
      {
        double system[3][3] = {
            {jtj[0][0], jtj[0][1], jtj[0][2]},
            {jtj[1][0], jtj[1][1], jtj[1][2]},
            {jtj[2][0], jtj[2][1], jtj[2][2]}};
        for (int idiag = 0; idiag < 3; ++idiag)
        {
          system[idiag][idiag] += lambda * std::max(jtj[idiag][idiag], 1.0);
        }

        double neg_jtr[3] = {-jtr[0], -jtr[1], -jtr[2]};
        if (!solve3x3(system, neg_jtr, delta))
        {
          lambda *= 10.0;
          continue;
        }

        const double a_try = a_cur + delta[0];
        const double b_try = b_cur + delta[1];
        const double kappa_try = kappa_cur + delta[2];
        std::vector<double> residuals_try;
        double chi2_try = 0.0;
        double xc_try = 0.0;
        double yc_try = 0.0;
        double r_try = 0.0;
        if (!residuals_from_local(a_try, b_try, kappa_try, residuals_try, chi2_try, &xc_try, &yc_try, &r_try) ||
            !(chi2_try < chi2_cur))
        {
          lambda *= 10.0;
          continue;
        }

        a_cur = a_try;
        b_cur = b_try;
        kappa_cur = kappa_try;
        residuals_cur.swap(residuals_try);
        chi2_cur = chi2_try;
        xc_cur = xc_try;
        yc_cur = yc_try;
        r_cur = r_try;
        lambda = std::max(1e-12, lambda * 0.3);
        accepted = true;
        break;
      }

      if (!accepted) break;

      const double rel_step = std::sqrt(
          std::pow(delta[0] / std::max({1.0, std::abs(a_cur), 0.01 * r_cur}), 2) +
          std::pow(delta[1] / std::max(1.0, std::abs(b_cur)), 2) +
          std::pow(delta[2] / std::max(std::abs(kappa_cur), 1e-12), 2));
      if (rel_step < 1e-10) break;
    }

    xc = xc_cur;
    yc = yc_cur;
    radius_cm = r_cur;
    npoints_used = static_cast<int>(points.size());
    return std::isfinite(xc) && std::isfinite(yc) && std::isfinite(radius_cm) && radius_cm > 0.0;
  }
}  // namespace G4UserDetail

class RecoPtTreeWriter : public SubsysReco
{
 public:
  RecoPtTreeWriter(const std::string& name,
                   const std::string& outputfile,
                   const std::string& track_map_name,
                   const std::string& cluster_container_name,
                   const bool write_tpc_circlefit,
                   const bool use_cluster_mover,
                   const int min_tpc_clusters,
                   const double circlefit_bfield_tesla)
    : SubsysReco(name)
    , m_outputfile(outputfile)
    , m_track_map_name(track_map_name)
    , m_cluster_container_name(cluster_container_name)
    , m_write_tpc_circlefit(write_tpc_circlefit)
    , m_use_cluster_mover(use_cluster_mover)
    , m_min_tpc_clusters(min_tpc_clusters)
    , m_circlefit_bfield_tesla(std::abs(circlefit_bfield_tesla))
  {
  }

  int Init(PHCompositeNode* /*topNode*/) override
  {
    m_outfile = TFile::Open(m_outputfile.c_str(), "RECREATE");
    if (!m_outfile)
    {
      std::cout << Name() << ": failed to open output file " << m_outputfile << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_track_tree = new TTree("reco_tracks", "Slim reconstructed track tree");
    m_track_tree->Branch("event", &m_track_event, "event/I");
    m_track_tree->Branch("track_id", &m_track_id, "track_id/i");
    m_track_tree->Branch("charge", &m_charge, "charge/I");
    m_track_tree->Branch("crossing", &m_crossing, "crossing/I");
    m_track_tree->Branch("nclusters", &m_nclusters, "nclusters/I");
    m_track_tree->Branch("quality", &m_quality, "quality/F");
    m_track_tree->Branch("pt", &m_pt, "pt/F");
    m_track_tree->Branch("eta", &m_eta, "eta/F");
    m_track_tree->Branch("phi", &m_phi, "phi/F");
    m_track_tree->Branch("best_truth_match_ok", &m_best_truth_match_ok, "best_truth_match_ok/I");
    m_track_tree->Branch("best_truth_track_id", &m_best_truth_track_id, "best_truth_track_id/I");
    m_track_tree->Branch("best_truth_pid", &m_best_truth_pid, "best_truth_pid/I");
    m_track_tree->Branch("best_truth_weight", &m_best_truth_weight, "best_truth_weight/F");

    if (m_write_tpc_circlefit)
    {
      m_track_tree->Branch("ntpc_seed_clusters", &m_ntpc_seed_clusters, "ntpc_seed_clusters/I");
      m_track_tree->Branch("ntpc_clusters", &m_ntpc_clusters, "ntpc_clusters/I");
      m_track_tree->Branch("ntpc_circle_points", &m_ntpc_circle_points, "ntpc_circle_points/I");
      m_track_tree->Branch("tpc_circle_ok", &m_tpc_circle_ok, "tpc_circle_ok/I");
      m_track_tree->Branch("tpc_circle_radius_cm", &m_tpc_circle_radius_cm, "tpc_circle_radius_cm/F");
      m_track_tree->Branch("tpc_circle_bfield_tesla", &m_tpc_circle_bfield_tesla_branch, "tpc_circle_bfield_tesla/F");
      m_track_tree->Branch("pt_tpc_circle", &m_pt_tpc_circle, "pt_tpc_circle/F");
    }

    m_selected_track_tree = new TTree("selected_reco_tracks", "One selected reconstructed track per event");
    m_selected_track_tree->Branch("event", &m_selected_event, "event/I");
    m_selected_track_tree->Branch("selected_match_ok", &m_selected_match_ok, "selected_match_ok/I");
    m_selected_track_tree->Branch("reco_track_id", &m_selected_reco_track_id, "reco_track_id/I");
    m_selected_track_tree->Branch("charge", &m_selected_charge, "charge/I");
    m_selected_track_tree->Branch("crossing", &m_selected_crossing, "crossing/I");
    m_selected_track_tree->Branch("nclusters", &m_selected_nclusters, "nclusters/I");
    m_selected_track_tree->Branch("quality", &m_selected_quality, "quality/F");
    m_selected_track_tree->Branch("pt", &m_selected_pt, "pt/F");
    m_selected_track_tree->Branch("eta", &m_selected_eta, "eta/F");
    m_selected_track_tree->Branch("phi", &m_selected_phi, "phi/F");
    m_selected_track_tree->Branch("best_truth_track_id", &m_selected_best_truth_track_id, "best_truth_track_id/I");
    m_selected_track_tree->Branch("best_truth_pid", &m_selected_best_truth_pid, "best_truth_pid/I");
    m_selected_track_tree->Branch("best_truth_weight", &m_selected_best_truth_weight, "best_truth_weight/F");
    if (m_write_tpc_circlefit)
    {
      m_selected_track_tree->Branch("ntpc_seed_clusters", &m_selected_ntpc_seed_clusters, "ntpc_seed_clusters/I");
      m_selected_track_tree->Branch("ntpc_clusters", &m_selected_ntpc_clusters, "ntpc_clusters/I");
      m_selected_track_tree->Branch("ntpc_circle_points", &m_selected_ntpc_circle_points, "ntpc_circle_points/I");
      m_selected_track_tree->Branch("tpc_circle_ok", &m_selected_tpc_circle_ok, "tpc_circle_ok/I");
      m_selected_track_tree->Branch("tpc_circle_radius_cm", &m_selected_tpc_circle_radius_cm, "tpc_circle_radius_cm/F");
      m_selected_track_tree->Branch("tpc_circle_bfield_tesla", &m_selected_tpc_circle_bfield_tesla, "tpc_circle_bfield_tesla/F");
      m_selected_track_tree->Branch("pt_tpc_circle", &m_selected_pt_tpc_circle, "pt_tpc_circle/F");
    }

    m_event_tree = new TTree("event_info", "Per-event reconstructed track counts");
    m_event_tree->Branch("event", &m_event, "event/I");
    m_event_tree->Branch("nreco", &m_event_nreco, "nreco/I");

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int InitRun(PHCompositeNode* topNode) override
  {
    if (!m_write_tpc_circlefit)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    m_cluster_container = findNode::getClass<TrkrClusterContainer>(topNode, m_cluster_container_name.c_str());
    if (!m_cluster_container)
    {
      std::cout << Name() << ": missing cluster container " << m_cluster_container_name
                << ", TPC circle fit output will be disabled" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    auto* tpcgeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
    if (m_use_cluster_mover && !tpcgeom)
    {
      std::cout << Name() << ": missing TPCGEOMCONTAINER, TPC circle fit output will be disabled" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    m_globalPositionWrapper.loadNodes(topNode);
    m_globalPositionWrapper.set_suppressCrossing(true);

    if (m_use_cluster_mover)
    {
      m_clusterMover.initialize_geometry(tpcgeom);
    }

    m_circlefit_ready = true;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    auto* trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name.c_str());
    if (!trackmap && !m_warned_missing_track_map)
    {
      std::cout << Name() << ": missing track map " << m_track_map_name << std::endl;
      m_warned_missing_track_map = true;
    }

    update_truth_matching(topNode);

    m_event_nreco = 0;
    reset_selected_track_branches();
    m_selected_event = m_event;
    if (trackmap)
    {
      for (auto it = trackmap->begin(); it != trackmap->end(); ++it)
      {
        auto* track = it->second;
        if (!track)
        {
          continue;
        }

        m_pt = track->get_pt();
        if (!std::isfinite(m_pt))
        {
          continue;
        }

        reset_circlefit_branches();
        reset_truth_match_branches();
        if (m_write_tpc_circlefit && m_circlefit_ready)
        {
          fill_tpc_circlefit(track);
        }
        fill_best_truth_match(track);

        m_track_event = m_event;
        m_track_id = track->get_id();
        m_charge = track->get_charge();
        m_crossing = track->get_crossing();
        m_nclusters = static_cast<int>(track->size_cluster_keys());
        if (m_nclusters == 0)
        {
          m_nclusters = static_cast<int>(track->size_clusters());
        }
        m_quality = track->get_quality();
        m_eta = track->get_eta();
        m_phi = track->get_phi();

        consider_selected_track();
        m_track_tree->Fill();
        ++m_event_nreco;
        ++m_total_tracks;
      }
    }

    m_selected_track_tree->Fill();
    if (m_selected_match_ok)
    {
      ++m_total_selected_matches;
    }
    m_event_tree->Fill();
    ++m_total_events;
    ++m_event;

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int End(PHCompositeNode* /*topNode*/) override
  {
    if (!m_outfile)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    m_outfile->cd();
    if (m_track_tree)
    {
      m_track_tree->Write();
    }
    if (m_selected_track_tree)
    {
      m_selected_track_tree->Write();
    }
    if (m_event_tree)
    {
      m_event_tree->Write();
    }
    m_outfile->Close();

    std::cout << Name() << ": wrote " << m_total_tracks
              << " reconstructed tracks from " << m_total_events
              << " events to " << m_outputfile
              << " with " << m_total_selected_matches
              << " selected truth-matched tracks" << std::endl;

    delete m_outfile;
    m_outfile = nullptr;
    m_track_tree = nullptr;
    m_selected_track_tree = nullptr;
    m_event_tree = nullptr;

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  void reset_circlefit_branches()
  {
    m_ntpc_seed_clusters = 0;
    m_ntpc_clusters = 0;
    m_ntpc_circle_points = 0;
    m_tpc_circle_ok = 0;
    m_tpc_circle_radius_cm = NAN;
    m_tpc_circle_bfield_tesla_branch = static_cast<float>(m_circlefit_bfield_tesla);
    m_pt_tpc_circle = NAN;
  }

  void reset_truth_match_branches()
  {
    m_best_truth_match_ok = 0;
    m_best_truth_track_id = std::numeric_limits<int>::min();
    m_best_truth_pid = std::numeric_limits<int>::min();
    m_best_truth_weight = NAN;
  }

  void reset_selected_track_branches()
  {
    m_selected_match_ok = 0;
    m_selected_reco_track_id = -1;
    m_selected_charge = 0;
    m_selected_crossing = 0;
    m_selected_nclusters = 0;
    m_selected_quality = NAN;
    m_selected_pt = NAN;
    m_selected_eta = NAN;
    m_selected_phi = NAN;
    m_selected_best_truth_track_id = std::numeric_limits<int>::min();
    m_selected_best_truth_pid = std::numeric_limits<int>::min();
    m_selected_best_truth_weight = NAN;
    m_selected_ntpc_seed_clusters = 0;
    m_selected_ntpc_clusters = 0;
    m_selected_ntpc_circle_points = 0;
    m_selected_tpc_circle_ok = 0;
    m_selected_tpc_circle_radius_cm = NAN;
    m_selected_tpc_circle_bfield_tesla = static_cast<float>(m_circlefit_bfield_tesla);
    m_selected_pt_tpc_circle = NAN;
  }

  void consider_selected_track()
  {
    if (!m_best_truth_match_ok)
    {
      return;
    }
    if (m_best_truth_track_id != m_selected_truth_track_id)
    {
      return;
    }
    if (m_require_selected_truth_pid && m_best_truth_pid != m_selected_truth_pid)
    {
      return;
    }
    if (m_selected_match_ok && std::isfinite(m_selected_best_truth_weight) &&
        std::isfinite(m_best_truth_weight) &&
        m_best_truth_weight <= m_selected_best_truth_weight)
    {
      return;
    }

    m_selected_match_ok = 1;
    m_selected_reco_track_id = static_cast<int>(m_track_id);
    m_selected_charge = m_charge;
    m_selected_crossing = m_crossing;
    m_selected_nclusters = m_nclusters;
    m_selected_quality = m_quality;
    m_selected_pt = m_pt;
    m_selected_eta = m_eta;
    m_selected_phi = m_phi;
    m_selected_best_truth_track_id = m_best_truth_track_id;
    m_selected_best_truth_pid = m_best_truth_pid;
    m_selected_best_truth_weight = m_best_truth_weight;
    m_selected_ntpc_seed_clusters = m_ntpc_seed_clusters;
    m_selected_ntpc_clusters = m_ntpc_clusters;
    m_selected_ntpc_circle_points = m_ntpc_circle_points;
    m_selected_tpc_circle_ok = m_tpc_circle_ok;
    m_selected_tpc_circle_radius_cm = m_tpc_circle_radius_cm;
    m_selected_tpc_circle_bfield_tesla = m_tpc_circle_bfield_tesla_branch;
    m_selected_pt_tpc_circle = m_pt_tpc_circle;
  }

  void update_truth_matching(PHCompositeNode* topNode)
  {
    if (!topNode)
    {
      return;
    }

    if (!m_svtxevalstack)
    {
      m_svtxevalstack.reset(new SvtxEvalStack(topNode));
      m_svtxevalstack->set_strict(false);
      m_svtxevalstack->set_use_initial_vertex(true);
      m_svtxevalstack->set_use_genfit_vertex(false);
      m_svtxevalstack->set_track_nodename(m_track_map_name);
    }

    m_svtxevalstack->set_verbosity(Verbosity());
    m_svtxevalstack->next_event(topNode);
  }

  template <class KeyIterator>
  void collect_tpc_positions(KeyIterator begin,
                             KeyIterator end,
                             const short int crossing,
                             std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>>& corrected_positions)
  {
    for (auto iter = begin; iter != end; ++iter)
    {
      const TrkrDefs::cluskey key = *iter;
      if (TrkrDefs::getTrkrId(key) != TrkrDefs::TrkrId::tpcId)
      {
        continue;
      }

      auto* cluster = m_cluster_container->findCluster(key);
      if (!cluster)
      {
        continue;
      }

      const Acts::Vector3 global =
          m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, crossing);
      corrected_positions.emplace_back(key, global);
      ++m_ntpc_clusters;
    }
  }

  void fill_tpc_circlefit(SvtxTrack* track)
  {
    if (!track || !m_cluster_container)
    {
      return;
    }

    const short int crossing = track->get_crossing();
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> corrected_positions;

    if (auto* tpc_seed = track->get_tpc_seed())
    {
      m_ntpc_seed_clusters = static_cast<int>(tpc_seed->size_cluster_keys());
      collect_tpc_positions(tpc_seed->begin_cluster_keys(), tpc_seed->end_cluster_keys(), crossing, corrected_positions);
    }
    else
    {
      collect_tpc_positions(track->begin_cluster_keys(), track->end_cluster_keys(), crossing, corrected_positions);
    }

    if (m_ntpc_clusters < m_min_tpc_clusters)
    {
      return;
    }

    const auto positions_for_fit =
        m_use_cluster_mover ? m_clusterMover.processTrack(corrected_positions) : corrected_positions;

    std::vector<G4UserDetail::FitPoint> fit_points;
    fit_points.reserve(positions_for_fit.size());
    for (const auto& [key, global] : positions_for_fit)
    {
      if (!std::isfinite(global.x()) || !std::isfinite(global.y()))
      {
        continue;
      }

      fit_points.push_back(G4UserDetail::FitPoint{
          global.x(),
          global.y(),
          static_cast<unsigned int>(TrkrDefs::getLayer(key)),
          static_cast<unsigned int>(crossing)});
    }

    double xc = NAN;
    double yc = NAN;
    double radius_cm = NAN;
    int npoints_used = 0;
    if (!G4UserDetail::fit_circle_xy(fit_points, xc, yc, radius_cm, npoints_used))
    {
      return;
    }

    if (npoints_used < m_min_tpc_clusters || !std::isfinite(radius_cm) || !(radius_cm > 0.0))
    {
      return;
    }

    const double pt_circle = 0.003 * m_circlefit_bfield_tesla * radius_cm;
    if (!std::isfinite(pt_circle) || !(pt_circle > 0.0))
    {
      return;
    }

    m_ntpc_circle_points = npoints_used;
    m_tpc_circle_ok = 1;
    m_tpc_circle_radius_cm = static_cast<float>(radius_cm);
    m_pt_tpc_circle = static_cast<float>(pt_circle);
  }

  void fill_best_truth_match(SvtxTrack* track)
  {
    if (!track || !m_svtxevalstack)
    {
      return;
    }

    auto* trackeval = m_svtxevalstack->get_track_eval();
    if (!trackeval)
    {
      return;
    }

    PHG4Particle* truth_particle = trackeval->max_truth_particle_by_nclusters(track);
    if (!truth_particle)
    {
      return;
    }

    m_best_truth_match_ok = 1;
    m_best_truth_track_id = truth_particle->get_track_id();
    m_best_truth_pid = truth_particle->get_pid();
    m_best_truth_weight = static_cast<float>(trackeval->get_nclusters_contribution(track, truth_particle));
  }

  std::string m_outputfile;
  std::string m_track_map_name;
  std::string m_cluster_container_name;

  bool m_write_tpc_circlefit{false};
  bool m_use_cluster_mover{true};
  int m_min_tpc_clusters{8};
  double m_circlefit_bfield_tesla{1.4};
  int m_selected_truth_track_id{G4USER::SELECTED_TRUTH_TRACK_ID};
  int m_selected_truth_pid{G4USER::SELECTED_TRUTH_PID};
  bool m_require_selected_truth_pid{G4USER::REQUIRE_SELECTED_TRUTH_PID};

  TFile* m_outfile{nullptr};
  TTree* m_track_tree{nullptr};
  TTree* m_selected_track_tree{nullptr};
  TTree* m_event_tree{nullptr};
  std::unique_ptr<SvtxEvalStack> m_svtxevalstack;

  TrkrClusterContainer* m_cluster_container{nullptr};
  TpcGlobalPositionWrapper m_globalPositionWrapper;
  TpcClusterMover m_clusterMover;
  bool m_circlefit_ready{false};

  int m_event{0};
  int m_track_event{0};
  int m_selected_event{0};
  unsigned int m_track_id{0};
  int m_charge{0};
  int m_crossing{0};
  int m_nclusters{0};
  int m_event_nreco{0};
  int m_best_truth_match_ok{0};
  int m_best_truth_track_id{std::numeric_limits<int>::min()};
  int m_best_truth_pid{std::numeric_limits<int>::min()};

  int m_ntpc_seed_clusters{0};
  int m_ntpc_clusters{0};
  int m_ntpc_circle_points{0};
  int m_tpc_circle_ok{0};
  int m_selected_match_ok{0};
  int m_selected_reco_track_id{-1};
  int m_selected_charge{0};
  int m_selected_crossing{0};
  int m_selected_nclusters{0};
  int m_selected_best_truth_track_id{std::numeric_limits<int>::min()};
  int m_selected_best_truth_pid{std::numeric_limits<int>::min()};
  int m_selected_ntpc_seed_clusters{0};
  int m_selected_ntpc_clusters{0};
  int m_selected_ntpc_circle_points{0};
  int m_selected_tpc_circle_ok{0};

  float m_quality{NAN};
  float m_pt{NAN};
  float m_eta{NAN};
  float m_phi{NAN};
  float m_tpc_circle_radius_cm{NAN};
  float m_tpc_circle_bfield_tesla_branch{1.4};
  float m_pt_tpc_circle{NAN};
  float m_best_truth_weight{NAN};
  float m_selected_quality{NAN};
  float m_selected_pt{NAN};
  float m_selected_eta{NAN};
  float m_selected_phi{NAN};
  float m_selected_best_truth_weight{NAN};
  float m_selected_tpc_circle_radius_cm{NAN};
  float m_selected_tpc_circle_bfield_tesla{1.4};
  float m_selected_pt_tpc_circle{NAN};

  std::size_t m_total_tracks{0};
  std::size_t m_total_events{0};
  std::size_t m_total_selected_matches{0};
  bool m_warned_missing_track_map{false};
};

void UserInit()
{
}

void UserDetector(PHG4Reco* /*g4Reco*/)
{
}

void UserAnalysisInit()
{
  if (!G4USER::WRITE_RECO_PT_TREE)
  {
    return;
  }

  if (G4USER::RECO_PT_TREE_OUTPUT.empty())
  {
    std::cout << "UserAnalysisInit: empty RECO_PT_TREE_OUTPUT, skipping reco pt tree writer" << std::endl;
    return;
  }

  Fun4AllServer* se = Fun4AllServer::instance();
  auto* reco_pt_tree_writer = new RecoPtTreeWriter(
      "RecoPtTreeWriter",
      G4USER::RECO_PT_TREE_OUTPUT,
      G4USER::TRACK_MAP_NAME,
      G4USER::CLUSTER_CONTAINER_NAME,
      G4USER::WRITE_TPC_CIRCLEFIT,
      G4USER::TPC_CIRCLEFIT_USE_CLUSTER_MOVER,
      G4USER::TPC_CIRCLEFIT_MIN_CLUSTERS,
      G4USER::TPC_CIRCLEFIT_BFIELD_TESLA);
  reco_pt_tree_writer->Verbosity(Enable::USER_VERBOSITY);
  se->registerSubsystem(reco_pt_tree_writer);
}

#endif
