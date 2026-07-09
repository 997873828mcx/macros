// Compare TPC-only transverse momentum resolution before/after DNL correction
// using a hybrid circle fit to the x-y projection of each track.
//
// Example:
// root -l -b -q '/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/TpcLaserDNL_StudyPtResolution_Hybrid.C("raw_dnl.root","corrected_dnl.root","pt_resolution_hybrid.root",1.4,8,20)'
// root -l -b -q '/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/TpcLaserDNL_StudyPtResolution_Hybrid.C("raw_dnl.root","corrected_dnl.root","pt_resolution_hybrid.root",1.4,8,20,5.0,20.0)'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{
  struct TrackKey
  {
    int job_id{0};
    int event{0};
    int trkid{0};

    bool operator<(const TrackKey& other) const
    {
      if (job_id != other.job_id) return job_id < other.job_id;
      if (event != other.event) return event < other.event;
      return trkid < other.trkid;
    }
  };

  struct FitPoint
  {
    double x{std::numeric_limits<double>::quiet_NaN()};
    double y{std::numeric_limits<double>::quiet_NaN()};
    unsigned int layer{0};
    unsigned int crossing{0};
  };

  struct TrackSample
  {
    double pt_true{std::numeric_limits<double>::quiet_NaN()};
    std::vector<FitPoint> points;
    int n_corrected_points{0};
  };

  struct TrackFitSummary
  {
    TrackKey key;
    double pt_true{std::numeric_limits<double>::quiet_NaN()};
    int npoints_before{0};
    int npoints_after{0};
    int ncorrected_points_after{0};
    double radius_before_cm{std::numeric_limits<double>::quiet_NaN()};
    double radius_after_cm{std::numeric_limits<double>::quiet_NaN()};
    double pt_before{std::numeric_limits<double>::quiet_NaN()};
    double pt_after{std::numeric_limits<double>::quiet_NaN()};
    double invpt_resid_before{std::numeric_limits<double>::quiet_NaN()};
    double invpt_resid_after{std::numeric_limits<double>::quiet_NaN()};
    double relpt_resid_before{std::numeric_limits<double>::quiet_NaN()};
    double relpt_resid_after{std::numeric_limits<double>::quiet_NaN()};
  };

  double safe_abs_quantile(std::vector<double> values, const double quantile)
  {
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](const double v) { return !std::isfinite(v); }),
                 values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    for (double& value : values) value = std::abs(value);
    const double q = std::clamp(quantile, 0.0, 1.0);
    const std::size_t index = static_cast<std::size_t>(q * static_cast<double>(values.size() - 1));
    std::nth_element(values.begin(), values.begin() + static_cast<long>(index), values.end());
    return values[index];
  }

  double safe_quantile(std::vector<double> values, const double quantile)
  {
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](const double v) { return !std::isfinite(v); }),
                 values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double q = std::clamp(quantile, 0.0, 1.0);
    const std::size_t index = static_cast<std::size_t>(q * static_cast<double>(values.size() - 1));
    std::nth_element(values.begin(), values.begin() + static_cast<long>(index), values.end());
    return values[index];
  }

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

    // --- Step 1: center the data ---
    // Subtracting the centroid keeps all subsequent sums O(arc-size²) rather
    // than O(R²), eliminating catastrophic cancellation for large-R tracks.
    double mean_x = 0.0, mean_y = 0.0;
    for (const auto& p : points) { mean_x += p.x; mean_y += p.y; }
    mean_x /= dN;
    mean_y /= dN;

    // --- Step 2: track direction phi from the principal axis ---
    // The eigenvector with the larger eigenvalue of the 2×2 spatial covariance
    // gives the chord direction, which equals the midpoint tangent for a
    // symmetrically sampled arc — a good phi estimate at any radius.
    double Cuu = 0.0, Cvv = 0.0, Cuv = 0.0;
    for (const auto& p : points)
    {
      const double u = p.x - mean_x;
      const double v = p.y - mean_y;
      Cuu += u * u;  Cvv += v * v;  Cuv += u * v;
    }
    double phi = 0.5 * std::atan2(2.0 * Cuv, Cuu - Cvv);

    // --- Step 3: local quadratic seed in the (s, n) track frame ---
    //
    // The previous 2-parameter parabola fit forced the linear term in n(s)
    // to vanish. That only holds if the local origin is exactly at the
    // tangent point. Here the origin is the centroid of the sampled arc, so
    // the correct local model is:
    //
    //   n(s) = a + b s + c s²
    //
    // We fit all three coefficients, then turn the local slope/curvature at
    // s=0 into an osculating-circle seed and refine that seed geometrically.
    const double cphi  = std::cos(phi);
    const double sphi  = std::sin(phi);
    double sum_s = 0.0, sum_s2 = 0.0, sum_s3 = 0.0, sum_s4 = 0.0;
    double sum_n = 0.0, sum_sn = 0.0, sum_s2n = 0.0;
    for (const auto& p : points)
    {
      const double u = p.x - mean_x;
      const double v = p.y - mean_y;
      const double s =  cphi * u + sphi * v;   // along-track coordinate
      const double n = -sphi * u + cphi * v;   // transverse (sagitta) coordinate
      sum_s   += s;
      sum_s2  += s * s;
      sum_s3  += s * s * s;
      sum_s4  += s * s * s * s;
      sum_n   += n;
      sum_sn  += s * n;
      sum_s2n += s * s * n;
    }
    double matrix[3][3] = {
        {dN,    sum_s,  sum_s2},
        {sum_s, sum_s2, sum_s3},
        {sum_s2,sum_s3, sum_s4}};
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

    // --- Step 4: convert the local quadratic to an osculating-circle seed ---
    // At s=0, the fitted point is (0,a) in the local frame and the local
    // tangent slope is b. The unit normal there is (-b,1)/sqrt(1+b²), and
    // the circle center sits one signed radius 1/kappa along that normal.
    const double slope_norm = std::sqrt(1.0 + b * b);
    const double point_s = 0.0;
    const double point_n = a;
    const double normal_s = -b / slope_norm;
    const double normal_n =  1.0 / slope_norm;
    double center_s = point_s + normal_s / kappa;
    double center_n = point_n + normal_n / kappa;
    double xc_cur = mean_x + center_s * cphi - center_n * sphi;
    double yc_cur = mean_y + center_s * sphi + center_n * cphi;
    double r_cur = std::abs(1.0 / kappa);
    if (!std::isfinite(xc_cur) || !std::isfinite(yc_cur) || !(r_cur > 0.0)) return false;

    // --- Step 5: geometric circle refinement in (xc, yc, R) ---
    // This is the previous hybrid method: curvature is used only to build the
    // seed, then the final nonlinear fit is done directly in the circle center
    // and radius parameters.
    for (int iter = 0; iter < 20; ++iter)
    {
      double jtj[3][3] = {};
      double jtr[3] = {};
      for (const auto& p : points)
      {
        const double dx = p.x - xc_cur;
        const double dy = p.y - yc_cur;
        const double d = std::sqrt(dx * dx + dy * dy);
        if (d < 1e-12) continue;

        const double resid = d - r_cur;
        const double jac[3] = {-dx / d, -dy / d, -1.0};
        for (int i = 0; i < 3; ++i)
        {
          jtr[i] += jac[i] * resid;
          for (int j = 0; j < 3; ++j) jtj[i][j] += jac[i] * jac[j];
        }
      }

      double neg_jtr[3] = {-jtr[0], -jtr[1], -jtr[2]};
      double delta[3] = {0.0, 0.0, 0.0};
      if (!solve3x3(jtj, neg_jtr, delta)) break;

      xc_cur += delta[0];
      yc_cur += delta[1];
      r_cur += delta[2];
      if (!(r_cur > 0.0) || !std::isfinite(r_cur)) return false;

      const double step = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
      if (step < 1e-10 * r_cur) break;
    }

    xc        = xc_cur;
    yc        = yc_cur;
    radius_cm = r_cur;
    npoints_used = static_cast<int>(points.size());
    return std::isfinite(xc) && std::isfinite(yc) && std::isfinite(radius_cm) && radius_cm > 0.0;
  }

  bool load_track_samples(TTree* tree,
                          const bool use_corrected_xy,
                          std::map<TrackKey, TrackSample>& out_samples,
                          std::string& x_branch_used,
                          std::string& y_branch_used)
  {
    if (!tree) return false;

    const bool has_corrected_x = (tree->GetBranch("xreco_corr") != nullptr);
    const bool has_corrected_y = (tree->GetBranch("yreco_corr") != nullptr);
    const bool has_crossing = (tree->GetBranch("crossing") != nullptr);
    const bool has_correction_applied = (tree->GetBranch("correction_applied") != nullptr);
    const bool has_job_id = (tree->GetBranch("job_id") != nullptr);
    const char* event_branch_name = (tree->GetBranch("event_local") != nullptr) ? "event_local" : "event";

    x_branch_used = (use_corrected_xy && has_corrected_x) ? "xreco_corr" : "xreco";
    y_branch_used = (use_corrected_xy && has_corrected_y) ? "yreco_corr" : "yreco";

    int job_id = 0;
    int event = 0;
    int trkid = 0;
    double pt_true = std::numeric_limits<double>::quiet_NaN();
    unsigned int layer = 0;
    unsigned int crossing = 0;
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    int correction_applied = 0;

    tree->SetBranchStatus("*", 0);
    if (has_job_id) tree->SetBranchStatus("job_id", 1);
    tree->SetBranchStatus(event_branch_name, 1);
    tree->SetBranchStatus("trkid", 1);
    tree->SetBranchStatus("pt", 1);
    tree->SetBranchStatus("layer", 1);
    tree->SetBranchStatus(x_branch_used.c_str(), 1);
    tree->SetBranchStatus(y_branch_used.c_str(), 1);
    if (has_crossing) tree->SetBranchStatus("crossing", 1);
    if (has_correction_applied) tree->SetBranchStatus("correction_applied", 1);

    if (has_job_id) tree->SetBranchAddress("job_id", &job_id);
    tree->SetBranchAddress(event_branch_name, &event);
    tree->SetBranchAddress("trkid", &trkid);
    tree->SetBranchAddress("pt", &pt_true);
    tree->SetBranchAddress("layer", &layer);
    tree->SetBranchAddress(x_branch_used.c_str(), &x);
    tree->SetBranchAddress(y_branch_used.c_str(), &y);
    if (has_crossing) tree->SetBranchAddress("crossing", &crossing);
    if (has_correction_applied) tree->SetBranchAddress("correction_applied", &correction_applied);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
      tree->GetEntry(i);
      if (!std::isfinite(pt_true) || !(pt_true > 0.0)) continue;
      if (!std::isfinite(x) || !std::isfinite(y)) continue;

      const TrackKey key{has_job_id ? job_id : 0, event, trkid};
      auto& sample = out_samples[key];
      if (!std::isfinite(sample.pt_true)) sample.pt_true = pt_true;
      sample.points.push_back(FitPoint{x, y, layer, has_crossing ? crossing : 0u});
      if (has_correction_applied && correction_applied == 1) ++sample.n_corrected_points;
    }

    return true;
  }

  TH1D* normalized_clone(const TH1D* source, const char* name)
  {
    if (!source) return nullptr;
    auto* clone = static_cast<TH1D*>(source->Clone(name));
    clone->SetDirectory(nullptr);
    const double integral = clone->Integral("width");
    if (integral > 0.0) clone->Scale(1.0 / integral);
    return clone;
  }

  TGraphErrors* make_sigma_graph(const std::vector<TrackFitSummary>& summaries,
                                 const bool use_after,
                                 const bool use_invpt,
                                 const int n_pt_bins,
                                 const double pt_min,
                                 const double pt_max,
                                 const char* name,
                                 const char* title)
  {
    std::vector<double> xvals;
    std::vector<double> yvals;
    std::vector<double> xerrs;
    std::vector<double> yerrs;

    const double width = (pt_max - pt_min) / static_cast<double>(n_pt_bins);
    for (int ibin = 0; ibin < n_pt_bins; ++ibin)
    {
      const double lo = pt_min + ibin * width;
      const double hi = lo + width;
      std::vector<double> residuals;
      for (const auto& summary : summaries)
      {
        if (!(summary.pt_true >= lo && summary.pt_true < hi)) continue;
        const double residual = use_invpt
                                    ? (use_after ? summary.invpt_resid_after : summary.invpt_resid_before)
                                    : (use_after ? summary.relpt_resid_after : summary.relpt_resid_before);
        if (!std::isfinite(residual)) continue;
        residuals.push_back(residual);
      }

      if (residuals.size() < 3) continue;

      double mean = 0.0;
      for (const double value : residuals) mean += value;
      mean /= static_cast<double>(residuals.size());

      double variance = 0.0;
      for (const double value : residuals) variance += (value - mean) * (value - mean);
      variance /= static_cast<double>(residuals.size());
      const double sigma = std::sqrt(std::max(0.0, variance));
      const double sigma_err = (residuals.size() > 1)
                                   ? sigma / std::sqrt(2.0 * static_cast<double>(residuals.size() - 1))
                                   : 0.0;

      xvals.push_back(0.5 * (lo + hi));
      xerrs.push_back(0.5 * width);
      yvals.push_back(sigma);
      yerrs.push_back(sigma_err);
    }

    auto* graph = new TGraphErrors(static_cast<int>(xvals.size()));
    graph->SetName(name);
    graph->SetTitle(title);
    for (int i = 0; i < static_cast<int>(xvals.size()); ++i)
    {
      graph->SetPoint(i, xvals[i], yvals[i]);
      graph->SetPointError(i, xerrs[i], yerrs[i]);
    }
    return graph;
  }
}  // namespace

void TpcLaserDNL_StudyPtResolution_Hybrid(const char* before_file,
                                          const char* after_file,
                                          const char* outfile = "tpc_dnl_pt_resolution_hybrid.root",
                                          double bfield_tesla = 1.4,
                                          int min_points = 8,
                                          int n_pt_bins = 20,
                                          double pt_plot_min = -1.0,
                                          double pt_plot_max = -1.0)
{
  if (!before_file || std::string(before_file).empty())
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: before_file is required\n");
    return;
  }
  if (!outfile || std::string(outfile).empty())
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: outfile is required\n");
    return;
  }
  const std::string after_path = (after_file && std::string(after_file).size() > 0) ? std::string(after_file) : std::string(before_file);
  if (!(bfield_tesla > 0.0) || !std::isfinite(bfield_tesla))
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: bfield_tesla must be positive\n");
    return;
  }
  if (min_points < 3) min_points = 3;
  if (n_pt_bins < 1) n_pt_bins = 1;

  std::unique_ptr<TFile> f_before(TFile::Open(before_file, "READ"));
  std::unique_ptr<TFile> f_after(TFile::Open(after_path.c_str(), "READ"));
  if (!f_before || f_before->IsZombie())
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: cannot open before file %s\n", before_file);
    return;
  }
  if (!f_after || f_after->IsZombie())
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: cannot open after file %s\n", after_path.c_str());
    return;
  }

  auto* tt_before = dynamic_cast<TTree*>(f_before->Get("dnl"));
  auto* tt_after = dynamic_cast<TTree*>(f_after->Get("dnl"));
  if (!tt_before)
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: missing dnl tree in %s\n", before_file);
    return;
  }
  if (!tt_after)
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: missing dnl tree in %s\n", after_path.c_str());
    return;
  }

  std::map<TrackKey, TrackSample> before_samples;
  std::map<TrackKey, TrackSample> after_samples;
  std::string before_x_branch;
  std::string before_y_branch;
  std::string after_x_branch;
  std::string after_y_branch;
  if (!load_track_samples(tt_before, false, before_samples, before_x_branch, before_y_branch))
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: failed to load before samples\n");
    return;
  }
  if (!load_track_samples(tt_after, true, after_samples, after_x_branch, after_y_branch))
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: failed to load after samples\n");
    return;
  }

  std::vector<TrackFitSummary> summaries;
  summaries.reserve(std::min(before_samples.size(), after_samples.size()));

  for (const auto& entry : before_samples)
  {
    const auto after_it = after_samples.find(entry.first);
    if (after_it == after_samples.end()) continue;

    double xc_before = 0.0;
    double yc_before = 0.0;
    double radius_before = std::numeric_limits<double>::quiet_NaN();
    int npoints_before = 0;
    if (!fit_circle_xy(entry.second.points, xc_before, yc_before, radius_before, npoints_before)) continue;

    double xc_after = 0.0;
    double yc_after = 0.0;
    double radius_after = std::numeric_limits<double>::quiet_NaN();
    int npoints_after = 0;
    if (!fit_circle_xy(after_it->second.points, xc_after, yc_after, radius_after, npoints_after)) continue;

    if (npoints_before < min_points || npoints_after < min_points) continue;

    const double pt_true = std::isfinite(entry.second.pt_true) ? entry.second.pt_true : after_it->second.pt_true;
    if (!std::isfinite(pt_true) || !(pt_true > 0.0)) continue;

    const double pt_before = 0.003 * bfield_tesla * radius_before;
    const double pt_after = 0.003 * bfield_tesla * radius_after;
    if (!(pt_before > 0.0) || !(pt_after > 0.0)) continue;

    TrackFitSummary summary;
    summary.key = entry.first;
    summary.pt_true = pt_true;
    summary.npoints_before = npoints_before;
    summary.npoints_after = npoints_after;
    summary.ncorrected_points_after = after_it->second.n_corrected_points;
    summary.radius_before_cm = radius_before;
    summary.radius_after_cm = radius_after;
    summary.pt_before = pt_before;
    summary.pt_after = pt_after;
    summary.invpt_resid_before = 1.0 / pt_before - 1.0 / pt_true;
    summary.invpt_resid_after = 1.0 / pt_after - 1.0 / pt_true;
    summary.relpt_resid_before = (pt_before - pt_true) / pt_true;
    summary.relpt_resid_after = (pt_after - pt_true) / pt_true;
    summaries.push_back(summary);
  }

  if (summaries.empty())
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: no matched tracks passed the fit and min_points=%d requirement\n", min_points);
    return;
  }

  double pt_min = std::numeric_limits<double>::infinity();
  double pt_max = -std::numeric_limits<double>::infinity();
  std::vector<double> invpt_before_vals;
  std::vector<double> invpt_after_vals;
  std::vector<double> relpt_before_vals;
  std::vector<double> relpt_after_vals;
  std::vector<double> pt_before_vals;
  std::vector<double> pt_after_vals;
  invpt_before_vals.reserve(summaries.size());
  invpt_after_vals.reserve(summaries.size());
  relpt_before_vals.reserve(summaries.size());
  relpt_after_vals.reserve(summaries.size());
  pt_before_vals.reserve(summaries.size());
  pt_after_vals.reserve(summaries.size());

  for (const auto& summary : summaries)
  {
    pt_min = std::min(pt_min, summary.pt_true);
    pt_max = std::max(pt_max, summary.pt_true);
    invpt_before_vals.push_back(summary.invpt_resid_before);
    invpt_after_vals.push_back(summary.invpt_resid_after);
    relpt_before_vals.push_back(summary.relpt_resid_before);
    relpt_after_vals.push_back(summary.relpt_resid_after);
    pt_before_vals.push_back(summary.pt_before);
    pt_after_vals.push_back(summary.pt_after);
  }

  if (!(pt_max > pt_min))
  {
    const double center = pt_min;
    pt_min = std::max(0.0, center * 0.9);
    pt_max = std::max(center * 1.1, center + 0.1);
  }

  const double invpt_range_raw = safe_abs_quantile(invpt_before_vals, 0.99);
  const double invpt_range_corr = safe_abs_quantile(invpt_after_vals, 0.99);
  const double relpt_range_raw = safe_abs_quantile(relpt_before_vals, 0.99);
  const double relpt_range_corr = safe_abs_quantile(relpt_after_vals, 0.99);

  double invpt_range = std::max(invpt_range_raw, invpt_range_corr);
  double relpt_range = std::max(relpt_range_raw, relpt_range_corr);
  if (!(invpt_range > 0.0) || !std::isfinite(invpt_range)) invpt_range = 0.1;
  if (!(relpt_range > 0.0) || !std::isfinite(relpt_range)) relpt_range = 0.5;
  invpt_range *= 1.2;
  relpt_range *= 1.2;

  if (!(pt_plot_max > pt_plot_min))
  {
    const double pt_lo_before = safe_quantile(pt_before_vals, 0.01);
    const double pt_hi_before = safe_quantile(pt_before_vals, 0.99);
    const double pt_lo_after = safe_quantile(pt_after_vals, 0.01);
    const double pt_hi_after = safe_quantile(pt_after_vals, 0.99);
    pt_plot_min = std::min(pt_lo_before, pt_lo_after);
    pt_plot_max = std::max(pt_hi_before, pt_hi_after);
    if (!(pt_plot_max > pt_plot_min) || !std::isfinite(pt_plot_min) || !std::isfinite(pt_plot_max))
    {
      pt_plot_min = std::min(safe_quantile(pt_before_vals, 0.0), safe_quantile(pt_after_vals, 0.0));
      pt_plot_max = std::max(safe_quantile(pt_before_vals, 1.0), safe_quantile(pt_after_vals, 1.0));
    }
    if (!(pt_plot_max > pt_plot_min) || !std::isfinite(pt_plot_min) || !std::isfinite(pt_plot_max))
    {
      pt_plot_min = std::max(0.0, pt_min * 0.5);
      pt_plot_max = std::max(pt_min * 1.5, pt_min + 1.0);
    }
    const double pt_span = pt_plot_max - pt_plot_min;
    const double pt_pad = std::max(0.1, 0.05 * pt_span);
    pt_plot_min = std::max(0.0, pt_plot_min - pt_pad);
    pt_plot_max += pt_pad;
  }

  auto fout = std::unique_ptr<TFile>(TFile::Open(outfile, "RECREATE"));
  if (!fout || fout->IsZombie())
  {
    printf("TpcLaserDNL_StudyPtResolution_Hybrid: cannot create output file %s\n", outfile);
    return;
  }

  auto* meta_before = new TNamed("before_file", before_file);
  auto* meta_after = new TNamed("after_file", after_path.c_str());
  const TString before_xy_title = TString::Format("%s,%s", before_x_branch.c_str(), before_y_branch.c_str());
  const TString after_xy_title = TString::Format("%s,%s", after_x_branch.c_str(), after_y_branch.c_str());
  const TString fit_config_title = TString::Format("fit_method=hybrid_seed_kappa_refine_R; bfield_tesla=%.6g; min_points=%d; n_pt_bins=%d; pt_plot_min=%.6g; pt_plot_max=%.6g; pt_fit=0.003*B[T]*R[cm]",
                                                   bfield_tesla,
                                                   min_points,
                                                   n_pt_bins,
                                                   pt_plot_min,
                                                   pt_plot_max);
  auto* meta_before_xy = new TNamed("before_xy_branches", before_xy_title.Data());
  auto* meta_after_xy = new TNamed("after_xy_branches", after_xy_title.Data());
  auto* meta_config = new TNamed("fit_config",
                                 fit_config_title.Data());
  meta_before->Write();
  meta_after->Write();
  meta_before_xy->Write();
  meta_after_xy->Write();
  meta_config->Write();

  auto* fit_tree = new TTree("ptfits", "TPC-only circle-fit pt comparison before/after DNL correction");
  int out_job_id = 0;
  int out_event = 0;
  int out_trkid = 0;
  double out_pt_true = std::numeric_limits<double>::quiet_NaN();
  int out_npoints_before = 0;
  int out_npoints_after = 0;
  int out_ncorrected_points_after = 0;
  double out_radius_before_cm = std::numeric_limits<double>::quiet_NaN();
  double out_radius_after_cm = std::numeric_limits<double>::quiet_NaN();
  double out_pt_before = std::numeric_limits<double>::quiet_NaN();
  double out_pt_after = std::numeric_limits<double>::quiet_NaN();
  double out_invpt_resid_before = std::numeric_limits<double>::quiet_NaN();
  double out_invpt_resid_after = std::numeric_limits<double>::quiet_NaN();
  double out_relpt_resid_before = std::numeric_limits<double>::quiet_NaN();
  double out_relpt_resid_after = std::numeric_limits<double>::quiet_NaN();
  double out_abs_invpt_improvement = std::numeric_limits<double>::quiet_NaN();
  double out_abs_relpt_improvement = std::numeric_limits<double>::quiet_NaN();

  fit_tree->Branch("job_id", &out_job_id, "job_id/I");
  fit_tree->Branch("event", &out_event, "event/I");
  fit_tree->Branch("trkid", &out_trkid, "trkid/I");
  fit_tree->Branch("pt_true", &out_pt_true, "pt_true/D");
  fit_tree->Branch("npoints_before", &out_npoints_before, "npoints_before/I");
  fit_tree->Branch("npoints_after", &out_npoints_after, "npoints_after/I");
  fit_tree->Branch("ncorrected_points_after", &out_ncorrected_points_after, "ncorrected_points_after/I");
  fit_tree->Branch("radius_before_cm", &out_radius_before_cm, "radius_before_cm/D");
  fit_tree->Branch("radius_after_cm", &out_radius_after_cm, "radius_after_cm/D");
  fit_tree->Branch("pt_before", &out_pt_before, "pt_before/D");
  fit_tree->Branch("pt_after", &out_pt_after, "pt_after/D");
  fit_tree->Branch("invpt_resid_before", &out_invpt_resid_before, "invpt_resid_before/D");
  fit_tree->Branch("invpt_resid_after", &out_invpt_resid_after, "invpt_resid_after/D");
  fit_tree->Branch("relpt_resid_before", &out_relpt_resid_before, "relpt_resid_before/D");
  fit_tree->Branch("relpt_resid_after", &out_relpt_resid_after, "relpt_resid_after/D");
  fit_tree->Branch("abs_invpt_improvement", &out_abs_invpt_improvement, "abs_invpt_improvement/D");
  fit_tree->Branch("abs_relpt_improvement", &out_abs_relpt_improvement, "abs_relpt_improvement/D");

  auto* h_invpt_before = new TH1D("h_invpt_resid_before", ";1/p_{T}^{fit} - 1/p_{T}^{true} [GeV^{-1}];Tracks", 200, -invpt_range, invpt_range);
  auto* h_invpt_after = new TH1D("h_invpt_resid_after", ";1/p_{T}^{fit} - 1/p_{T}^{true} [GeV^{-1}];Tracks", 200, -invpt_range, invpt_range);
  auto* h_relpt_before = new TH1D("h_relpt_resid_before", ";(p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true};Tracks", 200, -relpt_range, relpt_range);
  auto* h_relpt_after = new TH1D("h_relpt_resid_after", ";(p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true};Tracks", 200, -relpt_range, relpt_range);
  auto* h_pt_before = new TH1D("h_pt_fit_before", ";p_{T}^{fit} before correction [GeV/c];Tracks", 200, pt_plot_min, pt_plot_max);
  auto* h_pt_after = new TH1D("h_pt_fit_after", ";p_{T}^{fit} after correction [GeV/c];Tracks", 200, pt_plot_min, pt_plot_max);
  auto* h_abs_invpt_improvement = new TH1D("h_abs_invpt_improvement", ";|#Delta(1/p_{T})|_{before} - |#Delta(1/p_{T})|_{after} [GeV^{-1}];Tracks", 200, -invpt_range, invpt_range);
  auto* h_abs_relpt_improvement = new TH1D("h_abs_relpt_improvement", ";|#Delta p_{T}/p_{T}|_{before} - |#Delta p_{T}/p_{T}|_{after};Tracks", 200, -relpt_range, relpt_range);
  auto* h_npoints_before = new TH1D("h_npoints_before", ";TPC points used in before fit;Tracks", 60, -0.5, 59.5);
  auto* h_npoints_after = new TH1D("h_npoints_after", ";TPC points used in after fit;Tracks", 60, -0.5, 59.5);
  auto* h_ncorrected_points_after = new TH1D("h_ncorrected_points_after", ";Corrected TPC points contributing after fit;Tracks", 60, -0.5, 59.5);
  auto* h2_invpt_before = new TH2D("h2_invpt_resid_vs_pt_before", ";p_{T}^{true} [GeV/c];1/p_{T}^{fit} - 1/p_{T}^{true} [GeV^{-1}]", n_pt_bins, pt_min, pt_max, 160, -invpt_range, invpt_range);
  auto* h2_invpt_after = new TH2D("h2_invpt_resid_vs_pt_after", ";p_{T}^{true} [GeV/c];1/p_{T}^{fit} - 1/p_{T}^{true} [GeV^{-1}]", n_pt_bins, pt_min, pt_max, 160, -invpt_range, invpt_range);
  auto* h2_relpt_before = new TH2D("h2_relpt_resid_vs_pt_before", ";p_{T}^{true} [GeV/c];(p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true}", n_pt_bins, pt_min, pt_max, 160, -relpt_range, relpt_range);
  auto* h2_relpt_after = new TH2D("h2_relpt_resid_vs_pt_after", ";p_{T}^{true} [GeV/c];(p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true}", n_pt_bins, pt_min, pt_max, 160, -relpt_range, relpt_range);

  for (const auto& summary : summaries)
  {
    out_job_id = summary.key.job_id;
    out_event = summary.key.event;
    out_trkid = summary.key.trkid;
    out_pt_true = summary.pt_true;
    out_npoints_before = summary.npoints_before;
    out_npoints_after = summary.npoints_after;
    out_ncorrected_points_after = summary.ncorrected_points_after;
    out_radius_before_cm = summary.radius_before_cm;
    out_radius_after_cm = summary.radius_after_cm;
    out_pt_before = summary.pt_before;
    out_pt_after = summary.pt_after;
    out_invpt_resid_before = summary.invpt_resid_before;
    out_invpt_resid_after = summary.invpt_resid_after;
    out_relpt_resid_before = summary.relpt_resid_before;
    out_relpt_resid_after = summary.relpt_resid_after;
    out_abs_invpt_improvement = std::abs(summary.invpt_resid_before) - std::abs(summary.invpt_resid_after);
    out_abs_relpt_improvement = std::abs(summary.relpt_resid_before) - std::abs(summary.relpt_resid_after);
    fit_tree->Fill();

    h_invpt_before->Fill(summary.invpt_resid_before);
    h_invpt_after->Fill(summary.invpt_resid_after);
    h_relpt_before->Fill(summary.relpt_resid_before);
    h_relpt_after->Fill(summary.relpt_resid_after);
    h_pt_before->Fill(summary.pt_before);
    h_pt_after->Fill(summary.pt_after);
    h_abs_invpt_improvement->Fill(out_abs_invpt_improvement);
    h_abs_relpt_improvement->Fill(out_abs_relpt_improvement);
    h_npoints_before->Fill(summary.npoints_before);
    h_npoints_after->Fill(summary.npoints_after);
    h_ncorrected_points_after->Fill(summary.ncorrected_points_after);
    h2_invpt_before->Fill(summary.pt_true, summary.invpt_resid_before);
    h2_invpt_after->Fill(summary.pt_true, summary.invpt_resid_after);
    h2_relpt_before->Fill(summary.pt_true, summary.relpt_resid_before);
    h2_relpt_after->Fill(summary.pt_true, summary.relpt_resid_after);
  }

  auto* gr_sigma_invpt_before = make_sigma_graph(summaries, false, true, n_pt_bins, pt_min, pt_max,
                                                 "gr_sigma_invpt_before",
                                                 ";p_{T}^{true} [GeV/c];#sigma(1/p_{T}^{fit} - 1/p_{T}^{true}) [GeV^{-1}]");
  auto* gr_sigma_invpt_after = make_sigma_graph(summaries, true, true, n_pt_bins, pt_min, pt_max,
                                                "gr_sigma_invpt_after",
                                                ";p_{T}^{true} [GeV/c];#sigma(1/p_{T}^{fit} - 1/p_{T}^{true}) [GeV^{-1}]");
  auto* gr_sigma_relpt_before = make_sigma_graph(summaries, false, false, n_pt_bins, pt_min, pt_max,
                                                 "gr_sigma_relpt_before",
                                                 ";p_{T}^{true} [GeV/c];#sigma((p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true})");
  auto* gr_sigma_relpt_after = make_sigma_graph(summaries, true, false, n_pt_bins, pt_min, pt_max,
                                                "gr_sigma_relpt_after",
                                                ";p_{T}^{true} [GeV/c];#sigma((p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true})");

  auto* h_invpt_before_norm = normalized_clone(h_invpt_before, "h_invpt_resid_before_norm");
  auto* h_invpt_after_norm = normalized_clone(h_invpt_after, "h_invpt_resid_after_norm");
  auto* h_relpt_before_norm = normalized_clone(h_relpt_before, "h_relpt_resid_before_norm");
  auto* h_relpt_after_norm = normalized_clone(h_relpt_after, "h_relpt_resid_after_norm");
  auto* h_pt_before_norm = normalized_clone(h_pt_before, "h_pt_fit_before_norm");
  auto* h_pt_after_norm = normalized_clone(h_pt_after, "h_pt_fit_after_norm");

  if (h_invpt_before_norm)
  {
    h_invpt_before_norm->SetLineColor(kBlack);
    h_invpt_before_norm->SetLineWidth(2);
    h_invpt_before_norm->SetStats(0);
  }
  if (h_invpt_after_norm)
  {
    h_invpt_after_norm->SetLineColor(kRed + 1);
    h_invpt_after_norm->SetLineWidth(2);
    h_invpt_after_norm->SetStats(0);
  }
  if (h_relpt_before_norm)
  {
    h_relpt_before_norm->SetLineColor(kBlack);
    h_relpt_before_norm->SetLineWidth(2);
    h_relpt_before_norm->SetStats(0);
  }
  if (h_relpt_after_norm)
  {
    h_relpt_after_norm->SetLineColor(kRed + 1);
    h_relpt_after_norm->SetLineWidth(2);
    h_relpt_after_norm->SetStats(0);
  }
  if (h_pt_before_norm)
  {
    h_pt_before_norm->SetLineColor(kBlack);
    h_pt_before_norm->SetLineWidth(2);
    h_pt_before_norm->SetStats(0);
  }
  if (h_pt_after_norm)
  {
    h_pt_after_norm->SetLineColor(kRed + 1);
    h_pt_after_norm->SetLineWidth(2);
    h_pt_after_norm->SetStats(0);
  }

  auto* c_invpt = new TCanvas("c_invpt_residual_compare", "invpt residual comparison", 900, 700);
  if (h_invpt_before_norm) h_invpt_before_norm->Draw("hist");
  if (h_invpt_after_norm) h_invpt_after_norm->Draw(h_invpt_before_norm ? "hist same" : "hist");
  auto* leg_invpt = new TLegend(0.58, 0.74, 0.88, 0.88);
  leg_invpt->SetBorderSize(0);
  leg_invpt->AddEntry(h_invpt_before_norm, "Before correction", "l");
  leg_invpt->AddEntry(h_invpt_after_norm, "After correction", "l");
  leg_invpt->Draw();

  auto* c_relpt = new TCanvas("c_relpt_residual_compare", "relpt residual comparison", 900, 700);
  if (h_relpt_before_norm) h_relpt_before_norm->Draw("hist");
  if (h_relpt_after_norm) h_relpt_after_norm->Draw(h_relpt_before_norm ? "hist same" : "hist");
  auto* leg_relpt = new TLegend(0.58, 0.74, 0.88, 0.88);
  leg_relpt->SetBorderSize(0);
  leg_relpt->AddEntry(h_relpt_before_norm, "Before correction", "l");
  leg_relpt->AddEntry(h_relpt_after_norm, "After correction", "l");
  leg_relpt->Draw();

  auto* c_pt = new TCanvas("c_pt_fit_compare", "fitted pt comparison", 900, 700);
  if (h_pt_before_norm) h_pt_before_norm->Draw("hist");
  if (h_pt_after_norm) h_pt_after_norm->Draw(h_pt_before_norm ? "hist same" : "hist");
  auto* leg_pt = new TLegend(0.58, 0.74, 0.88, 0.88);
  leg_pt->SetBorderSize(0);
  leg_pt->AddEntry(h_pt_before_norm, "Before correction", "l");
  leg_pt->AddEntry(h_pt_after_norm, "After correction", "l");
  leg_pt->Draw();

  auto* c_sigma_invpt = new TCanvas("c_sigma_invpt_vs_pt", "sigma invpt vs pt", 900, 700);
  if (gr_sigma_invpt_before)
  {
    gr_sigma_invpt_before->SetLineColor(kBlack);
    gr_sigma_invpt_before->SetMarkerColor(kBlack);
    gr_sigma_invpt_before->SetMarkerStyle(20);
    gr_sigma_invpt_before->SetTitle(";p_{T}^{true} [GeV/c];#sigma(1/p_{T}^{fit} - 1/p_{T}^{true}) [GeV^{-1}]");
    gr_sigma_invpt_before->Draw("AP");
  }
  if (gr_sigma_invpt_after)
  {
    gr_sigma_invpt_after->SetLineColor(kRed + 1);
    gr_sigma_invpt_after->SetMarkerColor(kRed + 1);
    gr_sigma_invpt_after->SetMarkerStyle(21);
    gr_sigma_invpt_after->Draw(gr_sigma_invpt_before ? "P same" : "AP");
  }
  auto* leg_sigma_invpt = new TLegend(0.58, 0.74, 0.88, 0.88);
  leg_sigma_invpt->SetBorderSize(0);
  leg_sigma_invpt->AddEntry(gr_sigma_invpt_before, "Before correction", "lp");
  leg_sigma_invpt->AddEntry(gr_sigma_invpt_after, "After correction", "lp");
  leg_sigma_invpt->Draw();

  auto* c_sigma_relpt = new TCanvas("c_sigma_relpt_vs_pt", "sigma relpt vs pt", 900, 700);
  if (gr_sigma_relpt_before)
  {
    gr_sigma_relpt_before->SetLineColor(kBlack);
    gr_sigma_relpt_before->SetMarkerColor(kBlack);
    gr_sigma_relpt_before->SetMarkerStyle(20);
    gr_sigma_relpt_before->SetTitle(";p_{T}^{true} [GeV/c];#sigma((p_{T}^{fit} - p_{T}^{true}) / p_{T}^{true})");
    gr_sigma_relpt_before->Draw("AP");
  }
  if (gr_sigma_relpt_after)
  {
    gr_sigma_relpt_after->SetLineColor(kRed + 1);
    gr_sigma_relpt_after->SetMarkerColor(kRed + 1);
    gr_sigma_relpt_after->SetMarkerStyle(21);
    gr_sigma_relpt_after->Draw(gr_sigma_relpt_before ? "P same" : "AP");
  }
  auto* leg_sigma_relpt = new TLegend(0.58, 0.74, 0.88, 0.88);
  leg_sigma_relpt->SetBorderSize(0);
  leg_sigma_relpt->AddEntry(gr_sigma_relpt_before, "Before correction", "lp");
  leg_sigma_relpt->AddEntry(gr_sigma_relpt_after, "After correction", "lp");
  leg_sigma_relpt->Draw();

  fout->cd();
  fit_tree->Write();
  h_invpt_before->Write();
  h_invpt_after->Write();
  h_relpt_before->Write();
  h_relpt_after->Write();
  h_pt_before->Write();
  h_pt_after->Write();
  h_abs_invpt_improvement->Write();
  h_abs_relpt_improvement->Write();
  h_npoints_before->Write();
  h_npoints_after->Write();
  h_ncorrected_points_after->Write();
  h2_invpt_before->Write();
  h2_invpt_after->Write();
  h2_relpt_before->Write();
  h2_relpt_after->Write();
  if (h_invpt_before_norm) h_invpt_before_norm->Write();
  if (h_invpt_after_norm) h_invpt_after_norm->Write();
  if (h_relpt_before_norm) h_relpt_before_norm->Write();
  if (h_relpt_after_norm) h_relpt_after_norm->Write();
  if (h_pt_before_norm) h_pt_before_norm->Write();
  if (h_pt_after_norm) h_pt_after_norm->Write();
  if (gr_sigma_invpt_before) gr_sigma_invpt_before->Write();
  if (gr_sigma_invpt_after) gr_sigma_invpt_after->Write();
  if (gr_sigma_relpt_before) gr_sigma_relpt_before->Write();
  if (gr_sigma_relpt_after) gr_sigma_relpt_after->Write();
  c_invpt->Write();
  c_relpt->Write();
  c_pt->Write();
  c_sigma_invpt->Write();
  c_sigma_relpt->Write();
  fout->Write();
  fout->Close();

  printf("TpcLaserDNL_StudyPtResolution_Hybrid: matched %zu tracks\n", summaries.size());
  printf("  before branches: x=%s y=%s\n", before_x_branch.c_str(), before_y_branch.c_str());
  printf("  after branches:  x=%s y=%s\n", after_x_branch.c_str(), after_y_branch.c_str());
  printf("  output: %s\n", outfile);
}
