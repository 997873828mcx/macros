// Draw presentation-style phase CDF comparisons from
// TpcLaserDNL_ApplyPhaseCDFCorrection output.

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <vector>

namespace
{
  bool apply_sphenix_style()
  {
    if (gROOT && gROOT->GetStyle("sPHENIX"))
    {
      gROOT->SetStyle("sPHENIX");
      gROOT->ForceStyle();
      return true;
    }

    std::vector<TString> style_paths;
    const TString macro_file = __FILE__;
    const Ssiz_t last_slash = macro_file.Last('/');
    if (last_slash >= 0)
    {
      const TString macro_dir = macro_file(0, last_slash);
      style_paths.push_back(macro_dir + "/../../macros/sPHENIXStyle/sPhenixStyle.C");
    }
    style_paths.push_back("/sphenix/user/dcxchenxi/develope/macros/macros/sPHENIXStyle/sPhenixStyle.C");
    style_paths.push_back("develope/macros/macros/sPHENIXStyle/sPhenixStyle.C");

    for (const auto& style_path : style_paths)
    {
      if (style_path.IsNull()) continue;
      if (gSystem && gSystem->AccessPathName(style_path, kReadPermission)) continue;
      if (gROOT && gROOT->LoadMacro(style_path) == 0)
      {
        gROOT->ProcessLine("SetsPhenixStyle();");
        return (gROOT->GetStyle("sPHENIX") != nullptr);
      }
    }

    printf("TpcLaserDNL_PlotPhaseCDF: could not load sPHENIX style, using fallback settings\n");
    if (gStyle)
    {
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
    }
    return false;
  }

  TString make_sample_label(int layer, int side, int npads)
  {
    TString label;
    if (side == 0)
    {
      label += "south side";
    }
    else if (side == 1)
    {
      label += "north side";
    }
    else if (side >= 0)
    {
      label += TString::Format("side %d", side);
    }

    if (!label.IsNull()) label += ", ";
    label += TString::Format("layer %d", layer);
    if (npads > 0) label += TString::Format(", %d pads cluster", npads);
    return label;
  }

  std::unique_ptr<TH1D> clone_hist(const TH1D* source, const TString& name)
  {
    if (!source) return nullptr;
    std::unique_ptr<TH1D> hist(static_cast<TH1D*>(source->Clone(name)));
    hist->SetDirectory(nullptr);
    hist->SetStats(0);
    return hist;
  }

  std::unique_ptr<TH1D> make_cdf_from_counts(const TH1D* source, const TString& name)
  {
    auto cdf = clone_hist(source, name);
    if (!cdf) return nullptr;

    cdf->Reset("ICES");
    const double total = source->Integral();
    if (total <= 0.0) return cdf;

    double cumulative = 0.0;
    for (int bin = 1; bin <= source->GetNbinsX(); ++bin)
    {
      cumulative += source->GetBinContent(bin);
      cdf->SetBinContent(bin, std::clamp(cumulative / total, 0.0, 1.0));
      cdf->SetBinError(bin, 0.0);
    }
    return cdf;
  }

  std::unique_ptr<TH1D> get_cdf_hist(TFile* fin,
                                     const TString& cdf_path,
                                     const TString& source_path,
                                     const TString& draw_name)
  {
    if (auto* h_cdf = dynamic_cast<TH1D*>(fin->Get(cdf_path)))
    {
      return clone_hist(h_cdf, draw_name);
    }

    if (auto* h_counts = dynamic_cast<TH1D*>(fin->Get(source_path)))
    {
      printf("TpcLaserDNL_PlotPhaseCDF: missing %s; building CDF from %s\n",
             cdf_path.Data(),
             source_path.Data());
      return make_cdf_from_counts(h_counts, draw_name);
    }

    printf("TpcLaserDNL_PlotPhaseCDF: missing %s and %s\n",
           cdf_path.Data(),
           source_path.Data());
    return nullptr;
  }
}  // namespace

void TpcLaserDNL_PlotPhaseCDF(const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/mean1400_sig004_noalign/completed/pion_plus_phasecdf_2pad3pad_corrected.root",
                              const char* outdir = "",
                              int layer = 7,
                              int side = 0,
                              int npads = -1,
                              bool save_root_canvas = true)
{
  apply_sphenix_style();

  std::unique_ptr<TFile> fin(TFile::Open(infile, "READ"));
  if (!fin || fin->IsZombie())
  {
    printf("TpcLaserDNL_PlotPhaseCDF: cannot open input %s\n", infile);
    return;
  }

  TString output_dir = outdir;
  if (output_dir.IsNull()) output_dir = gSystem->DirName(infile);
  if (gSystem) gSystem->mkdir(output_dir, kTRUE);

  std::vector<int> pads_to_plot;
  if (npads > 0)
  {
    pads_to_plot.push_back(npads);
  }
  else
  {
    pads_to_plot.push_back(2);
    pads_to_plot.push_back(3);
  }

  for (const int pads : pads_to_plot)
  {
    const TString tag = TString::Format("layer%02d_side%d_pads%02d", layer, side, pads);
    const TString true_cdf_path = TString::Format("cdf_maps/phase_true_cdf/h_cdf_phase_true_%s", tag.Data());
    const TString reco_cdf_path = TString::Format("cdf_maps/phase_reco_cdf/h_cdf_phase_reco_%s", tag.Data());
    const TString true_source_path = TString::Format("cdf_maps/phase_true/h_phase_true_%s", tag.Data());
    const TString reco_source_path = TString::Format("cdf_maps/phase_reco/h_phase_reco_%s", tag.Data());

    auto h_true = get_cdf_hist(fin.get(),
                               true_cdf_path,
                               true_source_path,
                               TString("h_phase_true_cdf_draw_") + tag);
    auto h_reco = get_cdf_hist(fin.get(),
                               reco_cdf_path,
                               reco_source_path,
                               TString("h_phase_reco_cdf_draw_") + tag);
    if (!h_true || !h_reco) continue;

    h_true->UseCurrentStyle();
    h_reco->UseCurrentStyle();
    h_true->SetTitle(";phase;CDF");
    h_reco->SetTitle(";phase;CDF");
    h_true->SetMinimum(0.0);
    h_true->SetMaximum(1.08);
    h_true->SetLineColor(kBlack);
    h_true->SetMarkerColor(kBlack);
    h_true->SetLineWidth(2);
    h_true->SetMarkerStyle(20);
    h_true->SetMarkerSize(0.0);
    h_reco->SetLineColor(kRed + 1);
    h_reco->SetMarkerColor(kRed + 1);
    h_reco->SetLineWidth(2);
    h_reco->SetMarkerStyle(20);
    h_reco->SetMarkerSize(0.0);

    const TString canvas_name = "c_phase_cdf_" + tag;
    const TString output_base = TString::Format("%s/phase_cdf_%s",
                                                output_dir.Data(),
                                                tag.Data());

    std::unique_ptr<TCanvas> c(new TCanvas(canvas_name, canvas_name, 900, 700));
    c->UseCurrentStyle();
    c->SetTicks(1, 1);
    c->SetTopMargin(0.05);

    h_true->Draw("HIST");
    h_reco->Draw("HIST SAME");

    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    if (gStyle)
    {
      label.SetTextFont(gStyle->GetTextFont());
      label.SetTextSize(gStyle->GetTextSize());
    }
    label.DrawLatex(0.20, 0.88, "#bf{#it{sPHENIX}} Simulation");
    label.DrawLatex(0.20, 0.82, "Single #pi^{+}, p_{T}=10 GeV");
    label.DrawLatex(0.20, 0.76, make_sample_label(layer, side, pads));

    std::unique_ptr<TLegend> legend(new TLegend(0.58, 0.76, 0.88, 0.88));
    if (gStyle)
    {
      legend->SetBorderSize(gStyle->GetLegendBorderSize());
      legend->SetFillColor(gStyle->GetLegendFillColor());
      legend->SetTextFont(gStyle->GetLegendFont());
      legend->SetTextSize(gStyle->GetLegendTextSize());
    }
    legend->SetFillStyle(0);
    legend->AddEntry(h_true.get(), "True phase CDF", "l");
    legend->AddEntry(h_reco.get(), "Reconstructed phase CDF", "l");
    legend->Draw();

    c->SaveAs(output_base + ".pdf");
    c->SaveAs(output_base + ".png");

    if (save_root_canvas)
    {
      std::unique_ptr<TFile> fout(TFile::Open(output_base + ".root", "RECREATE"));
      if (fout && !fout->IsZombie())
      {
        fout->cd();
        c->Write();
        fout->Close();
      }
    }

    printf("TpcLaserDNL_PlotPhaseCDF: wrote %s.[pdf,png,root]\n", output_base.Data());
  }
}
