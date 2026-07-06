// Draw presentation-style dRphi vs phase_reco histograms from a
// TpcLaserDNL_Analyze summary ROOT file.

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

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

    printf("TpcLaserDNL_PlotPhaseReco2D: could not load sPHENIX style, using fallback settings\n");
    if (gStyle)
    {
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
    }
    return false;
  }
}  // namespace

void TpcLaserDNL_PlotPhaseReco2D(const char* infile = "/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/mean1400_sig004_noalign/completed/pion_plus_summary.root",
                                 const char* outdir = "",
                                 int layer = 7,
                                 bool logz = false,
                                 bool save_root_canvas = true,
                                 int npads = -1,
                                 int side = -1,
                                 const char* extra_label = "")
{
  apply_sphenix_style();

  std::unique_ptr<TFile> fin(TFile::Open(infile, "READ"));
  if (!fin || fin->IsZombie())
  {
    printf("TpcLaserDNL_PlotPhaseReco2D: cannot open input %s\n", infile);
    return;
  }

  TString hist_path;
  if (side >= 0 && npads > 0)
  {
    hist_path = TString::Format("by_layer_by_side_by_npads/side%d/dRphi_vs_phase_reco/h2_dRphi_vs_phaseReco_layer%02d_side%d_pads%02d",
                                side,
                                layer,
                                side,
                                npads);
  }
  else if (side >= 0)
  {
    hist_path = TString::Format("by_layer_by_side/side%d/dRphi_vs_phase_reco/h2_dRphi_vs_phaseReco_layer%02d_side%d",
                                side,
                                layer,
                                side);
  }
  else if (npads > 0)
  {
    hist_path = TString::Format("by_layer_by_npads/dRphi_vs_phase_reco/h2_dRphi_vs_phaseReco_layer%02d_pads%02d",
                                layer,
                                npads);
  }
  else
  {
    hist_path = TString::Format("by_layer_2d/dRphi_vs_phase_reco/h2_dRphi_vs_phaseReco_layer%02d",
                                layer);
  }
  auto* h_in = dynamic_cast<TH2D*>(fin->Get(hist_path));
  if (!h_in)
  {
    printf("TpcLaserDNL_PlotPhaseReco2D: missing histogram %s in %s\n", hist_path.Data(), infile);
    return;
  }

  TString output_dir = outdir;
  if (output_dir.IsNull())
  {
    output_dir = gSystem->DirName(infile);
  }
  if (gSystem) gSystem->mkdir(output_dir, kTRUE);

  TString output_tag = TString::Format("layer%02d", layer);
  if (side >= 0) output_tag += TString::Format("_side%d", side);
  if (npads > 0) output_tag += TString::Format("_pads%02d", npads);

  const TString canvas_name = "c_dRphi_vs_phase_reco_" + output_tag;
  const TString output_base = TString::Format("%s/dRphi_vs_phase_reco_%s",
                                              output_dir.Data(),
                                              output_tag.Data());

  std::unique_ptr<TH2D> h(static_cast<TH2D*>(h_in->Clone(TString::Format("%s_draw", h_in->GetName()))));
  h->SetDirectory(nullptr);
  h->UseCurrentStyle();
  h->SetStats(0);
  h->SetTitle(";phase;R#Delta#phi [cm]");
  h->GetZaxis()->SetTitle("Counts");

  std::unique_ptr<TCanvas> c(new TCanvas(canvas_name, canvas_name, 900, 700));
  c->UseCurrentStyle();
  c->SetTicks(1, 1);
  c->SetRightMargin(0.15);
  if (logz) c->SetLogz(1);

  h->Draw("COLZ");

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextFont(42);
  label.SetTextSize(0.040);
  label.DrawLatex(0.20, 0.88, "#bf{#it{sPHENIX}} Simulation");

  TString sample_label = "Single #pi^{+}, p_{T}=10 GeV";
  TString detail_label;
  if (side == 0)
  {
    detail_label = "south side";
  }
  else if (side == 1)
  {
    detail_label = "north side";
  }
  else if (side >= 0)
  {
    detail_label = TString::Format("side %d", side);
  }

  if (!detail_label.IsNull())
  {
    detail_label += TString::Format(", layer %d", layer);
  }
  else
  {
    detail_label = TString::Format("layer %d", layer);
  }

  std::vector<TString> detail_lines;
  detail_lines.push_back(sample_label);
  detail_lines.push_back(detail_label);
  if (npads > 0)
  {
    detail_lines.push_back(TString::Format("%d pads cluster", npads));
  }
  if (extra_label && extra_label[0] != '\0')
  {
    detail_lines.push_back(extra_label);
  }

  label.SetTextSize(0.032);
  double y_label = 0.82;
  for (const auto& line : detail_lines)
  {
    label.DrawLatex(0.20, y_label, line);
    y_label -= 0.055;
  }

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

  printf("TpcLaserDNL_PlotPhaseReco2D: wrote %s.[pdf,png,root]\n", output_base.Data());
}
