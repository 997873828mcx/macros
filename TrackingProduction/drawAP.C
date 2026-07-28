// =============================================================================
// drawAP.C
// -----------------------------------------------------------------------------
//   root -l -b -q 'drawAP.C("yourFile.root")'
//   root -l -b -q 'drawAP.C("yourFile.root",false)'  // no plotting cuts
// =============================================================================
void drawAP(const char* infile = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/myKShortReco/ap_79516.root",
            bool applyCuts = true,
            const char* outtag = "AP_plot",
            const char* customCut = "")
{
  //gROOT->LoadMacro("/sphenix/user/dcxchenxi/develope/macros/TrackingProduction/sPhenixStyle.C");
  //gROOT->ProcessLine("SetsPhenixStyle();");

  // ----------------------------------------------
  // open file and fetch the tree
  // ----------------------------------------------

  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "[drawAP] cannot open " << infile << std::endl;
    return;
  }

  TTree* t = nullptr;
  f->GetObject("pairTree", t);  // <-- alpha and qT live here
  if (!t)
  {
    std::cerr << "[drawAP] pairTree not found in " << infile << std::endl;
    return;
  }

  // Default selection for the current 79513_Skip production. Campaign-level
  // cuts are already applied before pairTree is written, but repeating them
  // here makes this macro safe for a looser input. The opposite-sign cut is
  // needed because the production intentionally stores same-sign pairs too.
  const char* requiredBranches[] = {
      "alpha", "qT", "cross1", "cross2", "charge1", "charge2",
      "px1", "py1", "px2", "py2", "pca_x", "pca_y", "pca_z",
      "pca1_z", "pca2_z", "pairDCA", "cosThetaReco", "npoints1", "npoints2"};
  if (applyCuts && (!customCut || std::string(customCut).empty()))
  {
    for (const char* branch : requiredBranches)
    {
      if (!t->GetBranch(branch))
      {
        std::cerr << "[drawAP] required branch missing: " << branch << std::endl;
        return;
      }
    }
  }

  TCut selectedCuts =
      //"cross1 == cross2"
      " charge1*charge2 < 0"
      " && abs(pca_z) < 20.0"
      //" && abs(pca1_z-pca2_z) < 1.0"
      " && sqrt(px1*px1+py1*py1) > 0.2"
      " && sqrt(px2*px2+py2*py2) > 0.2"
      " && sqrt(pca_x*pca_x+pca_y*pca_y) > 2.0"
      //" && abs(alpha) < 0.90"
      " && abs(pairDCA) < 0.5"
      " && cosThetaReco > 0.70"
      " && npoints1 > 20"
      " && npoints2 > 20";
  TCut allCuts = selectedCuts;
  if (customCut && std::string(customCut).size() > 0)
  {
    allCuts = TCut(customCut);
  }
  if (!applyCuts)
  {
    allCuts = "";
  }

  std::cout << "[drawAP] plotting with ";
  if (!applyCuts)
  {
    std::cout << "no plotting cuts";
  }
  else if (customCut && std::string(customCut).size() > 0)
  {
    std::cout << "custom cuts: " << customCut;
  }
  else
  {
    std::cout << "analysis cuts";
  }
  std::cout << std::endl;

  // ----------------------------------------------
  // make the AP histogram and fill it
  // ----------------------------------------------
  TCanvas* c1 = new TCanvas("c1", "AP", 1100, 800);
  c1->cd();

  TH2F* hAP = new TH2F("hAP",
                       "Armenteros-Podolanski; #alpha; q_{T}^{+}  (GeV/c)",
                       200, -1, 1,       // α-axis
                       200, 0.0, 0.25);  // qT-axis

  t->Draw("qT:alpha>>hAP", allCuts, "colz");

  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  TString pdfName = Form("%s_%s.pdf", outtag, applyCuts ? "withCuts" : "noCuts");
  c1->SaveAs(pdfName);

  std::cout << "[drawAP] histogram filled with " << hAP->GetEntries()
            << " entries" << std::endl;

  TString rootName = Form("%s_%s.root", outtag, applyCuts ? "withCuts" : "noCuts");
  TFile fout(rootName, "RECREATE");
  hAP->SetName(applyCuts ? "h_ArmPod_WithCut" : "h_ArmPod_NoCut");
  hAP->Write();
  fout.Close();
}
