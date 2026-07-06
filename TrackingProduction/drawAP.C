// =============================================================================
// drawAP.C
// -----------------------------------------------------------------------------
//   root -l drawAP.C\(\"yourFile.root\"\)
//   root -l drawAP.C\(\"yourFile.root\",false\)  // no plotting cuts
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

  // ----------------------------------------------
  // define the cuts  (variable names match pairTree)
  // ----------------------------------------------
  /* TCut crosscut    = "(cross1 == cross2)";
  TCut ptcut       = "sqrt(px1*px1 + py1*py1) > 0.1 && "
                     "sqrt(px2*px2 + py2*py2) > 0.1";

  // - Lproj   = |projected_pathlength⃗|    (branch you stored)
  // - pairDCA = projected_pair_dca         (branch you stored)
  TCut paircut     = "cosThetaReco > 0.94 && "

                     "abs(pairDCA) < 0.15";

  TCut chargecut   = "charge1 != charge2";

  // dca_xy?   = |dca3dxy?| that you wrote to the tree
  TCut trackdcacut = "dca_xy1 > 0.01 && dca_xy2 > 0.01";

  //TCut allCuts = crosscut && ptcut && paircut && chargecut && trackdcacut;
  TCut allCuts =  paircut && ptcut && crosscut && chargecut; */



  // - Lproj   = |projected_pathlength⃗|    (branch you stored)
  // - pairDCA = projected_pair_dca         (branch you stored)
  /*   TCut paircut     = "cosThetaReco > 0.94 && "
                       "Lproj > 0.20         && "
                       "abs(pairDCA) < 0.15"; */

  TCut crosscut = "(cross1 == cross2)";
  

  TCut DIRAcut = "cosThetaReco > 0.94";
  
  TCut pair_cut = "abs(pairDCA) < 0.15";
  TCut proj_cut = "Lproj > 0.2";
  TCut chargecut = "charge1 != charge2";

  // dca_xy?   = |dca3dxy?| that you wrote to the tree
  TCut trackdcacut = "dca_xy1 > 0.03 && dca_xy2 > 0.03";
  TCut ptcut =
      "sqrt(px1*px1 + py1*py1) > 0.2 && "
      "sqrt(px2*px2 + py2*py2) > 0.2";

  TCut selectedCuts = crosscut && DIRAcut && pair_cut && proj_cut && chargecut && trackdcacut && ptcut;
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
