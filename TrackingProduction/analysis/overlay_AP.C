//  overlay_AP.C --------------------------------------------------------
//  root -l overlay_AP.C

void overlay_AP()
{
  //--------------------------------------------------------------------
  // 1) open both files and fetch the histograms
  //--------------------------------------------------------------------
  TFile fData("../kshort_analysis.root");        // the file you create at the end of kshort_analysis()
  TFile fMC  ("Ks0_AP_sim.root");            // written by the generator macro

  if (!fData.IsOpen() || !fMC.IsOpen()) {
    std::cerr << "Could not open one of the input ROOT files!\n";
    return;
  }

  TH2D *hData = (TH2D*) fData.Get("h_ArmPod_WithCut");   
  TH2D *hMC   = (TH2D*) fMC  .Get("hArm_WithCut");       // <-- MC histo saved in Ks0_AP_sim.root

  if (!hData || !hMC) {
    std::cerr << "Histogram(s) not found!\n";
    return;
  }

  //--------------------------------------------------------------------
  // 2) optional: normalise MC so that total number of counts matches
  //--------------------------------------------------------------------
  const double scale = hData->Integral() / hMC->Integral();
  hMC->Scale(scale);                         

  //--------------------------------------------------------------------
  // 3) drawing
  //--------------------------------------------------------------------
  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c","A-P data vs. MC",900,700);
  c->SetLogz();                              // log Z helps a lot with wide dynamic range

  // 3a) draw the data as a coloured density map
  hData->SetTitle("Armenteros-Podolanski   (data : colour,  MC : contours)");
  hData->GetZaxis()->SetTitleOffset(1.25);
  hData->Draw("COLZ");                       // coloured pixels for the real distribution

  // 3b) overlay MC – two alternative styles
  //  ---- contour lines (nice for log‐Z)
  hMC->SetLineColor(kRed+1);
  hMC->SetLineWidth(2);
  hMC->SetContour(10);                       // number of contour levels
  hMC->Draw("CONT3 SAME");                  // “CONT3” = smooth closed contours

  /*  // ---- OR scatter points / markers  ---------------------------
  hMC->SetMarkerStyle(24);
  hMC->SetMarkerColor(kAzure+1);
  hMC->SetMarkerSize(0.4);
  hMC->Draw("SAME P");
  */

  //--------------------------------------------------------------------
  // 4) legend
  //--------------------------------------------------------------------
  auto leg = new TLegend(0.15,0.80,0.45,0.90);
  leg->AddEntry(hData,"real data","f");
  leg->AddEntry(hMC ,"MC (scaled)","l");
  leg->Draw();

  c->SaveAs("AP_data_vs_MC.pdf");
}