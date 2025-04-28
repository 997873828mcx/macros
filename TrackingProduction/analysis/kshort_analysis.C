#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>

//______________________________________________________________________________
void kshort_analysis()
{
  // 1. Create a TChain for the "DecayTree" in multiple ROOT files
  TChain chain("DecayTree");
  chain.Add("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/myKShortReco/*.root");

  // 2. Create a TTreeReader to read the TChain
  TTreeReader reader(&chain);

  // 3. Define TTreeReaderValue objects for each branch used

  // track 1
  TTreeReaderValue<Int_t> track_1_bunch_crossing(reader, "track_1_bunch_crossing");
  TTreeReaderValue<Float_t> track_1_px(reader, "track_1_px");
  TTreeReaderValue<Float_t> track_1_py(reader, "track_1_py");
  TTreeReaderValue<Float_t> track_1_pz(reader, "track_1_pz");
  TTreeReaderValue<Float_t> track_1_pE(reader, "track_1_pE");
  TTreeReaderValue<Float_t> track_1_pT(reader, "track_1_pT");
  TTreeReaderValue<Float_t> track_1_ip_xy(reader, "track_1_IP_xy");
  // Fix this if your data has a proper branch for charge, e.g. "track_1_charge"
  TTreeReaderValue<Char_t> track_1_charge(reader, "track_1_charge");
  TTreeReaderValue<UInt_t> track_1_mvtx_nhits(reader, "track_1_MVTX_nHits");

  // track 2
  TTreeReaderValue<Int_t> track_2_bunch_crossing(reader, "track_2_bunch_crossing");
  TTreeReaderValue<Float_t> track_2_px(reader, "track_2_px");
  TTreeReaderValue<Float_t> track_2_py(reader, "track_2_py");
  TTreeReaderValue<Float_t> track_2_pz(reader, "track_2_pz");
  TTreeReaderValue<Float_t> track_2_pE(reader, "track_2_pE");
  TTreeReaderValue<Float_t> track_2_pT(reader, "track_2_pT");
  TTreeReaderValue<Float_t> track_2_ip_xy(reader, "track_2_IP_xy");
  // Likewise, ensure this points to the real charge branch
  TTreeReaderValue<Char_t> track_2_charge(reader, "track_2_charge");
  TTreeReaderValue<UInt_t> track_2_mvtx_nhits(reader, "track_2_MVTX_nHits");

  // V0 / Kshort-related

  TTreeReaderValue<Float_t> kshort_px(reader, "K_S0_px"); // or "k_S0_m", etc.
  TTreeReaderValue<Float_t> kshort_py(reader, "K_S0_py");
  TTreeReaderValue<Float_t> kshort_pz(reader, "K_S0_pz");
  TTreeReaderValue<Float_t> kshort_pE(reader, "K_S0_pE");
  TTreeReaderValue<Float_t> kshort_pT(reader, "K_S0_pT");
  TTreeReaderValue<Float_t> kshort_mass(reader, "K_S0_mass"); // or "k_S0_m", etc.
  TTreeReaderValue<Float_t> kshort_decay_length(reader, "K_S0_decayLength");
  TTreeReaderValue<Float_t> kshort_dira(reader, "K_S0_DIRA");
  TTreeReaderValue<Float_t> track_1_track_2_DCA_xy(reader, "track_1_track_2_DCA_xy");

  // 4. Define histograms
  // Adjust binning and ranges as desired
  TH1D *h_kshort_mass = new TH1D("h_kshort_mass", "K^{0}_{S} Mass;M(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts", 100, 0.48, 0.51);
  TH1D *h_kshort_dlength = new TH1D("h_kshort_dlength", "K^{0}_{S} Decay Length;L_{xy} [cm];Counts", 100, 0.0, 30.0);
  TH1D *h_kshort_dira = new TH1D("h_kshort_dira", "K^{0}_{S} DIRA;cos#theta;Counts", 100, 0.9, 1.0);
  TH1D *h_pairDCA = new TH1D("h_pairDCA", "Track Pair DCA_{xy};DCA_{xy} [cm];Counts", 100, 0.0, 0.3);

  // (Optional) track pT distributions for reference
  TH1D *h_track1_pT = new TH1D("h_track1_pT", "Track1 p_{T};p_{T} [GeV/c];Counts", 100, 0.0, 5.0);
  TH1D *h_track2_pT = new TH1D("h_track2_pT", "Track2 p_{T};p_{T} [GeV/c];Counts", 100, 0.0, 5.0);

  TH1D *h_kshort_mass_PE =
      new TH1D("h_kshort_mass_PE",
               "K^{0}_{S} mass from 4-vector in tree;"
               "M_{4vec}(K^{0}_{S})  [GeV/c^{2}];Counts",
               120, 0.4, 0.6);
  TH2D *h_ArmPod_NoCut = new TH2D("h_ArmPod_NoCut",
                                  "Armenteros-Podolanski (Before Cut);#alpha = (p^{+}_{L} - p^{-}_{L}) / (p^{+}_{L} + p^{-}_{L});p_{T}^{+} (GeV/c)",
                                  200, -2.0, 2.0, 200, 0.0, 0.5);

  TH2D *h_ArmPod_WithCut = new TH2D("h_ArmPod_WithCut",
                                    "Armenteros-Podolanski (After Cut);#alpha = (p^{+}_{L} - p^{-}_{L}) / (p^{+}_{L} + p^{-}_{L});p_{T}^{+} (GeV/c)",
                                    200, -2.0, 2.0, 200, 0.0, 0.5);

  TH1D *h_cosTheta = new TH1D("h_cosTheta", "cos#theta in K_{S}^{0} rest frame; cos#theta; Counts", 100, -1.0, 1.0);

  // 5. Main event loop
  Long64_t nEntries = chain.GetEntries();
  std::cout << "Total entries: " << nEntries << std::endl;

  const double massKs = 0.497611; // PDG Kshort mass ~0.4976 GeV
  const double massPi = 0.13957;  // PDG pion mass ~0.13957 GeV
  // 2-body formula: pCM = sqrt( [M^2 - (m1+m2)^2]*[M^2 - (m1-m2)^2 ] ) / (2M )
  double M2 = massKs * massKs;
  double sumM = (massPi + massPi);
  double sumM2 = sumM * sumM; // (2*mPi)^2
  double difM = (massPi - massPi);
  double difM2 = difM * difM; // (0)^2
  double pCM = TMath::Sqrt((M2 - sumM2) * (M2 - difM2)) / (2. * massKs);
  // ~0.206 GeV/c for Kshort->pi+ pi-

  for (Long64_t i = 0; i < nEntries; i++)
  {
    // Read next entry
    reader.Next();

    // Retrieve track momenta, etc.
    Float_t t1pT = *track_1_pT;
    Float_t t2pT = *track_2_pT;

    Float_t t1px = *track_1_px;
    Float_t t1py = *track_1_py;
    Float_t t1pz = *track_1_pz;
    Float_t t2px = *track_2_px;
    Float_t t2py = *track_2_py;
    Float_t t2pz = *track_2_pz;

    Char_t t1q = *track_1_charge;
    Char_t t2q = *track_2_charge;

    Float_t ksPx = *kshort_px;
    Float_t ksPy = *kshort_py;
    Float_t ksPz = *kshort_pz;
    // Float_t ksE = *kshort_pE;
    Float_t ksMass = *kshort_mass;
    Float_t ksE = TMath::Sqrt(ksPx * ksPx + ksPy * ksPy + ksPz * ksPz + massKs * massKs);

    // double M2_meas = ksMass * ksMass;
    // double pCM_meas = TMath::Sqrt((M2_meas - sumM2) * (M2_meas - difM2)) / (2. * ksMass);

    TLorentzVector vKs(*kshort_px, *kshort_py, *kshort_pz, *kshort_pE);
    TLorentzVector vPi1(*track_1_px, *track_1_py, *track_1_pz, *track_1_pE);
    TLorentzVector vPi2(*track_2_px, *track_2_py, *track_2_pz, *track_2_pE);
    TVector3 ksDirLab = vKs.Vect().Unit(); // unit direction in lab
    TLorentzVector dirKsLab(ksDirLab.x(), ksDirLab.y(), ksDirLab.z(), 0.0);
    TVector3 b = -(vKs.BoostVector());

    TLorentzVector vPi1_rest = vPi1; // make copies so we don't overwrite
    TLorentzVector vPi2_rest = vPi2;
    TLorentzVector dirKs_rest = dirKsLab;

    vPi1_rest.Boost(b);
    vPi2_rest.Boost(b);
    dirKs_rest.Boost(b);

    TLorentzVector vPiPlus_rest;
    if (t1q > 0)
      vPiPlus_rest = vPi1_rest; // track1 is +
    else
      vPiPlus_rest = vPi2_rest; // track2 is +

    double dp = vPiPlus_rest.Vect().Dot(dirKs_rest.Vect());
    double magPi = vPiPlus_rest.Vect().Mag();
    double magDir = dirKs_rest.Vect().Mag();
    // double cosTheta = (magPi > 0. && magDir > 0.) ? (dp / (magPi * magDir)) : 0.0;

    // ---- Example selection cuts ----

    // 1) same bunch crossing
    if (*track_1_bunch_crossing != *track_2_bunch_crossing)
      continue;

    // 2) minimal pT on each track
    if (t1pT < 0.2 || t2pT < 0.2)
      continue;

    // 3) opposite charge for pi+ pi-
    if (t1q == t2q)
      continue;

    // 4) pointing angle cut (DIRA)
    if (*kshort_dira < 0.94)
      continue;

    // 5) minimal decay length in xy-plane
    // if (*kshort_decay_length < 0.2) continue;

    // 6) pair DCA cut
    // if (TMath::Abs(*track_1_track_2_DCA_xy) > 0.04) continue;

    // 7) each track's IP_xy > 0.03 ? (example)
    // if (*track_1_ip_xy < 0.03 || *track_2_ip_xy < 0.03) continue;

    // if (!((*track_1_mvtx_nhits > 0) && (*track_2_mvtx_nhits > 0))) continue;

    // 1) Magnitude of Kshort momentum
    double ksMag = TMath::Sqrt(ksPx * ksPx + ksPy * ksPy + ksPz * ksPz);

    // invariant mass from the candidate 4–vector
    double m2_PE = ksE * ksE - ksMag * ksMag; // E^2 - |p|^2
    if (m2_PE > 0)
    {
      double mInv_PE = std::sqrt(m2_PE);
      h_kshort_mass_PE->Fill(mInv_PE);
    }

    // 2) Longitudinal components (pL) for each track along Kshort flight direction
    double p1DotKs = (t1px * ksPx + t1py * ksPy + t1pz * ksPz);
    double pL1 = (ksMag > 0.) ? (p1DotKs / ksMag) : 0.;
    double p2DotKs = (t2px * ksPx + t2py * ksPy + t2pz * ksPz);
    double pL2 = (ksMag > 0.) ? (p2DotKs / ksMag) : 0.;

    double pL_plus = 0.;
    double pL_minus = 0.;

    if (t1q > 0)
    {
      pL_plus = pL1;
      pL_minus = pL2;
    }
    else
    {
      pL_plus = pL2;
      pL_minus = pL1;
    }

    // pL_plus = TMath::Abs(pL_plus);
    // pL_minus = TMath::Abs(pL_minus);

    double alpha = 0.;
    double denom = (pL_plus + pL_minus);
    if (TMath::Abs(denom) > 1e-10)
    {
      alpha = (pL_plus - pL_minus) / denom;
    }

    // if (alpha < -1.0 || alpha > 1.0) continue;

    double pT_pos = 0.;
    double pL_plus_vecx = (ksMag > 0.) ? (pL_plus / ksMag) * ksPx : 0.;
    double pL_plus_vecy = (ksMag > 0.) ? (pL_plus / ksMag) * ksPy : 0.;
    double pL_plus_vecz = (ksMag > 0.) ? (pL_plus / ksMag) * ksPz : 0.;

    double px_perp, py_perp, pz_perp;
    if (t1q > 0)
    {
      px_perp = t1px - pL_plus_vecx;
      py_perp = t1py - pL_plus_vecy;
      pz_perp = t1pz - pL_plus_vecz;
      pT_pos = TMath::Sqrt(px_perp * px_perp + py_perp * py_perp + pz_perp * pz_perp);
    }
    else
    {
      px_perp = t2px - pL_plus_vecx;
      py_perp = t2py - pL_plus_vecy;
      pz_perp = t2pz - pL_plus_vecz;
      pT_pos = TMath::Sqrt(px_perp * px_perp + py_perp * py_perp + pz_perp * pz_perp);
    }

    double pT_meas = TMath::Sqrt(px_perp * px_perp + py_perp * py_perp + pz_perp * pz_perp);

    double beta = 0.;
    if (ksE > 0.)
      beta = ksMag / ksE; // p/E
    if (beta < 1e-9)
      continue; // extremely slow Kshort => rare or unphysical
    double rAlpha = (2. * pCM) / (beta * massKs);

    double factor = (beta * massKs) / (2. * pCM);
    double cosTheta = factor * alpha;

    /*    double beta = 0.;
       if (ksE > 0.) beta = ksMag / ksE;  // p/E
       if (beta < 1e-9) continue;         // extremely slow Kshort => rare or unphysical
       double rAlpha = (2. * pCM_meas) / (beta * ksMass);
       if (alpha < -rAlpha || alpha > rAlpha) continue;
       double factor = (beta * ksMass) / (2. * pCM_meas);
       double cosTheta = factor * alpha; */

    //  Guard against rounding if cosTheta^2>1 => sinTheta^2<0
    double cos2 = cosTheta * cosTheta;
    if (cos2 > 1.)
      cos2 = 1.;
    double sinTheta = TMath::Sqrt(1. - cos2);

    double pT_exp = pCM * sinTheta; // "expected" transverse momentum
    // double pT_exp = pCM_meas * sinTheta;
    double pT2_meas = pT_meas * pT_meas;
    double pT2_exp = pT_exp * pT_exp;

    h_ArmPod_NoCut->Fill(alpha, pT_pos);

    if (alpha < -rAlpha || alpha > rAlpha)
      continue;
    if (pT2_meas < 0.95 * pT2_exp)
      continue;
    if (pT2_meas > 1.05 * pT2_exp)
      continue;

    h_ArmPod_WithCut->Fill(alpha, pT_pos);

    //    ---- Fill histograms ----
    h_kshort_mass->Fill(*kshort_mass);
    h_kshort_dlength->Fill(*kshort_decay_length);
    h_kshort_dira->Fill(*kshort_dira);
    h_pairDCA->Fill(*track_1_track_2_DCA_xy);

    h_track1_pT->Fill(t1pT);
    h_track2_pT->Fill(t2pT);

   
  }

  TCanvas *c1 = new TCanvas("c1", "Kshort Analysis", 1400, 1000);
  c1->Divide(3, 3);

  c1->cd(1);
  h_kshort_mass->Draw();

  c1->cd(2);
  h_kshort_dlength->Draw();

  c1->cd(3);
  h_kshort_dira->Draw();

  c1->cd(4);
  h_pairDCA->Draw();

  c1->cd(5);
  h_track1_pT->Draw();

  c1->cd(6);
  h_track2_pT->Draw();
  c1->cd(7);
  h_kshort_mass_PE->Draw();

  c1->SaveAs("kshort_analysis_plots.pdf");

  /* TCanvas* c2 = new TCanvas("c2", "Armenteros-Podolanski", 800, 600);
  c2->Divide(1, 2);

  c2->cd(1);
  h_ArmPod_NoCut->Draw("COLZ");
  gPad->SetLogz();  // Optional: log scale helps see the distribution better

  c2->cd(2);
  // Draw the uncut data first
  h_ArmPod_NoCut->Draw("COLZ");
  // Overlay the cut data with different color and transparency
  h_ArmPod_WithCut->SetMarkerColor(kRed);
  h_ArmPod_WithCut->Draw("SAME");

  c2->SaveAs("ArmenterosPodolanski_Comparison.pdf"); */

  // Alternative approach: Draw both on same plot with different colors
  TCanvas *c3 = new TCanvas("c3", "Armenteros-Podolanski Overlay", 800, 600);
  c3->SetLogz(); // Optional: log scale helps see the distribution better
  h_ArmPod_NoCut->Draw("COLZ");
  /* h_ArmPod_WithCut->SetMarkerStyle(24); // Open circle marker style (less obtrusive)
  h_ArmPod_WithCut->SetMarkerColor(kRed);
  h_ArmPod_WithCut->SetMarkerSize(0.5);
  Color_t transparentRed = TColor::GetColorTransparent(kRed, 0.3);
  h_ArmPod_WithCut->SetMarkerColor(transparentRed);
  h_ArmPod_WithCut->Draw("SAME P"); */
  TFile fout("kshort_analysis.root","RECREATE");   

  h_kshort_mass        ->Write();
  h_kshort_dlength     ->Write();
  h_kshort_dira        ->Write();
  h_pairDCA            ->Write();

  h_track1_pT          ->Write();
  h_track2_pT          ->Write();
  h_kshort_mass_PE     ->Write();

  h_ArmPod_NoCut       ->Write();   // <-- needed by overlay
  h_ArmPod_WithCut     ->Write();   // <--           "

  fout.Close();                     // flush to disk
  std::cout << "Saved histograms in  kshort_analysis.root\n";


  c3->SaveAs("ArmenterosPodolanski_Overlay.pdf");
}