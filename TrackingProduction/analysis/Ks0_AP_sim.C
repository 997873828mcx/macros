/*  Ks0_AP_sim.C  -----------------------------------------------------------
 *  Build K0s momentum pdf from STAR invariant yield
 *  shoot K0s, decay to π+π-, make Armenteros–Podolanski plot
 *  ------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <vector>
#include "TMath.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

void Ks0_AP_sim()
{
  /* ----------------------------------------------------------------- *
   * 1.  read CSV table                                                 *
   * ----------------------------------------------------------------- */
  std::ifstream fin("HEPData-ins722757-v1-Figure_8_left_up.csv");  
  if (!fin)
  {
    std::cerr << "CSV not found!\n";
    return;
  }

  std::vector<double> vPt, vYield, vErr;
  std::string line;
  while (std::getline(fin, line))
  {
    if (line.size() == 0 || line[0] == '#') continue;
    std::stringstream ss(line);
    double pt, y, eplus, eminus;
    char comma;
    ss >> pt >> comma >> y >> comma >> eplus >> comma >> eminus;
    vPt.push_back(pt);
    vYield.push_back(y);
    vErr.push_back(eplus);  // symmetrised already
  }
  fin.close();

  const int N = vPt.size();
  TGraphErrors gInv(N);
  for (int i = 0; i < N; ++i)
  {
    gInv.SetPoint(i, vPt[i], vYield[i]);
    gInv.SetPointError(i, 0., vErr[i]);
  }
  gInv.SetName("gInv");
  gInv.SetTitle("K^{0}_{S} Invariant Yield; p_{T} [GeV/c]; 1/(2#pi p_{T}) d^{2}N/dp_{T}dy");
  

  /* ----------------------------------------------------------------- *
   * 2.  Lévy–Tsallis fit  (eq. 11 in STAR paper)                       *
   * ----------------------------------------------------------------- */
  const double massKs = 0.497611;  // GeV
  TF1* levy = new TF1("levy",
    Form(
      "[0]*( ([2]-1.)*([2]-2.) )"
      " / ( 2.*TMath::Pi()*[2]*[1]*( [2]*[1] + %g*([2]-2.) ) )"
      " * pow( 1. + ( sqrt(x*x + %g ) - %g ) / ([2]*[1]) , -[2] )"
      , massKs                      // replaces first %g  → m0 in front‐factor
      , massKs*massKs                  // replaces second %g → under the sqrt: m0^2
      , massKs                      // replaces third %g  → the “- m0” inside
    ),
    0.0, 5.0
  );

  levy->SetParNames("dN/dy", "C", "n");
  levy->SetParameters(0.13, 0.2, 8.0);
  gInv.Fit(levy, "R0Q");  // R=use range, 0=quiet, Q=quiet fit print

  TCanvas* cYield = new TCanvas("cYield","K^{0}_{S} Invariant Yield",800,600);

  // 2) Draw the data points with errors:
  gInv.SetMarkerStyle(20);
  gInv.SetMarkerSize(1.0);
  gInv.SetLineColor(kBlack);
  gInv.Draw("AP");  // A=draw axes, P=draw points with error‐bars

  // 3) Overlay the fit curve:
  levy->SetLineColor(kBlack);
  levy->SetLineWidth(2);
  levy->Draw("same");  // redraw TF1 on top of the graph

  // 4) Add a legend:
  auto leg = new TLegend(0.60, 0.65, 0.88, 0.88);
  leg->AddEntry(&gInv, "STAR data", "p");
  leg->AddEntry(levy, "Levy-Tsallis fit", "l");
  leg->Draw();

  // 5) Save to file:
  cYield->SaveAs("Figure8_Ks0_Invariant.pdf");

  // 3) Build the normalized pT pdf: f(pt) = 2π pt I(pt)


   TH1D hPtPdf("hPtPdf","K^{0}_{S} p_{T} pdf; p_{T} [GeV/c]; prob. density",
    100, 0.0, 5.0);
    for (int bin=1; bin<=hPtPdf.GetNbinsX(); ++bin) {
      double pt   = hPtPdf.GetBinCenter(bin);
      double dpt  = hPtPdf.GetBinWidth(bin);
      double dens = levy->Eval(pt);          // I(p_T) itself
      hPtPdf.AddBinContent(bin, dens * dpt); // area piece
  }
  hPtPdf.Scale( 1.0 / hPtPdf.Integral("width") );   // normalise to 1

  /* ----------------------------------------------------------------- *
   * 4.  Monte‑Carlo generation                                         *
   * ----------------------------------------------------------------- */

  const int Nevents = 1e6;
  
  const double massPi = 0.13957;  // PDG pion mass ~0.13957 GeV
  // 2-body formula: pCM = sqrt( [M^2 - (m1+m2)^2]*[M^2 - (m1-m2)^2 ] ) / (2M )
  double M2 = massKs * massKs;
  double sumM = (massPi + massPi);
  double sumM2 = sumM * sumM; // (2*masspi)^2
  double difM = (massPi - massPi);
  double difM2 = difM * difM; // (0)^2
  double pCM = TMath::Sqrt((M2 - sumM2) * (M2 - difM2)) / (2. * massKs);

  TRandom3 R(0);
  TH2D hArm("hArm", "Armenteros-Podolanski; #alpha; p_{T}^{+}  [GeV/c]",
            200, -2, 2, 200, 0, 0.5);//set bin 200?

  TLorentzVector vKs, vPiPlus, vPiMinus;

  for (int ie = 0; ie < Nevents; ++ie)
  {
    /* --- sample parent momentum magnitude and direction -------- */

    double pt   = hPtPdf.GetRandom(&R);        
    double  y   = -0.5 + R.Rndm();             
    double  phi = 2.*TMath::Pi()*R.Rndm();            
    
    double pz   =  pt * sinh(y);               
    double  p   =  pt * cosh(y);               
    TVector3 pVec( pt*cos(phi), pt*sin(phi), pz );
    
    double eK = TMath::Sqrt(p*p + massKs*massKs);
    TLorentzVector vKs;  vKs.SetVectM( pVec, massKs );

    /* --- decay in rest frame  isotropically -------------------- */
    double cosd = 2. * R.Rndm() - 1.;
    double sind = TMath::Sqrt(1. - cosd * cosd);
    double phid = 2. * TMath::Pi() * R.Rndm();
    // π+ momentum in CM
    TVector3 pStar(
      pCM*sind*TMath::Cos(phid),
      pCM*sind*TMath::Sin(phid),
      pCM*cosd
    );
    TLorentzVector piStarPlus (pStar, TMath::Sqrt(pCM*pCM + massPi*massPi));
    TLorentzVector piStarMinus(-pStar,  TMath::Sqrt(pCM*pCM + massPi*massPi));

    /* --- boost daughters to lab -------------------------------- */
    TVector3 beta = vKs.BoostVector();
    vPiPlus = piStarPlus;
    vPiPlus.Boost(beta);
    vPiMinus = piStarMinus;
    vPiMinus.Boost(beta);

    /* --- Armenteros variables ---------------------------------- */
    double pLplus = vPiPlus.Vect().Dot(pVec.Unit());
    double pLminus = vPiMinus.Vect().Dot(pVec.Unit());
    double alpha = (pLplus - pLminus) / (pLplus + pLminus);

    // transverse momentum of the positive daughter w.r.t. parent
    TVector3 pParUnit = pVec.Unit();
    TVector3 pTplus = vPiPlus.Vect() - pParUnit * pLplus;
    double pTpos = pTplus.Mag();

    hArm.Fill(alpha, pTpos);
  }

  /* ----------------------------------------------------------------- *
   * 5.  save everything                                               *
   * ----------------------------------------------------------------- */
  TFile fout("Ks0_AP_sim.root", "RECREATE");
  gInv.Write();
  levy->Write();
  hPtPdf.Write();
  hArm.Write();
  fout.Close();

  std::cout << "Finished.  Results in Ks0_AP_sim.root\n";
}
