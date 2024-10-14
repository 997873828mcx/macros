#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <vector>


void event_display(const int eve=0){

    // Open the ROOT file
    TFile *file = TFile::Open("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_51576_122_resid.root");


    // Check if the file was opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    TFile *field = TFile::Open("magnetic_field_histograms.root");


     // Check if the file was opened successfully
    if (!field || field->IsZombie()) {
        std::cerr << "Error opening magnetic field file" << std::endl;
        file->Close();
        return;
    }

    TH3F *hist_bx = (TH3F *) field->Get("hist_bx");
    TH3F *hist_by = (TH3F *) field->Get("hist_by");
    TH3F *hist_bz = (TH3F *) field->Get("hist_bz");

    if (!hist_bx || !hist_by || !hist_bz) {
        std::cerr << "Error retrieving magnetic field histograms from file" << std::endl;
        file->Close();
        field->Close();
        return;
    }

    // Get the tree from the file
    TTree *tree = (TTree *) file->Get("residualtree");

    // Check if the tree was retrieved successfully
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        field->Close();
        return;
    }



    // Declare variables to hold the data
    std::vector<float> *clusgx = nullptr;
    std::vector<float> *clusgy = nullptr;
    std::vector<float> *clusgz = nullptr;
    float pt, eta, dedx, quality, phi, pcax, pcay, pcaz, phi0, the0, px, py, pz;
    int nhits, charge, ntpc, nintt, nmaps, event;
    std::string m_magField = "FIELDMAP_TRACKING";

    const int numColors = 11;
    int color[numColors] = {kOrange+1, kOrange-1, kRed-1, kRed-7, kPink-6, kMagenta+2, kViolet, kViolet+4, kBlue+2, kAzure+1, kCyan+4};

    // Set branch addresses
    tree->SetBranchAddress("clusgx", &clusgx);
    tree->SetBranchAddress("clusgy", &clusgy);
    tree->SetBranchAddress("clusgz", &clusgz);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("dedx", &dedx);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("nhits", &nhits);
    tree->SetBranchAddress("ntpc", &ntpc);
    tree->SetBranchAddress("nintt", &nintt);
    tree->SetBranchAddress("nmaps", &nmaps);
    tree->SetBranchAddress("pcax", &pcax);
    tree->SetBranchAddress("pcay", &pcay);
    tree->SetBranchAddress("pcaz", &pcaz);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("event", &event);

    Long64_t nentries = tree->GetEntries();



    // Variables to store the data for the track to be processed
    std::vector<double> x, y, z;
    std::vector<double> tx, ty, tz;

    tx.reserve(3000);
    ty.reserve(3000);
    tz.reserve(3000);
    //float tx[5000], ty[5000], tz[5000];
    int nEntriesInAnEvent=0;



    size_t npoints = 0;
    Long64_t entryNumber = -1; // To keep track of the entry number

    double test_x=0;
    double test_y=0;
    double test_z=0;

    int N;
    /*if (n<nentries)
        {
            N=n;
        }
    else{
        N=nentries;
    }*/


    TGraph2D *initial_point = new TGraph2D(1);

    initial_point->SetPoint(0,0,0,0);

    initial_point->SetMarkerStyle(20);
    initial_point->SetMarkerSize(0.5);
    initial_point->SetMarkerColor(kBlack);

    for (Long64_t i=0; i < nentries; i++)
    {
        tree->GetEntry(i);
        if(event == eve)
        {
            nEntriesInAnEvent += 1;
        }
    }
    //nEntriesInAnEvent=10;
    nentries=nEntriesInAnEvent;
    int nEntriesWeWant = 10;
    int nthEntry=0;

    std::vector<TGraph2D*> gr1;
    std::vector<TGraph2D*> gr2;


    //gr2->SetTitle(Form("Entry %lld;X [cm];Y [cm];Z [cm]", entryNumber));
    /*gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(0.5);
    gr2->SetMarkerColor(kRed);*/


    // Loop over entries until the Nth entry to find the first track without hits in INTT and MAPS
    for (Long64_t eEntry = 0; eEntry < nentries; eEntry++) {
        // Get the specified entry (track)
        tree->GetEntry(eEntry);
        if(event == eve && nthEntry < nEntriesWeWant) {

            //for (Long64_t eEntry=0; eEntry < nEntriesInAnEvent; eEntry++) {
            if (nintt == 0 && nmaps == 0) {
                //if (1) {
                TGraph2D* graph1 = new TGraph2D();
                TGraph2D* graph2 = new TGraph2D();
                // Get the number of hit points
                npoints = clusgx->size();
                // Check if there are any hits to plot
                if (npoints == 0) {
                    std::cerr << "No hits to plot for entry " << eEntry << std::endl;
                    continue; // Skip this entry and continue the loop
                }

                tx.clear();
                ty.clear();
                tz.clear();
                x.clear();
                y.clear();
                z.clear();

                tx.push_back(pcax);
                ty.push_back(pcay);
                tz.push_back(pcaz);
                //std::cout << "origin" << "x:" << pcax << "y:" << pcay << "z:" << pcaz << std::endl;




                // Convert std::vector<float> to std::vector<double> for TGraph2D
                x.assign(clusgx->begin(), clusgx->end());
                y.assign(clusgy->begin(), clusgy->end());
                z.assign(clusgz->begin(), clusgz->end());

                // Prepare the track propagation
                //TVector3 track(px, py, pz);
                TVector3 mom(px, py, pz);
                TVector3 dPosition(0, 0, 0);
                // std::cout << "momentum" << "px:" << px << "py:" << py << "pz:" << pz << std::endl;
                double phi = mom.Phi();
                double theta = mom.Theta();


                int nsteps = 10000;

                double p = std::sqrt(px*px+py*py+pz*pz);

                for (int j = 1; j < nsteps; j++) {
                    //double length = j * length_step;
                    double length = 0.01;
                    /*const int xBin = hist_bx->GetXaxis()->FindBin(tx[j-1]);
                    const int yBin = hist_bx->GetYaxis()->FindBin(ty[j-1]);
                    const int zBin = hist_bx->GetZaxis()->FindBin(tz[j-1]);*/

                    mom.SetMagThetaPhi(p, theta, phi);
                    /*dPosition.SetMagThetaPhi(length, theta, phi);
                    tx.push_back(tx[j - 1] + dPosition.X());
                    ty.push_back(ty[j - 1] + dPosition.Y());
                    tz.push_back(tz[j - 1] + dPosition.Z());*/


                    /*const float bx = -1.0 * hist_bx->GetBinContent(xBin, yBin, zBin);
                    const float by = -1.0 * hist_by->GetBinContent(xBin, yBin, zBin);
                    const float bz = -1.0 * hist_bz->GetBinContent(xBin, yBin, zBin);*/
                    const float bx = -1.0 * hist_bx->Interpolate(tx[j-1], ty[j-1], tz[j-1]);
                    const float by = -1.0 * hist_by->Interpolate(tx[j-1], ty[j-1], tz[j-1]);
                    const float bz = -1.0 * hist_bz->Interpolate(tx[j-1], ty[j-1], tz[j-1]);
                    //phi += (charge * bz * 0.3 * 0.01 * length) / (std::sqrt(mom.X() * mom.X() + mom.Y() * mom.Y()));
                    pt = mom.Pt();
                    double dTheta = -charge*length*0.01*0.3*(mom.X()*by-mom.Y()*bx)/p/pt;
                    double dThetaMiddle = dTheta/2.0;
                    double thetaMiddle = theta+dThetaMiddle;
                    theta+=dTheta;
                    //1st expression for dPhi
                    //double dPhi = (charge*length*0.01*0.3*(mom.Z()*bx-mom.X()*bz)-mom.Z()*p*std::sin(phi)*dTheta)/mom.Pt()/std::cos(phi)/p;

                    //2nd expressioni for dPhi
                    double dPhi = (-charge*length*0.01*0.3*(mom.Y()*bz-mom.Z()*by)+mom.Z()*p*std::cos(phi)*dTheta)/mom.Pt()/std::sin(phi)/p;
                    double dPhiMiddle = dPhi/2.0;
                    double phiMiddle = phi+dPhiMiddle;
                    phi+=dPhi;

                    dPosition.SetMagThetaPhi(length, thetaMiddle, phiMiddle);
                    //dPosition.SetMagThetaPhi(length, theta, phi);
                    tx.push_back(tx[j - 1] + dPosition.X());
                    ty.push_back(ty[j - 1] + dPosition.Y());
                    tz.push_back(tz[j - 1] + dPosition.Z());



                }

                //gr2[i] = new TGraph2D(3001, x, y, z);
                for (int k = 0; k < x.size(); k++) {
                    //std::cout<<"position"<<tx[i]<<ty[i]<<tz[i]<<std::endl;
                    graph1->SetPoint(k, x[k], y[k], z[k]);
                }

                for (int k = 0; k < tx.size(); k++) {
                    //std::cout<<"position"<<tx[i]<<ty[i]<<tz[i]<<std::endl;
                    graph2->SetPoint(k, tx[k], ty[k], tz[k]);
                }

              /*  graph1->SetMarkerStyle(20);
                graph1->SetMarkerSize(0.2);
                graph1->SetMarkerColor(kBlue);
                graph2->SetMarkerStyle(20);
                graph2->SetMarkerSize(0.2);
                graph2->SetMarkerColor(kRed);*/

                gr1.push_back(graph1);
                gr2.push_back(graph2);

                entryNumber = eEntry; // Store the entry number
                //break; // Exit the loop after processing the first matching track

                nthEntry++;
                if (nthEntry >= nEntriesWeWant) {
                    break; // Exit the loop after processing the desired number of entries
                }
            }

        }

    }

   /* // Calculate range of the track
    double xmin = *std::min_element(tx.begin(), tx.end());
    double xmax = *std::max_element(tx.begin(), tx.end());
    double ymin = *std::min_element(ty.begin(), ty.end());
    double ymax = *std::max_element(ty.begin(), ty.end());
    double zmin = *std::min_element(tz.begin(), tz.end());
    double zmax = *std::max_element(tz.begin(), tz.end());*/

/*    // Add some padding to the ranges
    double padding = 0.1 * std::max({xmax-xmin, ymax-ymin, zmax-zmin});
    xmin -= padding; xmax += padding;
    ymin -= padding; ymax += padding;
    zmin -= padding; zmax += padding;*/

    // Create TGraph2D for the propagated track


    // Set the ranges before adding points
/*    gr2->GetXaxis()->SetLimits(-10, 10);
    gr2->GetYaxis()->SetLimits(ymin, ymax);
    gr2->GetZaxis()->SetLimits(zmin, zmax);*/

/*    for (int i=0; i<tx.size();i++)
    {
        std::cout<<"position"<<tx[i]<<ty[i]<<tz[i]<<std::endl;
        gr2->SetPoint(i, tx[i], ty[i], tz[i]);
    }*/



    // Add a test point at (0,0,0) to ensure something is plotted
    /*gr2->SetPoint(tx.size(), 0, 0, 0);
    std::cout << "Added test point at (0,0,0)" << std::endl;*/


/*    // Check if a track was found
    if (entryNumber == -1) {
        std::cerr << "No track found without hits in INTT and MAPS" << std::endl;
        file->Close();
        return;
    }*/

    // Create a canvas to draw the graphs
    TCanvas *c = new TCanvas("c", "3D Event Display", 800, 800);
    /*c->SetFillColor(kWhite);
    c->SetFrameFillColor(kWhite);*/
    c->cd();
    TH3F* htemp = new TH3F("htemp"," ", 10, -100, 100, 10, -100, 100, 10, -100 , 100);
    htemp->GetXaxis()->SetTitle("X");
    htemp->GetYaxis()->SetTitle("Y");
    htemp->GetZaxis()->SetTitle("Z");
    htemp->Draw();

    /*// Create TGraph2D for the hits
    TGraph2D *gr1 = new TGraph2D(npoints, &x[0], &y[0], &z[0]);
    gr1->SetTitle(Form("Entry %lld;X [cm];Y [cm];Z [cm]", entryNumber));
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.5);
    gr1->SetMarkerColor(kBlue);*/

    for(size_t i=0; i < gr1.size(); i++){
        if(gr1[i])
            {
            gr1[i]->SetMarkerStyle(20);
                gr1[i]->SetMarkerSize(0.2);
                gr1[i]->SetMarkerColor(color[i % numColors]);
            gr1[i]->Draw("PSAME");
            }
        else
            {
            std::cout<<"the graph does not exist"<<std::endl;
            }
        if(gr2[i])
            {
            gr2[i]->SetMarkerStyle(20);
                gr2[i]->SetMarkerSize(0.2);
                gr2[i]->SetMarkerColor(color[i % numColors]);
            gr2[i]->Draw("PSAME");

            }
    }






/*    // Force 3D view
    TView3D *view = new TView3D();
    view->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
    c->SetView(view);*/


    /*// Draw the first graph
    gr2->Draw("P0");
    gr2->GetHistogram()->GetXaxis()->SetLimits(-60, 60);
    gr2->GetHistogram()->GetYaxis()->SetLimits(-60, 60);
    gr2->GetHistogram()->GetZaxis()->SetRangeUser(-100, 100);

    // Draw the second graph on the same canvas
    //initial_point->Draw("P0 SAME");

    gr1->Draw("P0 SAME");*/
    //gr1->Draw("P0");

/*    // Add a legend
    TLegend *legend = new TLegend(0.75, 0.85, 0.95, 0.95);
    //legend->AddEntry(gr1, "Hits", "p");
    legend->AddEntry(gr2, "Propagated Track", "p");
    legend->Draw();*/

    // Update the canvas
    //c->Update();

    // Optionally, save the canvas to a file
    //c->SaveAs(Form("event_display_%lld.png", entryNumber));

    /*// Clean up
    delete gr1;
    delete gr2;
    delete c;
    file->Close();*/
}





