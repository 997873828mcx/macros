/*
#include <fstream>


void event_display(int n){

    TFile *file = TFile::Open("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_51576_0_resid.root");


    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    TTree *tree = (TTree *) file->Get("residualtree");

    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }

    std::vector<float> *clusgx = nullptr;
    //std::vector<float> *stategx = nullptr;
    std::vector<float> *clusgy = nullptr;
    //std::vector<float> *stategy = nullptr;
    std::vector<float> *clusgz = nullptr;
    float pt, eta, dedx, quality, phi, pcax, pcay, pcaz, phi0, the0, px, py, pz;
    int nhits, charge, ntpc, nintt, nmaps;
    float tx[5000], ty[5000], tz[5000];

    TVector3 track;

    tree->SetBranchAddress("clusgx", &clusgx);
    //tree->SetBranchAddress("stategx", &stategx);
    tree->SetBranchAddress("clusgy", &clusgy);
    //tree->SetBranchAddress("stategy", &stategy);
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

    Long64_t nentries = tree->GetEntries();




    if (n < 0 || n >= nentries) {
        std::cerr << "Event number out of range" << std::endl;
        file->Close();
        return;
    }

    tree->GetEntry(n);
    size_t npoints = clusgx->size();
    //if (npoints == 0)

    std::vector<double> x(clusgx->begin(), clusgx->end());
    std::vector<double> y(clusgy->begin(), clusgy->end());
    std::vector<double> z(clusgz->begin(), clusgz->end());

    track.SetXYZ(px, py, pz);
    phi0 = track.Phi();
    the0 = track.Theta();


    for(int j = 1; j < 3000; j++){

        float length = j*0.01;
        float phi = phi0 + (0.5 * 0.001 * 0.3 * 1 * length * charge)/(sqrt(px*px + py*py));

        track.SetMagThetaPhi(length, the0, phi);

        tx[j] = pcax + track.X();
        ty[j] = pcay + track.Y();
        tz[j] = pcaz + track.Z();
    }

    TGraph2D *gr2 = new TGraph2D(3000, tx, ty, tz);

    TGraph2D *gr1 = new TGraph2D(npoints, &x[0], &y[0], &z[0]);
    gr1->SetTitle(Form("Entry %d;X [cm];Y [cm];Z [cm]", n));
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.5);
    gr1->SetMarkerColor(kBlue);
    gr1->Draw();


    gr2->SetTitle(Form("Entry %d;X [cm];Y [cm];Z [cm]", n));
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(0.5);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("SAME");




}*/


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

void event_display(){

    // Open the ROOT file
    TFile *file = TFile::Open("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_51576_11_resid.root");

    // Check if the file was opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Get the tree from the file
    TTree *tree = (TTree *) file->Get("residualtree;20");

    // Check if the tree was retrieved successfully
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }

    // Declare variables to hold the data
    std::vector<float> *clusgx = nullptr;
    std::vector<float> *clusgy = nullptr;
    std::vector<float> *clusgz = nullptr;
    float pt, eta, dedx, quality, phi, pcax, pcay, pcaz, phi0, the0, px, py, pz;
    int nhits, charge, ntpc, nintt, nmaps;

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

    Long64_t nentries = tree->GetEntries();

    // Variables to store the data for the track to be processed
    std::vector<double> x, y, z;
    std::vector<double> tx, ty, tz;
    std::vector<double> vx, vy, vz;
    size_t npoints = 0;
    Long64_t entryNumber = -1; // To keep track of the entry number

    // Loop over all entries to find the first track without hits in INTT and MAPS
    for (Long64_t i = 0; i < nentries; i++) {
        // Get the specified entry (track)
        tree->GetEntry(i);

        // Use logical AND '&&' instead of bitwise AND '&'
        if(nintt == 0 && nmaps == 0) {
        //if(1) {
            // Get the number of hit points
            tx.push_back(pcax + track.X());
            ty.push_back(pcay + track.Y());
            tz.push_back(pcaz + track.Z());
            npoints = clusgx->size();

            // Check if there are any hits to plot
            if (npoints == 0) {
                std::cerr << "No hits to plot for entry " << i << std::endl;
                continue; // Skip this entry and continue the loop
            }

            // Convert std::vector<float> to std::vector<double> for TGraph2D
            x.assign(clusgx->begin(), clusgx->end());
            y.assign(clusgy->begin(), clusgy->end());
            z.assign(clusgz->begin(), clusgz->end());

            // Prepare the track propagation
            TVector3 track(px, py, pz);
            phi0 = track.Phi();
            the0 = track.Theta();

            // Use std::vector<double> for track points
            double length_step = 0.01;
            int nsteps = 10000;

            for (int j = 0; j < nsteps; j++) {
                double length = j * length_step;
                double phi = phi0 - (charge * 0.001 * 1.4 * length) /(sqrt(px*px + py*py)); // Corrected formula with B=1.4T
                track.SetMagThetaPhi(length, the0, phi);

                tx.push_back(pcax + track.X());
                ty.push_back(pcay + track.Y());
                tz.push_back(pcaz + track.Z());
            }

            entryNumber = i; // Store the entry number
            break; // Exit the loop after processing the first matching track
        }
    }

    // Check if a track was found
    if (entryNumber == -1) {
        std::cerr << "No track found without hits in INTT and MAPS" << std::endl;
        file->Close();
        return;
    }

    // Create TGraph2D for the hits
    TGraph2D *gr1 = new TGraph2D(npoints, &x[0], &y[0], &z[0]);
    gr1->SetTitle(Form("Entry %lld;X [cm];Y [cm];Z [cm]", entryNumber));
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.5);
    gr1->SetMarkerColor(kBlue);

    // Create TGraph2D for the propagated track
    TGraph2D *gr2 = new TGraph2D(tx.size(), &tx[0], &ty[0], &tz[0]);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(0.5);
    gr2->SetMarkerColor(kRed);

    // Create a canvas to draw the graphs
    TCanvas *c = new TCanvas("c", "3D Track Display", 800, 600);

    // Draw the first graph
    //gr2->Draw("P0");

    // Draw the second graph on the same canvas
   // gr1->Draw("P0 SAME");
    gr1->Draw("P0");

    // Add a legend
    TLegend *legend = new TLegend(0.75, 0.85, 0.95, 0.95);
    legend->AddEntry(gr1, "Hits", "p");
    legend->AddEntry(gr2, "Propagated Track", "p");
    legend->Draw();

    // Update the canvas
    c->Update();

    // Optionally, save the canvas to a file
    c->SaveAs(Form("event_display_%lld.png", entryNumber));

    /*// Clean up
    delete gr1;
    delete gr2;
    delete c;
    file->Close();*/
}

