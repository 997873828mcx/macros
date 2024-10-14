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
#include <TSystem.h>
#include <phool/recoConsts.h>

void extract_mag_data() {

    // Open the ROOT file
    std::string filename = "/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/release/release_new/new.3/share/calibrations/Field/Map/sphenix3dtrackingmapxyz.root";
    TFile *file = TFile::Open(filename.c_str());

    // Check if the file was opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }



    // Get the tree from the file
    TTree *tree = (TTree *) file->Get("fieldmap;1");

    // Check if the tree was retrieved successfully
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }

    // Declare variables to hold the data from the branches
    float x, y, z, bx, by, bz;


    // Set branch addresses to read the data from the tree
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("bx", &bx);
    tree->SetBranchAddress("by", &by);
    tree->SetBranchAddress("bz", &bz);

    int nbins_x, nbins_y, nbins_z;
    float x_min, x_max, y_min, y_max, z_min, z_max;

// Create 3D histograms for Bx, By, Bz, and B (total magnitude)
if(filename == "/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/release/release_new/new.3/share/calibrations/Field/Map/sphenix3dbigmapxyz_steel_rebuild.root") {
     nbins_x = 271;  // Number of bins in x
     nbins_y = 271;  // Number of bins in y
     nbins_z = 331;  // Number of bins in z
     x_min = -271, x_max = 271;  // Set appropriate min and max for x
     y_min = -271, y_max = 271;  // Set appropriate min and max for y
     z_min = -331, z_max = 331;  // Set appropriate min and max for z
}

else if(filename=="/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/release/release_new/new.3/share/calibrations/Field/Map/sphenix3dtrackingmapxyz.root")
{
     nbins_x = 111;  // Number of bins in x
     nbins_y = 111;  // Number of bins in y
     nbins_z = 111;  // Number of bins in z
     x_min = -111, x_max = 111;  // Set appropriate min and max for x
     y_min = -111, y_max = 111;  // Set appropriate min and max for y
     z_min = -111, z_max = 111;  // Set appropriate min and max for z
} else {
    std::cerr << "Unknown file. Cannot determine histogram settings." << std::endl;
    file->Close();
    return;
}
    // 3D histograms for Bx, By, Bz, and the total field magnitude
    TH3F *hist_bx = new TH3F("hist_bx", "Bx (3D)", nbins_x, x_min, x_max, nbins_y, y_min, y_max, nbins_z, z_min, z_max);
    TH3F *hist_by = new TH3F("hist_by", "By (3D)", nbins_x, x_min, x_max, nbins_y, y_min, y_max, nbins_z, z_min, z_max);
    TH3F *hist_bz = new TH3F("hist_bz", "Bz (3D)", nbins_x, x_min, x_max, nbins_y, y_min, y_max, nbins_z, z_min, z_max);
    // TH3F *hist_b = new TH3F("hist_b", "B (Total Field Magnitude)", nbins_x, x_min, x_max, nbins_y, y_min, y_max, nbins_z, z_min, z_max);


    // Open output text file to save coordinates
    std::ofstream outfile("coordinates_output_track.txt");
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file" << std::endl;
        file->Close();
        return;
    }

    // Get the total number of entries in the tree
    Long64_t nentries = tree->GetEntries();

    std::set<std::tuple<int, int, int>> visitedBins;
    // Loop over each entry in the tree
    for (Long64_t i = 0; i < nentries; i++) {
        // Load the data for the i-th entry
        tree->GetEntry(i);
        // Save the coordinates to the output file
        outfile << "Entry " << i << ": x = " << x << ", y = " << y << ", z = " << z << std::endl;


// Find the appropriate bin for the given coordinates
        int bin_x = hist_bx->GetXaxis()->FindBin(x);
        int bin_y = hist_bx->GetYaxis()->FindBin(y);
        int bin_z = hist_bx->GetZaxis()->FindBin(z);
        if((hist_bx->GetXaxis()->GetBinCenter(bin_x)!=x)||hist_bx->GetYaxis()->GetBinCenter(bin_y)!=y||hist_bx->GetZaxis()->GetBinCenter(bin_z)!=z)
        {
            std::cerr<<"it's not in the center"<<bin_x<<" "<<x<<std::endl;
        }

        std::tuple<int, int, int> currentBin(bin_x, bin_y, bin_z);

        // Check if the bin has already been filled
        if (visitedBins.find(currentBin) == visitedBins.end()) {
            // If not filled yet, set the content and add to the visited set
            hist_bx->SetBinContent(bin_x, bin_y, bin_z, bx);
            hist_by->SetBinContent(bin_x, bin_y, bin_z, by);
            hist_bz->SetBinContent(bin_x, bin_y, bin_z, bz);
            visitedBins.insert(currentBin);
        } else {
            // If already filled, handle as needed (e.g., skip, or log warning)
            std::cerr << "Warning: Overlapping point at bin (" << bin_x << ", " << bin_y << ", " << bin_z << ")" << std::endl;
        }

    /*    hist_by->SetBinContent(bin_x, bin_y, bin_z, by);
        hist_bz->SetBinContent(bin_x, bin_y, bin_z, bz);*/

        if (i % 100000 == 0) {
            std::cout << "Processed " << i << " / " << nentries << " entries." << std::endl;
            std::cout<<hist_bx->GetBinContent(bin_x, bin_y, bin_z)<<std::endl;
        }

    }
    outfile.close();

    // Create a canvas to draw the histograms
    /*TCanvas *c = new TCanvas("c", "Magnetic Field 3D Components", 1200, 900);
    c->Divide(2, 2);  // Divide the canvas into 2x2 subplots*/

   /* // Draw the 3D histograms
    c->cd(1);
    hist_bx->GetXaxis()->SetTitle("X");
    hist_bx->GetYaxis()->SetTitle("Y");
    hist_bx->GetZaxis()->SetTitle("Bx");
    hist_bx->Draw("LEGO2");*/

    /*c->cd(2);
    hist_by->GetXaxis()->SetTitle("X");
    hist_by->GetYaxis()->SetTitle("Y");
    hist_by->GetZaxis()->SetTitle("By");
    hist_by->Draw("BOX2");

    c->cd(3);
    hist_bz->GetXaxis()->SetTitle("X");
    hist_bz->GetYaxis()->SetTitle("Y");
    hist_bz->GetZaxis()->SetTitle("Bz");
    hist_bz->Draw("BOX2");*/



    // Save the canvas to a file
    //c->SaveAs("magnetic_field_3d.png");

    // Save the histograms to a new ROOT file
    TFile *outfile_root = new TFile("magnetic_field_histograms_track.root", "RECREATE");
    if (!outfile_root || outfile_root->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        return;
    }

    hist_bx->Write();
    hist_by->Write();
    hist_bz->Write();
    file->Close();
    outfile_root->Close();



    std::cout << "Coordinates saved to coordinates_output_track.txt" << std::endl;

}
