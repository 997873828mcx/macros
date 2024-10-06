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

void extract_mag_data(){

// Open the ROOT file
    TFile *file = TFile::Open("sphenix3dbigmapxyz_steel_rebuild.root");  // Replace with your ROOT file

    // Check if the file was opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }


    FILE *ffield = fopen("field_map.dat","w");
    // Get the tree from the file
    TTree *tree = (TTree *) file->Get("fieldmap;7");  // Replace with your tree name

    // Check if the tree was retrieved successfully
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }

    // Declare variables to hold the data from the branches
    float x, y, z, bx, by, bz;
    std::set<float> coords[3];
    std::vector<float> bfield[3];

    // Set branch addresses to read the data from the tree
    tree->SetBranchAddress("x", &x);   // Replace "x" with your actual branch name
    tree->SetBranchAddress("y", &y);   // Replace "y" with your actual branch name
    tree->SetBranchAddress("z", &z);   // Replace "z" with your actual branch name
    tree->SetBranchAddress("bx", &bx); // Replace "bx" with your actual branch name
    tree->SetBranchAddress("by", &by); // Replace "by" with your actual branch name
    tree->SetBranchAddress("bz", &bz); // Replace "bz" with your actual branch name
/*
    // Open a text file to save the coordinates
    std::ofstream outfile("coordinates_output.txt");  // This will save to coordinates_output.txt

    // Check if the file opened successfully
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file" << std::endl;
        file->Close();
        return;
    }*/

    // Get the total number of entries in the tree
    Long64_t nentries = tree->GetEntries();

    // Loop over each entry in the tree
    for (Long64_t i = 0; i < nentries; i++) {
        // Load the data for the i-th entry
        tree->GetEntry(i);

        fprintf(ffield,"OUT %g %g %g %g %g %g\n", x, y, z, bx, by, bz);
        coords[0].insert(x);
        coords[1].insert(y);
        coords[2].insert(z);
        bfield[0].push_back(bx); //field intensity in x direction
        bfield[1].push_back(by);
        bfield[2].push_back(bz);
    }

    for (int i = 0; i < 3; i++)
    {
        min[i] = *coords[i].begin(); //minimal coordinate in each direction
        max[i] = *coords[i].rbegin();// maximal coordinate in each direction
        steps[i] = coords[i].size();//number of grids in each direction
        std::cout << "We have " << coords[i].size() << " from " << *(coords[i].begin()) << " to " << *(coords[i].rbegin()) << std::endl;
    }
    // Close the output file
    outfile.close();

    // Close the ROOT file
    file->Close();

    std::cout << "Coordinates saved to coordinates_output.txt" << std::endl;
}
