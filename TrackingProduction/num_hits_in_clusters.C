// ROOT macro to plot the 'num_hits' distribution from 'tpc_clustertree'


void num_hits_in_clusters() {
   
    const char* inputFileName = "outputFile_kso_52077_0_clusterizer_edgeOn_staticOff_1219_0.root";
    
    
    const char* outputFileName = "output_num_hits.root";
    
    
    
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file '" << inputFileName << "'." << std::endl;
        return;
    }
    std::cout << "Opened input file: " << inputFileName << std::endl;
    
 
    
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("tpc_clustertree"));
    if (!tree) {
        std::cerr << "Error: TTree 'tpc_clustertree' not found in file '" << inputFileName << "'." << std::endl;
        inputFile->Close();
        return;
    }
    std::cout << "Retrieved TTree: tpc_clustertree" << std::endl;
    
    
    int nbins = 40;
    double xmin = 0;
    double xmax = 40;
    
    // Create a 1D histogram for 'num_hits'
    TH1F* h_num_hits = new TH1F("h_num_hits", "Number of Hits per Cluster;num_hits;Entries", nbins, xmin, xmax);
    
    
    tree->Draw("num_hits>>h_num_hits", "", "goff");
    std::cout << "Filled histogram 'h_num_hits' with 'num_hits' data." << std::endl;
    
    
    h_num_hits->SetLineColor(kBlue);
    h_num_hits->SetLineWidth(2);
    h_num_hits->SetFillColor(kGreen-10);
    
    
    TCanvas* c1 = new TCanvas("c1", "Number of Hits per Cluster", 800, 600);
    h_num_hits->Draw("HIST");
    
    // Optionally, add statistics box
    gStyle->SetOptStat(1110); // Display entries, mean, and RMS
    c1->Update();

    
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file '" << outputFileName << "'." << std::endl;
        inputFile->Close();
        delete c1;
        delete h_num_hits;
        return;
    }
    std::cout << "Created output file: " << outputFileName << std::endl;
    
    // Write the histogram to the output file
    h_num_hits->Write();
    std::cout << "Wrote histogram 'h_num_hits' to '" << outputFileName << "'." << std::endl;
    
    // -------------------------------
    // 9. Save the Canvas as an Image (Optional)
    // -------------------------------
    
    // Save the canvas as a PNG image
    c1->SaveAs("num_hits_histogram.png");
    std::cout << "Saved histogram plot as 'num_hits_histogram.png'." << std::endl;
    
    // -------------------------------
    // 10. Clean Up and Close Files
    // -------------------------------
    
    // Close the output file
    outputFile->Close();
    std::cout << "Closed output file." << std::endl;
    
    // Close the input file
    //inputFile->Close();
    std::cout << "Closed input file." << std::endl;
    
    // Delete dynamically allocated objects to free memory
    //delete c1;
    //delete h_num_hits;
    
    std::cout << "Histogram plotting and saving completed successfully." << std::endl;
}
