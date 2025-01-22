{
  // Open the ROOT file
  TFile *f = TFile::Open("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_5_resid_edgeOff_staticOff_acts_0120_1_5.root");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file" << std::endl;
    return;
  }

  // Retrieve the "clustertree" tree from the file
  TTree *tree = (TTree*)f->Get("clustertree");
  if (!tree) {
    std::cerr << "Error: clustertree not found" << std::endl;
    f->Close();
    return;
  }

  // Prepare pointers for the branches
  std::vector<short> *clust_crossings = nullptr;
  unsigned long m_scluskey = 0;  // variable to hold cluster key

  // Set branch addresses
  tree->SetBranchAddress("clust_crossings", &clust_crossings);
  tree->SetBranchAddress("cluskey", &m_scluskey);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    // Skip if no crossing data available
    if (!clust_crossings) continue;

    std::set<short> uniqueCrossings;
    // Collect unique crossing numbers excluding -1
    for (auto crossing : *clust_crossings) {
      if (1) {
        uniqueCrossings.insert(crossing);
      }
    }

    // If there are more than one unique crossing number, print the cluster key and crossings
    if (uniqueCrossings.size()) {
      std::cout << "Cluster key " << m_scluskey << " has distinct crossings: ";
      for (auto val : uniqueCrossings) {
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }
  }

  f->Close();
}