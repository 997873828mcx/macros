void plot_unused_clusters() {
    //this script is used to plot the clusters which were not used in tracks with different color
    //for clusters used in seeds and not used in seeds seperately
  // Open files that contain cluster information
  TFile* f_tracking = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff_1110_0.root");  //PHCASeeding output
  TFile* f_residual = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_resid_edgeOn_staticOff.root");   // TrackResiduals output
  TFile* outfile = new TFile("unused_cluster_plots.root", "RECREATE");
  // Get trees
  TTree* t_tracking = (TTree*)f_tracking->Get("tracking_clusters");
  TTree* t_residual = (TTree*)f_residual->Get("residualtree");
  TTree* t_residual_clus = (TTree*)f_residual->Get("clustertree");
  
  unsigned long residual_cluskey;
  t_residual_clus->SetBranchAddress("cluskey", &residual_cluskey);
  // Variables to read
  std::set<unsigned long long> all_clusters;
  unsigned long long cluskey;   
  float x, y, z;
  int used_in_seed;
  t_tracking->SetBranchAddress("cluskey", &cluskey);
  t_tracking->SetBranchAddress("x", &x);
  t_tracking->SetBranchAddress("y", &y);
  t_tracking->SetBranchAddress("z", &z);
  t_tracking->SetBranchAddress("used_in_seed", &used_in_seed);

  for(int i = 0; i < t_tracking->GetEntries(); i++) {
    t_tracking->GetEntry(i);
        t_tracking->GetEntry(i);
    if(i < 5) {  // Print first 5 entries
        std::cout << "Entry " << i << ": cluskey = " << cluskey 
                 << " (0x" << std::hex << cluskey << std::dec << ")" 
                 << " x = " << x << " y = " << y << " z = " << z 
                 << " used_in_seed = " << used_in_seed << std::endl;
    }
    all_clusters.insert(cluskey);
  }

  for(int i = 0; i < 5 && i < t_residual_clus->GetEntries(); i++) {
    t_residual_clus->GetEntry(i);
    std::cout << "Clustertree Entry " << i << ": cluskey = " << residual_cluskey 
              << " (0x" << std::hex << residual_cluskey << std::dec << ")" << std::endl;
}

  std::cout << "Size of all_clusters after insertion: " << all_clusters.size() << std::endl;
  // Get clusters used in tracks first
  std::set<unsigned long long> clusters_in_tracks;
  std::vector<unsigned long long>* track_cluskeys = nullptr;
  int n_tpc=0;
  t_residual->SetBranchAddress("cluskeys", &track_cluskeys);
  t_residual->SetBranchAddress("ntpc", &n_tpc);


  
  int missing_clusters = 0;
  int num_tpc=0;
  int num_in_tpc_clusters=0;
   int total_residual_clusters = 0;
 
    int found_in_all_clusters = 0;
  for(int i = 0; i < t_residual->GetEntries(); i++) {
    t_residual->GetEntry(i);
    num_tpc+=n_tpc;
    for(auto key : *track_cluskeys) {
        
      clusters_in_tracks.insert(key);
      // Check if this cluster exists in tracking_clusters tree
      if(all_clusters.find(key) == all_clusters.end()) {
        //std::cout << "WARNING: Cluster " << key << " used in track but not found in tracking_clusters tree!" << std::endl;
        missing_clusters++;
      }
      else
      {
        num_in_tpc_clusters++;
      }
    }
  }

     for(int i = 0; i < t_residual_clus->GetEntries(); i++) {
        t_residual_clus->GetEntry(i);
        total_residual_clusters++;
        
        if(all_clusters.find((unsigned long long)residual_cluskey) != all_clusters.end()) {
    found_in_all_clusters++;
}
    }

  // Create a single TGraph2D for unused clusters
  TGraph2D* g_seed_clusters = new TGraph2D(); // For clusters used in seeds
  TGraph2D* g_noseed_clusters = new TGraph2D(); // For clusters not used in seeds



  // Fill graphs only with clusters NOT in tracks
  int n_seed = 0, n_noseed = 0;
  for(int i = 0; i < t_tracking->GetEntries(); i++) {
    t_tracking->GetEntry(i);
    
    // Only process if cluster is NOT in tracks
    if(clusters_in_tracks.find(cluskey) == clusters_in_tracks.end()) {
        if(used_in_seed) {
            g_seed_clusters->SetPoint(n_seed++, x, y, z);
        } else {
            g_noseed_clusters->SetPoint(n_noseed++, x, y, z);
        }
    }
  }

  // Create canvas and frame
  TCanvas* c = new TCanvas("c", "Unused Clusters 3D", 1200, 1200);
  c->cd();

  TH3D* frame = new TH3D("frame", "Clusters Not Used in Tracks;X;Y;Z",
                         100, -110, 110,
                         100, -110, 110,
                         100, -110, 110);
  frame->GetXaxis()->SetTitle("X");
  frame->GetYaxis()->SetTitle("Y");
  frame->GetZaxis()->SetTitle("Z");
  frame->SetStats(0);

  // Style graphs
  g_seed_clusters->SetMarkerColor(kBlue);
  g_seed_clusters->SetMarkerStyle(20);
  g_seed_clusters->SetMarkerSize(0.2);
  
  g_noseed_clusters->SetMarkerColor(kRed);
  g_noseed_clusters->SetMarkerStyle(20);
  g_noseed_clusters->SetMarkerSize(0.2);

  // Draw
  frame->Draw();
  g_seed_clusters->Draw("P");
  g_noseed_clusters->Draw("P SAME");

  g_seed_clusters->GetHistogram()->GetXaxis()->SetLimits(-80,80);
  g_seed_clusters->GetHistogram()->GetYaxis()->SetLimits(-80,80);
  g_seed_clusters->GetHistogram()->GetZaxis()->SetRangeUser(-250,250);

/*   // Add legend
  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(g_seed_clusters, "Used in Seeds", "p");
  leg->AddEntry(g_noseed_clusters, "Not Used in Seeds", "p");
  leg->Draw(); */

  c->SetTheta(30);
  c->SetPhi(30);

  // Print counts
  std::cout << "Total clusters not in tracks: " << n_seed + n_noseed << std::endl;
  std::cout << "   Used in seeds: " << n_seed << std::endl;
  std::cout << "   Not used in seeds: " << n_noseed << std::endl;
  std::cout << "   Number of clusters in tpc in all tracks: " << num_tpc << std::endl;
  std::cout << "   Number of clusters in all available tpc clusters from PHCASeeding module in all tracks: " << num_in_tpc_clusters << std::endl;

  std::cout << "Total clusters in all_clusters: " << all_clusters.size() << std::endl;
    std::cout << "Total clusters in TrackResiduals clustertree: " << total_residual_clusters << std::endl;
    std::cout << "Clusters found in both trees: " << found_in_all_clusters << std::endl;

  outfile->cd();
  c->Write();
  outfile->Close();
  f_tracking->Close();
  f_residual->Close();
}