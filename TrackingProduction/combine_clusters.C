void separate_cluster_and_hits() {
  // Open input files
  TFile* f_resid = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_resid_edgeOn_staticOff.root", "READ");
  TFile* f_seeding = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff_1110_0.root", "READ");

  // Get trees from residual file
  TTree* t_residual_clus = (TTree*)f_resid->Get("clustertree");
  TTree* t_residual = (TTree*)f_resid->Get("residualtree");
  TTree* t_hits = (TTree*)f_resid->Get("hittree"); // If available

  // Get tracking_clusters tree from the tracking file
  TTree* t_seeding = (TTree*)f_seeding->Get("tracking_clusters");

  // --------------------------
  // Read tracking_clusters tree
  // --------------------------
  ULong64_t tr_cluskey;
  float tr_x, tr_y, tr_z;
  int tr_used_in_seed;

  t_tracking->SetBranchAddress("cluskey", &tr_cluskey);
  t_tracking->SetBranchAddress("x", &tr_x);
  t_tracking->SetBranchAddress("y", &tr_y);
  t_tracking->SetBranchAddress("z", &tr_z);
  t_tracking->SetBranchAddress("used_in_seed", &tr_used_in_seed);

  struct TrackingInfo {
    float x;
    float y;
    float z;
    int   used_in_seed;
  };

  std::map<ULong64_t, TrackingInfo> tracking_map;

  for (int i = 0; i < t_tracking->GetEntries(); i++) {
    t_tracking->GetEntry(i);
    TrackingInfo info = {tr_x, tr_y, tr_z, tr_used_in_seed};
    tracking_map[tr_cluskey] = info;
  }

  // --------------------------
  // Determine which clusters are used in tracks from residualtree
  // --------------------------
  std::vector<ULong64_t>* track_cluskeys = nullptr;
  t_residual->SetBranchAddress("cluskeys", &track_cluskeys);

  std::set<ULong64_t> clusters_in_tracks;
  for (int i = 0; i < t_residual->GetEntries(); i++) {
    t_residual->GetEntry(i);
    for (auto key : *track_cluskeys) {
      clusters_in_tracks.insert(key);
    }
  }

  // --------------------------
  // Read hittree - we will just copy its info into a new tree
  // --------------------------
  ULong64_t hit_cluskey;
  int nhits = 0;
  float avg_adc = 0;
  if (t_hits) {
    t_hits->SetBranchAddress("cluskey", &hit_cluskey);
    if (t_hits->GetBranch("nhits")) t_hits->SetBranchAddress("nhits", &nhits);
    if (t_hits->GetBranch("avg_adc")) t_hits->SetBranchAddress("avg_adc", &avg_adc);
  }

  struct HitInfo {
    int nhits;
    float avg_adc;
  };

  std::map<ULong64_t, HitInfo> hit_map;
  if (t_hits) {
    for (int i = 0; i < t_hits->GetEntries(); i++) {
      t_hits->GetEntry(i);
      HitInfo hinfo = {nhits, avg_adc};
      hit_map[hit_cluskey] = hinfo;
    }
  }

  // --------------------------
  // Read clustertree from residual file (which has all clusters)
  // --------------------------
  ULong64_t res_cluskey;
  int layer = -1;
  float adc = 0.0;
  float charge = 0.0;

  t_residual_clus->SetBranchAddress("cluskey", &res_cluskey);
  if (t_residual_clus->GetBranch("layer")) t_residual_clus->SetBranchAddress("layer", &layer);
  if (t_residual_clus->GetBranch("adc")) t_residual_clus->SetBranchAddress("adc", &adc);
  if (t_residual_clus->GetBranch("charge")) t_residual_clus->SetBranchAddress("charge", &charge);

  // --------------------------
  // Create the output file and two separate trees:
  // 1) combined_clusters: cluster-level info
  // 2) combined_hits: hit-level info
  // --------------------------
  TFile* outFile = new TFile("separate_cluster_and_hits.root", "RECREATE");

  // Cluster tree
  TTree* t_combined_clusters = new TTree("combined_clusters", "Combined cluster info");
  ULong64_t out_cluskey;
  float out_x=0, out_y=0, out_z=0;
  int out_used_in_seed=0;
  int out_used_in_track=0;
  int out_layer=-1;
  float out_adc=0, out_charge=0;
  t_combined_clusters->Branch("cluskey", &out_cluskey, "cluskey/l");
  t_combined_clusters->Branch("x", &out_x, "x/F");
  t_combined_clusters->Branch("y", &out_y, "y/F");
  t_combined_clusters->Branch("z", &out_z, "z/F");
  t_combined_clusters->Branch("used_in_seed", &out_used_in_seed, "used_in_seed/I");
  t_combined_clusters->Branch("used_in_track", &out_used_in_track, "used_in_track/I");
  t_combined_clusters->Branch("layer", &out_layer, "layer/I");
  t_combined_clusters->Branch("adc", &out_adc, "adc/F");
  t_combined_clusters->Branch("charge", &out_charge, "charge/F");

  // Hit tree (separated out)
  TTree* t_combined_hits = new TTree("combined_hits", "Hit-level info");
  ULong64_t out_hit_cluskey;
  int out_nhits=0;
  float out_avg_adc=0;
  t_combined_hits->Branch("cluskey", &out_hit_cluskey, "cluskey/l");
  t_combined_hits->Branch("nhits", &out_nhits, "nhits/I");
  t_combined_hits->Branch("avg_adc", &out_avg_adc, "avg_adc/F");

  // Fill cluster tree
  for (int i = 0; i < t_residual_clus->GetEntries(); i++) {
    t_residual_clus->GetEntry(i);
    out_cluskey = res_cluskey;
    out_layer = layer;
    out_adc = adc;
    out_charge = charge;

    if (tracking_map.find(res_cluskey) != tracking_map.end()) {
      out_x = tracking_map[res_cluskey].x;
      out_y = tracking_map[res_cluskey].y;
      out_z = tracking_map[res_cluskey].z;
      out_used_in_seed = tracking_map[res_cluskey].used_in_seed;
    } else {
      out_x = out_y = out_z = 0;
      out_used_in_seed = 0;
    }

    out_used_in_track = (clusters_in_tracks.find(res_cluskey) != clusters_in_tracks.end()) ? 1 : 0;

    t_combined_clusters->Fill();
  }

  // Fill hit tree (just copy from hit_map)
  for (auto &hit_entry : hit_map) {
    out_hit_cluskey = hit_entry.first;
    out_nhits = hit_entry.second.nhits;
    out_avg_adc = hit_entry.second.avg_adc;
    t_combined_hits->Fill();
  }

  // Write out both trees
  outFile->cd();
  t_combined_clusters->Write();
  t_combined_hits->Write();
  outFile->Close();

  f_resid->Close();
  f_track->Close();

  std::cout << "Wrote separate_cluster_and_hits.root with combined_clusters and combined_hits trees." << std::endl;
}
