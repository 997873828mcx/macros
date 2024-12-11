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
  int tr_used_in_seed, tr_passed_straight;


  t_seeding->SetBranchAddress("cluskey", &tr_cluskey);
  t_seeding->SetBranchAddress("x", &tr_x);
  t_seeding->SetBranchAddress("y", &tr_y);
  t_seeding->SetBranchAddress("z", &tr_z);
  t_seeding->SetBranchAddress("passed_straight", &tr_passed_straight);
  t_seeding->SetBranchAddress("used_in_seed", &tr_used_in_seed);

  struct TrackingInfo {
    float x;
    float y;
    float z;
    int   used_in_seed;
    int  passed_straight;
  };

  std::map<ULong64_t, TrackingInfo> tracking_map;

  for (int i = 0; i < t_seeding->GetEntries(); i++) {
    t_seeding->GetEntry(i);
    TrackingInfo info = {tr_x, tr_y, tr_z, tr_passed_straight, tr_used_in_seed};
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
  int m_runnumber, m_segment, m_job, m_event;
  ULong64_t m_bco, m_bcotr;
  UInt_t m_hitsetkey;
  Float_t m_hitgx, m_hitgy, m_hitgz;
  Int_t m_hitlayer, m_sector, m_side, m_staveid, m_chipid;
  Int_t m_strobeid, m_ladderzid, m_ladderphiid, m_timebucket, m_hitpad, m_hittbin;
  Int_t m_col, m_row, m_segtype, m_tileid, m_strip;
  Float_t m_adc, m_zdriftlength;

  t_hits->SetBranchAddress("m_runnumber", &m_runnumber);
  t_hits->SetBranchAddress("m_segment", &m_segment);
  t_hits->SetBranchAddress("m_job", &m_job);
  t_hits->SetBranchAddress("m_event", &m_event);
  t_hits->SetBranchAddress("m_bco", &m_bco);
  t_hits->SetBranchAddress("m_bcotr", &m_bcotr);
  t_hits->SetBranchAddress("m_hitsetkey", &m_hitsetkey);
  t_hits->SetBranchAddress("m_hitgx", &m_hitgx);
  t_hits->SetBranchAddress("m_hitgy", &m_hitgy);
  t_hits->SetBranchAddress("m_hitgz", &m_hitgz);
  t_hits->SetBranchAddress("m_hitlayer", &m_hitlayer);
  t_hits->SetBranchAddress("m_sector", &m_sector);
  t_hits->SetBranchAddress("m_side", &m_side);
  t_hits->SetBranchAddress("m_staveid", &m_staveid);
  t_hits->SetBranchAddress("m_chipid", &m_chipid);
  t_hits->SetBranchAddress("m_strobeid", &m_strobeid);
  t_hits->SetBranchAddress("m_ladderzid", &m_ladderzid);
  t_hits->SetBranchAddress("m_ladderphiid", &m_ladderphiid);
  t_hits->SetBranchAddress("m_timebucket", &m_timebucket);
  t_hits->SetBranchAddress("m_hitpad", &m_hitpad);
  t_hits->SetBranchAddress("m_hittbin", &m_hittbin);
  t_hits->SetBranchAddress("m_col", &m_col);
  t_hits->SetBranchAddress("m_row", &m_row);
  t_hits->SetBranchAddress("m_segtype", &m_segtype);
  t_hits->SetBranchAddress("m_tileid", &m_tileid);
  t_hits->SetBranchAddress("m_strip", &m_strip);
  t_hits->SetBranchAddress("m_adc", &m_adc);
  t_hits->SetBranchAddress("m_zdriftlength", &m_zdriftlength);

/*   struct HitInfo {
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
  } */

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

  // Combined hits tree: copy from hittree input
  TTree* t_combined_hits = new TTree("combined_hits", "Hit-level info");
  // The branches will be the same as the input hittree
  t_combined_hits->Branch("runnumber", &m_runnumber, "runnumber/I");
  t_combined_hits->Branch("segment", &m_segment, "segment/I");
  t_combined_hits->Branch("job", &m_job, "job/I");
  t_combined_hits->Branch("event", &m_event, "event/I");
  t_combined_hits->Branch("gl1bco", &m_bco, "gl1bco/l");
  t_combined_hits->Branch("trbco", &m_bcotr, "trbco/l");
  t_combined_hits->Branch("hitsetkey", &m_hitsetkey, "hitsetkey/i");
  t_combined_hits->Branch("gx", &m_hitgx, "gx/F");
  t_combined_hits->Branch("gy", &m_hitgy, "gy/F");
  t_combined_hits->Branch("gz", &m_hitgz, "gz/F");
  t_combined_hits->Branch("layer", &m_hitlayer, "layer/I");
  t_combined_hits->Branch("sector", &m_sector, "sector/I");
  t_combined_hits->Branch("side", &m_side, "side/I");
  t_combined_hits->Branch("stave", &m_staveid, "stave/I");
  t_combined_hits->Branch("chip", &m_chipid, "chip/I");
  t_combined_hits->Branch("strobe", &m_strobeid, "strobe/I");
  t_combined_hits->Branch("ladderz", &m_ladderzid, "ladderz/I");
  t_combined_hits->Branch("ladderphi", &m_ladderphiid, "ladderphi/I");
  t_combined_hits->Branch("timebucket", &m_timebucket, "timebucket/I");
  t_combined_hits->Branch("pad", &m_hitpad, "pad/I");
  t_combined_hits->Branch("tbin", &m_hittbin, "tbin/I");
  t_combined_hits->Branch("col", &m_col, "col/I");
  t_combined_hits->Branch("row", &m_row, "row/I");
  t_combined_hits->Branch("segtype", &m_segtype, "segtype/I");
  t_combined_hits->Branch("tile", &m_tileid, "tile/I");
  t_combined_hits->Branch("strip", &m_strip, "strip/I");
  t_combined_hits->Branch("adc", &m_adc, "adc/F");
  t_combined_hits->Branch("zdriftlength", &m_zdriftlength, "zdriftlength/F");

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

/*   // Fill hit tree (just copy from hit_map)
  for (auto &hit_entry : hit_map) {
    out_hit_cluskey = hit_entry.first;
    out_nhits = hit_entry.second.nhits;
    out_avg_adc = hit_entry.second.avg_adc;
    t_combined_hits->Fill();
  } */

 for (int i = 0; i < t_hits->GetEntries(); i++) {
    t_hits->GetEntry(i);
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
