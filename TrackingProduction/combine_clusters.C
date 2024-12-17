void combine_clusters() {
  // Open input files
  TFile* f_resid = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_52077_0_resid_edgeOn_staticOff_1216_0.root", "READ");
  TFile* f_seed = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_52077_0_clusters_edgeOn_staticOff_1216_0.root", "READ");

  // Get trees from residual file
  TTree* t_residual_clus = (TTree*)f_resid->Get("clustertree");
  TTree* t_residual = (TTree*)f_resid->Get("residualtree");
  TTree* t_hits = (TTree*)f_resid->Get("hittree");

  // Get tracking_clusters tree from the tracking file
  TTree* t_seeding = (TTree*)f_seed->Get("tracking_clusters");

  // --------------------------
  // Read tracking_clusters tree
  // --------------------------
  unsigned long long tr_cluskey;
  float tr_x, tr_y, tr_z;
  int tr_used_in_seed;

  t_seeding->SetBranchAddress("cluskey", &tr_cluskey);
  t_seeding->SetBranchAddress("x", &tr_x);
  t_seeding->SetBranchAddress("y", &tr_y);
  t_seeding->SetBranchAddress("z", &tr_z);
  t_seeding->SetBranchAddress("used_in_seed", &tr_used_in_seed);

  struct TrackingInfo {
    float x;
    float y;
    float z;
    int used_in_seed;
  };

  std::map<unsigned long long, TrackingInfo> seeding_map;
  for (int i = 0; i < t_seeding->GetEntries(); i++) {
    t_seeding->GetEntry(i);
    TrackingInfo info = {tr_x, tr_y, tr_z, tr_used_in_seed};
    seeding_map[tr_cluskey] = info;
  }

  // --------------------------
  // Determine which clusters are used in tracks from residualtree
  // --------------------------
  std::vector<unsigned long long>* track_cluskeys = nullptr;
  t_residual->SetBranchAddress("cluskeys", &track_cluskeys);

  std::set<unsigned long long> clusters_in_tracks;
  for (int i = 0; i < t_residual->GetEntries(); i++) {
    t_residual->GetEntry(i);
    for (auto key : *track_cluskeys) {
      clusters_in_tracks.insert(key);
    }
  }

  // --------------------------
  // Read hittree from residual file
  // Copy hits into a separate tree unchanged
  // --------------------------
  int m_runnumber, m_segment, m_job, m_event;
  ULong64_t m_bco, m_bcotr;
  UInt_t m_hitsetkey;
  Float_t m_hitgx, m_hitgy, m_hitgz;
  Int_t m_hitlayer, m_sector, m_side, m_staveid, m_chipid;
  Int_t m_strobeid, m_ladderzid, m_ladderphiid, m_timebucket, m_hitpad, m_hittbin;
  Int_t m_col, m_row, m_segtype, m_tileid, m_strip;
  Float_t m_adc_hit, m_zdriftlength;

t_hits->SetBranchAddress("run", &m_runnumber);
t_hits->SetBranchAddress("segment", &m_segment);
t_hits->SetBranchAddress("job", &m_job);
t_hits->SetBranchAddress("event", &m_event);
t_hits->SetBranchAddress("gl1bco", &m_bco);
t_hits->SetBranchAddress("trbco", &m_bcotr);
t_hits->SetBranchAddress("hitsetkey", &m_hitsetkey);
t_hits->SetBranchAddress("gx", &m_hitgx);
t_hits->SetBranchAddress("gy", &m_hitgy);
t_hits->SetBranchAddress("gz", &m_hitgz);
t_hits->SetBranchAddress("layer", &m_hitlayer);
t_hits->SetBranchAddress("sector", &m_sector);
t_hits->SetBranchAddress("side", &m_side);
t_hits->SetBranchAddress("stave", &m_staveid);
t_hits->SetBranchAddress("chip", &m_chipid);
t_hits->SetBranchAddress("strobe", &m_strobeid);
t_hits->SetBranchAddress("ladderz", &m_ladderzid);
t_hits->SetBranchAddress("ladderphi", &m_ladderphiid);
t_hits->SetBranchAddress("timebucket", &m_timebucket);
t_hits->SetBranchAddress("pad", &m_hitpad);
t_hits->SetBranchAddress("tbin", &m_hittbin);
t_hits->SetBranchAddress("col", &m_col);
t_hits->SetBranchAddress("row", &m_row);
t_hits->SetBranchAddress("segtype", &m_segtype);
t_hits->SetBranchAddress("tile", &m_tileid);
t_hits->SetBranchAddress("strip", &m_strip);
t_hits->SetBranchAddress("adc", &m_adc_hit);
t_hits->SetBranchAddress("zdriftlength", &m_zdriftlength);

  // --------------------------
  // Read clustertree from residual file
  // This has all the cluster info you want to keep
  // --------------------------
  unsigned long m_scluskey;
  Float_t m_scluslx, m_scluslz, m_sclusgx, m_sclusgy, m_sclusgz, m_sclusgr;
  Float_t m_sclusphi, m_scluseta, m_adc_clus, m_scluselx, m_scluselz, m_clusmaxadc;
  Int_t m_scluslayer, m_phisize, m_zsize, m_clussector;
  
t_residual_clus->SetBranchAddress("cluskey", &m_scluskey);
t_residual_clus->SetBranchAddress("run", &m_runnumber);
t_residual_clus->SetBranchAddress("segment", &m_segment);
t_residual_clus->SetBranchAddress("job", &m_job);
t_residual_clus->SetBranchAddress("event", &m_event);
t_residual_clus->SetBranchAddress("gl1bco", &m_bco);
t_residual_clus->SetBranchAddress("trbco", &m_bcotr);
t_residual_clus->SetBranchAddress("lx", &m_scluslx);
t_residual_clus->SetBranchAddress("lz", &m_scluslz);
t_residual_clus->SetBranchAddress("gx", &m_sclusgx);
t_residual_clus->SetBranchAddress("gy", &m_sclusgy);
t_residual_clus->SetBranchAddress("gz", &m_sclusgz);
t_residual_clus->SetBranchAddress("r", &m_sclusgr);
t_residual_clus->SetBranchAddress("phi", &m_sclusphi);
t_residual_clus->SetBranchAddress("eta", &m_scluseta);
t_residual_clus->SetBranchAddress("adc", &m_adc_clus);
t_residual_clus->SetBranchAddress("phisize", &m_phisize);
t_residual_clus->SetBranchAddress("zsize", &m_zsize);
t_residual_clus->SetBranchAddress("layer", &m_scluslayer);
t_residual_clus->SetBranchAddress("erphi", &m_scluselx);
t_residual_clus->SetBranchAddress("ez", &m_scluselz);
t_residual_clus->SetBranchAddress("maxadc", &m_clusmaxadc);
t_residual_clus->SetBranchAddress("sector", &m_clussector);
t_residual_clus->SetBranchAddress("side", &m_side);
t_residual_clus->SetBranchAddress("stave", &m_staveid);
t_residual_clus->SetBranchAddress("chip", &m_chipid);
t_residual_clus->SetBranchAddress("strobe", &m_strobeid);
t_residual_clus->SetBranchAddress("ladderz", &m_ladderzid);
t_residual_clus->SetBranchAddress("ladderphi", &m_ladderphiid);
t_residual_clus->SetBranchAddress("timebucket", &m_timebucket);
t_residual_clus->SetBranchAddress("segtype", &m_segtype);
t_residual_clus->SetBranchAddress("tile", &m_tileid);

  // --------------------------
  // Create output file and trees
  // --------------------------
  TFile* outFile = new TFile("separate_cluster_and_hits_1216.root", "RECREATE");

  // Combined cluster tree
  // Keep all original cluster info plus used_in_seed and used_in_track
  unsigned long out_cluskey;
  int out_used_in_seed=0;
  int out_used_in_track=0;
  float out_x=0, out_y=0, out_z=0; // (optional) from seeding_clusters

  TTree* t_combined_clusters = new TTree("combined_clusters", "Combined cluster info");
  t_combined_clusters->Branch("cluskey", &out_cluskey, "cluskey/l");
  t_combined_clusters->Branch("run", &m_runnumber, "run/I");
  t_combined_clusters->Branch("segment", &m_segment, "segment/I");
  t_combined_clusters->Branch("job", &m_job, "job/I");
  t_combined_clusters->Branch("event", &m_event, "event/I");
  t_combined_clusters->Branch("gl1bco", &m_bco, "gl1bco/l");
  t_combined_clusters->Branch("trbco", &m_bcotr, "trbco/l");
  t_combined_clusters->Branch("lx", &m_scluslx, "lx/F");
  t_combined_clusters->Branch("lz", &m_scluslz, "lz/F");
  t_combined_clusters->Branch("gx", &m_sclusgx, "gx/F");
  t_combined_clusters->Branch("gy", &m_sclusgy, "gy/F");
  t_combined_clusters->Branch("gz", &m_sclusgz, "gz/F");
  t_combined_clusters->Branch("r", &m_sclusgr, "r/F");
  t_combined_clusters->Branch("phi", &m_sclusphi, "phi/F");
  t_combined_clusters->Branch("eta", &m_scluseta, "eta/F");
  t_combined_clusters->Branch("adc", &m_adc_clus, "adc/F");
  t_combined_clusters->Branch("phisize", &m_phisize, "phisize/I");
  t_combined_clusters->Branch("zsize", &m_zsize, "zsize/I");
  t_combined_clusters->Branch("layer", &m_scluslayer, "layer/I");
  t_combined_clusters->Branch("erphi", &m_scluselx, "erphi/F");
  t_combined_clusters->Branch("ez", &m_scluselz, "ez/F");
  t_combined_clusters->Branch("maxadc", &m_clusmaxadc, "maxadc/F");
  t_combined_clusters->Branch("sector", &m_clussector, "sector/I");
  t_combined_clusters->Branch("side", &m_side, "side/I");
  t_combined_clusters->Branch("stave", &m_staveid, "stave/I");
  t_combined_clusters->Branch("chip", &m_chipid, "chip/I");
  t_combined_clusters->Branch("strobe", &m_strobeid, "strobe/I");
  t_combined_clusters->Branch("ladderz", &m_ladderzid, "ladderz/I");
  t_combined_clusters->Branch("ladderphi", &m_ladderphiid, "ladderphi/I");
  t_combined_clusters->Branch("timebucket", &m_timebucket, "timebucket/I");
  t_combined_clusters->Branch("segtype", &m_segtype, "segtype/I");
  t_combined_clusters->Branch("tile", &m_tileid, "tile/I");

  // Additional branches from tracking_clusters
  t_combined_clusters->Branch("used_in_seed", &out_used_in_seed, "used_in_seed/I");
  t_combined_clusters->Branch("used_in_track", &out_used_in_track, "used_in_track/I");
  t_combined_clusters->Branch("seed_x", &out_x, "seed_x/F");
  t_combined_clusters->Branch("seed_y", &out_y, "seed_y/F");
  t_combined_clusters->Branch("seed_z", &out_z, "seed_z/F");

  // Combined hits tree
  TTree* t_combined_hits = new TTree("combined_hits", "Hit-level info");
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
  t_combined_hits->Branch("adc", &m_adc_hit, "adc/F");
  t_combined_hits->Branch("zdriftlength", &m_zdriftlength, "zdriftlength/F");

  // Fill combined_clusters
  for (int i = 0; i < t_residual_clus->GetEntries(); i++) {
    t_residual_clus->GetEntry(i);
    out_cluskey = m_scluskey;

    // Check tracking info
    if (seeding_map.find((unsigned long long)m_scluskey) != seeding_map.end()) {
    out_used_in_seed = seeding_map[(unsigned long long)m_scluskey].used_in_seed;
    out_x = seeding_map[(unsigned long long)m_scluskey].x;
    out_y = seeding_map[(unsigned long long)m_scluskey].y;
    out_z = seeding_map[(unsigned long long)m_scluskey].z;
} else {
      out_used_in_seed = 0;
      out_x = out_y = out_z = 0;
    }

    // Check if used in track
    out_used_in_track = (clusters_in_tracks.find((unsigned long long)m_scluskey) != clusters_in_tracks.end()) ? 1 : 0;

    t_combined_clusters->Fill();
  }

  // Fill combined_hits by looping over original hittree again
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
  f_seed->Close();

  std::cout << "Wrote separate_cluster_and_hits.root with combined_clusters and combined_hits trees." << std::endl;
}
