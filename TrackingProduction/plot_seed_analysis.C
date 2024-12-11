void plot_seed_analysis()
{
    TFile* f = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff_1104_1.root");
    TTree* t = (TTree*)f->Get("seed_analysis");

    TFile* outfile = new TFile("seed_analysis_plots_1104_30_1.root", "RECREATE");





    float x, y, z;
    int layer, is_rejected, seed_id;
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("z", &z);
    t->SetBranchAddress("layer", &layer);
    t->SetBranchAddress("is_rejected", &is_rejected);
    t->SetBranchAddress("seed_id", &seed_id);

        int max_seed_id = -1;
    for(int i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        if(seed_id > max_seed_id) max_seed_id = seed_id;
    }


    const int MAX_SEEDS = 30;
    //const int MAX_SEEDS = max_seed_id+1;
        // First pass: count how many used clusters in each seed
    std::vector<int> clusters_per_seed(MAX_SEEDS, 0);

    for(int i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        if(seed_id >= MAX_SEEDS) continue;
        if(!is_rejected) {
            clusters_per_seed[seed_id]++;
        }
    }

    // Create graphs for long and short seeds
    TGraph2D* g_long_seeds = new TGraph2D();     // seeds with 6+ clusters (red)
    TGraph2D* g_short_seeds = new TGraph2D();    // seeds with <6 clusters (blue)
    TGraph2D* g_rejected = new TGraph2D();       // rejected clusters (purple)

   /* std::vector<TGraph2D*> g_seeds(MAX_SEEDS);     // used clusters for each seed
    std::vector<TGraph2D*> g_rejected(MAX_SEEDS);  // rejected clusters for each seed*/

    /*// Initialize graphs
    for(int i = 0; i < MAX_SEEDS; i++) {
        g_seeds[i] = new TGraph2D();
        g_rejected[i] = new TGraph2D();

        // Different colors for each seed
        g_seeds[i]->SetMarkerColor(kBlue );
        g_seeds[i]->SetMarkerStyle(20);
        g_seeds[i]->SetMarkerSize(0.2);

        g_rejected[i]->SetMarkerColor(kRed );
        g_rejected[i]->SetMarkerStyle(20);
        g_rejected[i]->SetMarkerSize(0.5);
    }*/

   /* // Get the data
    float x, y, z;
    int layer, is_rejected, seed_id;
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("z", &z);
    t->SetBranchAddress("layer", &layer);
    t->SetBranchAddress("is_rejected", &is_rejected);
    t->SetBranchAddress("seed_id", &seed_id);

    // Arrays to count points for each seed
    std::vector<int> n_seeds(MAX_SEEDS, 0);
    std::vector<int> n_rejected(MAX_SEEDS, 0);*/

/*    // Fill the graphs
    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        if(seed_id >= MAX_SEEDS) continue;  // Skip seeds beyond our limit

        if(is_rejected)
        {
            g_rejected[seed_id]->SetPoint(n_rejected[seed_id]++, x, y, z);
        }
        else
        {
            g_seeds[seed_id]->SetPoint(n_seeds[seed_id]++, x, y, z);
        }
    }*/

    // Style the graphs
    g_long_seeds->SetMarkerColor(kBlue);
    g_long_seeds->SetMarkerStyle(20);
    g_long_seeds->SetMarkerSize(0.2);

    g_short_seeds->SetMarkerColor(kGreen+2);
    g_short_seeds->SetMarkerStyle(20);
    g_short_seeds->SetMarkerSize(0.2);

    g_rejected->SetMarkerColor(kRed);
    g_rejected->SetMarkerStyle(20);
    g_rejected->SetMarkerSize(0.5);

    // Second pass: fill graphs
    int n_long = 0, n_short = 0, n_rejected = 0;
    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        if(seed_id >= MAX_SEEDS) continue;

        if(is_rejected) {
            g_rejected->SetPoint(n_rejected++, x, y, z);
        }
        else {
            if(clusters_per_seed[seed_id] >= 6) {
                g_long_seeds->SetPoint(n_long++, x, y, z);
            } else {
                g_short_seeds->SetPoint(n_short++, x, y, z);
            }
        }
    }


/*
  TH2D* h = new TH2D("h", "Cluster Positions;X;Y;Z",
                   100, -80, 80,   // X bins and range
                   100, -80, 80);  // Y bins and range
h->SetMinimum(-250);  // Set Z minimum
h->SetMaximum(250);   // Set Z maximum
*/

/*// Set histogram for each graph
g_short_seeds->SetHistogram((TH2D*)h->Clone());
g_long_seeds->SetHistogram((TH2D*)h->Clone());
g_rejected->SetHistogram((TH2D*)h->Clone());*/
    // Create 3D canvas
    TCanvas* c = new TCanvas("c", "Seeds and Rejected Clusters 3D", 1200, 1200);
    c->cd();

    TH3D* frame = new TH3D("frame", "Cluster Positions;X;Y;Z",
                           100, -80, 80,
                           100, -80, 80,
                           100, -250, 250);



       /* frame->GetXaxis()->SetRangeUser(-80, 80);
    frame->GetYaxis()->SetRangeUser(-80, 80);
    frame->GetZaxis()->SetRangeUser(-250, 250);*/
    frame->GetXaxis()->SetTitle("X");
    frame->GetYaxis()->SetTitle("Y");
    frame->GetZaxis()->SetTitle("Z");
    frame->SetStats(0);


    // Draw
    frame->Draw();


    /*frame->GetXaxis()->SetLimits(-80,80);
    frame->GetYaxis()->SetLimits(-80,80);
    frame->GetZaxis()->SetRangeUser(-250,250);*/

    /*g_short_seeds->GetHistogram()->GetXaxis()->SetLimits(-80,80);
    g_short_seeds->GetHistogram()->GetYaxis()->SetLimits(-80,80);
    g_short_seeds->GetHistogram()->GetZaxis()->SetRangeUser(-250,250);*/

     g_short_seeds->GetHistogram()->GetXaxis()->SetLimits(-80,80);
    g_short_seeds->GetHistogram()->GetYaxis()->SetLimits(-80,80);
    g_short_seeds->GetHistogram()->GetZaxis()->SetRangeUser(-250,250);
        g_short_seeds->Draw("P");




    //g_short_seeds->Draw("P");
if(n_long > 0) g_long_seeds->Draw("P SAME");
    if(n_rejected > 0) g_rejected->Draw("P SAME");




/*g_seeds[0]->Draw("P");
    g_seeds[0]->GetHistogram()->GetXaxis()->SetLimits(-80,80);
    g_seeds[0]->GetHistogram()->GetYaxis()->SetLimits(-80,80);
    g_seeds[0]->GetHistogram()->GetZaxis()->SetRangeUser(-250,250);

    for(int i = 1; i < MAX_SEEDS; i++) {
        if(n_seeds[i] > 0) {
            g_seeds[i]->Draw("P SAME");
        }
    }

    // Draw rejected clusters
    for(int i = 0; i < MAX_SEEDS; i++) {
        if(n_rejected[i] > 0) {
            g_rejected[i]->Draw("P SAME");
        }
    }*/


 /*   // Add legend
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    for(int i = 0; i < MAX_SEEDS; i++) {
        if(n_seeds[i] > 0 || n_rejected[i] > 0) {
            leg->AddEntry(g_seeds[i], Form("Seed %d Used", i), "p");
            leg->AddEntry(g_rejected[i], Form("Seed %d Rejected", i), "p");
        }
    }
    leg->Draw();*/

/* TLegend* leg = new TLegend(0.7, 0.1, 0.9, 0.3);  // Parameters are x1,y1,x2,y2 in NDC coordinates
leg->AddEntry(g_long_seeds, "Blue: Seeds #geq 6 clusters", "p");
leg->AddEntry(g_short_seeds, "Green: Seeds < 6 clusters", "p");
leg->AddEntry(g_rejected, "Red: Rejected Clusters", "p");
leg->SetBorderSize(0);  // Optional: removes the border
leg->SetFillStyle(0);   // Optional: makes background transparent
leg->Draw();*/

    c->SetTheta(30);
    c->SetPhi(30);

    outfile->cd();
    c->Write();
    outfile->Close();
    f->Close();
}