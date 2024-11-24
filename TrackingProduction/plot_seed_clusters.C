void plot_seed_clusters()
{
    TFile *f = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff_1120_0.root");
    TTree *t = (TTree *)f->Get("tracking_clusters");

    TFile *outfile = new TFile("failed_straight_line_cluster_plots_1120.root", "RECREATE");

    // Create graphs
    TGraph2D *g_passed_straight = new TGraph2D(); // For clusters that passed straight line condition
    TGraph2D *g_failed_straight = new TGraph2D(); // For clusters that failed straight line condition

    // Get the data
    float x, y, z;
    // int layer, is_rejected;
    int is_passed_straight;
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("z", &z);
    // t->SetBranchAddress("layer", &layer);
    t->SetBranchAddress("passed_straight", &is_passed_straight);

    // Fill the graphs
    int n_passed = 0, n_failed = 0;
    for (int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        if (is_passed_straight)
        {
            g_passed_straight->SetPoint(n_passed++, x, y, z);
        }
        else
        {
            g_failed_straight->SetPoint(n_failed++, x, y, z);
        }
    }

    // Create 3D canvas
    TCanvas *c = new TCanvas("c", "Straight Line Condition Clusters 3D", 1200, 1200);
    c->cd();

    TH3D *frame = new TH3D("frame", "Cluster Positions;X;Y;Z",
                           100, -150, 150,
                           100, -150, 150,
                           100, -300, 300);

    frame->GetXaxis()->SetTitle("X");
    frame->GetYaxis()->SetTitle("Y");
    frame->GetZaxis()->SetTitle("Z");

    g_passed_straight->SetMarkerColor(kBlue);
    g_passed_straight->SetMarkerStyle(20);
    g_passed_straight->SetMarkerSize(0.2);

    g_failed_straight->SetMarkerColor(kRed);
    g_failed_straight->SetMarkerStyle(20);
    g_failed_straight->SetMarkerSize(0.2);

    // Draw
    frame->Draw();

    frame->GetXaxis()->SetLimits(-90, 90);
    frame->GetYaxis()->SetLimits(-90, 90);
    frame->GetZaxis()->SetLimits(-280, 280);

    g_failed_straight->Draw("P SAME");

    //g_passed_straight->Draw("P SAME");

    c->SetTheta(30);
    c->SetPhi(30);

    outfile->cd();
    c->Write();
    f->Close();
}