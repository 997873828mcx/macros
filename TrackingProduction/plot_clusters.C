void plot_clusters()
{
    TFile* f = new TFile("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff.root");
    TTree* t = (TTree*)f->Get("clusters");

    TFile* outfile = new TFile("cluster_plots.root", "RECREATE");

    // Create separate graphs for used and unused clusters
    TGraph2D* g_used = new TGraph2D();
    TGraph2D* g_unused = new TGraph2D();

    // Get the data
    float x, y, z;
    int is_used;
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("z", &z);
    t->SetBranchAddress("is_used", &is_used);

    // Fill the graphs
    int n_used = 0, n_unused = 0;
    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        if(is_used)
        {
            g_used->SetPoint(n_used++, x, y, z);
        }
        else
        {
            g_unused->SetPoint(n_unused++, x, y, z);
        }
    }

    // Create 3D canvas
    TCanvas* c = new TCanvas("c", "Clusters 3D", 1200, 1200);
c->cd();

    /*TView* view = TView::CreateView(1);
    view->SetRange(-110, -110, -110, 110, 110, 110);*/
    TH3D* frame = new TH3D("frame", "Cluster Positions;X;Y;Z",
                           100, -110, 110,  // x bins and range
                           100, -110, 110,  // y bins and range
                           100, -110, 110); // z bins and range

    frame->GetXaxis()->SetTitle("X");
    frame->GetYaxis()->SetTitle("Y");
    frame->GetZaxis()->SetTitle("Z");

    //frame->SetStats(0);

    // Style the graphs
    g_used->SetMarkerColor(kBlue);
    g_used->SetMarkerStyle(20);
    g_used->SetMarkerSize(0.2);

    //g_used->SetTitle("Cluster Positions;X;Y;Z");

    g_unused->SetMarkerColor(kRed);
    g_unused->SetMarkerStyle(20);
    g_unused->SetMarkerSize(0.2);

    // Draw

    frame->Draw();
    g_used->Draw("P");
    g_unused->Draw("PSAME");

    /*// Add legend
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(g_used, "Used Clusters", "p");
    leg->AddEntry(g_unused, "Unused Clusters", "p");
    leg->Draw();*/
    c->SetTheta(30);
    c->SetPhi(30);



   /* c->SetTheta(30);
    c->SetPhi(30);*/
    outfile->cd();
   c->Write();
    outfile->Close();
    //f->Close();
}