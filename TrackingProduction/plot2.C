void plot2() {
    //TFile *file = TFile::Open("/sphenix/user/dcxchenxi/Beam_Residuals_2/output/outputFile_kso_41992_9_resid.root");
    //TFile *file = TFile::Open("/sphenix/user/dcxchenxi/Beam_Residuals_2/output/merged_2.root");
    TFile *file = TFile::Open("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/myKShortReco/inReconstruction/merged_51576.root");


    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    //TTree *tree = (TTree *) file->Get("residualtree;6");
    TTree *tree = (TTree *) file->Get("ntp_reco_info");

    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }
    //tree->Print();

    std::vector<float> *clusgx = nullptr;
    std::vector<float> *stategx = nullptr;
    float invariant_pt, eta1, eta2, dedx, quality, phi, invariant_mass, dca3dxy1, dca3dxy2, pt1, pt2;
    int nhits, charge;

    tree->SetBranchAddress("invariant_mass", &invariant_mass);
    //tree->SetBranchAddress("stategx", &stategx);
    tree->SetBranchAddress("invariant_pt", &invariant_pt);
    /*tree->SetBranchAddress("eta1", &eta1);
    tree->SetBranchAddress("dedx", &dedx);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("nhits", &nhits);
*/
    tree->SetBranchAddress("eta1", &eta1);
    tree->SetBranchAddress("eta2", &eta2);
    tree->SetBranchAddress("pt1", &pt1);
    tree->SetBranchAddress("pt2", &pt2);
    tree->SetBranchAddress("dca3dxy1", &dca3dxy1);
    tree->SetBranchAddress("dca3dxy2", &dca3dxy2);



/*        tree->SetBranchAddress("rzint", &rzint);
    tree->SetBranchAddress("rzslope", &rzslope);*/
    TH1F *hist1 = new TH1F("h_invariant_mass", "invariant_mass distribution", 200, 0, 3);//0-60
    tree->Draw("invariant_mass>>h_invariant_mass");

    TH1F *hist2 = new TH1F("h_invariant_pt", "invariant_pt distribution", 200, 0, 5);
    tree->Draw("invariant_pt>>h_invariant_pt");

        TH1F *hist3 = new TH1F("h_dca3dxy1", "dca3dxy1", 200, 0, 20);
    tree->Draw("dca3dxy1>>h_dca3dxy1");

          TH1F *hist4 = new TH1F("h_dca3dxy2", "dca3dxy2", 200, 0, 20);
    tree->Draw("dca3dxy2>>h_dca3dxy2");

    TH1F *hist5 = new TH1F("h_pt1", "pt1", 200, 0, 20);
    tree->Draw("pt1>>h_pt1");

    TH1F *hist6 = new TH1F("h_pt2", "pt2", 200, 0, 20);
    tree->Draw("pt2>>h_pt2");

    TH1F *hist7 = new TH1F("h_eta1", "eta1", 200, -5, 5);
    tree->Draw("eta1>>h_eta1");

    TH1F *hist8 = new TH1F("h_eta2", "eta2", 200, -5, 5);
    tree->Draw("eta2>>h_eta2");

    /*TH1F *hist3 = new TH1F("h_dedx", "dedx_distribution", 200, 0, 12000);
    //tree->Draw("dedx>>h_dedx");

    TH1F *hist4 = new TH1F("h_charge", "charge_distribution", 200, -2, 2);
    //tree->Draw("charge>>h_charge");

    TH1F *hist5 = new TH1F("h_quality", "quality_distribution", 200, 0, 150);
    //tree->Draw("quality>>h_quality");

    TH1F *hist6 = new TH1F("h_nhits", "nhits_distribution", 200, 30, 60);
    //tree->Draw("nhits>>h_nhits");

    TH1F *hist7 = new TH1F("h_phi", "phi_distribution", 200, -3.5, 3.5);*/
    //tree->Draw("phi>>h_phi");

   /* TGraph *graph_dedx_pt = new TGraph();
    graph_dedx_pt->SetTitle("dE/dx vs. pT;pT (GeV/c);dE/dx");

    TGraph *graph_dedx_pt_charge_pos = new TGraph();
    //graph_dedx_pt_charge_pos->SetTitle("dE/dx vs. pT*charge;pT*charge (GeV/c);dE/dx");
    graph_dedx_pt_charge_pos->SetMarkerColor(kRed);
    //graph_dedx_pt_charge_pos->SetMarkerStyle(20);

    TGraph *graph_dedx_pt_charge_neg = new TGraph();
    //graph_dedx_pt_charge_neg->SetTitle("dE/dx vs. pT*charge;pT*charge (GeV/c);dE/dx");
    graph_dedx_pt_charge_neg->SetMarkerColor(kBlue);
    //graph_dedx_pt_charge_neg->SetMarkerStyle(20);*/


   /* Long64_t nentries = tree->GetEntries();
    std::cout<<"Number of entries:"<<nentries<<std::endl;*/

    /*for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        bool all_points_satisfy = true;
        for (size_t j = 0; j < clusgx->size(); j++) {
            float xresidual = clusgx->at(j) - stategx->at(j);
            if (fabs(xresidual) >= 0.015) {
                all_points_satisfy = false;
                break;
            }

            // && abs(-1*rzint) / sqrt((rzslope * rzslope) + 1) < 5
        }
        if (all_points_satisfy && nhits > 30 && quality < 100) {
        //if (1) {
            // std::cout<<"value of pt:"<<pt<<std::endl;
            hist1->Fill(pt);
            hist2->Fill(eta);
            hist3->Fill(dedx);
            hist4->Fill(charge);
            hist5->Fill(quality);
            hist6->Fill(nhits);
            hist7->Fill(phi);
            //graph_dedx_pt->SetPoint(graph_dedx_pt->GetN(), pt, dedx);
        }
        if (all_points_satisfy && nhits > 30 && quality < 100 && pt < 2) {
            graph_dedx_pt->SetPoint(graph_dedx_pt->GetN(), pt, dedx);
            float pt_charge = pt * charge;
            if (charge > 0) {
                graph_dedx_pt_charge_pos->SetPoint(graph_dedx_pt_charge_pos->GetN(), pt_charge, dedx);
            } else {
                graph_dedx_pt_charge_neg->SetPoint(graph_dedx_pt_charge_neg->GetN(), pt_charge, dedx);
            }
        }


        //std::cout<<rzint<<" * "<<rzslope<<std::endl;

        //hist6->Fill(abs(-1*rzint) / sqrt((rzslope * rzslope) + 1));
    }*/

   /* TFile *outputFile = new TFile("dedx_vs_pt_3.root", "RECREATE");
    graph_dedx_pt_charge_pos->Write("dedx_vs_pt_charge_pos");
    graph_dedx_pt_charge_neg->Write("dedx_vs_pt_charge_neg");
    outputFile->Close();*/

    //const char* cuts = "nhits > 30 && abs(clusgx - stategx) < 0.015 && abs(rzint)/sqrt((rzslope*rzslope) + 1) < 5";


    TCanvas *c1 = new TCanvas("c1");


    c1->Divide(1, 2);

    c1->cd(1);
    gPad->SetLogy();
    hist1->Draw();

    c1->cd(2);
    gPad->SetLogy();
    hist2->Draw();
    c1->SaveAs("/sphenix/user/dcxchenxi/files/dist_pt_mass_09_16.pdf");


    /*c1->cd(3);
    gPad->SetLogy();
    hist3->Draw();

    c1->cd(4);
    //gPad->SetLogy();
    hist4->Draw();

    c1->cd(5);
    gPad->SetLogy();
    hist5->Draw();

    c1->cd(6);
    gPad->SetLogy();
    hist6->Draw();*/


    /*TCanvas *c2 = new TCanvas("c2");
    //c2->SetLogy();

    // 创建TMultiGraph来叠加两个TGraph
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(graph_dedx_pt_charge_pos);
    mg->Add(graph_dedx_pt_charge_neg);

    // 绘制叠加的图形
    mg->SetTitle("dE/dx vs. pT*charge;pT*charge (GeV/c);dE/dx");
    mg->Draw("AP");
    //c2->BuildLegend();

    // 设置mg的x和y轴范围
    mg->GetXaxis()->SetLimits(-0.5, 0.5); // 设置x轴范围
    mg->GetYaxis()->SetRangeUser(0, 5000); // 设置y轴范围

    c2->SaveAs("~/Downloads/dedx_vs_pt_charge5.pdf");*/
    //c2->SaveAs("~/Downloads/distribution1.pdf");


    TCanvas *c2 = new TCanvas("c2");


    c2->Divide(2, 3);

    c2->cd(1);
    gPad->SetLogy();
    hist3->Draw();

    c2->cd(2);
    gPad->SetLogy();
    hist4->Draw();

    c2->cd(3);
    gPad->SetLogy();
    hist5->Draw();

    c2->cd(4);
    gPad->SetLogy();
    hist6->Draw();

    c2->cd(5);
    gPad->SetLogy();
    hist7->Draw();

    c2->cd(6);
    gPad->SetLogy();
    hist8->Draw();
    c2->SaveAs("/sphenix/user/dcxchenxi/files/dist_09_16.pdf");


  /*  TCanvas *c3 = new TCanvas("c3");
    c3->SetLogy();
    graph_dedx_pt->Draw("AP");
    c3->SaveAs("~/Downloads/dedx_vs_pt_total_after_cuts5.pdf");*/

    //file->Close();
    /* h_pt->GetXaxis()->SetRangeUser(0,900);
     h_pt->GetYaxis()->SetRangeUser(0,1);*/

}





