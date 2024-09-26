void plot() {

    TFile *file = TFile::Open("/sphenix/user/dcxchenxi/Beam_Residuals_2/output_frawley/selected_data_frawley_ww.root");


    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    TTree *tree = (TTree *) file->Get("residualtree;3");

    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        file->Close();
        return;
    }
    //tree->Print();

    std::vector<float> *clusgx = nullptr;
    std::vector<float> *stategx = nullptr;
    std::vector<float> *clusgy = nullptr;
    std::vector<float> *stategy = nullptr;
    float pt, eta, dedx, quality, phi, pz;
    int nhits, charge, ntpc, nintt, nmaps;

    tree->SetBranchAddress("clusgx", &clusgx);
    tree->SetBranchAddress("stategx", &stategx);
    tree->SetBranchAddress("clusgy", &clusgy);
    tree->SetBranchAddress("stategy", &stategy);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("dedx", &dedx);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("nhits", &nhits);
    tree->SetBranchAddress("ntpc", &ntpc);
    tree->SetBranchAddress("nintt", &nintt);
    tree->SetBranchAddress("nmaps", &nmaps);

    double x[121];
    double y[121];
    int n = 0;
    x[n] = 0.01929;	y[n] = 85.13;	n++;
    x[n] = 0.01964;	y[n] = 83.26;	n++;
    x[n] = 0.02024;	y[n] = 79.2;	n++;
    x[n] = 0.02062;	y[n] = 77.03;	n++;
    x[n] = 0.02087;	y[n] = 74.92;	n++;
    x[n] = 0.02164;	y[n] = 70.86;	n++;
    x[n] = 0.02244;	y[n] = 67.03;	n++;
    x[n] = 0.02313;	y[n] = 64.12;	n++;
    x[n] = 0.02399;	y[n] = 60.31;	n++;
    x[n] = 0.02503;	y[n] = 56.74;	n++;
    x[n] = 0.02596;	y[n] = 53.37;	n++;
    x[n] = 0.02676;	y[n] = 50.77;	n++;
    x[n] = 0.02725;	y[n] = 49.38;	n++;
    x[n] = 0.02809;	y[n] = 47.23;	n++;
    x[n] = 0.02931;	y[n] = 44.18;	n++;
    x[n] = 0.03021;	y[n] = 42.02;	n++;
    x[n] = 0.03152;	y[n] = 39.53;	n++;
    x[n] = 0.03249;	y[n] = 37.6;	n++;
    x[n] = 0.0337;	y[n] = 35.57;	n++;
    x[n] = 0.03494;	y[n] = 33.46;	n++;
    x[n] = 0.03602;	y[n] = 32;	n++;
    x[n] = 0.03713;	y[n] = 30.27;	n++;
    x[n] = 0.03874;	y[n] = 28.32;	n++;
    x[n] = 0.03969;	y[n] = 27.24;	n++;
    x[n] = 0.04091;	y[n] = 26.05;	n++;
    x[n] = 0.04243;	y[n] = 24.65;	n++;
    x[n] = 0.04321;	y[n] = 23.97;	n++;
    x[n] = 0.04508;	y[n] = 22.55;	n++;
    x[n] = 0.04675;	y[n] = 21.21;	n++;
    x[n] = 0.04849;	y[n] = 19.95;	n++;
    x[n] = 0.04998;	y[n] = 19.08;	n++;
    x[n] = 0.04998;	y[n] = 19.08;	n++;
    x[n] = 0.05152;	y[n] = 18.25;	n++;
    x[n] = 0.05343;	y[n] = 17.27;	n++;
    x[n] = 0.05541;	y[n] = 16.42;	n++;
    x[n] = 0.05711;	y[n] = 15.71;	n++;
    x[n] = 0.05995;	y[n] = 14.61;	n++;
    x[n] = 0.06255;	y[n] = 13.82;	n++;
    x[n] = 0.06409;	y[n] = 13.3;	n++;
    x[n] = 0.06606;	y[n] = 12.72;	n++;
    x[n] = 0.0681;	y[n] = 12.23;	n++;
    x[n] = 0.07062;	y[n] = 11.64;	n++;
    x[n] = 0.0728;	y[n] = 11.13;	n++;
    x[n] = 0.07642;	y[n] = 10.41;	n++;
    x[n] = 0.07973;	y[n] = 9.849;	n++;
    x[n] = 0.08169;	y[n] = 9.526;	n++;
    x[n] = 0.0842;	y[n] = 9.213;	n++;
    x[n] = 0.08627;	y[n] = 8.911;	n++;
    x[n] = 0.08839;	y[n] = 8.667;	n++;
    x[n] = 0.09167;	y[n] = 8.336;	n++;
    x[n] = 0.09279;	y[n] = 8.153;	n++;
    x[n] = 0.09564;	y[n] = 7.885;	n++;
    x[n] = 0.09681;	y[n] = 7.712;	n++;
    x[n] = 0.1004;	y[n] = 7.459;	n++;
    x[n] = 0.1035;	y[n] = 7.214;	n++;
    x[n] = 0.108;	y[n] = 6.9;	n++;
    x[n] = 0.1113;	y[n] = 6.637;	n++;
    x[n] = 0.1168;	y[n] = 6.348;	n++;
    x[n] = 0.1234;	y[n] = 6.039;	n++;
    x[n] = 0.1287;	y[n] = 5.808;	n++;
    x[n] = 0.1343;	y[n] = 5.586;	n++;
    x[n] = 0.1385;	y[n] = 5.433;	n++;
    x[n] = 0.1427;	y[n] = 5.314;	n++;
    x[n] = 0.1471;	y[n] = 5.197;	n++;
    x[n] = 0.1535;	y[n] = 5.027;	n++;
    x[n] = 0.1592;	y[n] = 4.889;	n++;
    x[n] = 0.1651;	y[n] = 4.755;	n++;
    x[n] = 0.1691;	y[n] = 4.702;	n++;
    x[n] = 0.1765;	y[n] = 4.573;	n++;
    x[n] = 0.1841;	y[n] = 4.473;	n++;
    x[n] = 0.1945;	y[n] = 4.326;	n++;
    x[n] = 0.2079;	y[n] = 4.207;	n++;
    x[n] = 0.2209;	y[n] = 4.092;	n++;
    x[n] = 0.2333;	y[n] = 4.002;	n++;
    x[n] = 0.2494;	y[n] = 3.914;	n++;
    x[n] = 0.2698;	y[n] = 3.828;	n++;
    x[n] = 0.2973;	y[n] = 3.744;	n++;
    x[n] = 0.3198;	y[n] = 3.702;	n++;
    x[n] = 0.3481;	y[n] = 3.682;	n++;
    x[n] = 0.3859;	y[n] = 3.641;	n++;
    x[n] = 0.4151;	y[n] = 3.641;	n++;
    x[n] = 0.4464;	y[n] = 3.621;	n++;
    x[n] = 0.4772;	y[n] = 3.621;	n++;
    x[n] = 0.5164;	y[n] = 3.641;	n++;
    x[n] = 0.5724;	y[n] = 3.662;	n++;
    x[n] = 0.627;	y[n] = 3.682;	n++;
    x[n] = 0.7036;	y[n] = 3.723;	n++;
    x[n] = 0.7895;	y[n] = 3.765;	n++;
    x[n] = 0.8699;	y[n] = 3.786;	n++;
    x[n] = 1.116;	y[n] = 3.914;	n++;
    x[n] = 1.229;	y[n] = 3.958;	n++;
    x[n] = 1.379;	y[n] = 4.024;	n++;
    x[n] = 1.737;	y[n] = 4.138;	n++;
    x[n] = 2.122;	y[n] = 4.254;	n++;
    x[n] = 2.484;	y[n] = 4.35;	n++;
    x[n] = 2.839;	y[n] = 4.423;	n++;
    x[n] = 3.109;	y[n] = 4.473;	n++;
    x[n] = 3.447;	y[n] = 4.523;	n++;
    x[n] = 3.915;	y[n] = 4.599;	n++;
    x[n] = 4.501;	y[n] = 4.676;	n++;
    x[n] = 5.302;	y[n] = 4.781;	n++;
    x[n] = 6.322;	y[n] = 4.862;	n++;
    x[n] = 7.357;	y[n] = 4.943;	n++;
    x[n] = 8.408;	y[n] = 5.027;	n++;
    x[n] = 9.904;	y[n] = 5.111;	n++;
    x[n] = 11.88;	y[n] = 5.197;	n++;
    x[n] = 13.74;	y[n] = 5.255;	n++;
    x[n] = 16.38;	y[n] = 5.314;	n++;
    x[n] = 18.72;	y[n] = 5.373;	n++;
    x[n] = 21.27;	y[n] = 5.433;	n++;
    x[n] = 24.01;	y[n] = 5.464;	n++;
    x[n] = 26.46;	y[n] = 5.525;	n++;
    x[n] = 28.8;	y[n] = 5.556;	n++;
    x[n] = 45.95;	y[n] = 5.68;	n++;
    x[n] = 54.12;	y[n] = 5.744;	n++;
    x[n] = 27.44;	y[n] = 5.556;	n++;
    x[n] = 29.51;	y[n] = 5.556;	n++;
    x[n] = 33.52;	y[n] = 5.618;	n++;
    x[n] = 38.54;	y[n] = 5.649;	n++;
    x[n] = 43.24;	y[n] = 5.68;	n++;
    x[n] = 47.65;	y[n] = 5.712;	n++;
    for (int i = 0; i < n; i++) {
        //y[i] *= 1000;
    }

    TGraph *pion_graph = new TGraph(n, x, y);
    pion_graph->SetMarkerStyle(20);
    pion_graph->SetMarkerColor(kGreen);
    //xyscanGraph->SetLineColor(kRed);


    double kaonx[110];
    double kaony[110];
    n = 0;
    kaonx[n] = 0.06722; kaony[n] = 87.53; n++;
    kaonx[n] = 0.06928; kaony[n] = 82.8; n++;
    kaonx[n] = 0.07142; kaony[n] = 79.2; n++;
    kaonx[n] = 0.07361; kaony[n] = 75.75; n++;
    kaonx[n] = 0.07542; kaony[n] = 72.86; n++;
    kaonx[n] = 0.07821; kaony[n] = 68.92; n++;
    kaonx[n] = 0.08062; kaony[n] = 65.92; n++;
    kaonx[n] = 0.0836; kaony[n] = 62.01; n++;
    kaonx[n] = 0.08618; kaony[n] = 58.66; n++;
    kaonx[n] = 0.08937; kaony[n] = 55.49; n++;
    kaonx[n] = 0.09212; kaony[n] = 53.08; n++;
    kaonx[n] = 0.09495; kaony[n] = 50.49; n++;
    kaonx[n] = 0.09847; kaony[n] = 47.49; n++;
    kaonx[n] = 0.1027; kaony[n] = 45.17; n++;
    kaonx[n] = 0.1059; kaony[n] = 42.73; n++;
    kaonx[n] = 0.1098; kaony[n] = 40.2; n++;
    kaonx[n] = 0.1139; kaony[n] = 38.23; n++;
    kaonx[n] = 0.1174; kaony[n] = 36.17; n++;
    kaonx[n] = 0.121; kaony[n] = 34.4; n++;
    kaonx[n] = 0.124; kaony[n] = 33.09; n++;
    kaonx[n] = 0.1293; kaony[n] = 30.95; n++;
    kaonx[n] = 0.1349; kaony[n] = 28.96; n++;
    kaonx[n] = 0.1417; kaony[n] = 26.94; n++;
    kaonx[n] = 0.1514; kaony[n] = 24.37; n++;
    kaonx[n] = 0.157; kaony[n] = 22.8; n++;
    kaonx[n] = 0.1648; kaony[n] = 21.21; n++;
    kaonx[n] = 0.172; kaony[n] = 19.84; n++;
    kaonx[n] = 0.1805; kaony[n] = 18.46; n++;
    kaonx[n] = 0.1918; kaony[n] = 16.89; n++;
    kaonx[n] = 0.2001; kaony[n] = 15.89; n++;
    kaonx[n] = 0.2063; kaony[n] = 15.19; n++;
    kaonx[n] = 0.2139; kaony[n] = 14.53; n++;
    kaonx[n] = 0.2205; kaony[n] = 13.9; n++;
    kaonx[n] = 0.2328; kaony[n] = 12.79; n++;
    kaonx[n] = 0.2415; kaony[n] = 12.1; n++;
    kaonx[n] = 0.2519; kaony[n] = 11.32; n++;
    kaonx[n] = 0.2629; kaony[n] = 10.89; n++;
    kaonx[n] = 0.2742; kaony[n] = 10.3; n++;
    kaonx[n] = 0.2844; kaony[n] = 9.74; n++;
    kaonx[n] = 0.2967; kaony[n] = 9.265; n++;
    kaonx[n] = 0.325; kaony[n] = 8.244; n++;
    kaonx[n] = 0.339; kaony[n] = 7.842; n++;
    kaonx[n] = 0.3516; kaony[n] = 7.584; n++;
    kaonx[n] = 0.3646; kaony[n] = 7.254; n++;
    kaonx[n] = 0.3781; kaony[n] = 6.978; n++;
    kaonx[n] = 0.3945; kaony[n] = 6.637; n++;
    kaonx[n] = 0.4141; kaony[n] = 6.313; n++;
    kaonx[n] = 0.4373; kaony[n] = 6.005; n++;
    kaonx[n] = 0.4535; kaony[n] = 5.808; n++;
    kaonx[n] = 0.4732; kaony[n] = 5.586; n++;
    kaonx[n] = 0.5028; kaony[n] = 5.284; n++;
    kaonx[n] = 0.5407; kaony[n] = 5.027; n++;
    kaonx[n] = 0.5886; kaony[n] = 4.755; n++;
    kaonx[n] = 0.6446; kaony[n] = 4.523; n++;
    kaonx[n] = 0.7017; kaony[n] = 4.302; n++;
    kaonx[n] = 0.741; kaony[n] = 4.161; n++;
    kaonx[n] = 0.7825; kaony[n] = 4.092; n++;
    kaonx[n] = 0.8365; kaony[n] = 3.98; n++;
    kaonx[n] = 0.8365; kaony[n] = 3.98; n++;
    kaonx[n] = 0.8942; kaony[n] = 3.914; n++;
    kaonx[n] = 0.9733; kaony[n] = 3.807; n++;
    kaonx[n] = 1.053; kaony[n] = 3.765; n++;
    kaonx[n] = 1.139; kaony[n] = 3.723; n++;
    kaonx[n] = 1.218; kaony[n] = 3.682; n++;
    kaonx[n] = 1.318; kaony[n] = 3.662; n++;
    kaonx[n] = 1.443; kaony[n] = 3.641; n++;
    kaonx[n] = 1.6; kaony[n] = 3.621; n++;
    kaonx[n] = 1.752; kaony[n] = 3.641; n++;
    kaonx[n] = 1.884; kaony[n] = 3.641; n++;
    kaonx[n] = 2.076; kaony[n] = 3.641; n++;
    kaonx[n] = 2.219; kaony[n] = 3.662; n++;
    kaonx[n] = 2.387; kaony[n] = 3.682; n++;
    kaonx[n] = 2.63; kaony[n] = 3.744; n++;
    kaonx[n] = 2.863; kaony[n] = 3.786; n++;
    kaonx[n] = 3.116; kaony[n] = 3.807; n++;
    kaonx[n] = 3.351; kaony[n] = 3.828; n++;
    kaonx[n] = 3.715; kaony[n] = 3.871; n++;
    kaonx[n] = 4.093; kaony[n] = 3.914; n++;
    kaonx[n] = 4.483; kaony[n] = 3.98; n++;
    kaonx[n] = 4.821; kaony[n] = 4.024; n++;
    kaonx[n] = 5.312; kaony[n] = 4.069; n++;
    kaonx[n] = 5.783; kaony[n] = 4.115; n++;
    kaonx[n] = 6.333; kaony[n] = 4.184; n++;
    kaonx[n] = 7.106; kaony[n] = 4.254; n++;
    kaonx[n] = 7.735; kaony[n] = 4.302; n++;
    kaonx[n] = 8.627; kaony[n] = 4.35; n++;
    kaonx[n] = 9.563; kaony[n] = 4.399; n++;
    kaonx[n] = 10.8; kaony[n] = 4.473; n++;
    kaonx[n] = 11.82; kaony[n] = 4.548; n++;
    kaonx[n] = 12.71; kaony[n] = 4.573; n++;
    kaonx[n] = 13.59; kaony[n] = 4.624; n++;
    kaonx[n] = 14.88; kaony[n] = 4.676; n++;
    kaonx[n] = 16.6; kaony[n] = 4.728; n++;
    kaonx[n] = 17.96; kaony[n] = 4.781; n++;
    kaonx[n] = 19.91; kaony[n] = 4.835; n++;
    kaonx[n] = 21.54; kaony[n] = 4.889; n++;
    kaonx[n] = 23.17; kaony[n] = 4.943; n++;
    kaonx[n] = 24.77; kaony[n] = 4.971; n++;
    kaonx[n] = 26.31; kaony[n] = 5.027; n++;
    kaonx[n] = 28.82; kaony[n] = 5.055; n++;
    kaonx[n] = 31.18; kaony[n] = 5.083; n++;
    kaonx[n] = 33.53; kaony[n] = 5.14; n++;
    kaonx[n] = 35.63; kaony[n] = 5.168; n++;
    kaonx[n] = 38.32; kaony[n] = 5.197; n++;
    kaonx[n] = 40.96; kaony[n] = 5.226; n++;
    kaonx[n] = 44.58; kaony[n] = 5.284; n++;
    kaonx[n] = 47.66; kaony[n] = 5.343; n++;
    kaonx[n] = 50.63; kaony[n] = 5.343; n++;
    kaonx[n] = 54.13; kaony[n] = 5.373; n++;
    kaonx[n] = 57.51; kaony[n] = 5.403; n++;

    for (int i = 0; i < n; i++) {
       // kaony[i] *= 1000;
    }

    TGraph *kaon_graph = new TGraph(n, kaonx, kaony);
    kaon_graph->SetMarkerStyle(20);
    kaon_graph->SetMarkerColor(kRed);
    //00000000000

    double protonx[108];
double protony[108];
n = 0;
protonx[n] = 0.127;	protony[n] = 87.05;	n++;
protonx[n] = 0.1317;	protony[n] = 83.26;	n++;
protonx[n] = 0.1349;	protony[n] = 80.08;	n++;
protonx[n] = 0.1383;	protony[n] = 76.18;	n++;
protonx[n] = 0.1425;	protony[n] = 72.86;	n++;
protonx[n] = 0.1478;	protony[n] = 70.08;	n++;
protonx[n] = 0.1523;	protony[n] = 66.29;	n++;
protonx[n] = 0.158;	protony[n] = 62.71;	n++;
protonx[n] = 0.1628;	protony[n] = 58.99;	n++;
protonx[n] = 0.1679;	protony[n] = 56.42;	n++;
protonx[n] = 0.1751;	protony[n] = 53.08;	n++;
protonx[n] = 0.1805;	protony[n] = 50.49;	n++;
protonx[n] = 0.1883;	protony[n] = 47.49;	n++;
protonx[n] = 0.1953;	protony[n] = 44.67;	n++;
protonx[n] = 0.2025;	protony[n] = 42.02;	n++;
protonx[n] = 0.2113;	protony[n] = 39.53;	n++;
protonx[n] = 0.2218;	protony[n] = 36.57;	n++;
protonx[n] = 0.23;	protony[n] = 34.59;	n++;
protonx[n] = 0.2415;	protony[n] = 31.65;	n++;
protonx[n] = 0.2535;	protony[n] = 29.28;	n++;
protonx[n] = 0.2661;	protony[n] = 27.7;	n++;
protonx[n] = 0.2742;	protony[n] = 26.05;	n++;
protonx[n] = 0.2844;	protony[n] = 24.65;	n++;
protonx[n] = 0.2931;	protony[n] = 23.31;	n++;
protonx[n] = 0.3077;	protony[n] = 21.93;	n++;
protonx[n] = 0.3172;	protony[n] = 20.86;	n++;
protonx[n] = 0.3289;	protony[n] = 19.62;	n++;
protonx[n] = 0.3411;	protony[n] = 18.66;	n++;
protonx[n] = 0.3516;	protony[n] = 17.75;	n++;
protonx[n] = 0.3646;	protony[n] = 16.79;	n++;
protonx[n] = 0.3804;	protony[n] = 15.8;	n++;
protonx[n] = 0.3945;	protony[n] = 15.11;	n++;
protonx[n] = 0.4116;	protony[n] = 14.21;	n++;
protonx[n] = 0.4321;	protony[n] = 13.37;	n++;
protonx[n] = 0.4563;	protony[n] = 12.3;	n++;
protonx[n] = 0.4703;	protony[n] = 11.77;	n++;
protonx[n] = 0.4877;	protony[n] = 11.13;	n++;
protonx[n] = 0.5028;	protony[n] = 10.71;	n++;
protonx[n] = 0.5182;	protony[n] = 10.24;	n++;
protonx[n] = 0.5374;	protony[n] = 9.849;	n++;
protonx[n] = 0.5539;	protony[n] = 9.421;	n++;
protonx[n] = 0.571;	protony[n] = 9.061;	n++;
protonx[n] = 0.585;	protony[n] = 8.764;	n++;
protonx[n] = 0.6067;	protony[n] = 8.476;	n++;
protonx[n] = 0.633;	protony[n] = 8.018;	n++;
protonx[n] = 0.6644;	protony[n] = 7.542;	n++;
protonx[n] = 0.6974;	protony[n] = 7.214;	n++;
protonx[n] = 0.7276;	protony[n] = 6.9;	n++;
protonx[n] = 0.7546;	protony[n] = 6.637;	n++;
protonx[n] = 0.7825;	protony[n] = 6.384;	n++;
protonx[n] = 0.8164;	protony[n] = 6.106;	n++;
protonx[n] = 0.8467;	protony[n] = 5.906;	n++;
protonx[n] = 0.8942;	protony[n] = 5.649;	n++;
protonx[n] = 0.9386;	protony[n] = 5.403;	n++;
protonx[n] = 0.9972;	protony[n] = 5.14;	n++;
protonx[n] = 1.053;	protony[n] = 4.971;	n++;
protonx[n] = 1.099;	protony[n] = 4.781;	n++;
protonx[n] = 1.167;	protony[n] = 4.599;	n++;
protonx[n] = 1.248;	protony[n] = 4.448;	n++;
protonx[n] = 1.31;	protony[n] = 4.35;	n++;
protonx[n] = 1.318;	protony[n] = 4.326;	n++;
protonx[n] = 1.375;	protony[n] = 4.231;	n++;
protonx[n] = 1.443;	protony[n] = 4.138;	n++;
protonx[n] = 1.506;	protony[n] = 4.069;	n++;
protonx[n] = 1.619;	protony[n] = 3.958;	n++;
protonx[n] = 1.71;	protony[n] = 3.914;	n++;
protonx[n] = 1.817;	protony[n] = 3.828;	n++;
protonx[n] = 1.931;	protony[n] = 3.765;	n++;
protonx[n] = 2.051;	protony[n] = 3.723;	n++;
protonx[n] = 2.26;	protony[n] = 3.682;	n++;
protonx[n] = 2.445;	protony[n] = 3.662;	n++;
protonx[n] = 2.646;	protony[n] = 3.641;	n++;
protonx[n] = 2.863;	protony[n] = 3.641;	n++;
protonx[n] = 3.023;	protony[n] = 3.641;	n++;
protonx[n] = 3.251;	protony[n] = 3.641;	n++;
protonx[n] = 3.434;	protony[n] = 3.641;	n++;
protonx[n] = 3.648;	protony[n] = 3.641;	n++;
protonx[n] = 3.876;	protony[n] = 3.662;	n++;
protonx[n] = 4.118;	protony[n] = 3.682;	n++;
protonx[n] = 4.565;	protony[n] = 3.702;	n++;
protonx[n] = 5.154;	protony[n] = 3.765;	n++;
protonx[n] = 5.818;	protony[n] = 3.786;	n++;
protonx[n] = 6.371;	protony[n] = 3.828;	n++;
protonx[n] = 7.02;	protony[n] = 3.892;	n++;
protonx[n] = 7.782;	protony[n] = 3.936;	n++;
protonx[n] = 8.679;	protony[n] = 4.002;	n++;
protonx[n] = 9.621;	protony[n] = 4.069;	n++;
protonx[n] = 10.73;	protony[n] = 4.115;	n++;
protonx[n] = 12.11;	protony[n] = 4.184;	n++;
protonx[n] = 13.43;	protony[n] = 4.231;	n++;
protonx[n] = 15.07;	protony[n] = 4.302;	n++;
protonx[n] = 16.91;	protony[n] = 4.374;	n++;
protonx[n] = 17.32;	protony[n] = 4.399;	n++;
protonx[n] = 18.4;	protony[n] = 4.423;	n++;
protonx[n] = 20.03;	protony[n] = 4.473;	n++;
protonx[n] = 21.16;	protony[n] = 4.498;	n++;
protonx[n] = 22.21;	protony[n] = 4.548;	n++;
protonx[n] = 23.88;	protony[n] = 4.573;	n++;
protonx[n] = 26;	protony[n] = 4.624;	n++;
protonx[n] = 28.82;	protony[n] = 4.702;	n++;
protonx[n] = 31.56;	protony[n] = 4.755;	n++;
protonx[n] = 34.99;	protony[n] = 4.808;	n++;
protonx[n] = 38.08;	protony[n] = 4.862;	n++;
protonx[n] = 41.21;	protony[n] = 4.889;	n++;
protonx[n] = 44.85;	protony[n] = 4.943;	n++;
protonx[n] = 48.83;	protony[n] = 4.999;	n++;
protonx[n] = 52.83;	protony[n] = 5.055;	n++;
protonx[n] = 56.47;	protony[n] = 5.083;	n++;
    for (int i = 0; i < n; i++) {
       // protony[i] *= 1000;
    }

    TGraph *proton_graph = new TGraph(n, protonx, protony);
    proton_graph->SetMarkerStyle(20);
    proton_graph->SetMarkerColor(kBlue);









/*        tree->SetBranchAddress("rzint", &rzint);
    tree->SetBranchAddress("rzslope", &rzslope);*/
    /*TH1F *hist_pt_without_cuts = new TH1F("h_pt_without_cuts", "Pt distribution without cuts", 200, 0, 250);//0-60
    TH1F *hist_pt_quality_cuts = new TH1F("h_pt_quality_cuts", "Pt distribution with good quality", 200, 0, 60);
    TH1F *hist_pt_quality_acceptance_cuts = new TH1F("h_pt_quality_acceptance_cuts",
                                                     "Pt distribution with good quality and eta acceptance", 200, 0,
                                                     60);
    TH1F *hist_pt_quality_acceptance_nhits_cuts = new TH1F("h_pt_quality_acceptance_nhits_cuts",
                                                           "Pt distribution with good quality, eta acceptance, and nhits",
                                                           200, 0, 8);
    TH1F *hist_pt_quality_acceptance_nhits_residual_cuts = new TH1F("h_pt_quality_acceptance_nhits_residual_cuts",
                                                                    "Pt distribution with good quality, eta acceptance, nhits, and residual",
                                                                    200, 0, 8);*/

    TH1F *hist_pt_without_cuts = new TH1F("h_pt_without_cuts", "Pt distribution without cuts", 200, 0, 50);//0-60
    TH1F *hist_pt_pt_cuts = new TH1F("h_pt_pt_cuts", "Pt distribution with pt > 0.2", 200, 0, 25);
    TH1F *hist_pt_pt_quality_cuts = new TH1F("h_pt_pt_quality_cuts",
                                                     "Pt distribution with pt > 0.2, quality < 100", 200, 0,
                                                     20);
    TH1F *hist_pt_pt_quality_ntpc_cuts = new TH1F("h_pt_pt_quality_ntpc_cut",
                                                           "Pt distribution with pt > 0.2, quality < 100, ntpc > 20",
                                                           200, 0, 6);
    TH1F *hist_pt_pt_quality_ntpc_nmaps_cuts = new TH1F("h_pt_pt_quality_ntpc_nmaps_cuts",
                                                                    "Pt distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps>2",
                                                                    200, 0, 6);

    TH1F *hist_pt_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_pt_pt_quality_ntpc_nmaps_nintt_cuts",
                                                        "Pt distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps>2, nintt > 0",
                                                        200, 0, 4);



    //tree->Draw("pt>>h_pt");

    TH1F *hist_eta_without_cuts = new TH1F("h_eta_without_cuts", "Eta distribution without cuts", 200, -6, 6);
    TH1F *hist_eta_pt_cuts = new TH1F("h_eta_pt_cuts", "Eta distribution with pt > 0.2", 200, -3, 3);
    TH1F *hist_eta_pt_quality_cuts = new TH1F("h_eta_pt_quality_cuts", "Eta distribution with pt > 0.2, quality < 100", 200, -3, 3);
    TH1F *hist_eta_pt_quality_ntpc_cuts = new TH1F("h_eta_pt_quality_ntpc_cuts", "Eta distribution with pt > 0.2, quality < 100, ntpc > 20", 200, -3, 3);
    TH1F *hist_eta_pt_quality_ntpc_nmaps_cuts = new TH1F("h_eta_pt_quality_ntpc_nmaps_cuts", "Eta distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2", 200, -3, 3);
    TH1F *hist_eta_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_eta_pt_quality_ntpc_nmaps_nintt_cuts", "Eta distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0", 200, -3, 3);
    //tree->Draw("eta>>h_eta");

    /*TH1F *hist_dedx_without_cuts = new TH1F("h_dedx_without_cuts", "dedx distribution without cuts", 200, 0, 20000);
    TH1F *hist_dedx_quality_cuts = new TH1F("h_dedx_quality_cuts", "dedx distribution with good quality", 200, 0,
                                            12000);
    TH1F *hist_dedx_quality_acceptance_cuts = new TH1F("h_dedx_quality_acceptance_cuts",
                                                       "dedx distribution with good quality in TPC acceptance", 200, 0,
                                                       12000);
    TH1F *hist_dedx_quality_acceptance_nhits_cuts = new TH1F("h_dedx_quality_acceptance_nhits_cuts",
                                                             "dedx distribution with good quality, eta acceptance, and nhits",
                                                             200, 0, 7000);
    TH1F *hist_dedx_quality_acceptance_nhits_residual_cuts = new TH1F("h_dedx_quality_acceptance_nhits_residual_cuts",
                                                                      "dedx distribution with good quality, eta acceptance, nhits, and residual",
                                                                      200, 0, 7000);*/

    // dedx histograms
TH1F *hist_dedx_without_cuts = new TH1F("h_dedx_without_cuts", "Dedx distribution without cuts", 200, 0, 10000);
TH1F *hist_dedx_pt_cuts = new TH1F("h_dedx_pt_cuts", "Dedx distribution with pt > 0.2", 200, 0, 10000);
TH1F *hist_dedx_pt_quality_cuts = new TH1F("h_dedx_pt_quality_cuts", "Dedx distribution with pt > 0.2, quality < 100", 200, 0, 10000);
TH1F *hist_dedx_pt_quality_ntpc_cuts = new TH1F("h_dedx_pt_quality_ntpc_cuts", "Dedx distribution with pt > 0.2, quality < 100, ntpc > 20", 200, 0, 10000);
TH1F *hist_dedx_pt_quality_ntpc_nmaps_cuts = new TH1F("h_dedx_pt_quality_ntpc_nmaps_cuts", "Dedx distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2", 200, 0, 9000);
TH1F *hist_dedx_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_dedx_pt_quality_ntpc_nmaps_nintt_cuts", "Dedx distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0", 200, 0, 7000);

    //tree->Draw("dedx>>h_dedx");

   /* TH1F *hist_charge_without_cuts = new TH1F("h_charge_without_cuts", "charge distribution without cuts", 200, -2, 2);
    TH1F *hist_charge_quality_cuts = new TH1F("h_charge_quality_cuts", "charge distribution with good quality", 200, -2,
                                              2);
    TH1F *hist_charge_quality_acceptance_cuts = new TH1F("h_charge_quality_acceptance_cuts",
                                                         "charge distribution with good quality in TPC acceptance", 200,
                                                         -2, 2);
    TH1F *hist_charge_quality_acceptance_nhits_cuts = new TH1F("h_charge_quality_acceptance_nhits_cuts",
                                                               "charge distribution with good quality, eta acceptance, and nhits",
                                                               200, -2, 2);
    TH1F *hist_charge_quality_acceptance_nhits_residual_cuts = new TH1F(
            "h_charge_quality_acceptance_nhits_residual_cuts",
            "charge distribution with good quality, eta acceptance, nhits, and residual", 200, -2, 2);*/
    //tree->Draw("charge>>h_charge");

    // charge histograms
TH1F *hist_charge_without_cuts = new TH1F("h_charge_without_cuts", "Charge distribution without cuts", 200, -2, 2);
TH1F *hist_charge_pt_cuts = new TH1F("h_charge_pt_cuts", "Charge distribution with pt > 0.2", 200, -2, 2);
TH1F *hist_charge_pt_quality_cuts = new TH1F("h_charge_pt_quality_cuts", "Charge distribution with pt > 0.2, quality < 100", 200, -2, 2);
TH1F *hist_charge_pt_quality_ntpc_cuts = new TH1F("h_charge_pt_quality_ntpc_cuts", "Charge distribution with pt > 0.2, quality < 100, ntpc > 20", 200, -2, 2);
TH1F *hist_charge_pt_quality_ntpc_nmaps_cuts = new TH1F("h_charge_pt_quality_ntpc_nmaps_cuts", "Charge distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2", 200, -2, 2);
TH1F *hist_charge_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_charge_pt_quality_ntpc_nmaps_nintt_cuts", "Charge distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0", 200, -2, 2);

    //TH1F *hist5 = new TH1F("h_quality", "quality_distribution", 200, 0, 150);
    /*TH1F *hist_quality_without_cuts = new TH1F("h_quality_without_cuts", "quality distribution without cuts", 200, 0,
                                               20000);
    TH1F *hist_quality_quality_cuts = new TH1F("h_quality_quality_cuts", "good quality distribution", 200, 0, 150);
    TH1F *hist_quality_quality_acceptance_cuts = new TH1F("h_quality_quality_acceptance_cuts",
                                                          "good quality distribution in TPC acceptance", 200, 0, 150);
    TH1F *hist_quality_quality_acceptance_nhits_cuts = new TH1F("h_quality_quality_acceptance_nhits_cuts",
                                                                "quality distribution with good quality, eta acceptance, and nhits",
                                                                200, 0, 150);
    TH1F *hist_quality_quality_acceptance_nhits_residual_cuts = new TH1F(
            "h_quality_quality_acceptance_nhits_residual_cuts",
            "quality distribution with good quality, eta acceptance, nhits, and residual", 200, 0, 150);*/
    //tree->Draw("quality>>h_quality");

    // quality histograms
TH1F *hist_quality_without_cuts = new TH1F("h_quality_without_cuts", "Quality distribution without cuts", 200, 0, 20000);
TH1F *hist_quality_pt_cuts = new TH1F("h_quality_pt_cuts", "Quality distribution with pt > 0.2", 200, 0, 150);
TH1F *hist_quality_pt_quality_cuts = new TH1F("h_quality_pt_quality_cuts", "Quality distribution with pt > 0.2, quality < 100", 200, 0, 150);
TH1F *hist_quality_pt_quality_ntpc_cuts = new TH1F("h_quality_pt_quality_ntpc_cuts", "Quality distribution with pt > 0.2, quality < 100, ntpc > 20", 200, 0, 150);
TH1F *hist_quality_pt_quality_ntpc_nmaps_cuts = new TH1F("h_quality_pt_quality_ntpc_nmaps_cuts", "Quality distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2", 200, 0, 150);
TH1F *hist_quality_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_quality_pt_quality_ntpc_nmaps_nintt_cuts", "Quality distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0", 200, 0, 150);



    /*TH1F *hist_hits_without_cuts = new TH1F("h_nhits_without_cuts", "nhits distribution without cuts", 200, 0, 60);
    TH1F *hist_hits_quality_cuts = new TH1F("h_nhits_quality_cuts", "nhits distribution with good quality", 200, 0, 60);
    TH1F *hist_hits_quality_acceptance_cuts = new TH1F("h_nhits_quality_acceptance_cuts",
                                                       "nhits distribution with good quality in TPC acceptance", 200, 0,
                                                       60);
    TH1F *hist_hits_quality_acceptance_nhits_cuts = new TH1F("h_nhits_quality_acceptance_nhits_cuts",
                                                             "nhits distribution with good quality, eta acceptance, and nhits",
                                                             200, 0, 60);
    TH1F *hist_hits_quality_acceptance_nhits_residual_cuts = new TH1F("h_nhits_quality_acceptance_nhits_residual_cuts",
                                                                      "nhits distribution with good quality, eta acceptance, nhits, and residual",
                                                                      200, 0, 60);*/

    // nhits histograms
TH1F *hist_hits_without_cuts = new TH1F("h_nhits_without_cuts", "Nhits distribution without cuts", 200, 0, 60);
TH1F *hist_hits_pt_cuts = new TH1F("h_nhits_pt_cuts", "Nhits distribution with pt > 0.2", 200, 0, 60);
TH1F *hist_hits_pt_quality_cuts = new TH1F("h_nhits_pt_quality_cuts", "Nhits distribution with pt > 0.2, quality < 100", 200, 0, 60);
TH1F *hist_hits_pt_quality_ntpc_cuts = new TH1F("h_nhits_pt_quality_ntpc_cuts", "Nhits distribution with pt > 0.2, quality < 100, ntpc > 20", 200, 0, 60);
TH1F *hist_hits_pt_quality_ntpc_nmaps_cuts = new TH1F("h_nhits_pt_quality_ntpc_nmaps_cuts", "Nhits distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2", 200, 0, 60);
TH1F *hist_hits_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_nhits_pt_quality_ntpc_nmaps_nintt_cuts", "Nhits distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0", 200, 0, 60);

    //tree->Draw("nhits>>h_nhits");

    /*TH1F *hist_phi_without_cuts = new TH1F("h_phi_without_cuts", "phi distribution without cuts", 200, -3.5, 3.5);
    TH1F *hist_phi_quality_cuts = new TH1F("h_phi_quality_cuts", "phi distribution with good quality", 200, -3.5, 3.5);
    TH1F *hist_phi_quality_acceptance_cuts = new TH1F("h_phi_quality_acceptance_cuts",
                                                      "phi distribution with good quality in TPC acceptance", 200,
                                                      -3.5, 3.5);
    TH1F *hist_phi_quality_acceptance_nhits_cuts = new TH1F("h_phi_quality_acceptance_nhits_cuts",
                                                            "phi distribution with good quality, eta acceptance, and nhits",
                                                            200, -3.5, 3.5);
    TH1F *hist_phi_quality_acceptance_nhits_residual_cuts = new TH1F("h_phi_quality_acceptance_nhits_residual_cuts",
                                                                     "phi distribution with good quality, eta acceptance, nhits, and residual",
                                                                     200, -3.5, 3.5);*/
    //tree->Draw("phi>>h_phi");

    // phi histograms
TH1F *hist_phi_without_cuts = new TH1F("h_phi_without_cuts", "Phi distribution without cuts", 200, -3.5, 3.5);
TH1F *hist_phi_pt_cuts = new TH1F("h_phi_pt_cuts", "Phi distribution with pt > 0.2", 200, -3.5, 3.5);
TH1F *hist_phi_pt_quality_cuts = new TH1F("h_phi_pt_quality_cuts", "Phi distribution with pt > 0.2, quality < 100", 200, -3.5, 3.5);
TH1F *hist_phi_pt_quality_ntpc_cuts = new TH1F("h_phi_pt_quality_ntpc_cuts", "Phi distribution with pt > 0.2, quality < 100, ntpc > 20", 200, -3.5, 3.5);
TH1F *hist_phi_pt_quality_ntpc_nmaps_cuts = new TH1F("h_phi_pt_quality_ntpc_nmaps_cuts", "Phi distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2", 200, -3.5, 3.5);
TH1F *hist_phi_pt_quality_ntpc_nmaps_nintt_cuts = new TH1F("h_phi_pt_quality_ntpc_nmaps_nintt_cuts", "Phi distribution with pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0", 200, -3.5, 3.5);

    /*// 定义histograms
    TH2F *hist_dedx_pt_without_cuts = new TH2F("hist_dedx_pt_without_cuts", "dE/dx vs. pT without cuts; pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_quality_cuts = new TH2F("hist_dedx_pt_quality_cuts", "dE/dx vs. pT (quality cuts); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_quality_acceptance_cuts = new TH2F("hist_dedx_pt_quality_acceptance_cuts", "dE/dx vs. pT (quality + acceptance cuts); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_quality_acceptance_nhits_cuts = new TH2F("hist_dedx_pt_quality_acceptance_nhits_cuts", "dE/dx vs. pT (quality + acceptance + nhits cuts); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_quality_acceptance_nhits_residual_cuts = new TH2F("hist_dedx_pt_quality_acceptance_nhits_residual_cuts", "dE/dx vs. pT (quality + acceptance + nhits + residual cuts); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);

    TH2F *hist_dedx_pt_charge_without_cuts = new TH2F("hist_dedx_pt_charge_without_cuts", "dE/dx vs. pT*charge (without cuts); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_quality_cuts = new TH2F("hist_dedx_pt_charge_quality_cuts", "dE/dx vs. pT*charge (quality cuts); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_quality_acceptance_cuts = new TH2F("hist_dedx_pt_charge_quality_acceptance_cuts", "dE/dx vs. pT*charge (quality + acceptance cuts); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_quality_acceptance_nhits_cuts = new TH2F("hist_dedx_pt_charge_quality_acceptance_nhits_cuts", "dE/dx vs. pT*charge (quality + acceptance + nhits cuts); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_quality_acceptance_nhits_residual_cuts = new TH2F("hist_dedx_pt_charge_quality_acceptance_nhits_residual_cuts", "dE/dx vs. pT*charge (quality + acceptance + nhits + residual cuts); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);*/

    // 定义histograms for dE/dx vs. pT with gradual cuts
    TH2F *hist_dedx_pt_without_cuts = new TH2F("hist_dedx_pt_without_cuts", "dE/dx vs. pT without cuts; pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_pt_cuts = new TH2F("hist_dedx_pt_pt_cuts", "dE/dx vs. pT (pt > 0.2); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_pt_quality_cuts = new TH2F("hist_dedx_pt_pt_quality_cuts", "dE/dx vs. pT (pt > 0.2, quality < 100); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_pt_quality_ntpc_cuts = new TH2F("hist_dedx_pt_pt_quality_ntpc_cuts", "dE/dx vs. pT (pt > 0.2, quality < 100, ntpc > 20); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_pt_quality_ntpc_nmaps_cuts = new TH2F("hist_dedx_pt_pt_quality_ntpc_nmaps_cuts", "dE/dx vs. pT (pt > 0.2, quality < 100, ntpc > 20, nmaps > 2); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_pt_quality_ntpc_nmaps_nintt_cuts = new TH2F("hist_dedx_pt_pt_quality_ntpc_nmaps_nintt_cuts", "dE/dx vs. pT (pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0); pT (GeV/c); dE/dx", 200, 0, 2, 200, 0, 4000);

// 定义histograms for dE/dx vs. pT*charge with gradual cuts
    TH2F *hist_dedx_pt_charge_without_cuts = new TH2F("hist_dedx_pt_charge_without_cuts", "dE/dx vs. pT*charge without cuts; pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_pt_cuts = new TH2F("hist_dedx_pt_charge_pt_cuts", "dE/dx vs. pT*charge (pt > 0.2); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_pt_quality_cuts = new TH2F("hist_dedx_pt_charge_pt_quality_cuts", "dE/dx vs. pT*charge (pt > 0.2, quality < 100); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_pt_quality_ntpc_cuts = new TH2F("hist_dedx_pt_charge_pt_quality_ntpc_cuts", "dE/dx vs. pT*charge (pt > 0.2, quality < 100, ntpc > 20); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_pt_quality_ntpc_nmaps_cuts = new TH2F("hist_dedx_pt_charge_pt_quality_ntpc_nmaps_cuts", "dE/dx vs. pT*charge (pt > 0.2, quality < 100, ntpc > 20, nmaps > 2); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);
    TH2F *hist_dedx_pt_charge_pt_quality_ntpc_nmaps_nintt_cuts = new TH2F("hist_dedx_pt_charge_pt_quality_ntpc_nmaps_nintt_cuts", "dE/dx vs. pT*charge (pt > 0.2, quality < 100, ntpc > 20, nmaps > 2, nintt > 0); pT*charge (GeV/c); dE/dx", 200, -2, 2, 200, 0, 4000);

    TGraph *graph_dedx_p_without_cuts = new TGraph();
    graph_dedx_p_without_cuts->SetTitle("dE/dx vs. p;p (GeV/c);dE/dx");

    TGraph *graph_dedx_p_with_all_cuts = new TGraph();
    graph_dedx_p_with_all_cuts->SetTitle("dE/dx vs. p;p (GeV/c);dE/dx");

   /* TGraph *graph_dedx_pt_without_residual_cut = new TGraph();
    graph_dedx_pt_without_residual_cut->SetTitle("dE/dx vs. pT without residual cut;pT (GeV/c);dE/dx");

    TGraph *graph_dedx_pt_with_residual_cut = new TGraph();
    graph_dedx_pt_with_residual_cut->SetTitle("dE/dx vs. pT with residual cut;pT (GeV/c);dE/dx");*/

   //old one
    /*TGraph *graph_dedx_pt_charge_pos = new TGraph();
    //graph_dedx_pt_charge_pos->SetTitle("dE/dx vs. pT*charge;pT*charge (GeV/c);dE/dx");
    graph_dedx_pt_charge_pos->SetMarkerColor(kRed);
    //graph_dedx_pt_charge_pos->SetMarkerStyle(20);*/

   /* TGraph *graph_dedx_pt_charge_pos_without_residual_cut = new TGraph();
    graph_dedx_pt_charge_pos_without_residual_cut->SetMarkerColor(kRed);

    TGraph *graph_dedx_pt_charge_pos_with_residual_cut = new TGraph();
    graph_dedx_pt_charge_pos_with_residual_cut->SetMarkerColor(kRed);*/

/*    TGraph *graph_dedx_pt_charge_neg = new TGraph();
    //graph_dedx_pt_charge_neg->SetTitle("dE/dx vs. pT*charge;pT*charge (GeV/c);dE/dx");
    graph_dedx_pt_charge_neg->SetMarkerColor(kBlue);
    //graph_dedx_pt_charge_neg->SetMarkerStyle(20);*/

   /* TGraph *graph_dedx_pt_charge_neg_without_residual_cut = new TGraph();
    graph_dedx_pt_charge_neg_without_residual_cut->SetMarkerColor(kBlue);

    TGraph *graph_dedx_pt_charge_neg_with_residual_cut = new TGraph();
    graph_dedx_pt_charge_neg_with_residual_cut->SetMarkerColor(kBlue);*/


    Long64_t nentries = tree->GetEntries();
    std::cout << "Number of entries:" << nentries << std::endl;

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


/*    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        hist_pt_without_cuts->Fill(pt);
        hist_eta_without_cuts->Fill(eta);
        hist_dedx_without_cuts->Fill(dedx);
        hist_quality_without_cuts->Fill(quality);
        hist_hits_without_cuts->Fill(nhits);
        hist_phi_without_cuts->Fill(phi);
        hist_charge_without_cuts->Fill(charge);

        hist_dedx_pt_without_cuts->Fill(pt, dedx);
        hist_dedx_pt_charge_without_cuts->Fill(pt * charge, dedx);



        if (quality < 100) {
            //if (1) {
            // std::cout<<"value of pt:"<<pt<<std::endl;
            hist_pt_quality_cuts->Fill(pt);
            hist_eta_quality_cuts->Fill(eta);
            hist_dedx_quality_cuts->Fill(dedx);
            hist_charge_quality_cuts->Fill(charge);
            hist_quality_quality_cuts->Fill(quality);
            hist_hits_quality_cuts->Fill(nhits);
            hist_phi_quality_cuts->Fill(phi);

            hist_dedx_pt_quality_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_quality_cuts->Fill(pt * charge, dedx);
        }
        if (quality < 100 && eta > -1 && eta < 1) {
            hist_pt_quality_acceptance_cuts->Fill(pt);
            hist_eta_quality_acceptance_cuts->Fill(eta);
            hist_dedx_quality_acceptance_cuts->Fill(dedx);
            hist_charge_quality_acceptance_cuts->Fill(charge);
            hist_quality_quality_acceptance_cuts->Fill(quality);
            hist_hits_quality_acceptance_cuts->Fill(nhits);
            hist_phi_quality_acceptance_cuts->Fill(phi);

            hist_dedx_pt_quality_acceptance_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_quality_acceptance_cuts->Fill(pt * charge, dedx);

        }

        if (quality < 100 && eta > -1 && eta < 1 && nhits > 30) {
            //graph_dedx_pt_without_residual_cut->SetPoint(graph_dedx_pt_without_residual_cut->GetN(), pt, dedx);

            hist_dedx_pt_quality_acceptance_nhits_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_quality_acceptance_nhits_cuts->Fill(pt * charge, dedx);

            hist_pt_quality_acceptance_nhits_cuts->Fill(pt);
            hist_eta_quality_acceptance_nhits_cuts->Fill(eta);
            hist_dedx_quality_acceptance_nhits_cuts->Fill(dedx);
            hist_charge_quality_acceptance_nhits_cuts->Fill(charge);
            hist_quality_quality_acceptance_nhits_cuts->Fill(quality);
            hist_hits_quality_acceptance_nhits_cuts->Fill(nhits);
            hist_phi_quality_acceptance_nhits_cuts->Fill(phi);
            *//*float pt_charge = pt * charge;
            if (charge > 0) {
                graph_dedx_pt_charge_pos_without_residual_cut->SetPoint(
                        graph_dedx_pt_charge_pos_without_residual_cut->GetN(), pt_charge, dedx);
            } else {
                graph_dedx_pt_charge_neg_without_residual_cut->SetPoint(
                        graph_dedx_pt_charge_neg_without_residual_cut->GetN(), pt_charge, dedx);
            }*//*
        }
        if (quality < 100 && eta > -1 && eta < 1 && nhits > 30) {
            bool all_points_satisfy = true;
            for (size_t j = 0; j < clusgx->size(); j++) {
                // 计算全局坐标的r和phi
                float r_clus = std::sqrt(clusgx->at(j) * clusgx->at(j) + clusgy->at(j) * clusgy->at(j));
                float phi_clus = std::atan2(clusgy->at(j), clusgx->at(j));
                float x_local_clus = r_clus * phi_clus;

                float r_state = std::sqrt(stategx->at(j) * stategx->at(j) + stategy->at(j) * stategy->at(j));
                float phi_state = std::atan2(stategy->at(j), stategx->at(j));
                float x_local_state = r_state * phi_state;

                // 计算local x residual
                float x_local_residual = x_local_clus - x_local_state;

                // 检查residual是否满足条件
                if (std::fabs(x_local_residual) >= 1) {
                    all_points_satisfy = false;
                    break;
                }
            }

            if (all_points_satisfy) {
                //graph_dedx_pt_with_residual_cut->SetPoint(graph_dedx_pt_with_residual_cut->GetN(), pt, dedx);
                hist_pt_quality_acceptance_nhits_residual_cuts->Fill(pt);
                hist_eta_quality_acceptance_nhits_residual_cuts->Fill(eta);
                hist_dedx_quality_acceptance_nhits_residual_cuts->Fill(dedx);
                hist_charge_quality_acceptance_nhits_residual_cuts->Fill(charge);
                hist_quality_quality_acceptance_nhits_residual_cuts->Fill(quality);
                hist_hits_quality_acceptance_nhits_residual_cuts->Fill(nhits);
                hist_phi_quality_acceptance_nhits_residual_cuts->Fill(phi);


                float p;
                p = std::sqrt(pt*pt+pz*pz);
                graph_dedx_p_with_all_cuts->SetPoint(graph_dedx_p_with_all_cuts->GetN(), p, dedx);

                *//*float pt_charge = pt * charge;
                if (charge > 0) {
                    graph_dedx_pt_charge_pos_with_residual_cut->SetPoint(
                            graph_dedx_pt_charge_pos_with_residual_cut->GetN(), pt_charge, dedx);
                } else {
                    graph_dedx_pt_charge_neg_with_residual_cut->SetPoint(
                            graph_dedx_pt_charge_neg_with_residual_cut->GetN(), pt_charge, dedx);
                }*//*

                hist_dedx_pt_quality_acceptance_nhits_residual_cuts->Fill(pt, dedx);
                hist_dedx_pt_charge_quality_acceptance_nhits_residual_cuts->Fill(pt * charge, dedx);
            }
        }


        //std::cout<<rzint<<" * "<<rzslope<<std::endl;

        //hist6->Fill(abs(-1*rzint) / sqrt((rzslope * rzslope) + 1));
    }*/

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Without cuts
        hist_pt_without_cuts->Fill(pt);
        hist_eta_without_cuts->Fill(eta);
        hist_dedx_without_cuts->Fill(dedx);
        hist_quality_without_cuts->Fill(quality);
        hist_hits_without_cuts->Fill(ntpc);
        hist_phi_without_cuts->Fill(phi);
        hist_charge_without_cuts->Fill(charge);

        hist_dedx_pt_without_cuts->Fill(pt, dedx);
        hist_dedx_pt_charge_without_cuts->Fill(pt * charge, dedx);

        float p;
        p = std::sqrt(pt*pt+pz*pz);
        graph_dedx_p_with_all_cuts->SetPoint(graph_dedx_p_with_all_cuts->GetN(), p, dedx);

        // Apply pt cut
        if (pt > 0.2) {
            hist_pt_pt_cuts->Fill(pt);
            hist_eta_pt_cuts->Fill(eta);
            hist_dedx_pt_cuts->Fill(dedx);
            hist_quality_pt_cuts->Fill(quality);
            hist_hits_pt_cuts->Fill(ntpc);
            hist_phi_pt_cuts->Fill(phi);
            hist_charge_pt_cuts->Fill(charge);

            hist_dedx_pt_pt_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_pt_cuts->Fill(pt * charge, dedx);
        }

        // Apply quality cut
        if (pt > 0.2 && quality < 100) {
            hist_pt_pt_quality_cuts->Fill(pt);
            hist_eta_pt_quality_cuts->Fill(eta);
            hist_dedx_pt_quality_cuts->Fill(dedx);
            hist_quality_pt_quality_cuts->Fill(quality);
            hist_hits_pt_quality_cuts->Fill(ntpc);
            hist_phi_pt_quality_cuts->Fill(phi);
            hist_charge_pt_quality_cuts->Fill(charge);

            hist_dedx_pt_pt_quality_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_pt_quality_cuts->Fill(pt * charge, dedx);
        }

        // Apply ntpc cut
        if (pt > 0.2 && quality < 100 && ntpc > 20) {
            hist_pt_pt_quality_ntpc_cuts->Fill(pt);
            hist_eta_pt_quality_ntpc_cuts->Fill(eta);
            hist_dedx_pt_quality_ntpc_cuts->Fill(dedx);
            hist_quality_pt_quality_ntpc_cuts->Fill(quality);
            hist_hits_pt_quality_ntpc_cuts->Fill(ntpc);
            hist_phi_pt_quality_ntpc_cuts->Fill(phi);
            hist_charge_pt_quality_ntpc_cuts->Fill(charge);

            hist_dedx_pt_pt_quality_ntpc_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_pt_quality_ntpc_cuts->Fill(pt * charge, dedx);
        }

        // Apply nmaps cut
        if (pt > 0.2 && quality < 100 && ntpc > 20 && nmaps > 2) {
            hist_pt_pt_quality_ntpc_nmaps_cuts->Fill(pt);
            hist_eta_pt_quality_ntpc_nmaps_cuts->Fill(eta);
            hist_dedx_pt_quality_ntpc_nmaps_cuts->Fill(dedx);
            hist_quality_pt_quality_ntpc_nmaps_cuts->Fill(quality);
            hist_hits_pt_quality_ntpc_nmaps_cuts->Fill(ntpc);
            hist_phi_pt_quality_ntpc_nmaps_cuts->Fill(phi);
            hist_charge_pt_quality_ntpc_nmaps_cuts->Fill(charge);

            hist_dedx_pt_pt_quality_ntpc_nmaps_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_pt_quality_ntpc_nmaps_cuts->Fill(pt * charge, dedx);
        }

        // Apply nintt cut
        if (pt > 0.2 && quality < 100 && ntpc > 20 && nmaps > 2 && nintt > 0) {
            hist_pt_pt_quality_ntpc_nmaps_nintt_cuts->Fill(pt);
            hist_eta_pt_quality_ntpc_nmaps_nintt_cuts->Fill(eta);
            hist_dedx_pt_quality_ntpc_nmaps_nintt_cuts->Fill(dedx);
            hist_quality_pt_quality_ntpc_nmaps_nintt_cuts->Fill(quality);
            hist_hits_pt_quality_ntpc_nmaps_nintt_cuts->Fill(ntpc);
            hist_phi_pt_quality_ntpc_nmaps_nintt_cuts->Fill(phi);
            hist_charge_pt_quality_ntpc_nmaps_nintt_cuts->Fill(charge);

            hist_dedx_pt_pt_quality_ntpc_nmaps_nintt_cuts->Fill(pt, dedx);
            hist_dedx_pt_charge_pt_quality_ntpc_nmaps_nintt_cuts->Fill(pt * charge, dedx);

            float p;
                p = std::sqrt(pt*pt+pz*pz);
                graph_dedx_p_with_all_cuts->SetPoint(graph_dedx_p_with_all_cuts->GetN(), p, dedx);
        }
    }

    // 获取所有直方图的最大值
double max_pt = std::max({hist_pt_without_cuts->GetMaximum(), hist_pt_pt_cuts->GetMaximum(),
                          hist_pt_pt_quality_cuts->GetMaximum(), hist_pt_pt_quality_ntpc_cuts->GetMaximum(),
                          hist_pt_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_pt_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});
double max_eta = std::max({hist_eta_without_cuts->GetMaximum(), hist_eta_pt_cuts->GetMaximum(),
                           hist_eta_pt_quality_cuts->GetMaximum(), hist_eta_pt_quality_ntpc_cuts->GetMaximum(),
                           hist_eta_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_eta_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});
double max_dedx = std::max({hist_dedx_without_cuts->GetMaximum(), hist_dedx_pt_cuts->GetMaximum(),
                            hist_dedx_pt_quality_cuts->GetMaximum(), hist_dedx_pt_quality_ntpc_cuts->GetMaximum(),
                            hist_dedx_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_dedx_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});
double max_charge = std::max({hist_charge_without_cuts->GetMaximum(), hist_charge_pt_cuts->GetMaximum(),
                              hist_charge_pt_quality_cuts->GetMaximum(), hist_charge_pt_quality_ntpc_cuts->GetMaximum(),
                              hist_charge_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_charge_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});
double max_quality = std::max({hist_quality_without_cuts->GetMaximum(), hist_quality_pt_cuts->GetMaximum(),
                               hist_quality_pt_quality_cuts->GetMaximum(), hist_quality_pt_quality_ntpc_cuts->GetMaximum(),
                               hist_quality_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_quality_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});
double max_hits = std::max({hist_hits_without_cuts->GetMaximum(), hist_hits_pt_cuts->GetMaximum(),
                            hist_hits_pt_quality_cuts->GetMaximum(), hist_hits_pt_quality_ntpc_cuts->GetMaximum(),
                            hist_hits_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_hits_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});
double max_phi = std::max({hist_phi_without_cuts->GetMaximum(), hist_phi_pt_cuts->GetMaximum(),
                           hist_phi_pt_quality_cuts->GetMaximum(), hist_phi_pt_quality_ntpc_cuts->GetMaximum(),
                           hist_phi_pt_quality_ntpc_nmaps_cuts->GetMaximum(), hist_phi_pt_quality_ntpc_nmaps_nintt_cuts->GetMaximum()});

double max_range_pt = max_pt * 1.1;
double max_range_eta = max_eta * 1.1;
double max_range_dedx = max_dedx * 1.1;
double max_range_charge = max_charge * 1.1;
double max_range_quality = max_quality * 1.1;
double max_range_hits = max_hits * 1.1;
double max_range_phi = max_phi * 1.1;

// 第三种 y 轴的设置

hist_pt_without_cuts->GetYaxis()->SetRangeUser(1, max_pt * 1.1);
hist_pt_pt_cuts->GetYaxis()->SetRangeUser(1, max_pt * 1.1);
hist_pt_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_pt * 1.1);
hist_pt_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_pt * 1.1);
hist_pt_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_pt * 1.1);
hist_pt_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_pt * 1.1);

hist_eta_without_cuts->GetYaxis()->SetRangeUser(1, max_eta * 1.1);
hist_eta_pt_cuts->GetYaxis()->SetRangeUser(1, max_eta * 1.1);
hist_eta_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_eta * 1.1);
hist_eta_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_eta * 1.1);
hist_eta_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_eta * 1.1);
hist_eta_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_eta * 1.1);

hist_dedx_without_cuts->GetYaxis()->SetRangeUser(1, max_dedx * 1.1);
hist_dedx_pt_cuts->GetYaxis()->SetRangeUser(1, max_dedx * 1.1);
hist_dedx_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_dedx * 1.1);
hist_dedx_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_dedx * 1.1);
hist_dedx_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_dedx * 1.1);
hist_dedx_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_dedx * 1.1);

hist_charge_without_cuts->GetYaxis()->SetRangeUser(1, max_charge * 1.1);
hist_charge_pt_cuts->GetYaxis()->SetRangeUser(1, max_charge * 1.1);
hist_charge_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_charge * 1.1);
hist_charge_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_charge * 1.1);
hist_charge_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_charge * 1.1);
hist_charge_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_charge * 1.1);

hist_quality_without_cuts->GetYaxis()->SetRangeUser(1, max_quality * 1.1);
hist_quality_pt_cuts->GetYaxis()->SetRangeUser(1, max_quality * 1.1);
hist_quality_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_quality * 1.1);
hist_quality_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_quality * 1.1);
hist_quality_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_quality * 1.1);
hist_quality_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_quality * 1.1);

hist_hits_without_cuts->GetYaxis()->SetRangeUser(1, max_hits * 1.1);
hist_hits_pt_cuts->GetYaxis()->SetRangeUser(1, max_hits * 1.1);
hist_hits_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_hits * 1.1);
hist_hits_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_hits * 1.1);
hist_hits_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_hits * 1.1);
hist_hits_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_hits * 1.1);

hist_phi_without_cuts->GetYaxis()->SetRangeUser(1, max_phi * 1.1);
hist_phi_pt_cuts->GetYaxis()->SetRangeUser(1, max_phi * 1.1);
hist_phi_pt_quality_cuts->GetYaxis()->SetRangeUser(1, max_phi * 1.1);
hist_phi_pt_quality_ntpc_cuts->GetYaxis()->SetRangeUser(1, max_phi * 1.1);
hist_phi_pt_quality_ntpc_nmaps_cuts->GetYaxis()->SetRangeUser(1, max_phi * 1.1);
hist_phi_pt_quality_ntpc_nmaps_nintt_cuts->GetYaxis()->SetRangeUser(1, max_phi * 1.1);


    // 绘制直方图 - c1
TCanvas *c1 = new TCanvas("c1");
c1->Divide(6, 3);

// 绘制quality相关的直方图
c1->cd(1);
gPad->SetLogy();
hist_quality_without_cuts->Draw();

c1->cd(2);
gPad->SetLogy();
hist_quality_pt_cuts->Draw();

c1->cd(3);
gPad->SetLogy();
hist_quality_pt_quality_cuts->Draw();

c1->cd(4);
gPad->SetLogy();
hist_quality_pt_quality_ntpc_cuts->Draw();

c1->cd(5);
gPad->SetLogy();
hist_quality_pt_quality_ntpc_nmaps_cuts->Draw();

c1->cd(6);
gPad->SetLogy();
hist_quality_pt_quality_ntpc_nmaps_nintt_cuts->Draw();

// 绘制pt相关的直方图
c1->cd(7);
gPad->SetLogy();
hist_pt_without_cuts->Draw();

c1->cd(8);
gPad->SetLogy();
hist_pt_pt_cuts->Draw();

c1->cd(9);
gPad->SetLogy();
hist_pt_pt_quality_cuts->Draw();

c1->cd(10);
gPad->SetLogy();
hist_pt_pt_quality_ntpc_cuts->Draw();

c1->cd(11);
gPad->SetLogy();
hist_pt_pt_quality_ntpc_nmaps_cuts->Draw();

c1->cd(12);
gPad->SetLogy();
hist_pt_pt_quality_ntpc_nmaps_nintt_cuts->Draw();

// 绘制dedx相关的直方图
c1->cd(13);
gPad->SetLogy();
hist_dedx_without_cuts->Draw();

c1->cd(14);
gPad->SetLogy();
hist_dedx_pt_cuts->Draw();

c1->cd(15);
gPad->SetLogy();
hist_dedx_pt_quality_cuts->Draw();

c1->cd(16);
gPad->SetLogy();
hist_dedx_pt_quality_ntpc_cuts->Draw();

c1->cd(17);
gPad->SetLogy();
hist_dedx_pt_quality_ntpc_nmaps_cuts->Draw();

c1->cd(18);
gPad->SetLogy();
hist_dedx_pt_quality_ntpc_nmaps_nintt_cuts->Draw();

c1->SaveAs("/sphenix/user/dcxchenxi/files/dist_2024_0811_0.pdf");

TCanvas *c2 = new TCanvas("c2");
c2->Divide(6, 3);

// 绘制eta相关的直方图
c2->cd(1);
gPad->SetLogy();
hist_eta_without_cuts->Draw();

c2->cd(2);
gPad->SetLogy();
hist_eta_pt_cuts->Draw();

c2->cd(3);
gPad->SetLogy();
hist_eta_pt_quality_cuts->Draw();

c2->cd(4);
gPad->SetLogy();
hist_eta_pt_quality_ntpc_cuts->Draw();

c2->cd(5);
gPad->SetLogy();
hist_eta_pt_quality_ntpc_nmaps_cuts->Draw();

c2->cd(6);
gPad->SetLogy();
hist_eta_pt_quality_ntpc_nmaps_nintt_cuts->Draw();

// 绘制phi相关的直方图
c2->cd(7);
gPad->SetLogy();
hist_phi_without_cuts->Draw();

c2->cd(8);
gPad->SetLogy();
hist_phi_pt_cuts->Draw();

c2->cd(9);
gPad->SetLogy();
hist_phi_pt_quality_cuts->Draw();

c2->cd(10);
gPad->SetLogy();
hist_phi_pt_quality_ntpc_cuts->Draw();

c2->cd(11);
gPad->SetLogy();
hist_phi_pt_quality_ntpc_nmaps_cuts->Draw();

c2->cd(12);
gPad->SetLogy();
hist_phi_pt_quality_ntpc_nmaps_nintt_cuts->Draw();

// 绘制nhits相关的直方图
c2->cd(13);
gPad->SetLogy();
hist_hits_without_cuts->Draw();

c2->cd(14);
gPad->SetLogy();
hist_hits_pt_cuts->Draw();

c2->cd(15);
gPad->SetLogy();
hist_hits_pt_quality_cuts->Draw();

c2->cd(16);
gPad->SetLogy();
hist_hits_pt_quality_ntpc_cuts->Draw();

c2->cd(17);
gPad->SetLogy();
hist_hits_pt_quality_ntpc_nmaps_cuts->Draw();

c2->cd(18);
gPad->SetLogy();
hist_hits_pt_quality_ntpc_nmaps_nintt_cuts->Draw();

c2->SaveAs("/sphenix/user/dcxchenxi/files/dist_2024_0811_1.pdf");




    /*// 绘制直方图 - c3
    TCanvas *c3 = new TCanvas("c3");
    c3->Divide(5, 1);

    c3->cd(1);
    gPad->SetLogy();
    hist_charge_without_cuts->Draw();

    c3->cd(2);
    gPad->SetLogy();
    hist_charge_quality_cuts->Draw();

    c3->cd(3);
    gPad->SetLogy();
    hist_charge_quality_acceptance_cuts->Draw();

    c3->cd(4);
    gPad->SetLogy();
    hist_charge_quality_acceptance_nhits_cuts->Draw();

    c3->cd(5);
    gPad->SetLogy();
    hist_charge_quality_acceptance_nhits_residual_cuts->Draw();*/
    //c3->SaveAs("~/Downloads/dist_2024_0721_2.pdf");

    //c3->SaveAs("~/Downloads/dedx_vs_pt_total_after_cuts5.pdf");

    //file->Close();
    /* h_pt->GetXaxis()->SetRangeUser(0,900);
     h_pt->GetYaxis()->SetRangeUser(0,1);*/

    /*TCanvas *c4 = new TCanvas("c4");
    c4->Divide(2, 1);
    c4->cd(1);
    //gPad->SetLogy();
    graph_dedx_pt_without_residual_cut->GetXaxis()->SetLimits(0, 2);
    graph_dedx_pt_without_residual_cut->Draw("AP");
    gPad->SetLogy();
    c4->cd(2);
    //gPad->SetLogy();
    graph_dedx_pt_with_residual_cut->GetXaxis()->SetLimits(0, 2);
    graph_dedx_pt_with_residual_cut->Draw("AP");
    gPad->SetLogy();
    c4->SaveAs("~/Downloads/dist_2024_0716_3.pdf");


    TCanvas *c5 = new TCanvas("c5", "Canvas", 800, 600);
    c5->Divide(2, 1);

    c5->cd(1);
// 创建TMultiGraph来叠加两个TGraph
    TMultiGraph *mg1 = new TMultiGraph();
    mg1->Add(graph_dedx_pt_charge_pos_without_residual_cut);
    mg1->Add(graph_dedx_pt_charge_neg_without_residual_cut);

// 绘制叠加的图形
    mg1->SetTitle("dE/dx vs. pT*charge without residual cut;pT*charge (GeV/c);dE/dx");
    mg1->Draw("AP");
// 设置mg的x和y轴范围
//mg1->GetXaxis()->SetLimits(-2, 2); // 设置x轴范围
//mg1->GetYaxis()->SetRangeUser(0.1, 5000); // 设置y轴范围，并确保最小值为正
    gPad->SetLogy();
    mg1->GetXaxis()->SetLimits(-2, 2); // 设置x轴范围

    c5->cd(2);
// 创建TMultiGraph来叠加两个TGraph
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(graph_dedx_pt_charge_pos_with_residual_cut);
    mg2->Add(graph_dedx_pt_charge_neg_with_residual_cut);

// 绘制叠加的图形
    mg2->SetTitle("dE/dx vs. pT*charge with residual cut;pT*charge (GeV/c);dE/dx");
    mg2->Draw("AP");
// 设置mg的x和y轴范围
//mg2->GetXaxis()->SetLimits(-2, 2); // 设置x轴范围
//mg2->GetYaxis()->SetRangeUser(0.1, 5000); // 设置y轴范围，并确保最小值为正
    gPad->SetLogy();
    mg2->GetXaxis()->SetLimits(-2, 2); // 设置x轴范围

    c5->SaveAs("~/Downloads/dist_2024_0716_4.pdf");*/


    /*// 绘制直方图 - c4
    TCanvas *c4 = new TCanvas("c4", "Canvas", 1600, 800);
    c4->Divide(5, 2);

    c4->cd(1);
    gPad->SetLogz();
    hist_dedx_pt_without_cuts->Draw("COLZ");

    c4->cd(2);
    gPad->SetLogz();
    hist_dedx_pt_quality_cuts->Draw("COLZ");

    c4->cd(3);
    gPad->SetLogz();
    hist_dedx_pt_quality_acceptance_cuts->Draw("COLZ");

    c4->cd(4);
    gPad->SetLogz();
    hist_dedx_pt_quality_acceptance_nhits_cuts->Draw("COLZ");

    c4->cd(5);
    gPad->SetLogz();
    hist_dedx_pt_quality_acceptance_nhits_residual_cuts->Draw("COLZ");

    c4->cd(6);
    gPad->SetLogz();
   // hist_dedx_pt_charge_pos_without_cuts->Draw("COLZ");
    hist_dedx_pt_charge_neg_without_cuts->Draw("COLZ");

    c4->cd(7);
    gPad->SetLogz();
    hist_dedx_pt_charge_pos_quality_cuts->Draw("COLZ");
    hist_dedx_pt_charge_neg_quality_cuts->Draw("COLZ SAME");

    c4->cd(8);
    gPad->SetLogz();
    hist_dedx_pt_charge_pos_quality_acceptance_cuts->Draw("COLZ");
    hist_dedx_pt_charge_neg_quality_acceptance_cuts->Draw("COLZ SAME");

    c4->cd(9);
    gPad->SetLogz();
    hist_dedx_pt_charge_pos_quality_acceptance_nhits_cuts->Draw("COLZ");
    hist_dedx_pt_charge_neg_quality_acceptance_nhits_cuts->Draw("COLZ SAME");

    c4->cd(10);
    gPad->SetLogz();
    hist_dedx_pt_charge_pos_quality_acceptance_nhits_residual_cuts->Draw("COLZ");
    hist_dedx_pt_charge_neg_quality_acceptance_nhits_residual_cuts->Draw("COLZ SAME");

    c4->SaveAs("~/Downloads/dist_2024_0721_5_colz.pdf");*/

  //  file->Close();
// 绘制直方图
  TCanvas *c6 = new TCanvas("c6", "Canvas", 1200, 800);
c6->Divide(6, 2);

c6->cd(1);
gPad->SetLogz();
hist_dedx_pt_without_cuts->Draw("COLZ");

c6->cd(2);
gPad->SetLogz();
hist_dedx_pt_pt_cuts->Draw("COLZ");

c6->cd(3);
gPad->SetLogz();
hist_dedx_pt_pt_quality_cuts->Draw("COLZ");

c6->cd(4);
gPad->SetLogz();
hist_dedx_pt_pt_quality_ntpc_cuts->Draw("COLZ");

c6->cd(5);
gPad->SetLogz();
hist_dedx_pt_pt_quality_ntpc_nmaps_cuts->Draw("COLZ");

c6->cd(6);
gPad->SetLogz();
hist_dedx_pt_pt_quality_ntpc_nmaps_nintt_cuts->Draw("COLZ");

c6->cd(7);
gPad->SetLogz();
hist_dedx_pt_charge_without_cuts->Draw("COLZ");

c6->cd(8);
gPad->SetLogz();
hist_dedx_pt_charge_pt_cuts->Draw("COLZ");

c6->cd(9);
gPad->SetLogz();
hist_dedx_pt_charge_pt_quality_cuts->Draw("COLZ");

c6->cd(10);
gPad->SetLogz();
hist_dedx_pt_charge_pt_quality_ntpc_cuts->Draw("COLZ");

c6->cd(11);
gPad->SetLogz();
hist_dedx_pt_charge_pt_quality_ntpc_nmaps_cuts->Draw("COLZ");

c6->cd(12);
gPad->SetLogz();
hist_dedx_pt_charge_pt_quality_ntpc_nmaps_nintt_cuts->Draw("COLZ");
    c6->SaveAs("/sphenix/user/dcxchenxi/files/dist_2024_0811_2.pdf");


    TCanvas *c7 = new TCanvas("c7", "Canvas", 1200, 800);

    c7->SetLogx();
    c7->SetLogy();

    graph_dedx_p_with_all_cuts->GetXaxis()->SetLimits(0.01, 100); // 设置 X 轴范围
    graph_dedx_p_with_all_cuts->GetYaxis()->SetRangeUser(2, 10000); // 设置 Y 轴范围
    graph_dedx_p_with_all_cuts -> Draw("AP");
    /*xyscanGraph->SetMinimum(2); // 设置 y 轴的最小值
    xyscanGraph->SetMaximum(10000); // 设置 y 轴的最大值*/
    pion_graph->Draw("P same");
    kaon_graph->Draw("P same");
    proton_graph->Draw("P same");
    c7->SaveAs("/sphenix/user/dcxchenxi/files/dist_compare_0811.pdf");

}





