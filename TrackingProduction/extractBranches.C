#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

void extractBranches() {
    // 打开原始文件
    TFile *inputFile = TFile::Open("/sphenix/user/dcxchenxi/Beam_Residuals_2/output_frawley/merged_frawley_ww.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // 获取树
    TTree *inputTree = (TTree *)inputFile->Get("residualtree");
    if (!inputTree) {
        std::cerr << "Error getting tree from input file" << std::endl;
        inputFile->Close();
        return;
    }

    // 定义变量
    std::vector<float> *clusgx = nullptr;
    std::vector<float> *stategx = nullptr;
    std::vector<float> *clusgy = nullptr;
    std::vector<float> *stategy = nullptr;
    float rzint, rzslope, pt, eta, dedx, quality, phi, pz;
    int nhits, charge, ntpc, nintt, nmaps;

    // 设置 branches 地址
    inputTree->SetBranchAddress("clusgx", &clusgx);
    inputTree->SetBranchAddress("stategx", &stategx);
    inputTree->SetBranchAddress("clusgy", &clusgy);
    inputTree->SetBranchAddress("stategy", &stategy);
    inputTree->SetBranchAddress("pt", &pt);
    inputTree->SetBranchAddress("pz", &pz);
    inputTree->SetBranchAddress("eta", &eta);
    inputTree->SetBranchAddress("dedx", &dedx);
    inputTree->SetBranchAddress("charge", &charge);
    inputTree->SetBranchAddress("quality", &quality);
    inputTree->SetBranchAddress("phi", &phi);
    inputTree->SetBranchAddress("nhits", &nhits);
    inputTree->SetBranchAddress("ntpc", &ntpc);
    inputTree->SetBranchAddress("nintt", &nintt);
    inputTree->SetBranchAddress("nmaps", &nmaps);

    // 创建输出文件
    TFile *outputFile = TFile::Open("/sphenix/user/dcxchenxi/Beam_Residuals_2/output_frawley/selected_data_frawley_ww.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        inputFile->Close();
        return;
    }

    // 创建新的树
    TTree *outputTree = new TTree("residualtree", "Extracted branches");

    // 克隆 branches 到新的树
    outputTree->Branch("clusgx", &clusgx);
    outputTree->Branch("stategx", &stategx);
    outputTree->Branch("clusgy", &clusgy);
    outputTree->Branch("stategy", &stategy);
    outputTree->Branch("pt", &pt);
    outputTree->Branch("pz", &pz);
    outputTree->Branch("eta", &eta);
    outputTree->Branch("dedx", &dedx);
    outputTree->Branch("charge", &charge);
    outputTree->Branch("quality", &quality);
    outputTree->Branch("phi", &phi);
    outputTree->Branch("nhits", &nhits);
    outputTree->Branch("ntpc", &ntpc);
    outputTree->Branch("nintt", &nintt);
    outputTree->Branch("nmaps", &nmaps);
    // 遍历事件并填充新树
    Long64_t nentries = inputTree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        inputTree->GetEntry(i);
        outputTree->Fill();
    }

    // 写入输出文件
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
    inputFile->Close();
}

