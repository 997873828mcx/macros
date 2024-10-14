#include <TFile.h>
#include <TH3F.h>
#include <iostream>

void extract_histogram_values() {
    // 打开 ROOT 文件
    TFile *file = TFile::Open("magnetic_field_histograms_track.root", "READ");

    // 检查文件是否成功打开
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // 从文件中获取直方图
    TH3F *hist_bx = (TH3F *)file->Get("hist_bx");
    TH3F *hist_by = (TH3F *)file->Get("hist_by");
    TH3F *hist_bz = (TH3F *)file->Get("hist_bz");

    // 检查直方图是否存在
    if (!hist_bx || !hist_by || !hist_bz) {
        std::cerr << "Error getting histogram(s) from file" << std::endl;
        file->Close();
        return;
    }

    // 输出每个 bin 的 bx, by, bz 值
    int nBinsX = hist_bx->GetNbinsX();
    int nBinsY = hist_bx->GetNbinsY();
    int nBinsZ = hist_bx->GetNbinsZ();

    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            for (int k = 1; k <= nBinsZ; ++k) {
                float bx = hist_bx->GetBinContent(i, j, k);
                float by = hist_by->GetBinContent(i, j, k);
                float bz = hist_bz->GetBinContent(i, j, k);
                std::cout << "Bin (" << i << ", " << j << ", " << k << "): bx = " << bx << ", by = " << by << ", bz = " << bz << std::endl;
            }
        }
    }

    // 关闭文件
    file->Close();
}