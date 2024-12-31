import ROOT


def test_root():
    with ROOT.TFile("test.root", "RECREATE") as f:
        h = ROOT.TH1F("h", "Test Histogram", 100, -4, 4)
        h.FillRandom("gaus")
        h.Write()


if __name__ == "__main__":
    test_root()
    print("ROOT operations completed successfully.")
