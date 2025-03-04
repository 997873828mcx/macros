#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <string>
#include <limits>

void check_hitkey_duplicates()
{
    // Input file path
    std::string input_path = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_52844_0_clusterizer_edgeOff_staticOff_acts_0205.root";

    // Open input file
    TFile *f_clusterizer = new TFile(input_path.c_str());
    if (!f_clusterizer || f_clusterizer->IsZombie())
    {
        std::cerr << "Error opening clusterizer file\n";
        return;
    }

    // Get the hit tree
    TTree *t_tpc_hit = (TTree *)f_clusterizer->Get("hitTree");
    if (!t_tpc_hit)
    {
        std::cerr << "Error: 'hitTree' not found in clusterizer file.\n";
        f_clusterizer->Close();
        delete f_clusterizer;
        return;
    }

    // Variables to read from tree
    ULong64_t tpc_hitkey = std::numeric_limits<uint64_t>::max();

    // Set branch address
    t_tpc_hit->SetBranchAddress("hitkey", &tpc_hitkey);

    // Map to store hitkey counts
    std::map<ULong64_t, int> hitkey_count;

    // First pass: count hitkeys
    std::cout << "Analyzing " << t_tpc_hit->GetEntries() << " entries..." << std::endl;

    for (Long64_t i = 0; i < t_tpc_hit->GetEntries(); i++)
    {
        t_tpc_hit->GetEntry(i);
        hitkey_count[tpc_hitkey]++;
    }

    // Analyze results
    int total_duplicates = 0;
    int max_duplicates = 0;
    ULong64_t most_duplicated_key = 0;

    for (const auto &[hitkey, count] : hitkey_count)
    {
        if (count > 1)
        {
            total_duplicates++;
            if (count > max_duplicates)
            {
                max_duplicates = count;
                most_duplicated_key = hitkey;
            }
            std::cout << "Hitkey " << hitkey << " appears " << count << " times" << std::endl;
        }
    }

    // Print summary
    std::cout << "\nSummary:" << std::endl;
    std::cout << "Total entries in tree: " << t_tpc_hit->GetEntries() << std::endl;
    std::cout << "Total unique hitkeys: " << hitkey_count.size() << std::endl;
    std::cout << "Number of duplicated hitkeys: " << total_duplicates << std::endl;
    if (total_duplicates > 0)
    {
        std::cout << "Most duplicated hitkey: " << most_duplicated_key
                  << " (appears " << max_duplicates << " times)" << std::endl;
    }

    // Clean up
    f_clusterizer->Close();
    delete f_clusterizer;
}