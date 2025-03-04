// check_duplicates.C
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <set>
#include <limits>

// Define a simple struct for the key.
struct HitIdentifier
{
    int event;
    int hitsetkey;
    int hitkey;

    // Define ordering for use in std::set.
    bool operator<(const HitIdentifier &other) const
    {
        if (event != other.event)
            return event < other.event;
        if (hitsetkey != other.hitsetkey)
            return hitsetkey < other.hitsetkey;
        return hitkey < other.hitkey;
    }
};

void check_duplicates()
{
    // Manually specify the file path.
    std::string filePath = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_53534_0_clusterizer_edgeOff_staticOff_acts_0218.root";

    // Open the file.
    TFile *f = TFile::Open(filePath.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return;
    }

    // Retrieve the hitTree.
    TTree *t = (TTree *)f->Get("hitTree");
    if (!t)
    {
        std::cerr << "Error: 'hitTree' not found in file " << filePath << std::endl;
        f->Close();
        return;
    }

    // Variables to hold branch data.
    int tpc_event = 0;
    int tpc_hitkey = 0;
    int tpc_hitsetkey = 0;

    // Set branch addresses (make sure these match your tree).
    t->SetBranchAddress("event", &tpc_event);
    t->SetBranchAddress("hitkey", &tpc_hitkey);
    t->SetBranchAddress("hitsetkey", &tpc_hitsetkey);

    int nEntries = t->GetEntries();
    std::set<HitIdentifier> uniqueKeys;
    int duplicateCount = 0;

    for (int i = 0; i < nEntries; i++)
    {
        t->GetEntry(i);
        HitIdentifier key = {tpc_event, tpc_hitsetkey, tpc_hitkey};

        // If insertion fails, the key is a duplicate.
        if (!uniqueKeys.insert(key).second)
        {
            duplicateCount++;
        }
    }

    std::cout << "Total entries in hitTree: " << nEntries << std::endl;
    std::cout << "Unique key combinations: " << uniqueKeys.size() << std::endl;
    std::cout << "Number of duplicate key combinations: " << duplicateCount << std::endl;

    f->Close();
}

// Call the function if the macro is run
void check_duplicates_macro()
{
    check_duplicates();
}