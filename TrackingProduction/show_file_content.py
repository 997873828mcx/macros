import uproot

# Replace with your ROOT file path
file_path = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff_1120_0.root"

# Open and check the structure
file = uproot.open(file_path)
print("Tree names in file:", file.keys())

# Get the tree you want to use
tree_name = "tracking_clusters"  # Replace with your actual tree name
tree = file[tree_name]
print("\nBranches in tree:", tree.keys())