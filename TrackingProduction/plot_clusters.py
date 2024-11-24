import uproot  
import numpy as np
import pandas as pd
import plotly.graph_objects as go  

# Read ROOT file 
file = uproot.open("/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_kso_53217_0_clusters_edgeOn_staticOff_1120_0.root")
tree = file["seed_clusters"]  
data = tree.arrays(library="pandas")  

# Create interactive 3D scatter plot
fig = go.Figure()

# Add used clusters (is_rejected == 0)
fig.add_trace(go.Scatter3d(
    x=data[data['is_rejected']==0]['x'],
    y=data[data['is_rejected']==0]['y'],
    z=data[data['is_rejected']==0]['z'],
    mode='markers',
    name='Used Clusters',
    marker=dict(
        size=2,
        color='blue',
        opacity=0.8
    )
))

# Add rejected clusters (is_rejected == 1)
fig.add_trace(go.Scatter3d(
    x=data[data['is_rejected']==1]['x'],
    y=data[data['is_rejected']==1]['y'],
    z=data[data['is_rejected']==1]['z'],
    mode='markers',
    name='Rejected Clusters',
    marker=dict(
        size=2,
        color='red',
        opacity=0.8
    )
))

# Update layout
fig.update_layout(
    title='Cluster Visualization',
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    )
)

# Show plot
fig.show()  # Opens in web browser with interactive controls