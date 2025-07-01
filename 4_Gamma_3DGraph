import pandas as pd
import plotly.express as px
import numpy as np
from tqdm import tqdm

# Load dataset
df = pd.read_csv("papa_gvr_output.csv")

# Remove rows with NaN GVR
df = df[df['GVR'].notna()]

# Round gamma values to avoid float mismatch
df['Gamma1'] = df['Gamma1'].round(3)
df['Gamma2'] = df['Gamma2'].round(3)
df['Gamma3'] = df['Gamma3'].round(3)
df['Gamma4'] = df['Gamma4'].round(3)

# Define Gamma1 values to iterate through
gamma1_values = np.round(np.arange(0.4, 2.1, 0.1), 2)

# Group once and compute median GVR for all combinations
summary = (
    df.groupby(['Gamma1', 'Gamma2', 'Gamma3', 'Gamma4'])['GVR']
    .median()
    .reset_index()
)

# Loop through Gamma1 values and plot
for g1 in tqdm(gamma1_values, desc="Plotting 3D scatter plots"):
    subset = summary[summary['Gamma1'] == g1]

    fig = px.scatter_3d(
        subset,
        x='Gamma2',
        y='Gamma3',
        z='Gamma4',
        color='GVR',
        color_continuous_scale='viridis',
        title=f'GVR Median at Gamma1 = {g1}',
    )

    fig.update_traces(marker=dict(size=5))
    fig.update_layout(
        scene=dict(
            xaxis_title='Gamma2',
            yaxis_title='Gamma3',
            zaxis_title='Gamma4',
        ),
        coloraxis_colorbar=dict(title='GVR')
    )

    fig.show()
