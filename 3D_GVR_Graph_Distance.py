import pandas as pd
import plotly.express as px

df = pd.read_csv("npm1_gvr_distance_output.csv")  # replace with your actual filename

subset = df[(df["Gamma1"] == 1.0) & (df["Gamma2"] == 0.5)]

fig = px.scatter_3d(
    subset,
    x='[NPM1]_Halo',
    y='[NPM1]_FAST',
    z='[NPM1]_WT',
    color='GVR',
    color_continuous_scale='viridis',
    title='GVR at Gamma1 = 1.0 and Gamma2 = 0.5',
)

fig.update_traces(marker=dict(size=5))
fig.update_layout(
    scene=dict(
        xaxis_title='[NPM1]_Halo',
        yaxis_title='[NPM1]_FAST',
        zaxis_title='[NPM1]_WT',
    ),
    coloraxis_colorbar=dict(title='GVR')
)

fig.show()