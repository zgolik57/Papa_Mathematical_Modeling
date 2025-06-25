import pandas as pd
import plotly.express as px

df = pd.read_csv("split_fast_gvr_with_gamma.csv")  # replace with your actual filename

subset = df[(df["Gamma"] == 1.0) & (df["[NPM1]_WT"] == 0.5)]

fig = px.scatter_3d(
    subset,
    x='[NPM1]_Halo',
    y='[NPM1]_NFAST',
    z='[NPM1]_CFAST',
    color='GVR',
    color_continuous_scale='viridis',
    title='GVR at Gamma = 1.0 and [WT] = 0.5',
)

fig.update_traces(marker=dict(size=5))
fig.update_layout(
    scene=dict(
        xaxis_title='[NPM1]_Halo',
        yaxis_title='[NPM1]_NFAST',
        zaxis_title='[NPM1]_CFAST',
    ),
    coloraxis_colorbar=dict(title='GVR')
)

fig.show()