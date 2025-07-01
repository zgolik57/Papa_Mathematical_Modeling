import pandas as pd
import plotly.express as px

df = pd.read_csv("papa_gvr_output.csv")

#unique identity for gamma set
df['GammaCombo'] = df.apply(
    lambda row: f"(γ1={row['Gamma1']}, γ2={row['Gamma2']}, γ3={row['Gamma3']}, γ4={row['Gamma4']})",
    axis=1
)

#group by gamma config
summary = df.groupby('GammaCombo').agg(
    GVR_median=('GVR', 'median'),
    GVR_std=('GVR', 'std'),
    Gamma1=('Gamma1', 'first'),
    Gamma2=('Gamma2', 'first'),
    Gamma3=('Gamma3', 'first'),
    Gamma4=('Gamma4', 'first'),
    Count=('GVR', 'count')
).reset_index()

#sort by median
summary = summary.sort_values(by='GVR_median').reset_index(drop=True)
summary['GammaIndex'] = summary.index

#plot with error bars
fig = px.scatter(
    summary,
    x='GammaIndex',
    y='GVR_median',
    #error_y='GVR_std',
    hover_name='GammaCombo',
    hover_data={
        'GVR_median': True,
        #"'GVR_std': True,
        'Gamma1': True,
        'Gamma2': True,
        'Gamma3': True,
        'Gamma4': True,
        'Count': True,
        'GammaIndex': False
    },
    labels={'GammaIndex': 'Gamma Configuration (sorted)', 'GVR_median': 'Median GVR'},
    title='Median GVR vs Gamma Config (All Concentrations)'
)

fig.update_traces(marker=dict(size=6))
fig.update_layout(hovermode='closest')
fig.show()
