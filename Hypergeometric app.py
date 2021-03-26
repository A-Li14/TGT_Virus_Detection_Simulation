import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from scipy.stats import hypergeom
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
import math
import numpy as np


# def hypergeo(N,n,K) :
#     """Return hypergeometric PMF. N=background size, n=# drawn, K=target size"""
#     def pmf(k) :
#         return math.comb(K,k)*math.comb(N,n-k)/math.comb(N+K,n)
#     return pmf


app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div([
    html.H1("Hypergeometric Model for Target Detection using Metagenomic Sequencing"),
    html.Br(),
    ###Input elements
    
    ###Update to look more like fragmentation app
    html.Div([
        dbc.Row([
            dbc.Col(html.Label(
                'Distribution Parameters:', 
                style= {'font-weight':'bold'}
            )),
            dbc.Col(),
            dbc.Col(html.Label(
                'Detection Criteria, Plot Controls:', 
                style= {'font-weight':'bold'}
            )),
            dbc.Col()
            
        ])
    ]),
    html.Div([
        dbc.Row([
            dbc.Col(html.Label(
                "Total Reads: "
            )),
            dbc.Col(
                dcc.Input(id='total_pop', value=1770410, type='number',min=1)
            ),
            dbc.Col(html.Label(
                "Target Threshold: "
            )),
            dbc.Col(
                dcc.Input(id='target_thresh', value=1, type='number',min=1)
            )
        ])
    ]),
    html.Div([
        dbc.Row([
            dbc.Col(html.Label(
                "Target Reads: "
            )),
            dbc.Col(
                dcc.Input(id="target_pop",value=284,type="number",min=1)
            ),
            dbc.Col(html.Label(
                "X-axis Range (10-folds): "
            )),
            dbc.Col(
                dcc.Input(id="range",value=1,type="number",min=1,max=100)
            )
        ])
    ]),
    html.Div([
        dbc.Row([
            dbc.Col(html.Label(
                "Reads Sequenced: "
            )),
            dbc.Col(
                dcc.Input(id="sample_size",value=3775,type="number",min=1)
            ),
            dbc.Col(html.Label(
                "Target Confidence: "
            )),
            dbc.Col(
                dcc.Input(id="target_conf",value=.99,type="number",min=0,max=1)
            )
        ])
    ]),


    # html.Div([html.H3("Distribution Parameters"),
    #           "Total Reads: ",
    #           dcc.Input(id='total_pop', value=1770410, type='number',min=1),
    #           "Target Reads: ",
    #           dcc.Input(id="target_pop",value=284,type="number",min=1),
    #           "Reads Sequenced: ",
    #           dcc.Input(id="sample_size",value=3775,type="number",min=1)]),
    # html.Br(),
    # html.Div([html.H3("Detection Criteria and Plot Controls"),
    #           "Target Threshold",
    #           dcc.Input(id="target_thresh",value=1,type="number",min=1),
    #           "X-axis Range (10-folds)",
    #           dcc.Input(id="range",value=1,type="number",min=1,max=100),
    #           "Target Confidence",
    #           dcc.Input(id="target_conf",value=.99,type="number",min=0,max=1)
    #           ]),
    html.Br(),
    ###Outputs
    dcc.Graph(id="confidence_plot"),
    dbc.Table(id="summary_table")
])


@app.callback(
    [Output("confidence_plot","figure"),
    Output("summary_table","children")],
    [Input("total_pop","value"),
    Input("target_pop","value"),
    Input("sample_size","value"),
    Input("range","value"),
    Input("target_thresh","value"),
    Input("target_conf","value")]
)
def update_hypergeo(total_population,target_population,sampled,x_range,target_threshold,target_confidence) :
    hyper = hypergeom(M=total_population,n=target_population,N=sampled)
    
    ###Create figure
    low = max(int(sampled/10**x_range),1)
    high = min(int(sampled*10**x_range),total_population)
    full_range = range(low,high,max(int((high-low)/100),1)) #Take 100 points between low and high
    conf = pd.DataFrame(
        # {"sample_size": [math.log2(i) for i in full_range],
        {"sample_size": [i for i in full_range],
        "unrounded_conf": [1-hypergeom.cdf(target_threshold-1,M=total_population,n=target_population,N=i) for i in full_range]
        }
    )
    conf["confidence"] = [round(i,3) for i in conf["unrounded_conf"]]
    fig = px.scatter(conf,x="sample_size",y="confidence",title="P(X>threshold) by Reads Sequenced",labels={'sample_size':'Sample Size', 'confidence':'Confidence'})
    
    ###Create summary table
    conf99 = conf[conf.unrounded_conf>target_confidence]
    min_reads = None
    if len(conf99) > 0 :
        min_reads = conf99.iloc[0,0]
    min_reads_label = "~Reads for "+str(target_confidence)+" Confidence"
    conf_level = round(1-hyper.cdf(target_threshold-1),3)
    target_ratio = round(target_population/total_population,3)
    prop_sampled = round(sampled/total_population,3)
    summary = pd.DataFrame({"Target_Read_Ratio":[target_ratio],
        "Sampled_Proportion":[prop_sampled],
        "Confidence_Level":[conf_level],
        min_reads_label:[min_reads]})
        # min_reads_label:[round(2**min_reads) if min_reads is not None else None]})
    summary_table = dbc.Table.from_dataframe(summary, striped=True, bordered=True, hover=True)
    return fig,summary_table



if __name__ == '__main__':
    app.run_server(debug=True)