import os
import statistics
import random
import numpy as np

import plotly.figure_factory as ff
import plotly.express as px


import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc

import dash_html_components as html
from dash.dependencies import Input, Output

#------------------------------
# read genome file
def readGenomeFile(fileName):
    os.chdir('C:/Users/Xandy/Documents/Python/MIT 2020/Run BLAST/Viral_Databases/refseq')
    genome= ''
    with open(fileName, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
#------------------------------

#------------------------------
# generate sequencing reads 
reads= []
def readsGenerator(genome, minLen, maxLen):  
    global reads
    readLen = random.randint(minLen, maxLen)
    startPos = random.randint(0,len(genome)-readLen)-1
    read1, read2, read3 = (genome[:startPos], 
                           genome[startPos:startPos+readLen], 
                           genome[startPos+readLen:])
    # if len(read2) !=0:
    #     reads.append(read2)
    for read in [read1, read2,read3]:
        if len(read) !=0:
            if len(read) > maxLen:
                readsGenerator(read, minLen, maxLen)
            else:
                reads.append(read)
#------------------------------

#------------------------------
# start the dash app
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])



app.layout = dbc.Container(

    html.Div(
        [
            html.Div([
                html.H1('Virus Fragmantation Simulation', style={'text-align':'center'})
            ]),

            dbc.Row(
                [
                    dbc.Col(
                        [
                        # 1st dive for vrius name
                        html.Div([
                            html.Hr(),

                            dbc.Row([

                            dbc.Col(html.Label(
                                'Virus genome:', 
                                style= {'font-weight':'bold'}
                            )),
                            dbc.Col(dbc.Select(
                                id='virusName',
                                options= [
                                    {'label':'FeLV', 'value':'Feline Leukemia Virus.fasta'},
                                    {'label':'MVM', 'value':'Minute virus of mice.fasta'},
                                    {'label':'PCV1', 'value':'Porcine circovirus 1.fasta'}
                                    ],
                                value= 'Feline Leukemia Virus.fasta',
                                style= {'width': '40%'}
                            ))

                            ])

                        ]),

                        # 2st input for virus count
                        html.Div([
                            html.Hr(),

                            dbc.Row([

                            dbc.Col(html.Label(
                                'Virus count: ', 
                                style={'font-weight':'bold'}
                            )),
                            dbc.Col(dbc.Input(
                                id= 'virusCounter',
                                type= 'number',
                                debounce=False,
                                min=1,
                                max=500,
                                step=1,
                                value=50,
                                style = {'width': '40%'}
                            ))])
                        ]),

                        # 3st input for min
                        html.Div([
                            html.Hr(),
                            dbc.Row(
                                [

                            dbc.Col(html.Label(
                                'Virus fragment length:', 
                                style={'font-weight':'bold'}
                            )),
                            dbc.Col(
                            dcc.RangeSlider(
                                id='readLengthSlider',
                                min=1,
                                max=10000,
                                step=10,
                                value=[500, 5000],
                                tooltip = { 'always_visible': True }
                            ), align="center")
                            
                            # dbc.Input(
                            #     id= 'min_length',
                            #     type= 'number',
                            #     debounce=False,
                            #     min=1,
                            #     max=5000,
                            #     step=1,
                            #     value=500,
                            #     style = {'width': '40%'}
                            # ))],
                            # dbc.Input(
                            #     id= 'max_length',
                            #     type= 'number',
                            #     debounce=False,
                            #     min=1,
                            #     max=50000,
                            #     step=1,
                            #     value=5000,
                            #     style = {'width': '40%'}
                            # ))]
                            
                            ]
                            ),
                            html.Hr()

                        ])
                        ], width=4, align="start"
                    ),


                    dbc.Col([

                            html.H5(id='virusName_output'),
                            dbc.ListGroupItem([
                            html.Div(id='virusSize_output'),
                            html.Div(id='virusCounter_output'),
                            html.Div(id='readLengthSlider_output1'),
                            html.Div(id='readLengthSlider_output2'),
                            html.Div(id='state1'),
                            html.Div(id='state2'),
                            html.Div(id='state3')
                                ])
                            ], width=3, align="center")
                            
                ]),



                    dbc.Col(dcc.Graph(id='thegraph'))

        ]
    )

,fluid=True)



#------------------------------
# Display Results:
#------------------------------
'''1: The "Input" "value" from the "callback" below is comning from the dcc.Dropdown and is past to 
   the function "virusName_func" below.
   
   2: The "Output" "children" from the "callback" below is comning from the function 
   "virusName_func" below and is past to the html.Div "virusName_output" above.
'''

@app.callback(Output('virusName_output', 'children'), [Input('virusName','value')])
def virusName_func(fileName):
    virus_name= fileName.split('.')[0]
    # return f'Virus name: {virus_name}', fileName
    return f'{virus_name}'

@app.callback(Output('virusSize_output', 'children'), [Input('virusName','value')])
def virusSize_func(fileName):
    genome = readGenomeFile(fileName)
    return f'Virus size (bases): {str(len(genome))}'

@app.callback(Output('virusCounter_output', 'children'), [Input('virusCounter','value')])
def virusCounter_func(virusCount):
     return f'Virus counts: {virusCount}'

@app.callback([Output('readLengthSlider_output1', 'children'), 
               Output('readLengthSlider_output2', 'children')], 
              [Input('readLengthSlider','value')])
def readLengthSlider_func(value):
    readLength = value[1]-value[0]
    return (f'Read range: {value}', f'Read length spread: {readLength}')

@app.callback([Output('thegraph','figure'),
               Output('state1', 'children'), 
               Output('state2', 'children'),
               Output('state3', 'children')], 
              [Input('virusName','value'),
               Input('virusCounter','value'),
               Input('readLengthSlider','value')])
def update_graph(fileName, count, minmax):
    global reads
    genome = readGenomeFile(fileName)
    genomes = [genome] * count

    minSize = minmax[0]
    maxSize = minmax[1]

    for each in genomes:
        readsGenerator (each,minSize,maxSize)
    
    reads2= reads
    reads=[]
    reads_len= [len(i) for i in reads2]
    mean= round(statistics.mean(reads_len),1)
    stdev= round(statistics.stdev(reads_len),1)
    readsAboveMin= [i for i in reads_len if i>=minSize]
    readsBelowMin= [i for i in reads_len if i<minSize]
    

    readsAboveMinPercent= round((len(readsAboveMin)/len(reads2)*100),1)
    readsBelowMinPercent= round((len(readsBelowMin)/len(reads2)*100),1)

    print('-----')
    arr1 = np.array(readsAboveMin)
    arr2 = np.array(readsBelowMin)
    hist_data = [arr1,arr2]

    print(readsAboveMin)
    print(readsBelowMin)
    print('-----')
    group_labels = [f"Read length > {minSize} (bases) = {readsAboveMinPercent} %",
                    f"Read length ≤ {minSize} (bases) = {readsBelowMinPercent} %"]
    # colors = ['#333F44', '#37AA9C']


    # fig = px.histogram(reads_len,nbins=100)

    
    # fig = ff.create_distplot(hist_data,
    #                          group_labels,
    #                          bin_size=.2, 
    #                          show_rug=True, 
    #                         curve_type='normal',show_hist=False,

    #                          histnorm='probability')
    fig = ff.create_distplot([np.concatenate(hist_data)],
                            ["Density"],
                            bin_size=.2, 
                            show_rug=True, 
                            curve_type='normal',show_hist=False,

                            histnorm='probability')

    fig.update_layout(showlegend=True,
                      height=600,
                      width=1000,
                      title_text='Density Plot of Read Lengths')

    return (fig, f'Total fragmants: {str(len(reads2))}' ,
            f'Read length ≤ {minSize} (bases) = {str(len(readsBelowMin))}',
            f'Read length > {minSize} (bases) = {str(len(readsAboveMin))}')
#------------------------------

if __name__ == '__main__':
    # type this in your browser: http://127.0.0.1:8050/
    app.run_server(debug=True, port='8050')
    



# I had to change the fragmentation algorithem so it fragments all 3 substrings. The only benefit of using Oxford Nanopore is the fact that it sequence in sears. Long reads are not benifciul in our purpeses as long as reads are not too short to maintain idenity.