import sys
import os
sys.path.append(os.getcwd())
import pandas as pd
import numpy as np
import plotly.express as px
import dash
from dash import dcc
from dash import html
from dash import ctx
from dash import ALL, MATCH
from dash.dependencies import Input, Output, State
import calc_utils as cu
import data_utils as du
import dash_bootstrap_components as dbc
import dash_daq as daq
import json
import base64
import io

import time

#----------------------- data init -------------------------------
cbio_df = pd.DataFrame()
included_studies = []   

features = ['age', 'age at diagnosis', 'sex', 'mutations', 'icds']
patient_data_df = pd.DataFrame(columns=features)

hover_data = {'age':True,
              'age at diagnosis':True,
              'sex':True,
              'mutations':True,
              'icds':True,
              'x':False,
              'y':False
              }

feature_weights = dict(zip(features,np.ones(len(patient_data_df.columns))))

studyIDs = du.getAllStudyIDsFromCBioPortal()

source_df = patient_data_df.copy()

#----------------- calculation inits --------------------

# load local gene sim (panther)
du.loadGeneSimLookUp('panther')
du._initOncoKBLookUp_()
du._initICDTree_()
#--------------------------------------------------------

if not source_df.empty:
    # Create initial distance matrix
    source_df = cu.getDistScores(source_df,feature_weights,
                                    calc_dem_dist=True,calc_icds_dist=True,calc_mut_dist=True)
else:
    source_df = pd.DataFrame(columns=features)
    
# Define the Plotly Dash app
app = dash.Dash(__name__,external_stylesheets=[dbc.themes.YETI])
#----------------------------------------------------------------

# Define the app layout
app.layout = html.Div([
    dbc.Row(
        [
        html.H1("Patient Similarity Visualization", style={'textAlign': 'center'}),    
        ],
    ),
    dbc.Row(
        [
        html.H2("IMPORTANT USER INFO", 
                style={'textAlign': 'center',
                       'color' : 'red'
                       }),
        html.H3('Please wait until the tab stopped \'Updating...\' !', 
                style={'textAlign': 'center',
                       'color' : 'red'
                       }),
        html.H4('Tweaking analysis parameters before the computation is over might crash the application.', 
                style={'textAlign': 'center',
                       'color' : 'red'
                       }),  
        ],
    ),
    dbc.Row([
        dbc.Col([
                html.Div([
                            dcc.Graph(id='scatter-plot'),
                            html.Div([
                                dcc.Graph(id='heatmap')
                            ], style={'textAlign': 'center'})
                ]),
        ], width=8),
        dbc.Col([
            html.Div([
                dcc.Upload(
                    id='upload-data',
                    children=html.Div([
                        dbc.Button('Upload File',color="info")
                        ], style={'textAlign': 'center'})
                ),
                html.P('Select coloring:'),
                dcc.Dropdown(options=[
                                {'label': 'Study/Cohort', 'value': 'cohort'},
                                {'label': 'Age', 'value': 'age'},
                                {'label': 'Sex', 'value': 'sex'},
                            ], 
                             value='cohort',id='color-dropdown'),
                html.P('Select studies to include in plot:'),
                dcc.Dropdown(options=studyIDs,id='cbiostudy-dropdown',multi=True),
                html.P('Feature weight: age'),
                dcc.Slider(
                    #id='age-weight',
                    id={'type':'feature-weight-slider','id':'age-weight','index':0},
                    min=0,
                    max=1,
                    step=0.1,
                    value=1        
                ),
                html.P('Feature weight: age at diagnoses'),
                dcc.Slider(
                    # id='age-at-diagnosis-weight',
                    id={'type':'feature-weight-slider','id':'age-at-diagnosis-weight','index':1},
                    min=0,
                    max=1,
                    step=0.1,
                    value=1        
                ),
                html.P('Feature weight: sex'),
                dcc.Slider(
                    #id='sex-weight',
                    id={'type':'feature-weight-slider','id':'sex-weight','index':2},
                    min=0,
                    max=1,
                    step=0.1,
                    value=1        
                ),
                html.P('Feature weight: mutation profile'),
                dcc.Slider(
                    #id='mutations-weight',
                    id={'type':'feature-weight-slider','id':'mutations-weight','index':3},
                    min=0,
                    max=1,
                    step=0.1,
                    value=1       
                ),
                html.P('Feature weight: ICD-code(s)'),
                dcc.Slider(
                    #id=('icd-weight'),
                    id={'type':'feature-weight-slider','id':'icds-weight','index':4},
                    min=0,
                    max=1,
                    step=0.1,
                    value=1        
                ),
                html.P('Select GO-annotation set:'),
                dcc.Dropdown(options=
                            [
                                {'label': 'GO Slim BP (recommended)', 'value': 'panther'},
                                {'label': 'GO BP (takes longer to compute - not filtered)', 'value': 'pygosemsim'},
                            ],                            
                             id='term-source-dropdown',
                             value='panther')
            ]),
            dbc.Row([
                dbc.Col([
                        html.P('Current feature weights:', style={'textAlign': 'center'}),
                        html.Pre(id='weights-paragraph'),
                        html.Div([
                            dcc.Markdown('''
                                        [Documentation & code base](https://github.com/janschneida/)
                                     ''')
                            ], style={'textAlign': 'center'}),                      
                ]),
                dbc.Col([                
                        html.P('Toggle OncoKB gene filter:', style={'textAlign': 'center'}),
                        daq.ToggleSwitch(
                            id = 'onco-kb-filter-switch',
                            label='Include all genes',
                            labelPosition='top',
                            value = False
                        ),
                        html.Div([
                            html.P(''),
                            dbc.Button('Save similarity scores & distance matrices',id='save-button',n_clicks=0,color="info"),
                            dcc.Download(id="download-scores")                
                        ], style={'textAlign': 'center'})
                ])
            ]),
        ],width=4)
    ])
])

# Define the app callback
@app.callback(
    Output('scatter-plot', 'figure'),
    Output('heatmap', 'figure'),
    Output('onco-kb-filter-switch','value'),
    Output('cbiostudy-dropdown','value',allow_duplicate=True),
    Input('weights-paragraph','children'),
    Input('color-dropdown','value'),
    Input('term-source-dropdown','value'),
    Input('onco-kb-filter-switch','value'),
    Input('cbiostudy-dropdown','value'),
    prevent_initial_call = True
)
def update_figure(weights,color,term_source,filter_genes,studies):
    triggered_id = ctx.triggered_id
    global patient_data_df
    global feature_weights
    global cbio_df
    global included_studies
    global source_df
    
    start = time.time()
    
    # recalc if new term source
    if triggered_id == 'term-source-dropdown':
        source_df=cu.getDistScores(source_df,feature_weights,calc_mut_dist=True,term_source=term_source)
    elif triggered_id == 'onco-kb-filter-switch':
        if filter_genes:
            filtered_patient_data = du.filterOncoGenesInPatientDataFrame(source_df)
            source_df=cu.getDistScores(filtered_patient_data,feature_weights,calc_mut_dist=True,term_source=term_source,filter_genes=True)
        else:
            if included_studies:
                cbio_studies_df = du.getDataFrameFromCBioPortalStudies(included_studies)
                source_df = pd.concat([cbio_studies_df,patient_data_df],ignore_index=True)
            source_df = cu.getDistScores(source_df,feature_weights,calc_mut_dist=True)
    elif triggered_id == 'cbiostudy-dropdown':
        included_studies = studies
        if len(studies) == 0:
            # no sstudies selected -> local data only
            source_df = patient_data_df.copy()
        else:
            cbio_studies_df = du.getDataFrameFromCBioPortalStudies(studies)
            # merge local patient data with external 
            source_df = pd.concat([cbio_studies_df,patient_data_df],ignore_index=True)
        source_df = cu.getDistScores(source_df,feature_weights,calc_dem_dist=True,calc_icds_dist=True,calc_mut_dist=True)
        # reset filter toggle to include all genes
        filter_genes = False
    if source_df.empty:
        return
    # set up scatter plot
    scatter_plot = px.scatter(data_frame=source_df,
                              x='x',
                              y='y',
                              color=color,
                              hover_name='patient_id',
                              hover_data=hover_data,
                              height=600
                              )
    scatter_plot.update_layout(transition_duration=700)
    scatter_plot.update_yaxes(title_text='',showticklabels=False)
    scatter_plot.update_xaxes(title_text='The displayed distances between the patient points are a '+
                     'representation of the relative similarity of the included patient cohorts.',
                     nticks=0, showticklabels=False)
    # set up heatmap
    x = list(source_df['patient_id'].astype(str).values)
    y = list(source_df['patient_id'].astype(str).values)
    heatmap = px.imshow(1 - cu.dist_scores,
          x=x,
          y=y,
          labels=dict(x="Patient ID", y="Patient ID", color="Similarity"),
          width=1000,
          height=1000
          )
    heatmap.update_xaxes(tickangle=-60)
    
    end = time.time()
    timee = end-start
    print('figure updating time for',triggered_id,timee)
    return scatter_plot, heatmap, filter_genes, included_studies

@app.callback(
    Output('weights-paragraph','children'),
    Input({'type': 'feature-weight-slider', 'id': ALL,'index': ALL}, 'value')
)
def update_weights_from_slider(value):
    # skip init call to callback
    global feature_weights
    start = time.time()
    if(not ctx.triggered_id):        
        return json.dumps(feature_weights,indent=4)
    global source_df
    global hover_data
    
    feature = str.replace(str(ctx.triggered_id['id']),'-',' ').replace(' weight','')
    weight = value[ctx.triggered_id['index']]
    feature_weights[feature] = weight
    
    if all(value == 0.0 for value in feature_weights.values()):
        raise dash.exceptions.PreventUpdate("All weights set to 0!")
    
    hover_data[feature] = weight > 0.0
    
    # no need to recalc anything, since only weights changed
    if feature == 'mutations' or feature == 'icds':
        source_df = cu.getDistScores(source_df,feature_weights)
    elif feature == 'age' or feature == 'age at diagnosis' or feature == 'sex':
        # only demographics left -> recalc dem dists
        source_df = cu.getDistScores(source_df,feature_weights,calc_dem_dist=True)
    end = time.time()
    feature_time=end-start
    print('time to adjust weights for',feature,':',feature_time)
    return json.dumps(feature_weights,indent=4)

@app.callback(
    Output('onco-kb-filter-switch','label'),
    Input('onco-kb-filter-switch', 'value')
)
def update_genes(value):
    if value:
        return 'Only include oncogenes/TSGs.'
    else:
        return 'Include all genes'

@app.callback(
    Output('download-scores','data'),
    Input('save-button','n_clicks'),
    prevent_initial_call = True
)
def save_similarities(n_clicks):

    path = os.getcwd()+'\\saved_files\\'
    writer_name = path+'saved_scores'+'.xlsx'
    with pd.ExcelWriter(writer_name) as writer:
        # Write each DataFrame to a separate sheet
        source_df.to_excel(writer, sheet_name='Patient data with coordinates', index=False)
        pd.DataFrame(cu.dist_scores,columns=source_df['patient_id'],index=source_df['patient_id']).to_excel(writer, sheet_name='Pairwise Distscores', index=True)
        pd.DataFrame(feature_weights.items(),columns=['feature','weight']).to_excel(writer, sheet_name='Used Weights', index=False)
        writer.save()
        return dcc.send_file(writer_name)

@app.callback(
    Output('cbiostudy-dropdown','options'),
    Output('cbiostudy-dropdown','value',allow_duplicate=True),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    #State('cbiostudy-dropdown','value'),
    prevent_initial_call = True
    )
def update_output(content,filename):#, list_of_names, included_studies):
    if content is not None:
        df = parse_content(content,filename)
    write_dir = 'resources/local_studies'
    cohort = df['cohort'][0]
    df.to_csv(write_dir+'/'+cohort)
    studyIDs.append(cohort)
    included_studies.append(cohort) 
    return studyIDs, included_studies

def parse_content(content, filename):
    content_type, content_string = content.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
        return df
    except Exception as e:
        print(e)


if __name__ == '__main__':
    # NOTE runs code twice with true
    app.run_server(debug=False)