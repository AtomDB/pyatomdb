#!/usr/bin/env python3

import os
import numpy
import pandas as pd
import plotly.graph_objects as go

from dash import Dash, html, dcc, callback, Output, Input, State, ctx, no_update
import dash_ag_grid as dag
import localsettings 
import astropy
import pyatomdb
import pickle
import dash_bootstrap_components as dbc
# ----------------------------------------------------------------------
# Constants and data
# ----------------------------------------------------------------------

Tlist = numpy.linspace(4, 9, 51)
Vellist = numpy.linspace(0, 1000, 11)

# Requires DATAFILEPATH env var to be set, e.g.:
# export DATAFILEPATH=/path/to/atomdb/data
responses = pd.read_csv(os.path.expandvars('$DATAFILEPATH/responses/manifest.txt'))

KBOLTZ = pyatomdb.const.KBOLTZ
HC_IN_KEV_A = pyatomdb.const.HC_IN_KEV_A

#linelist = astropy.io.fits.open(os.path.expandvars('$ATOMDB/apec_linelist.fits'))

linelist = pickle.load(open(os.path.expandvars('$ATOMDB/apec_linelist.pkl'), 'rb'))
CIE = pyatomdb.spectrum.CIESession()

verbose=False # for printing to the screen to debug
# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------

def Ztoelsymb(Z):
    """
    Returns element symbol of element with nuclear charge Z.

    PARAMETERS
    ----------
    Z  - nuclear charge of element (e.g. 6 for carbon)

    RETURNS
    -------
    element symbol (e.g. "C" for carbon)
    """
    elsymb = (
        'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca',
        'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr',
        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
        'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
        'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
        'Pa', 'U'
    )

    if Z < 1 or Z > 92:
        print("Z must be between 1 and 92. You have given Z= " + repr(Z))
        return -1
    return elsymb[Z - 1]


def int2roman(number):
    numerals = {
        1: "I", 4: "IV", 5: "V", 9: "IX",
        10: "X", 40: "XL", 50: "L", 90: "XC",
        100: "C", 400: "CD", 500: "D", 900: "CM",
        1000: "M"
    }
    result = ""
    for value, numeral in sorted(numerals.items(), reverse=True):
        while number >= value:
            result += numeral
            number -= value
    return result


def spectroscopic_name(Z, z1):
    """
    Converts Z,z1 to spectroscopic name, e.g. 6,5 to "C V"
    """
    elsymb = Ztoelsymb(Z)
    roman = int2roman(z1)
    return elsymb + ' ' + roman


def transform_value(value):
    return 10 ** value


# ----------------------------------------------------------------------
# Spectrum calculation helpers
# ----------------------------------------------------------------------

def calc_precalc_spectrum(telescope, response, temperature_value, broadening_value):
    """
    Use pre-calculated spectra from $DATAFILEPATH/precalc/<telescope>/<response>/
    master.npz and dat<iT>.npz files.
    """
    pre_calc_dir = os.path.expandvars(
        "$DATAFILEPATH/precalc/%s/%s" % (telescope, response)
    )

    masterfile = "%s/master.npz" % pre_calc_dir
    masterdata = numpy.load(masterfile)

    eedges = masterdata['ebins']

    # find the indices
    iT = numpy.argmin(numpy.abs(10 ** temperature_value - masterdata['Tlist']))
    iVel = numpy.argmin(numpy.abs(broadening_value - masterdata['Vellist']))

    # data file
    specdatafile = "%s/dat%i.npz" % (pre_calc_dir, iT)
    specdata = numpy.load(specdatafile)

    lines = specdata['lines']
    spec = specdata['spec'][iVel, :]
    ewidth = eedges[1:] - eedges[:-1]
    spec = spec / ewidth
    spec = numpy.append(spec[0], spec)

    return eedges, spec, lines


# ----------------------------------------------------------------------
# Dash app
# ----------------------------------------------------------------------

app = Dash(
    __name__,
    requests_pathname_prefix='/apps/newapp/ai/',
    routes_pathname_prefix='/apps/newapp/ai/',
)
server = app.server 

# Layout
app.layout = html.Div([
    
    # Header 
    html.H1(
            children='Interactive Plotting with AtomDB',
            className="app-header"
    ),
     
    # First Box: Plotting Spectra
    
    html.Div(className = 'atomdbinnerbox',
        children=[
            html.P("Here you can plot spectra for a range of instruments and temperatures, list the strong lines in a wavelength range, and by clicking on the points on the graph or on the table"),
            html.Details([
                html.Summary('Plot Interactive Spectra'),
                html.Div(id = 'interactive-plot-div',
                    className = 'atomdbinnerbox',
                    children=[
                        # Telescope and instrument selection
                        html.Div(className = 'atomdbbox',
                            id='telescope-div',
                            children=[
                                'Telescope:',
                                dcc.Dropdown(
                                    options=responses['Telescope'].unique(),
                                    value='Chandra cy28',
                                    id='telescope-dropdown'
                                ),
                                
                                'Instrument:',
                                dcc.Dropdown(
                                    options=list(
                                        responses['DisplayName'][responses['Telescope'] == 'Chandra cy28']
                                        ),
                                    value='ACIS-S HEG -1',
                                    id='response-dropdown'
                                ),
                            #],
                        #), # end DIV telescope-div
                       
                        #html.Div(className = 'atomdbbox',
                            #id='graph-settings-div',
                            #children=[
#                                html.Div(id='my-temperature',
                                    html.P(id='my-temperature'),
                                    dcc.Slider( # Temperature gauge
                                        min=Tlist[0],
                                        max=Tlist[-1],
                                        step=(Tlist[-1] - Tlist[0]) / (len(Tlist) - 1),
                                        id='slider-temperature',
                                        value=7,
                                        marks={i: f"10^{i}" for i in range(4, 10)},
                                    ), 
 #                                   ),
 #                    ],
                    
  #              ),

  #              html.Div(
  #                  [
   #                     html.Div(id='my-broadening'),
                                    html.P(id='my-broadening'),
                                    dcc.Slider(
                                        # Velocity slider
                                        min=Vellist[0],
                                        max=Vellist[-1],
                                        step=(Vellist[-1] - Vellist[0]) / (len(Vellist) - 1),
                                        id='slider-broadening',
                                        value=100,
                                        marks={i: f"{i}km/s" for i in range(0, 1001, 200)},
                                        ),
#                            ], # end settings div
                    
                        

#                      html.Div(
                    #[
                        #html.H2("Graph Settings:"),
                                    dcc.RadioItems(
                                        ['Y-linear', 'Y-log'],
                                        'Y-linear',
                                        inline=True,
                                        id='radio-yaxis'
                                    ),
                                    dcc.RadioItems(
                                        ['X-linear', 'X-log'],
                                        'X-linear',
                                        inline=True,
                                        id='radio-xaxis'
                                    ),
                                    dcc.RadioItems(
                                        ['keV', 'Å'],
                                        'keV',
                                        inline=True,
                                        id='radio-specunits'
                                    ), 
                                    html.Button(
                                        'Update Graph',
                                        id='update-button',
                                    ),
                            ], # end graph-settings-div children
                    
                        ), # end DIV graph-settings-div

                        html.Div(id='graph-div',
                            children=[
                                dcc.Graph(id='graph-content', responsive=False)
                            ]
                        ),
                ]
            ),]),


   

        html.Details([
            html.Summary('List Strong Lines in Region'),
            html.Div(id = 'linelist-controls-div',
                className='atomdbbox', # This is the Linelist Div
                children=\
                    [
                    html.H1('List Lines in a Wavelength Region',
                        className="app-header",
                            ),
                    
                    #id='linelist-controls-div',\
#                         children=\
#                    [
                   
                    'Spectral Range: Min:',
                    dcc.Input(0.0,id='linelist-min', type='number', style={'width':'60px'}),
                    ' Max:', 
                    dcc.Input(0.0,id='linelist-max', type='number', style={'width':'60px'}),
                   
                    dcc.RadioItems(
                        ['Å', 'keV'],
                        'Å',
                        inline=True,
                        id='linelist-units'
                    ),
                    'Min Emissivity to display:', 
                    dcc.Input(1e-18,id='linelist-mineps', type='number', style={'width':'80px'}),
                    html.Br(),

                    dcc.Checklist(options={
                                  'yes':'Show lines at specific temperature?'},
                                  id='linelist-temperature-checkbox'
                                 ),
                    html.Div(id='linelist-temperature-div',
                             children=[           
                               html.P(id='linelist-temperature-P'),
                               dcc.Slider( # Temperature gauge
                                 min=Tlist[0],
                                 max=Tlist[-1],
                                 step=(Tlist[-1] - Tlist[0]) / (len(Tlist) - 1),
                                 id='linelist-temperature-slider',
                                 value=7,
                                 marks={i: f"10^{i}" for i in range(4, 10)},
                               ),
                             ],
                             style={'display':'none'}
                           ), 


                    html.Button('Update Linelist',
                        id='linelist-button'
                    ),
                        
                    html.Div(id='linelist-container', children=[dag.AgGrid(id='linelist-table',
                           rowData = [],
                           columnDefs = [{'headerName':'Wavelength',"wrapText": True, "autoHeight": True, "flex":1, "cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Lambda', 'headerName':'Å',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5f')(params.value)"}}]},
                {'headerName':'Energy',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Energy', 'headerName':'keV',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5f')(params.value)"}}]},
                {'headerName':'Ion', 'field':'Spectroscopic',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}},
                {'field':'UpperLev', 'headerName':'Upper Level' ,"wrapText": True, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "flex":1},
                {'field':'LowerLev', 'headerName':'Lower Level' ,"wrapText": True, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "flex":1},
                {'headerName':'Peak Emissivity',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'PeakEmissivity', 'headerName':'ph cm^3 s^-1',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"} , "valueFormatter": {"function": "d3.format(',.5e')(params.value)"}}]},
                {'headerName':'Peak Temperature',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'headerName':'K', 'field':'PeakTemperature',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5e')(params.value)"}}]},
                 ],
                           dashGridOptions={
                           "enableCellTextSelection": True, 
                           "ensureDomOrder": True,
                           "theme": {
                               "function":"themeBalham.withParams({accentColor:'#772007',\
                                                                      headerTextColor:'#772007'})"}
                           },
                           csvExportParams={
                           "fileName": "atomdb_linelist_blank.csv"
                           },
                           columnSize="autoSize",
            
                           )]
                    ),
                    html.P("Found zero lines", id='linelist-output'),
                    html.Button("Download As CSV", id="linelist-csv-button", n_clicks=0),
                    html.Button("Show Selected Transition Info", id="linelist-selected-button", n_clicks=0)
                    ]
            )
                
        ]),
                      
        html.Details([
            html.Summary('Transition Details'),
            html.Div(id = 'transition-details-div',
                     className='atomdbbox',
                     children=[],
                     style={'width':'500px'}
            ),
            
        ],id="transition-details")             
                    
                   

        
     
    ])
])


                      
         

# ----------------------------------------------------------------------
# Callbacks
# ----------------------------------------------------------------------

@callback(
    Output('response-dropdown', 'options'),
    Output('response-dropdown', 'value'),
    Input('telescope-dropdown', 'value')
)
def set_telescope_options(telescope):
    new_options = list(responses['DisplayName'][responses['Telescope'] == telescope])
    return new_options, new_options[0]


@callback(
    Output('my-temperature', 'children'),
    Input('slider-temperature', 'value')
)
def update_temperature(temperature_value):
    return "Temperature: %.2e K, %.3f keV" % (
        10 ** temperature_value,
        10 ** temperature_value * KBOLTZ
    )


@callback(
    Output('my-broadening', 'children'),
    Input('slider-broadening', 'value')
)
def update_broadening(broadening_value):
    return "Velocity broadening: %i km/s" % (broadening_value)


@callback(
    Output('graph-content', 'figure'),
    State('telescope-dropdown', 'value'),
    State('response-dropdown', 'value'),
    State('slider-temperature', 'value'),
    State('slider-broadening', 'value'),
    Input('update-button', 'n_clicks'),
    State('radio-yaxis', 'value'),
    State('radio-xaxis', 'value'),
    State('radio-specunits', 'value')
)
def update_graph(telescope, response, temperature_value, broadening_value,
                 n_clicks, radioyaxis, radioxaxis, specunits):
    # Use pre-calculated spectra
    eedges, a, lines = calc_precalc_spectrum(
        telescope, response, temperature_value, broadening_value
    )

    # Renormalize the lines
    if len(lines) > 0:
        imaxeps = numpy.argmax(lines['Epsilon'])
        maxeps = lines['Epsilon'][imaxeps]
        ibin = numpy.max(
            numpy.where(eedges < HC_IN_KEV_A / lines['Lambda'][imaxeps])[0]
        )
        normfactor = a[ibin] / maxeps
    else:
        normfactor = 1.0

    lines['Epsilon'] *= normfactor

    if specunits.lower() == 'kev':
        linepositions = HC_IN_KEV_A / lines['Lambda']
        xtitle = "Energy (keV)"
    else:  # 'Å'
        # convert edges to Å, reverse order, etc.
        eedges_ang = HC_IN_KEV_A / eedges[::-1]
        a_ang = a[1:]
        a_ang = a_ang[::-1]
        a_ang = numpy.append(a_ang[0], a_ang)
        eedges = eedges_ang
        a = a_ang
        linepositions = lines['Lambda']
        xtitle = "Wavelength (Å)"

    ion_symbols = numpy.zeros(len(lines), dtype='<U250')
    for iline, l in enumerate(lines):
        ion_symbols[iline] = (
                        f'{spectroscopic_name(l["Element"], l["Ion"])} {l["UpperLev"]}->{l["LowerLev"]}'
        )

    fig = go.Figure()
    fig.update_layout(
        xaxis_title=xtitle,
        yaxis_title='Emissivity*Aeff (ph cm<sup>5</sup> s<sup>1</sup> kev<sup>-1</sup>)',
        plot_bgcolor='#FDEAE5',
        showlegend=False
    )

    fig.add_trace(
        go.Scatter(
            x=linepositions,
            y=lines['Epsilon'],
            mode='markers',
            hovertext=ion_symbols,
            marker_color='#80ADA0',
            error_y={
                'type': 'data',
                'symmetric': False,
                'arrayminus': lines['Epsilon'],
                'array': [0] * len(lines['Epsilon']),
                'thickness': 0.5,
                'width': 0,
            },
            hoverinfo='text',
            hoverlabel={'bgcolor': '#80ADA0',
                        'showarrow':True,
                        'font':{'color':'#5E2BFF'}},
        )
    )

    fig.update_yaxes(exponentformat='E')

    fig.add_trace(
        go.Scatter(
            x=eedges,
            y=a,
            line_shape='vh',
            name='spectrum',
            mode='lines',
            line_color='#772007'
        )
    )
    
    #fig.add_

    if radioyaxis == 'Y-log':
        fig.update_layout(yaxis_type="log")
    else:
        fig.update_layout(yaxis_type="linear")

    if radioxaxis == 'X-log':
        fig.update_layout(xaxis_type="log")
    else:
        fig.update_layout(xaxis_type="linear")

#    print("a",a)
    return fig

@callback(
     Output('linelist-container', 'children'),
     Output('linelist-output','children'),
     Input('linelist-button', 'n_clicks'),
     State('slider-temperature', 'value'),
     State('linelist-min', 'value'),
     State('linelist-max', 'value'),
     State('linelist-units', 'value'),
     State('linelist-mineps', 'value'),
     State('linelist-temperature-checkbox', 'value'),
     State('linelist-temperature-slider', 'value'),
     
 )
#def tmp(n_clicks, temperature_value, spec_min_in,
#                   spec_max_in, spec_units_in):
#  print('n_clicks', n_clicks)

def make_linelist(n_clicks, temperature_value, spec_min_in,
                   spec_max_in, spec_unit_in, mineps, do_single_temp_in, single_temp_in):

    T = 10** temperature_value
    spec_min_in=float(spec_min_in)
    spec_max_out=float(spec_max_in)
    
    
    if spec_unit_in=='Å':
        spec_units = 'Å'
        spec_min = spec_min_in
        spec_max = spec_max_in
        
    else:
        spec_units='keV'
        # protect against divide by zero
        spec_max_in=max([1e-10, spec_max_in])
        spec_min_in=max([1e-10, spec_min_in])
        spec_min = HC_IN_KEV_A/spec_max_in
        spec_max = HC_IN_KEV_A/spec_min_in
    
    
    do_single_temp=False
    if do_single_temp_in is not None:
      if "yes" in do_single_temp_in:
        do_single_temp=True
        single_temp = 10**single_temp_in # in K
   
    # This is the data for a single line
    
    if not do_single_temp:
    
      columnDefs = [{'headerName':'Wavelength',"wrapText": True, "autoHeight": True, "flex":1, "cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Lambda', 'headerName':'Å',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5f')(params.value)"}}]},
                {'headerName':'Energy',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Energy', 'headerName':'keV',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5f')(params.value)"}}]},
                {'headerName':'Ion', 'field':'Spectroscopic',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}},
                {'field':'UpperLev', 'headerName':'Upper Level' ,"wrapText": True, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "flex":1},
                {'field':'LowerLev', 'headerName':'Lower Level' ,"wrapText": True, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "flex":1},
                {'headerName':'Peak Emissivity',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'PeakEmissivity', 'headerName':'ph cm^3 s^-1',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"} , "valueFormatter": {"function": "d3.format(',.5e')(params.value)"}}]},
                {'headerName':'Peak Temperature',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'headerName':'K', 'field':'PeakTemperature',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5e')(params.value)"}}]},
                 ]
               
      s= linelist[(linelist['Lambda']>spec_min) & (linelist['Lambda']<spec_max) & (linelist['PeakEmissivity']>mineps)]
      grid=[dag.AgGrid(id='linelist-table',
                           rowData = s.to_dict('records'),
                           columnDefs = columnDefs,
                           columnSize="autoSize",
            
                           dashGridOptions={
                           "enableCellTextSelection": True, 
                           "ensureDomOrder": True,
                           "theme": {
                               "function":"themeBalham.withParams({accentColor:'#772007',\
                                                                      headerTextColor:'#772007'})"},
                           'rowSelection': {
                               'mode': 'singleRow',
                               'enableClickSelection': True,
                           },
    
                           },
                           csvExportParams={
                           "fileName": "atomdb_linelist_%f_%f_%e.csv"%(spec_min, spec_max, mineps)
                           },
                           )]

      textout="Found %i lines"%(len(s))
    
      return(grid, textout)


    else:
    
      columnDefs = [{'headerName':'Wavelength',"wrapText": True, "autoHeight": True, "flex":1, "cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Lambda', 'headerName':'Å',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5f')(params.value)"}}]},
                {'headerName':'Energy',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Energy', 'headerName':'keV',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "valueFormatter": {"function": "d3.format(',.5f')(params.value)"}}]},
                {'headerName':'Ion', 'field':'Spectroscopic',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}},
                {'field':'UpperLev', 'headerName':'Upper Level' ,"wrapText": True, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "flex":1},
                {'field':'LowerLev', 'headerName':'Lower Level' ,"wrapText": True, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"}, "flex":1},
                {'headerName':'Emissivity',"wrapText": True, "autoHeight": True, "flex":1,"cellStyle": {"whiteSpace": "pre-wrap"},\
                   'children':[{'field':'Epsilon', 'headerName':'ph cm^3 s^-1',"wrapText": True, "flex":1, "autoHeight": True,"cellStyle": {"whiteSpace": "pre-wrap"} , "valueFormatter": {"function": "d3.format(',.5e')(params.value)"}}]},
                 ]
      s = CIE.return_linelist(single_temp, [spec_min, spec_max], specunit='A', teunit='K', format_data='pandas', nearest=True)
      #s= linelist[(linelist['Lambda']>spec_min) & (linelist['Lambda']<spec_max) & (linelist['PeakEmissivity']>mineps)]
      grid=[dag.AgGrid(id='linelist-table',
                           rowData = s.to_dict('records'),
                           columnDefs = columnDefs,
                           dashGridOptions={
                           "enableCellTextSelection": True, 
                           "ensureDomOrder": True,
                           "theme": {
                               "function":"themeBalham.withParams({accentColor:'#772007',\
                                                                      headerTextColor:'#772007'})"},
#                           "rowSelection": {'mode': 'singleRow', 'enableClickSelection': True},
                           },
                           csvExportParams={
                           "fileName": "atomdb_linelist_%f_%f_%e.csv"%(spec_min, spec_max, mineps)
                           },
                           )]

      textout="Found %i lines"%(len(s))
    
      return(grid, textout)



                     
@callback(
    Output("linelist-table", "exportDataAsCsv"),
    Input("linelist-csv-button", "n_clicks"),
    )

def export_data_as_csv(n_clicks):
    if n_clicks:
        return True
    return False              


@callback(
    Output("linelist-temperature-div","style"),
    Input("linelist-temperature-checkbox","value")
)

def update_linelist_temperature_div_visibility(values):
    
    if values is None:
        return {'display':'none'}
    if "yes" in values:
        return {'display':'block'}
    return {'display':'none'}

# now detect clicks on the table
@callback(
  Output("transition-details-div","children"),
  Output("transition-details", "open"),
  #State("linelist-table","data"),
  Input("linelist-table","selectedRows"),
  Input("graph-content","clickData"),
)
def display_transition_data(selectedRows, clickData): 
  trigger_id = ctx.triggered_id
  if trigger_id=='linelist-table':
  
    if not selectedRows:
          return (False, False)
    # IN HERE MAKE A MAGIC FUNCTION....
    out = get_transition_information(selectedRows[0])
    
    
  elif trigger_id=='graph-content':
    if clickData['points'][0]['curveNumber']==0:
      rowInfo={}
    
      ht = clickData['points'][0]['hovertext']
      tmp = ht.split()
      elsymb = tmp[0]
      chgrm = tmp[1]
      upper, lower = tmp[2].split('->')
      rowInfo['Element'] = pyatomdb.atomic.elsymb_to_Z(elsymb)
      rowInfo['Ion'] = pyatomdb.atomic.roman_to_int(chgrm)
      rowInfo['UpperLev'] = int(upper)
      rowInfo['LowerLev'] = int(lower)
      out = get_transition_information(rowInfo)
    else:
      return (no_update,no_update)
    #out = False
  else:
    return (no_update,no_update)
  # now make a table!

  return(out, True)

    


def parse_reference_string(reference, txt=''):
  if txt=='':
    txt = reference
  if "NIST" in reference:
    web = 'https://www.nist.gov/pml/atomic-spectra-database'
  elif "afoster" in reference:
    web = 'http://adsabs.harvard.edu/abs/2012ApJ...756..128F'
  else:
    web= "http://adsabs.harvard.edu/abs/%s"%(reference)
    
  return "<a href=\"http://adsabs.harvard.edu/abs/%s\">"%(reference)+\
         "%s</a>"%(txt) 

def get_transition_information(rowInfo):
    """
    This is the one stop shop call to convert a row from the datatable in AtomDB
    Into a series of Dash AgGrid tables
    
    """
    
    # If no data, return nothing (prevents errors when a row has not been selected)
    if len(rowInfo)==0:
        return None
    
    # For storing open files, prevents repeated opening
    datacache={}
    
    # Get the raw information from the files
    dat, is_dr = pyatomdb.atomdb.get_transition_information(rowInfo['Element'], rowInfo['Ion'], rowInfo['UpperLev'], rowInfo['LowerLev'], datacache=datacache)
    
    # Add the appropriate formatting (e.g. replace data column labels with more verbose text)
    dat2 = pyatomdb.atomdb.make_transition_information_table(dat,rowInfo['Element'], rowInfo['Ion'], is_dr)
    
    if verbose:
      print('dat', dat)
      print('dat2', dat2)
    
    # Make into Dash AgGrid tables, which will be returned in 'ret'
    ret = []
    # add in a header
    ret.append(html.H1("Transition Information for %s %i→%s"%(pyatomdb.atomic.spectroscopic_name(rowInfo['Element'], rowInfo['Ion']), rowInfo['UpperLev'], rowInfo['LowerLev'])))
    # First print out the upper and lower levels
    for key in ['up','lo']:
      try:
        tmp = dat2[key]
        # Need to convert to the correct list of dicts format expected by AgGrid
        cdat = []
        for i in range(len(tmp['labels'])):
          cdat.append({'labels':tmp['labels'][i], 'data':tmp['data'][i]})
        
        # Define the columns
        columnDefs = [{'field': "labels", 'headerName': tmp['title']},
                      {'field': "data", 'headerName': '', 
                       'cellRenderer': 'markdown',
                       'linkTarget': '_blank' }]
        
        # add the graph to the returned list
        ret.append(dag.AgGrid(
            id="atomic-data-%s-table"%(key),
            rowData=cdat,
            #defaultColDef={"filter": True},
            columnDefs=columnDefs,
            columnSize="autoSize",
            dashGridOptions={"animateRows": False,
                             "domLayout": "autoHeight",
                             "theme": {
                               "function":"themeBalham.withParams({accentColor:'#772007',\
                                                                      headerTextColor:'#772007'})"},
            },
            style = {"height": None}
           ))
      except:
        raise()
        pass

    for key in dat2.keys():
      if key in ['up','lo']: continue
      try:
        tmp = dat2[key]
        cdat = []
        for i in range(len(tmp['descr'])):
          cdat.append({'labels':tmp['descr'][i], 'data':tmp['data'][i]})
          
        columnDefs = [{'field': "labels", 'headerName': tmp['title']},
                      {'field': "data", 'headerName': '', 
                       'cellRenderer': 'markdown',
                       'linkTarget': '_blank' }]
        
        ret.append(dag.AgGrid(
            id="atomic-data-%s-table"%(key),
            rowData=cdat,
            #defaultColDef={"filter": True},
            columnDefs=columnDefs,
#            columnSize="sizeToFit",\
            columnSize="autoSize",
            dashGridOptions={"animateRows": False,
                             "domLayout": "autoHeight",
                             "theme": {
                               "function":"themeBalham.withParams({accentColor:'#772007',\
                                                                      headerTextColor:'#772007'})"},
            },
            style = {"height": None}
           ))
      except:
        raise()
        pass


    return(ret)  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def make_transition_information_table(dat):
 
    return False
  # get the data

       #dag.AgGrid(id="line-table")
# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------

if __name__ == '__main__':
    # Access in browser at http://127.0.0.1:8051/
    app.run(host='127.0.0.1', port=8054, debug=True, use_reloader=True)
