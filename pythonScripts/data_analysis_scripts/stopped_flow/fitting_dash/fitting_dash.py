# dash app

import dash
import sys
from dash import dcc, html, Input, Output, State
sys.path.append("/Users/eleanordunietz/Desktop/Stanford/pythonScripts/data_analysis_scripts/stopped_flow")
from SF_analysis_fitting_ln import *
from SF_analysis_processing_ln import *

# Utility function to get concentration dropdowns based on substrate being fit
def get_concentration_dropdowns_for_substrate(processed_data, substrate):
    # Filter out rows with '-' in the concentration column and match the selected substrate
    processed_data_filtered = processed_data[
        (processed_data['substrate_concentration'] != '-') & 
        (processed_data['substrate'] == substrate)
    ]
    # Convert concentration values to float for proper sorting and get unique values
    unique_concentrations = processed_data_filtered['substrate_concentration'].astype(float).unique()
    # Sort the concentrations
    sorted_concentrations = sorted(unique_concentrations)
    return sorted_concentrations

# utility function to extract input parameters 
def select_parameters(substrate_concentration, fitting_params_dict):
    # Retrieve parameters for the specified concentration or the 'default' set
    params = fitting_params_dict.get(substrate_concentration, fitting_params_dict['default'])
    return params['initial_guesses'], params['time_cutoff'], params['bounds']

def run_dash(processed_data, substrate, input_parameters_dict, wavelength_to_fit, pH, solvent, date, kinetic_model):
    # Initialize your Dash app
    app = dash.Dash(__name__)

    options = get_concentration_dropdowns_for_substrate(processed_data, substrate)
    # Assuming 'options' is not empty, prepare dropdown options and set the default value
    if options:
        dropdown_options = [{'label': str(value), 'value': value} for value in options]
        default_value = options[0]  # Safely assign the first value as the default
    else:
        dropdown_options = []
        default_value = None  # Safely handle the case where 'options' is empty

    # App layout
    app.layout = html.Div([
        html.Label('Select substrate concentration (mM):'),
        dcc.Dropdown(
            id='concentration-dropdown',
            options=dropdown_options,
            value=default_value
        ),
        dcc.Graph(id='fitting-results-graph'),
        html.Button('Prev', id='prev-button', n_clicks=0),
        html.Button('Next', id='next-button', n_clicks=0),
        # Hidden div for storing current index of experiment
        html.Div(id='experiment-index', style={'display': 'none'}, children=0)
    ])

    @app.callback(
        [Output('fitting-results-graph', 'figure'),
        Output('experiment-index', 'children'),
        Output('prev-button', 'disabled'),
        Output('next-button', 'disabled')],
        [Input('concentration-dropdown', 'value'),
        Input('prev-button', 'n_clicks'),
        Input('next-button', 'n_clicks')],
        [State('experiment-index', 'children')]
    )
    def update_graph(concentration, prev_clicks, next_clicks, index):
        ctx = dash.callback_context

        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

        experiments = get_experiments_by_criteria(processed_data, pH=pH, substrate=substrate, solvent=solvent, date=date, substrate_concentration=str(concentration))

        # If the dropdown changes, reset index to 0
        if triggered_id == 'concentration-dropdown':
            index = 0
        elif triggered_id == 'prev-button' and index > 0:
            index -= 1
        elif triggered_id == 'next-button' and index < len(experiments) - 1:
            index += 1

        index = max(0, min(index, len(experiments) - 1))

        if not experiments.empty:
            row = experiments.iloc[index]
            # Call select_parameters once and unpack the input parameters defined
            initial_guesses, time_cutoff, bounds = select_parameters(row['substrate_concentration'], input_parameters_dict)
            _, figure = fit_kinetics_to_experiment(processed_data, row['push'], wavelength_to_fit,  
                                                time_cutoff=time_cutoff, 
                                                initial_guesses=initial_guesses,
                                                bounds=bounds,
                                                plot_fit=True,
                                                kinetic_model=kinetic_model)
        else:
            figure = go.Figure()

        # Determine whether to disable Prev or Next buttons
        disable_prev = index <= 0
        disable_next = index >= len(experiments) - 1

        return figure, index, disable_prev, disable_next

    app.run_server(debug=True, port=8054)
