import sys
sys.path.append("/Users/eleanordunietz/Desktop/Stanford/pythonScripts/data_analysis_scripts/stopped_flow")
from SF_analysis_processing_ln import *
import plotly.express as px
import plotly.graph_objs as go
import pandas as pd

# Function to fetch the time range for an experiment
def get_time_range_for_experiment(key, substrate, pH=None, substrate_concentration=None, solvent=None, index=None, time_step=20, time_range=None, wavelength_plotting_range=None):
    # Fetch experiments based on the criteria
    criteria = {'substrate': substrate, 'pH': pH, 'substrate_concentration': substrate_concentration, 'solvent': solvent}

    # Fetch experiments based on the criteria
    experiments = get_experiments_by_criteria(key, **criteria)

    if experiments.empty or index is None or index >= len(experiments):
        return None

    # Fetch the time range from the first experiment
    min_time = experiments.iloc[index]['data']['Time'].min()
    max_time = experiments.iloc[index]['data']['Time'].max()

    return min_time, max_time

def plot_wavelength_vs_intensity(key, substrate, pH=None, substrate_concentration=None, solvent=None, index=None, time_step=10, time_range=None, subtract_baseline_flag=False, wavelength_plotting_range=None):
    # Build the criteria dictionary with only non-None values
    criteria = {'substrate': substrate, 'pH': pH, 'substrate_concentration': substrate_concentration, 'solvent': solvent}

    # Fetch experiments based on the criteria
    experiments = get_experiments_by_criteria(key, **criteria)

    if experiments.empty or index is None or index >= len(experiments):
        return go.Figure()  # Return an empty figure if no valid experiment is found

    fig = go.Figure()

    # Select the specific experiment based on the index
    experiment = experiments.iloc[index]

    # Get push number and date for the title
    push_number = experiment.get('push', 'Unknown')
    experiment_date = experiment.get('date', 'Unknown Date')

    data = experiment['data']
    data['Time'] = pd.to_numeric(data['Time'], errors='coerce')
    data.dropna(subset=['Time'], inplace=True)

    # Apply time cutoff and filtering
    if time_range is not None:
        start_time, time_cutoff = time_range
        data = filter_by_time_cutoff(data, time_cutoff, start_time)

    # Baseline Subtraction if enabled
    if subtract_baseline_flag:
        baseline = find_baseline_for_push(key, experiment['push'])
        if baseline is not None:
            data = subtract_baseline(data, baseline)
    
    if not data.empty:
        time_min = data['Time'].min()
        time_max = data['Time'].max()
        time_points = data['Time'].unique()
        selected_time_points = time_points[::time_step]

        # Determine wavelength range
        wavelength_columns = [col for col in data.columns if col != 'Time']
        if wavelength_plotting_range is not None:
            first_wavelength, last_wavelength = wavelength_plotting_range
            wavelength_columns = [col for col in wavelength_columns if first_wavelength <= float(col) <= last_wavelength]
        else:
            first_wavelength = float(wavelength_columns[0])
            last_wavelength = float(wavelength_columns[-1])

        # Viridis color scale
        color_scale = px.colors.sequential.Sunsetdark

        # Plotting data for each time point
        for time_point in selected_time_points:
            time_point_data = data[data['Time'] == time_point][wavelength_columns]
            color_value = (time_point - time_min) / (time_max - time_min)  # Normalize time value for color scale
            color = color_scale[int(color_value * (len(color_scale) - 1))]  # Map to color scale
            fig.add_trace(go.Scatter(
                x=time_point_data.columns.astype(float),
                y=time_point_data.iloc[0].values,
                mode='lines',
                line=dict(color=color, width=2),
                showlegend=False  # Hide individual line legends
            ))

        # Add a separate scatter trace for the color bar
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                colorscale='Sunsetdark',
                cmin=time_min,
                cmax=time_max,
                colorbar=dict(title="Time"),
                size=10
            ),
            hoverinfo='none',  # Hide hover info
            showlegend=False  # Ensure this trace does not appear in the legend
        ))

    # Customize the layout of the Plotly figure
    fig.update_layout(
        title={
            'text': f"Experiment: {push_number} (Date: {experiment_date})",
        #     'font': {
        #         'family': 'Arial',
        #         'size': 18
        # }
    },
        xaxis_title="Wavelength (nm)",
        yaxis_title="Intensity",
        xaxis=dict(range=[first_wavelength, last_wavelength]),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(
            family="Arial",
            size=12,
            color="#E5ECF6",
        ),
        height=450,
        width=600

    )

    return fig


def plot_specified_wavelength_traces(key, push_number, wavelengths, time_cutoff=None, start_time=None, xaxis_type='linear'):
    """
    Plots specified wavelength traces from the dataset using Plotly.
    Args:
    - key: The key structure containing the experiment data.
    - push_number: The specific push number to plot data for.
    - wavelengths: A list of desired wavelengths to plot.
    """
    # Find the experiment with the given push number
    experiment = key.loc[(key['push'] == push_number) & key['data'].notnull()].iloc[0]
    if experiment is None:
        print(f"No data available for push number {push_number}.")
        return

    data = experiment['data']

    # Apply time cutoff if specified
    if time_cutoff is not None or start_time is not None:
        filtered_data = filter_by_time_cutoff(data, time_cutoff, start_time)
        if filtered_data is None:
            return  # Stop the function if the filtered data is None
        data = filtered_data

    # Create a Plotly figure
    fig = go.Figure()

    # Plot each specified wavelength
    for desired_wavelength in wavelengths:
        # Find the closest actual wavelength to the desired one
        closest_wavelength = find_closest_wavelength(data, desired_wavelength)

        # Add the trace to the figure
        fig.add_trace(go.Scatter(x=data['Time'], y=data[closest_wavelength], mode='lines', name=f"{closest_wavelength} nm"))

    # Update the layout with the x-axis type
    fig.update_layout(
        title={
            'text': f"Wavelength Traces for {push_number}",
            'font': {
                'family': 'Arial',
                'size': 18
            }
        },
        xaxis_title="Time",
        xaxis_type=xaxis_type,  # Set the x-axis type to 'log' or 'linear'
        yaxis_title="Intensity",
        legend_title="Wavelength",
        font=dict(
            family="Arial",
            size=12
        )
    )
    return fig