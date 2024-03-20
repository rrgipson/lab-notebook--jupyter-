# current, most recent code
import pandas as pd
import ipywidgets as widgets
from IPython.display import display, clear_output, HTML
import sys
sys.path.append("/Users/eleanordunietz/Desktop/Stanford/pythonScripts/data_analysis_scripts/stopped_flow")
from SF_plotting import *

param_map = {'substrate': 'Substrate', 'substrate_concentration': 'Substrate Concentration (mM)'}

# Global state for current index, plot arguments, and filtered data
current_index = 0
filtered_data = pd.DataFrame()  # Initialize as an empty DataFrame
# Initialize the arguments dictionary with default values for optional parameters
plot_args = {
    'substrate': None,
    'pH': None,
    'substrate_concentration': None,
    'solvent': None,
    'time_step': 20,
    'time_range': None,
    'subtract_baseline_flag': False,
    'wavelength_plotting_range': None
}
required_keys = ['substrate_concentration', 'solvent', 'pH', 'substrate'] 
# Create navigation buttons
prev_button = widgets.Button(description="Previous")
next_button = widgets.Button(description="Next")
navigation_info = widgets.Label()
# Output area for the plot
plot_output = widgets.Output()
toggle_baseline = widgets.Checkbox(
    value=False,
    description='Toggle Baseline',
    disabled=True)

def run_interactive_ui(processed_data):
    # Initialize UI components (e.g., create dropdowns, buttons)
    dropdowns = create_dropdowns_for_parameters(processed_data)
    
    # Setup event handlers
    for dropdown in dropdowns.values():
        dropdown.observe(on_dropdown_change, names='value')
    
    # Display UI components
    display(toggle_baseline)
    display(plot_output)
    display(widgets.HBox([prev_button, navigation_info, next_button]))
    update_navigation_buttons()
    update_plot()

def update_plot():
    global plot_args, current_index, filtered_data, required_keys

    with plot_output:
        clear_output(wait=True)     
        # Check if all required plot_args are specified before plotting
        if all(plot_args[key] is not None for key in required_keys):  
            temp_plot_args = plot_args.copy()
            temp_plot_args['index'] = current_index
            if not filtered_data.empty:  # Check if there's filtered data to plot
                fig = plot_wavelength_vs_intensity(filtered_data, **temp_plot_args)
                fig.show()
            else:
                print("No data to display. Please select specific values for all parameters.")
        else:
            print("Please select specific values for all parameters to display the plot.")

def update_navigation_buttons():
    global current_index, filtered_data
    total = len(filtered_data)
    prev_button.disabled = current_index <= 0
    next_button.disabled = current_index >= total - 1
    navigation_info.value = f"Spectrum {current_index + 1} of {total}" if total > 0 else "No spectra available"

def on_prev_clicked(b):
    global current_index
    if current_index > 0:
        current_index -= 1
        update_plot()
    update_navigation_buttons()

def on_next_clicked(b):
    global current_index, filtered_data
    if current_index < len(filtered_data) - 1:
        current_index += 1
        update_plot()
    update_navigation_buttons()

# Attach click event handlers to buttons
prev_button.on_click(on_prev_clicked)
next_button.on_click(on_next_clicked)


# Function to create and display dropdowns based on the diversity of options
def create_dropdowns_for_parameters(processed_data):
    parameters = ['substrate', 'pH', 'substrate_concentration', 'solvent']
    dropdowns = {}
    
    for parameter in parameters:
        unique_values = processed_data[parameter].dropna().unique()
        description = param_map.get(parameter, parameter) + ':'
        if len(unique_values) == 1:
            # If only one unique value, set it as the default value and don't disable the dropdown
            dropdowns[parameter] = widgets.Dropdown(
                options=[unique_values[0]],
                value=unique_values[0],
                description=description,
                disabled=True, 
            )
        else:
            # Provide a blank option for unselected state and list other unique values
            dropdowns[parameter] = widgets.Dropdown(
                options=[" "] + list(unique_values),
                description=description,
                disabled=False,
            )
        display(dropdowns[parameter])
    
    return dropdowns


def get_filtered_data(processed_data, selected_values):
    # Start with the full DataFrame
    filtered_df = processed_data
    
    # Iterate over the selected values and apply filters if the value is not blank or 'All'
    for param, value in selected_values.items():
        if value not in [' ', 'All']:
            filtered_df = filtered_df[filtered_df[param] == value]

    return filtered_df

# Updated to ensure plot is updated immediately when dropdown values change
def on_dropdown_change(change):
    global plot_args, current_index, filtered_data, dropdowns
    selected_values = {param: dropdowns[param].value for param in dropdowns if dropdowns[param].value not in [' ', 'All']}
    all_selected = all(value not in [' ', 'All'] for value in selected_values.values())
    
    # Update plot_args and filtered data only if all dropdowns have specific values
    if all_selected:
        plot_args.update(selected_values)
        filtered_data = get_filtered_data(processed_data, selected_values)
        current_index = 0
        update_plot()
    else:
        with plot_output:
            clear_output(wait=True)
            print("Please select specific values for all parameters to display the plot.")
    if all(plot_args[key] is not None for key in required_keys):
        update_navigation_buttons()




