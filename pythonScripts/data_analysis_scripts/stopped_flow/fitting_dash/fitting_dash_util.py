# utility functions


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
