import os
import pandas as pd

def find_baseline_for_push(key, push_number):
    """
    Find the baseline experiment for a given push number.
    
    Parameters:
    key (pd.DataFrame): The DataFrame containing the experiments.
    push_number (str): The push identifier to match.
    
    Returns:
    pd.Series: A pandas Series representing the averaged baseline data, or None if not found.
    """
    # Find the experiment with the given push number
    target_experiment = get_by_push(key, push_number)

    if target_experiment is None:
        # If the experiment with the given push number is not found
        print(f"No experiment found with push number {push_number}.")
        return None

    # Extracting the required properties from the target experiment
    target_properties = {
        'solvent': target_experiment['solvent'],
        'date': target_experiment['date'],
        'Ty': target_experiment['Ty'],
        'pH': target_experiment['pH']
    }

    # Find the baseline experiment
    baseline = key[
        (key['solvent'] == target_properties['solvent']) &
        (key['date'] == target_properties['date']) &
        (key['Ty'] == target_properties['Ty']) &
        (key['pH'] == target_properties['pH']) &
        (key['substrate'] == '-') &
        (key['substrate_concentration'] == '-')
    ]

    if baseline.empty:
        # If no baseline experiment is found
        print("No baseline found matching the criteria.")
        return None
    else:
        # print(f"Baseline: {baseline['push'].values[0]}")
        # Assuming the data is a DataFrame within the baseline dict under the key 'data'
        baseline_data = baseline['data'].values[0]
        # Averaging each column (each wavelength) across all rows (time points), excluding 'Time'
        averaged_baseline = baseline_data.drop(columns=['Time']).mean()
        return averaged_baseline

def subtract_baseline(experiment_data, baseline):    
    # Check if baseline is a Series and not None
    if not isinstance(baseline, pd.Series):
        print("Baseline data is not in the expected format.")
        return None

    # Subtract the baseline values from the experiment data
    adjusted_data = experiment_data.copy()
    for column in adjusted_data.columns:
        if column != 'Time':  # Skip the 'Time' column
            adjusted_data[column] = adjusted_data[column] - baseline.get(column, 0)
    return adjusted_data[:-1]

def detect_decimal_precision(data):
    # Sample a few wavelength headers to determine the decimal precision
    sample_headers = data.columns[1:10]  # Adjust the range as needed
    # Function to count decimal places
    def count_decimals(string):
        if '.' in string:
            return len(string.split('.')[1])
        return 0
    # Count the number of decimal places in the sample
    decimal_counts = [count_decimals(header) for header in sample_headers]
    # Find the most common decimal count
    most_common_precision = max(set(decimal_counts), key=decimal_counts.count)
    return most_common_precision

def find_closest_wavelength(data, desired_wavelength):
    # Detect the common decimal precision in the data
    precision = detect_decimal_precision(data)
    # Convert wavelength headers to floats
    wavelength_headers = [float(w) for w in data.columns[1:]]
    # Find closest wavelength
    closest_wavelength = min(wavelength_headers, key=lambda x: abs(x - desired_wavelength))

    return str(closest_wavelength)

def process_csv_file(file_path):
    """
    Process a single CSV file. Skip initial rows, transpose if necessary, and reshape.
    Adjusts the cutoff line based on the file name 
    """
    # List of files with a different cutoff line
    special_files = ['Pda00357.csv', 'Pda00327.csv', 'Pda00336.csv', 
                     'Pda00337.csv', 'Pda00338.csv', 'Pda00339.csv', 
                     'Pda00344.csv', 'Pda00348.csv', 'Pda00349.csv', 
                     'Pda00352.csv', 'Pda00355.csv', 'Pda00356.csv']

    data_lines = []
    # Extract the base name of the file from the path
    base_name = os.path.basename(file_path)
    cutoffLine = 27 if base_name in special_files else 26

    transposed = False
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i == cutoffLine -1:
                # Check if the line before cutoff indicates transposed data
                transposed = 'Wavelength' in line.split(',')[0]
                continue
            if i < cutoffLine:
                continue;
            if ('Count' in line):
                break;
            if i == cutoffLine and not transposed:
                line = 'Time' + line
            data_lines.append(line)

    # Extract headers from the first line
    headers = data_lines[0].split(',')[:-1]
    df = pd.DataFrame([sub.split(',') for sub in data_lines[1:]], columns=headers)
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    df = df.apply(pd.to_numeric, errors='coerce')

    # Transpose DataFrame if transposed
    if transposed:
        # Transpose the DataFrame
        df_transposed = df.transpose()   
        # Now the first row contains the new headers, set them and then drop the first row
        df_transposed.columns = df_transposed.iloc[0]
        df_transposed = df_transposed.drop(df_transposed.index[0])       
        # Reset the index to make sure 'Time' is a column and not the index
        df_transposed.reset_index(inplace=True)       
        # Rename 'index' to 'Time' to match your desired output
        df_transposed.rename(columns={'index': 'Time'}, inplace=True)
        df = df_transposed
        df = df.iloc[:,:-1]
    
    if not transposed: 
        df = df.iloc[:-1]

    # Standardize headers to strings
    df.columns = df.columns.astype(str)
    return df

def load_key_from_csv(key_csv_file_path):
    """
    Load the key data from a CSV file.
    """
    # return pd.read_csv(key_csv_file_path).to_dict(orient='records')
    return pd.read_csv(key_csv_file_path)

def get_experiments_by_criteria(key, **criteria):
    filtered_data = key
    for key, value in criteria.items():
        if value is not None:
            filtered_data = filtered_data[filtered_data[key] == value]

    return filtered_data

def get_by_push(key, push):
    """
    Get a single experiment from the key file that matches the specified push.
    """
    # Filter the DataFrame for the row with the matching 'push' value
    filtered_df = key[key['push'] == push]
    # Check if there is at least one match
    if not filtered_df.empty:
        # Return the first match as a Series
        return filtered_df.iloc[0]
    else:
        # Return None if no match is found
        return None

def process_all_csv_files(directory_path, key_file_path):
    """
    Process all CSV files in the specified directory and merge with key file.
    """
    
    # Read and merge with the key file
    key = load_key_from_csv(key_file_path)

    # Process each CSV file in the directory
    for file_name in os.listdir(directory_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(directory_path, file_name)
            # process the csv
            processed_data = process_csv_file(file_path)
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            print("Processed file:", file_name.replace('.csv', ''))
            # add the data to the experiment within the key array
            experiment = get_by_push(key, file_name.replace('.csv', ''))
            if experiment is not None:
                # Check if 'data' column exists, if not, create it
                if 'data' not in key.columns:
                    key['data'] = None

                # Assign the processed data to the 'data' column of the experiment
                key.at[experiment.name, 'data'] = processed_data

    df = pd.DataFrame(key)

    return df

def filter_by_time_cutoff(data, time_cutoff, start_time=None):
    """
    Filters the dataset to include only the data points where the time value 
    is less than or equal to the time_cutoff.
    """
    # Debugging: ensure the 'Time' column is of type float for comparison
    data['Time'] = data['Time'].astype(float)

    if start_time is not None:
        data = data[data['Time'] >= start_time]
        if time_cutoff is None:
            return data

    max_time = data['Time'].max()
    if time_cutoff > max_time:
        print(f"Time cutoff ({time_cutoff}) is beyond the maximum time ({max_time} s) value for the dataset.")
        return None
    # Filter the data
    filtered_data = data[data['Time'] <= time_cutoff]
    return filtered_data