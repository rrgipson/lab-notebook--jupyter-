import sys
sys.path.append("/Users/eleanordunietz/Desktop/Stanford/pythonScripts/data_analysis_scripts/stopped_flow")
from SF_analysis_processing_ln import get_by_push, find_closest_wavelength, filter_by_time_cutoff
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import seaborn as sns


# models:
def reversible_irreversible_model(t, k1, k_neg_1, k2, A0, B_offset=10):
    """
    Model for reversible formation and irreversible decay of an intermediate B during a single turnover,
    with an additional offset for B to account for time lag or baseline shift.
    
    :param t: Time data
    :param k1: Rate constant for the forward reaction A to B
    :param k_neg_1: Rate constant for the reverse reaction B to A
    :param k2: Rate constant for the decay of B to C
    :param A0: Initial concentration of A
    :param B_offset: Offset for the concentration of B
    :return: Concentration of B at time t with an offset
    """
    
    def rate_equations(y, t, k1, k_neg_1, k2):
        A, B, _ = y
        dA_dt = -k1 * A + k_neg_1 * B
        dB_dt = k1 * A - (k_neg_1 + k2) * B
        dC_dt = k2 * B
        return [dA_dt, dB_dt, dC_dt]

    y0 = [A0, 0, 0]  # Initial concentrations of A, B, and C
    solution = odeint(rate_equations, y0, t, args=(k1, k_neg_1, k2))
    return solution[:, 1] + B_offset  # Add the B_offset to the concentration of B

def reversible_first_order_model(t, k1, k_neg_1, A0, B_offset=0):
    """
    Model for a reversible first-order reaction A <-> B.
    
    :param t: Time data
    :param k1: Rate constant for the forward reaction A to B
    :param k_neg_1: Rate constant for the reverse reaction B to A
    :param A0: Initial concentration of A
    :param B_offset: Offset for the concentration of B, accounting for the baseline measurement
    :return: Concentration of B at time t
    """
    # Assuming that B is not being consumed in a secondary reaction, or that it's negligible
    B_eq = A0 * k1 / (k1 + k_neg_1)  # Equilibrium concentration of B
    # If k1 is not equal to k_neg_1, use the full expression
    if not np.isclose(k1, k_neg_1):
        B = B_eq + (A0 - B_eq) * np.exp(-(k1 + k_neg_1) * t)
    else:  # If k1 is very close to k_neg_1, avoid division by zero
        B = A0 * k1 * t * np.exp(-k1 * t)
    return B + B_offset

def consecutive_first_order_model(t, k1, k2, A0, B_offset=0):
    """
    Model for a consecutive first-order reaction where A converts to B, and then B converts to C.
    
    :param t: Time data
    :param k1: Rate constant for the conversion of A to B
    :param k2: Rate constant for the conversion of B to C
    :param A0: Initial concentration of A
    :param B_offset: Offset for the concentration of B, accounting for the baseline measurement
    :return: Concentration of B at time t
    """
    # Avoid division by zero in case k1 and k2 are very close
    if np.isclose(k1, k2):
        B = A0 * k1 * t * np.exp(-k1 * t)
    else:
        B = (k1 * A0) / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
    return B + B_offset

def first_order_formation_decay(time, k_formation, k_decay, reactant_initial=100, product_initial=0, product_offset=0):
    """
    First-order kinetic model for product formation and decay.
    :param time: Time in seconds (s).
    :param k_formation: Rate constant for product formation in per second (s^-1).
    :param k_decay: Rate constant for product decay in per second (s^-1).
    :param reactant_initial: Initial concentration of the reactant in micromolar (uM).
    :param product_initial: Initial concentration of the product in micromolar (uM) before decay.
    :param product_offset: Initial concentration of the product in micromolar (uM) before formation.
    :return: Concentration of product in micromolar (uM).
    """
    product_formed = reactant_initial * (1 - np.exp(-k_formation * time)) + product_offset
    product_decayed = product_initial * np.exp(-k_decay * time)
    return product_formed - product_decayed    

def first_order_product_formation_with_offset(time, k, reactant_initial=100, product_offset=0):
    """
    First-order kinetic model for product formation with an initial product offset.
    :param time: Time in seconds (s).
    :param k: Rate constant in per second (s^-1).
    :param reactant_initial: Initial concentration of the reactant in micromolar (uM).
    :param product_offset: Initial concentration of the product in micromolar (uM).
    :return: Concentration of product formed in micromolar (uM).
    """
    return reactant_initial * (1 - np.exp(-k * time)) + product_offset

def consecutive_first_order_corrected_model(t, k1, k2, I_A0, B_offset=0):
    """
    Model for a consecutive first-order reaction where A converts to B, and then B converts to C,
    corrected for the initial intensity contribution of A.

    :param t: Time data
    :param k1: Rate constant for the conversion of A to B
    :param k2: Rate constant for the conversion of B to C
    :param A0: Initial concentration of A
    :param I_A0: Initial intensity contribution of A
    :param B_offset: Offset for the concentration of B, accounting for the baseline measurement
    :return: Corrected intensity of B at time t
    """
    # Avoid division by zero in case k1 and k2 are very close
    # if np.isclose(k1, k2):
    #     B = A0 * k1 * t * np.exp(-k1 * t)

    B = (k1 * I_A0) / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
    
    # Subtract the initial intensity contribution of A from the calculated intensity of B
    B_intensity_corrected = B - I_A0 * np.exp(-k1 * t)
    
    return B_intensity_corrected + B_offset

# fitting functions: 
def fit_consecutive_first_order(kinetic_model, initial_guesses, time_data, intensity_data, **kwargs):
    expected_params = 4  # Number of expected parameters: k1, k2, A0, B_offset
    if initial_guesses is None:
        initial_guesses = [1, 0.1, 100, 0]  # Default initial guesses: k1, k2, A0, B_offset
    elif len(initial_guesses) != expected_params:
        print(f"Initial guesses should have {expected_params} elements: k1, k2, A0, B_offset")
        return None

    # Perform the curve fitting
    try:
        params, _ = curve_fit(lambda t, k1, k2, A0, B_offset: kinetic_model(t, k1, k2, A0, B_offset, **kwargs), time_data, intensity_data, p0=initial_guesses, bounds=([0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]))
        k1_fit, k2_fit, A0_fit, B_offset_fit = params
        print(f"Fitted rate constant for A to B (k1): {k1_fit} s^-1")
        print(f"Fitted rate constant for B to C (k2): {k2_fit} s^-1")
        print(f"Fitted initial concentration of A (A0): {A0_fit} (relative intensity)")
        print(f"Fitted baseline offset for B (B_offset): {B_offset_fit} (relative intensity)")

        # Calculate residuals and R-squared
        residuals = intensity_data - kinetic_model(time_data, *params)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((intensity_data - np.mean(intensity_data))**2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"R-squared (Goodness of Fit): {r_squared}")
        return params
    
    except RuntimeError as e:
        print(f"Curve fitting failed: {e}")
        return None
    
def fit_consecutive_first_order_corrected(kinetic_model, initial_guesses, time_data, intensity_data, bounds=None, **kwargs):
    expected_params = 4  # Number of expected parameters: k1, k2, A0, I_A0, B_offset
    if initial_guesses is None:
        # Default initial guesses: k1, k2, I_A0, B_offset
        initial_guesses = [10, 0.5, 0.030, 0]  # Assuming initial intensity contribution of A is 10
    elif len(initial_guesses) != expected_params:
        print(f"Initial guesses should have {expected_params} elements: k1, k2, I_A0, B_offset")
        return None
    
    # set default bounds if none are provided:
    if bounds is None: 
        bounds =([0,0,0,-np.inf], [np.inf, np.inf, np.inf, np.inf])

    # Perform the curve fitting
    try:
        params, _ = curve_fit(lambda t, k1, k2, I_A0, B_offset: kinetic_model(t, k1, k2, I_A0, B_offset, **kwargs),
                              time_data, intensity_data, p0=initial_guesses,
                              bounds=bounds)
        k1_fit, k2_fit, I_A0_fit, B_offset_fit = params
        # print(f"Fitted rate constant for A to B (k1): {k1_fit} s^-1")
        # print(f"Fitted rate constant for B to C (k2): {k2_fit} s^-1")
        # print(f"Fitted initial intensity of A (I_A0): {I_A0_fit} (relative intensity)")
        # print(f"Fitted baseline offset for B (B_offset): {B_offset_fit} (relative intensity)")

        # Calculate residuals and R-squared
        residuals = intensity_data - kinetic_model(time_data, *params)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((intensity_data - np.mean(intensity_data))**2)
        r_squared = 1 - (ss_res / ss_tot)
        # print(f"R-squared (Goodness of Fit): {r_squared}")
        return params, r_squared
    
    except RuntimeError as e:
        print(f"Curve fitting failed: {e}")
        return None

def fit_first_order_formation(kinetic_model, initial_guesses, time_data, intensity_data, **kwargs):
    expected_params = 3
    # Set default initial guesses if none are provided
    if initial_guesses is None:
        initial_guesses = [2, 50, 5]  # Default initial guesses including product offset 
    elif len(initial_guesses) != expected_params:
        print(f"Initial guesses should have {expected_params} elements: k (/s), reactant_initial (uM), product_offset (uM)")
        return None

    # Perform the curve fitting
    try:
    
        params, _ = curve_fit(lambda t, *p: kinetic_model(t, *p, **kwargs), time_data, intensity_data, p0=initial_guesses, bounds=([0, 0, 0], [150, 100, 100]))
        k_fit = params[0]  # Assuming the first parameter is always the rate constant
        print(f"Fitted rate constant (k): {k_fit} s^-1")  # Rate constant in per second

        if len(params) > 1:
            reactant_initial_fit = params[1]
            print(f"Fitted initial reactant concentration: {reactant_initial_fit} (relative intensity)")  # Concentration in micromolar
        # Calculate residuals and R-squared
        residuals = intensity_data - kinetic_model(time_data, *params)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((intensity_data - np.mean(intensity_data))**2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"R-squared (Goodness of Fit): {r_squared}")
        return params
    
    except RuntimeError as e:
        print(f"Curve fitting failed: {e}")
        return None
    
def fit_first_order_formation_decay(kinetic_model, initial_guesses, time_data, intensity_data, **kwargs):
    expected_params = 5  # Number of expected parameters for first_order_formation_decay
    if initial_guesses is None:
        initial_guesses = [1, 0.1, 50, 0, 10]  # Default initial guesses
    elif len(initial_guesses) != expected_params:
        print(f"Initial guesses should have {expected_params} elements: k_formation, k_decay, reactant_initial, product_initial, product_offset")
        return None
    
    # Perform the curve fitting
    try:
        params, _ = curve_fit(lambda t, *p: kinetic_model(t, *p, **kwargs), time_data, intensity_data, p0=initial_guesses, bounds=([0, 0.001, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf, np.inf]))
        k_formation_fit, k_decay_fit, reactant_initial_fit, product_initial_fit, product_offset_fit = params

        print(f"Fitted rate constant (k_formation): {k_formation_fit} s^-1")
        print(f"Fitted rate constant (k_decay): {k_decay_fit} s^-1")
        print(f"Fitted initial reactant concentration: {reactant_initial_fit} (relative intensity)")
        print(f"Fitted initial product concentration before decay: {product_initial_fit} (relative intensity)")

        # Calculate residuals and R-squared
        residuals = intensity_data - kinetic_model(time_data, *params)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((intensity_data - np.mean(intensity_data))**2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"R-squared (Goodness of Fit): {r_squared}")
        return params
    
    except RuntimeError as e:
        print(f"Curve fitting failed: {e}")
        return None

def fit_reversible_irreversible(kinetic_model, initial_guesses, time_data, intensity_data, **kwargs):
    expected_params = 5  # Now we have five parameters, including the offset
    if initial_guesses is None:
        initial_guesses = [0.01, 0.005, 0.001, 100, 0]  # Add a default guess for B_offset
    elif len(initial_guesses) != expected_params:
        print(f"Initial guesses should have {expected_params} elements: k1, k_neg_1, k2, A0, B_offset")
        return None
    
    # Perform the curve fitting
    try:
        # Set bounds with the additional parameter for B_offset
        bounds = ([0, 0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf])
        params, _ = curve_fit(lambda t, *p: kinetic_model(t, *p, **kwargs), time_data, intensity_data, p0=initial_guesses, bounds=bounds)
        k1_fit, k_neg_1_fit, k2_fit, A0_fit, B_offset_fit = params

        print(f"Fitted rate constant for A to B (k1): {k1_fit} s^-1")
        print(f"Fitted rate constant for B to A (k_neg_1): {k_neg_1_fit} s^-1")
        print(f"Fitted rate constant for B to C (k2): {k2_fit} s^-1")
        print(f"Fitted initial concentration of A (A0): {A0_fit} uM")
        print(f"Fitted B offset: {B_offset_fit} uM")

        # Calculate residuals and R-squared
        residuals = intensity_data - kinetic_model(time_data, *params, **kwargs)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((intensity_data - np.mean(intensity_data))**2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"R-squared (Goodness of Fit): {r_squared}")
        return params
    
    except RuntimeError as e:
        print(f"Curve fitting failed: {e}")
        return None


# main fitting function:
def fit_kinetics_to_experiment(key, push, wavelength, product_epsilon=None, time_cutoff=None, initial_guesses=None, bounds=None, plot_fit=True, kinetic_model='first_order_product_formation_with_offset', *args, **kwargs):
    # Map string names to actual function objects
    kinetic_model_functions = {
        'first_order_product_formation_with_offset': [first_order_product_formation_with_offset, fit_first_order_formation],
        'first_order_formation_decay': [first_order_formation_decay, fit_first_order_formation_decay],
        'consecutive_first_order_model': [consecutive_first_order_model, fit_consecutive_first_order],
        'reversible_irreversible_model': [reversible_irreversible_model, fit_reversible_irreversible],
        'consecutive_first_order_corrected_model': [consecutive_first_order_corrected_model, fit_consecutive_first_order_corrected]
    }

    # Get the actual function object based on the string
    kinetic_model_function = kinetic_model_functions.get(kinetic_model)[0]
    fit_function = kinetic_model_functions.get(kinetic_model)[1]

    if kinetic_model_function is None:
        print(f"Kinetic model '{kinetic_model}' is not recognized.")
        return None
    
    # Find the experiment
    experiment = get_by_push(key, push)
    if experiment is None:
        print(f"No experiment found with push {push}.")
        return None
    
    data = experiment['data']

    # Apply time cutoff if specified
    if time_cutoff is not None:
        data = filter_by_time_cutoff(data, time_cutoff)
        if data is None:
            return  # Stop the function if the filtered data is None
    
    # Extract time and intensity data after filtering
    time_data = data['Time']
    closest_wavelength = find_closest_wavelength(data, wavelength)
    intensity_data = data[closest_wavelength]
   

    if kinetic_model in kinetic_model_functions:
        if kinetic_model == 'consecutive_first_order_corrected_model':
            params, r_squared = fit_function(kinetic_model_function, initial_guesses, time_data, intensity_data,  bounds=bounds)
        else:params = fit_function(kinetic_model_function, initial_guesses, time_data, intensity_data)

    if params is None:
        print("Curve fitting failed")
        return None
    
    # Initialize an empty figure object that might be filled later
    fig = None
    
    # Plot the fitted curve if plot_fit is True
    if plot_fit:
        fig = go.Figure()
        # Add experimental data as a scatter plot
        fig.add_trace(go.Scatter(x=time_data, y=intensity_data, mode='lines', name=f"{closest_wavelength} nm (Data)"))
        
        time_range = np.linspace(time_data.min(), time_data.max(), 500)
        fit_curve = kinetic_model_function(time_range, *params)
        # Add fitted curve as a line plot
        fig.add_trace(go.Scatter(x=time_range, y=fit_curve, mode='lines', name='Fitted Curve', line=dict(dash='dash')))

        if kinetic_model == 'consecutive_first_order_corrected_model':
            # Add annotations for fitted parameters
            fig.add_annotation(xref="paper", yref="paper", 
                   x=1.02, y=0.3,  # Positions annotation at the bottom center outside the plotting area
                   xanchor='left', yanchor='top',
                   text=f'Fitted parameters: <br> k<sub>1</sub>: {params[0]:.4f} s<sup>-1</sup> <br> k<sub>2</sub>: {params[1]:.4f} s<sup>-1</sup> <br> Intensity A<sub>0</sub>: {params[2]:.4f} <br> B offset: {params[3]:.4f} <br> R<sup>2</sup>: {r_squared:.4f}',
                   showarrow=False, font=dict(size=12), align='left')

        # Update the layout
        fig.update_layout(title=f"Fitted Kinetics for Push {push} at {closest_wavelength} nm",
                          xaxis_title="Time (s)",
                          yaxis_title="Intensity",
                          template='plotly_dark',)


    return params, fig

