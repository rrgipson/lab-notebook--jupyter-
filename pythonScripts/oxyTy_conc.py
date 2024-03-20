
"""
Biquinoline Assay and Dithionite Volume Calculation Module

This script performs calculations for a biquinoline assay to determine the Tyrosinase (Ty) and Copper (Cu) concentrations.
It assumes that the Ty concentration is roughly half of the total Cu concentration due to the binuclear Cu site.
After determining the Cu concentration, the script then calculates the amount of dithionite needed to reduce the solution.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def biquinoline_assay_calculation(sample_absorbance, sample_added_uL):
    """
    Performs the biquinoline assay calculation to determine the Ty and Cu concentrations.

    The Ty and Cu concentrations are obtained from the absorbance measurements using a previously established calibration curve.
    The Cu concentration is then used to determine the volume of dithionite stock solution required to reduce the enzyme solution.

    Parameters:
    sample_absorbance (float): Absorbance of the sample at 547 nm.
    sample_added_uL (float): Volume in microliters added to the biquinoline assay.

    Returns:
    float: Adjusted Cu concentration.
    """
    # Calculate dilution factor 
    dilution_factor = 700 / sample_added_uL

    # Data from the calibration curve
    data = {
        "Cu_standard_uL": [40, 80, 150, 260, 280],
        "Cu_mM": [0.009, 0.018, 0.03375, 0.0585, 0.063],
        "abs_546": [0.0628, 0.11897, 0.21237, 0.3654, 0.39065]
    }

    # Convert to DataFrame
    calibration_df = pd.DataFrame(data)

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(calibration_df['Cu_mM'], calibration_df['abs_546'])

    # Function to calculate Cu concentration based on the absorbance value
    def calculate_cu_concentration(absorbance, slope, intercept):
        cu_concentration = (absorbance - intercept) / slope
        return cu_concentration

    # Function to plot the sample point along with the calibration curve
    def plot_with_sample(absorbance, cu_concentration, calibration_df, slope, intercept):
        plt.figure(figsize=(8, 6))
        # Plot the calibration points
        plt.plot(calibration_df['Cu_mM'], calibration_df['abs_546'], 'o', label='Calibration Data')
        # Plot the fitted line
        cu_values = np.linspace(min(calibration_df['Cu_mM']), max(calibration_df['Cu_mM']), 100)
        abs_values = slope * cu_values + intercept
        plt.plot(cu_values, abs_values, 'r', label='Fitted Line')
        # Plot the sample point
        plt.plot(cu_concentration, absorbance, 'g*', markersize=15, label='Sample Point')
        # Labels and title
        plt.xlabel('Cu Concentration (mM)')
        plt.ylabel('Absorbance at 546 nm')
        plt.title('Calibration Curve with Sample Point')
        plt.legend()
        plt.grid(True)
        # Adjust x-axis scale if necessary
        plt.xlim(0, max(calibration_df['Cu_mM']) + 0.02)
        plt.ylim(0, max(calibration_df['abs_546']) + 0.05)
        # Show the plot
        plt.show()

    # Calculate the Cu concentration for the sample
    sample_cu_concentration = calculate_cu_concentration(sample_absorbance, slope, intercept)

    # Plot the calibration curve with the sample point
    plot_with_sample(sample_absorbance, sample_cu_concentration, calibration_df, slope, intercept)

    # Adjust for dilution
    sample_cu_concentration_adjusted = sample_cu_concentration * dilution_factor
    # Display the calculated Cu concentration
    print('Calculated Ty concentration (adjusted for dilution):', sample_cu_concentration_adjusted / 2) # divide by 2 for the 2 Cu centers in Ty

    return sample_cu_concentration_adjusted

def dithionite_volume_calculation(concentration_of_Cu_mM, volume_of_enzyme_solution_mL, volume_of_dithionite_stock_mL, mass_of_dithionite_mg):
    """
    Calculates the necessary volume of dithionite solution to add to the enzyme solution.

    The calculated Cu concentration from the assay is used to compute the necessary volume of dithionite solution.
    The dithionite is prepared in a stock solution, and the volume to add is calculated to ensure proper reduction.

    Parameters:
    concentration_of_Cu_mM (float): Concentration of Cu, obtained from previous calculations.
    volume_of_enzyme_solution_mL (float): Volume of enzyme solution in mL.
    volume_of_dithionite_stock_mL (float): Volume of dithionite stock in mL.
    mass_of_dithionite_mg (float): Dithionite mass in mg.

    Returns:
    float: Volume to add of the dithionite stock solution in microliters.
    """
    # Adjustable parameters 
    excess_factor = 5  # Excess factor for dithionite

    # Constants
    molecular_weight_dithionite_g_per_mol = 174.107  # Molecular weight of dithionite in g/mol
    purity_of_dithionite = 0.85  # Purity of dithionite

    # Calculation for enzyme solution
    n_umol_enzyme = concentration_of_Cu_mM * volume_of_enzyme_solution_mL

    # Calculations for dithionite
    n_umol_dithionite_needed = n_umol_enzyme * purity_of_dithionite * excess_factor
    concentration_dithionite_stock_mM = (mass_of_dithionite_mg * purity_of_dithionite * 1000) / (molecular_weight_dithionite_g_per_mol * volume_of_dithionite_stock_mL)

    # Calculate the volume to add of the dithionite stock solution to the enzyme solution
    volume_to_add_uL = (n_umol_dithionite_needed * 1000) / concentration_dithionite_stock_mM

    # Print out the result
    print(f'Volume to add of the dithionite stock solution: {volume_to_add_uL:.2f} uL')


   
def oxy_calc_conc(abs, dilution_uL, final_volume):
    """
    Calculates the concentration of oxyTy in a solution and provides dilution instructions.

    Parameters:
    abs (float): Absorbance of the solution at 345 nm.
    dilution_uL (float): Volume of the Ty solution in microliters used for dilution in UV/Vis measurement (out of 200).
    final_volume: Final desired volume of the enzyme solution
    
    The function will output the concentration of Ty in micromolar (uM) and instructions for diluting to a specific concentration.
    """
    dilution = 200 / dilution_uL
    conc = (abs * dilution) / 16 # 16 is epsilon at 345 nm
    conc_uM = conc * 1000
    print(f"Ty concentration: {conc_uM:.2f} uM")

    if conc_uM > 30:
        dilution2_mL = (30 * final_volume) / conc_uM
        dilution2_uL = dilution2_mL * 1000
        stock_vol = final_volume - dilution2_mL
        print(f"To dilute to {final_volume} mL of 30 uM Ty, add {dilution2_uL:.2f} uL Ty stock to {stock_vol:.2f} mL buffer")

