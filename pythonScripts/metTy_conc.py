def met_calc_conc(A350_H2O2, A350_O2, dilution_uL_H2O2, dilution_uL_O2):
    """
    Calculate concentrations of deoxy, met, and total Tyrosinase (Ty).
    
    Parameters:
    A350_H2O2 (float): Absorbance of Ty in 1 mM H2O2 at 350 nm.
    A350_O2 (float): Absorbance of Ty in O2 at 350 nm.
    dilution_uL_H2O2 (int): Number of uL of Ty solution added to H2O2 solution (out of 200 total uL, range of 5-20).
    dilution_uL_O2 (int): Number of uL of Ty solution added to O2 solution (out of 200 total uL, range of 5-20).

    Returns:
    dict: Concentrations in mM and uM for total, met, and deoxy Tyrosinase.
    """
    # Constants
    epsilon_350 = 16  # Molar absorptivity
    dilution_factor_H2O2 = 200 / dilution_uL_H2O2
    dilution_factor_O2 = 200 / dilution_uL_O2
    
    # Calculate [Ty]total
    Ty_total = (A350_H2O2 * dilution_factor_H2O2) / epsilon_350
    
    # Estimate [metTy]
    metTy_estimated = 0.85 * Ty_total
    
    # Quantify [deoxy]
    deoxy = (A350_O2 * dilution_factor_O2) / epsilon_350
    
    # Final concentrations
    metTy_calculated = Ty_total - deoxy
    
    # Results dictionary
    results = {
        '[Ty] total (mM)': Ty_total,
        'Estimated [metTy](uM)': metTy_estimated * 1000,
        '[deoxyTy] (uM)': deoxy * 1000,
        'Calculated [metTy] (uM)': metTy_calculated * 1000
    }
   
    # Print results
    print("Concentration Results:")
    for key, value in results.items():
        print(f"{key}: {value:.2f}")
    return results

def dilute_to_target(metTy_calculated_uM, target_concentration_uM, final_volume_ml):
    """
    Calculates the volumes required for diluting to a target concentration.

    Parameters:
    metTy_calculated_uM (float): The calculated concentration of metTy (in uM) from 'met_calc_conc'.
    target_concentration_uM (float): The target concentration in uM.
    final_volume_ml (float): The final volume of the solution in mL.

    Returns:
    dict: Volumes in uL and mL required for the dilution.
    """
    # Calculate the volume needed for the dilution
    volume_to_dilute_ml = (target_concentration_uM * final_volume_ml) / metTy_calculated_uM
    # Calculate the volume of solvent needed for the dilution
    solvent_volume_ml = final_volume_ml - volume_to_dilute_ml
    
    # Results dictionary
    dilution_results = {
        'Volume of enzyme solution to dilute (uL)': volume_to_dilute_ml * 1000,
        'Solvent Volume (mL)': solvent_volume_ml
    }

    # Print short protocol
    print(f"Dilution Protocol:")
    print(f"1. Add {dilution_results['Volume of enzyme solution to dilute (uL)']:.2f} uL of enzyme solution.")
    print(f"2. Add {dilution_results['Solvent Volume (mL)']:.2f} mL of solvent.")
    print(f"3. Mix to reach a final volume of {final_volume_ml} mL.")

    return dilution_results
