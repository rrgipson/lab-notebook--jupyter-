def concentration_calculator(conc_str, Mw, final_volume_str, chemical_name):
    # Parse concentration string to extract value and units
    if 'mM' in conc_str:
        conc_units = 'mM'
        conc_value = float(conc_str.replace('mM', ''))
        conc_M = conc_value / 1000  # Convert mM to M
    elif 'uM' in conc_str:
        conc_units = 'uM'
        conc_value = float(conc_str.replace('uM', ''))
        conc_M = conc_value / 1000000  # Convert µM to M
    elif 'M' in conc_str:
        conc_units = 'M'
        conc_value = float(conc_str.replace('M', ''))
        conc_M = conc_value
    else:
        return "Invalid concentration string. Please include units mM, uM, or M."

    # Parse volume string to extract value and convert to liters if necessary
    if 'mL' in final_volume_str:
        final_volume_L = float(final_volume_str.replace('mL', '')) / 1000  # Convert mL to L
    elif 'L' in final_volume_str:
        final_volume_L = float(final_volume_str.replace('L', ''))
    else:
        return "Invalid volume string. Please include units mL or L."

    # Calculate the amount of compound needed in moles
    moles_needed = conc_M * final_volume_L
    
    # Convert moles to grams
    grams_needed = moles_needed * Mw
    
    # Determine if we should use grams or milligrams
    if grams_needed < 1:
        mass_needed_str = f"{grams_needed * 1000:.2f} mg"
    else:
        mass_needed_str = f"{grams_needed:.2f} g"
    
    # Generate protocol
    protocol = (f"To prepare {final_volume_str} of {conc_str} {chemical_name} solution:\n"
                f"1. Weigh out {mass_needed_str} of {chemical_name}.\n"
                "2. Dissolve the compound in a volume of deionized water slightly less than the final volume.\n"
                "3. Once fully dissolved, adjust the solution to the final desired volume with deionized water.")
    
    return protocol