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

def conc_from_uvvis(A, eps, v_add=2, v_tot=500, pathlength=1):
    '''Calculates concentrations from UV-vis Abs numbers.
    eps: takes numerical value or parses 'tph' Abs(280nm), 'besc' Abs(280nm), 'nonoate' Abs(250nm).
    Assumes diluted from stock solution by adding v_add of stock into cuvette to give v_tot uL of solution.'''
    
    #params
    eps_dict = {'tph' : 39310,
                'besc' : 34400,
                'nonoate' : 7250}

    if eps in eps_dict.keys():
        epsilon=eps_dict[eps]
    else:
        try:
            epsilon = float(eps)
        except:
            print('Please input eps as compound/protein name or numberical value in [1/cm*M]')
            return False

    #get concentration of whats in cuvette
    c1=A/(epsilon*pathlength) #M
    #concentration of the 2uL (the actual sample)
    #M1V1=M2V2
    v1=v_tot #uL
    v2=v_add #uL
    m2=(c1*v1)/v2
    m2_mM=m2*1000
    print(m2_mM,' mM')
    return m2_mM
