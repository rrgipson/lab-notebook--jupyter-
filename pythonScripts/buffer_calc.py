def buffer_calculator(buffer_type, pH, total_volume_ml, conc_buffer, temp_C, NaOH_concentration_M=None):
    # Dictionary to include a variety of buffer systems
    buffers = {
        'CHES': {'Mw': 207.27, 'pKa': 9.5},
        'Borate': {'Mw': 61.83, 'pKa': 9.24},
        'Phosphate': {
            'acid': {'name': 'Potassium phosphate monobasic', 'Mw': 136.09, 'pKa': 7.2},
            'base': {'name': 'Sodium phosphate dibasic', 'Mw': 141.96}
        },
        'Tris': {
            'acid': {'name': 'Tris hydrochloride', 'Mw': 157.60, 'pKa': 8.1},
            'base': {'name': 'Tris', 'Mw': 121.14}
        }
    }

    # Check if the buffer type is valid
    if buffer_type not in buffers:
        return f"Unknown buffer type: {buffer_type}"

    if buffer_type in ['CHES', 'Borate']:
        # Handle single-component buffers that require pH adjustment with NaOH
        buffer_properties = buffers[buffer_type]
        Mw = buffer_properties['Mw']
        pKa = buffer_properties['pKa']
        adjusted_pKa = pKa - 0.01 * (temp_C - 25)  # Temperature correction
        ratio_base_acid = 10 ** (pH - adjusted_pKa)
        total_amount_buffer_moles = total_volume_ml * conc_buffer / 1000
        amount_solid_g = total_amount_buffer_moles * Mw

        if NaOH_concentration_M:
            proportion_base = ratio_base_acid / (1 + ratio_base_acid)
            amount_base_moles = total_amount_buffer_moles * proportion_base
            NaOH_volume_ml = (amount_base_moles / NaOH_concentration_M) * 1000
            adjustment_str = f"2. Adjust the pH to {pH} by adding {NaOH_volume_ml:.2f} mL of {NaOH_concentration_M} M NaOH.\n"
        else:
            adjustment_str = ""

        protocol = f"To prepare {total_volume_ml} mL of {conc_buffer} M {buffer_type} buffer at pH {pH}:\n1. Add {amount_solid_g:.2f} g of {buffer_type} to the solution.\n" + adjustment_str

    elif buffer_type in ['Phosphate', 'Tris']:
        # Handle dual-component buffers directly by calculating the mix of acid and base
        acid_properties = buffers[buffer_type]['acid']
        base_properties = buffers[buffer_type]['base']
        pKa = acid_properties['pKa']
        adjusted_pKa = pKa - 0.01 * (temp_C - 25)  # Temperature correction
        ratio_base_acid = 10 ** (pH - adjusted_pKa)
        total_amount_buffer_moles = total_volume_ml * conc_buffer / 1000

        moles_acid = total_amount_buffer_moles / (1 + ratio_base_acid)
        moles_base = total_amount_buffer_moles - moles_acid

        amount_acid_g = moles_acid * acid_properties['Mw']
        amount_base_g = moles_base * base_properties['Mw']

        protocol = (f"To prepare {total_volume_ml} mL of {conc_buffer} M {buffer_type} buffer at pH {pH}:\n"
                    f"1. Add {amount_acid_g:.2f} g of {acid_properties['name']} and {amount_base_g:.2f} g of {base_properties['name']} to the solution.\n"
                    "2. Mix until fully dissolved and adjust the total volume with deionized water as needed.\n")

    else:
        protocol = "Buffer type not supported."

    print(protocol)
    return protocol

