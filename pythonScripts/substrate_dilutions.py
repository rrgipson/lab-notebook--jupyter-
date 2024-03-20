from sys import displayhook
import pandas as pd
from IPython.display import HTML

def dilutions(substrate_name, stock_conc, target_conc, target_vol, og_stock_vol, mw):
    """
    Substrate Dilution Calculator 
    This function calculates the necessary volumes for creating dilutions from a stock solution. 
    It simplifies the process of preparing solutions at various concentrations.

    Features:
    - Customizable Inputs: Set parameters like stock concentration, molecular weight, and desired target concentrations.
    - Automated Volume Calculation: Determines volumes of stock and buffer for each dilution.
    - Result Presentation: Outputs in an easy-to-read table format.
    
    Parameters:
    - stock_conc (float): The concentration of the stock solution in mM.
    - target_conc (list of float): A list of target concentrations in mM.
    - target_vol (float): The desired total volume of each dilution in mL.
    - og_stock_vol (float): The original volume of the stock solution in mL.
    - mw (float): The molecular weight of the substrate in g/mol.

    Returns:
    None: Prints a DataFrame with dilution information and other relevant data.
    """
    # Initialize variables
    v_total = []
    stock_vol = []
    buffer_vol = []

    # Calculate volumes and concentrations for each dilution step
    for i in range(len(target_conc)):
        if i == len(target_conc) - 1:
            concentration = stock_conc
        else:
            concentration = target_conc[i + 1]
        if i == 0:
            v_total_i = target_vol
        else:
            v_total_i = target_vol + stock_vol[i - 1]
        v_total.append(round(v_total_i, 2))

        stock_vol_i = (target_conc[i] * v_total[i]) / concentration
        stock_vol.append(round(stock_vol_i, 2))

        buffer_vol_i = v_total[i] - stock_vol[i]
        buffer_vol.append(round(buffer_vol_i, 2))

    # Calculate the mass of the stock solution
    stock_mass = stock_conc * og_stock_vol * mw / 1000

    # Create a pandas DataFrame to store the dilution step information
    data = {'Concentration (mM)': target_conc,
            'Buffer volume (mL)': buffer_vol,
            'Stock volume (mL)': stock_vol,
            'Total volume (mL)': v_total}
    df = pd.DataFrame(data)
    
    # Print the styled DataFrame
    html = df.to_html(index=False, border = 1)
    displayhook(HTML(html))

    # Print the stock solution information
    print(f"\nMolecular Weight (g/mol): {mw}")
    print(f"Stock Concentration (mM): {stock_conc}")
    print(f"Stock Volume (mL): {og_stock_vol}")
    print(f"Mass (mg): {round(stock_mass, 2)}")

    # Call the existing print_protocol function to print the protocol instructions
    print_protocol(df, substrate_name, stock_mass, stock_conc, target_vol, og_stock_vol, mw)


def print_protocol(df, substrate_name, stock_mass, stock_conc, target_vol, og_stock_vol, mw):
    print("\nDilution Protocol Instructions:\n")
    print(f"Step 1:")
    print(f"  - Prepare stock solution ({stock_conc} mM)")
    print(f"  - Weigh {round(stock_mass, 2)} mg of {substrate_name}.")
    print(f"  - Dissolve in  {og_stock_vol} mL.")
    print(f"  - Filter with 0.2 micron syringe filter.\n")
    for index, row in df.iloc[::-1].iterrows():
        print(f"Step {len(df) - index + 1}:")
        print(f"  - Prepare {row['Total volume (mL)']} mL of {row['Concentration (mM)']} mM solution.")
        print(f"  - Add {row['Stock volume (mL)']} mL of previous stock solution.")
        print(f"  - Add {row['Buffer volume (mL)']} mL of buffer to reach the total volume.")
        print("  - Mix thoroughly.\n")
