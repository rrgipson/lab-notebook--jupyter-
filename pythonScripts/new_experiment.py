import pandas as pd
from pathlib import Path
import datetime
import nbformat
from nbformat.v4 import new_notebook

def contains_space(s):
    return " " in s

def create_project_map(project_names):
    """Convert project names to a dictionary mapping original names to folder names."""
    # Remove any empty strings and strip leading/trailing whitespace
    project_names = [name.strip() for name in project_names if name]
    project_map =  {name: name.replace(' ', '_').lower() for name in project_names}
    return project_map

# Ensure the DataFrame CSV file exists
df_path = Path('experiments.csv')
if not df_path.exists():
    df = pd.DataFrame(columns=['Experiment ID', 'Description', 'Date', 'Type', 'Project', 'Notebook Link'])
    df.to_csv(df_path, index=False)
else:
    df = pd.read_csv(df_path)

def create_experiment(exp_id, exp_type, project, description, project_options):
    # Use a relative path from the experimentHub.ipynb location
    base_directory = Path('experiments')
    global df
    # Define a mapping for project names to folder names
    project_mapping = create_project_map(project_options)
    # Check if exp_id is already in the "Experiment ID" column
    if not df['Experiment ID'].str.contains(exp_id).any():
        if contains_space(exp_id):
            return
        
        # Map the project to its folder name
        folder_name = project_mapping.get(project, 'General')
        project_path = base_directory / folder_name

        # Create the project folder if it doesn't exist
        project_path.mkdir(parents=True, exist_ok=True)
        # Create the base experiment folder and subfolders within the project folder
        base_path = project_path / exp_id
        base_path.mkdir(parents=True, exist_ok=True)
        (base_path / 'rawData').mkdir(exist_ok=True)
        (base_path / 'templates').mkdir(exist_ok=True)
        
        # Adjust the relative path for the backlink to experimentHub.ipynb
        relative_path_to_hub = Path('..') / '..' / '..' / 'experimentHub.ipynb'

        # Create a new Jupyter notebook
        nb_path = base_path / f"{exp_id}.ipynb"
        nb = new_notebook()
        # Update the link to use the relative path
        nb.cells.append(nbformat.v4.new_markdown_cell(f"[Back to Experiment Hub]({relative_path_to_hub})\n# {description}\n**Hypothesis:**"))
        # Add other sections as before
        nb.cells.append(nbformat.v4.new_markdown_cell("## Materials\n\nMaterials | Cas # | Location | Molecular weight (g/mol)\n---------|----------|---------|---"))  # Materials section
        nb.cells.append(nbformat.v4.new_markdown_cell("## Method\nDescribe the method here."))  # Method section
        nb.cells.append(nbformat.v4.new_markdown_cell("## Analysis\nWrite your analysis plan here."))  # Analysis section
        nb.cells.append(nbformat.v4.new_markdown_cell("## Results\nDocument your results here."))  # Results section
        nb.cells.append(nbformat.v4.new_markdown_cell("## Conclusion\nSummarize the conclusions here."))  # Conclusion section

        with open(nb_path, 'w', encoding='utf-8') as f:
            nbformat.write(nb, f)

        
        # Prepare a new row to add to the DataFrame
        today = datetime.datetime.now().strftime("%y-%m-%d")
        notebook_link = f"[{exp_id}]({str(nb_path)})"
        new_row = pd.DataFrame([{'Experiment ID': exp_id, 'Description': description, 'Date': today, 'Type': exp_type, 'Project': project, 'Notebook Link': notebook_link}])
        
        # Use pd.concat to append the new row

        df = pd.concat([df, new_row], ignore_index=True)
        
        # Save the updated DataFrame
        df.to_csv(df_path, index=False)
    return df