# Digital Lab Notebook
This repository serves as my lab notebook, providing a clean design to manage my experiments efficiently. It combines the organizational framework of Notion with the ability to run executable code from Jupyter notebooks, optimized for VS Code's Flexoki theme.

## Structure
The notebook is organized into three main directories:

1. **experimentHub**: 
    
    Features the experimentHub.ipynb dashboard, a central hub inspired by Notion's table functionality, linking to individual experiment notebooks. It cleanly displays experiments.csv, which gets updated when new experiments are added either through the CSV or via an easy-to-use "Add Experiment" button. This button facilitates the creation of new experiments from an IPython notebook template. Simplified ipywidgets on the dashboard collect essential details like experiment ID, description, project, and type for easy organization and sorting. The dashboard also allows for filtering by project, facilitating easy access to relevant experiments.

2.  **experimentTemplates**: 
   
    Stores reusable templates for experimental protocols, streamlining the process of starting new experiments. (I also recently made a custom VS Code extension to make it even easier to copy templates for experimental prep directly into experiment notebooks while maintaining a readable format.)
   
3. **pythonScripts**: 
   
    Contains Python scripts that support the functionality of the notebook, including data analysis and automation scripts.

### Key Features
- **Dynamic Dashboard:** The experimentHub.ipynb dashboard dynamically updates to reflect new experiments, offering a clean display and simple navigation in VS Code.
- **Experiment Creation:** Easily create new experiments with predefined templates, automatically generating directories with necessary folders (templates and rawData) and a starter IPython notebook.
- **Optimized Workflow**: Designed with simplicity in mind, it's tailored to reduce visual clutter while maintaining detailed records of experiments. Integration with a custom VS Code extension for copying experimental templates enhances efficiency.


This project is a reflection of my ongoing efforts to optimize my experimental workflow. By leveraging automation and standardizing processes, it aims to thoroughly document each research step. This  approach not only enhances scientific communication but also ensures that results and conclusions are reproducible. 

As the Digital Lab Notebook evolves, future updates will include additional tools and templates to further support research needs, particularly in data analysis. This project is driven by the belief that a well-organized and systematically documented research process is key to advancing scientific knowledge.



