# UNIFACforACELFA

This repository contains the supplementary code for the article "A Molecular Thermodynamic Model for Improving Point-of-Care Diagnostics Incorporating the Lateral-Flow Immunoassay and Aqueous Two-Phase Systems" by Chawla *et al.*, *Ind. Eng. Chem. Res.* (2023).

The `TestLineAnalysis.ipynb` file contains the Python code used to analyze the test lines of LFA results.

The `MATLABFiles` directory contains relevant MATLAB files for the project. The majority were used to generate UNIFAC predictions, and `find_tie_line.m` contains the code used to determine tie lines for an LFA. Within the folder, `unifac.m`, `unifacSetUp.m`, and `unifacAijLLE` were adapted from Elliott and Lira (https://chethermo.net/, https://sourceforge.net/projects/chethermo/). The only major modifications to these files involved commenting explanations, and they were treated as a package to run the UNIFAC equations. The files `CalculatingK`, `comp_calc`, and `dlngammadDAB` were used to perform the actual mathematical computations to solve the governing equations of the system. To calculate the partition coefficient, run the `CalculatingK` file. 