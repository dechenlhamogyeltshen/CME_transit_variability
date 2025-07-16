# CME_Transit_Variability
Python scripts and data to investigate the variability of CME (Coronal Mass Ejection) transit times to Earth and arrival speeds at Earth. Uses the HUXt (Heliospheric Upwind Extrapolation with time dependence) solar wind model and 50 years of ambient solar wind solutions from the MAS (Magnethydrodynamic Algorithm around a Sphere) coronal model (https://www.predsci.com/mhdweb/home.php).

code/article.py contains script to simulate the CME propagations using HUXt.
code/analysis.ipynb contains script to generate plots and perform analysis tests.

## Citations

HUXt software: DOI:[10.5281/zenodo.4889326](https://doi.org/10.5281/zenodo.4889326) 

Barnard and Owens. (2022), *HUXt - An open source, computationally efficient reduced-physics solar wind model, written in Python*, Frontiers in Physics [10.3389/fphy.2022.1005621](https://doi.org/10.3389/fphy.2022.1005621)

Owens et al. (2020),  *A Computationally Efficient, Time-Dependent Model of the Solar Wind for Use as a Surrogate to Three-Dimensional Numerical Magnetohydrodynamic Simulations*,  Sol Phys, DOI:[10.1007/s11207-020-01605-3.svg](https://doi.org/10.1007/s11207-020-01605-3)
