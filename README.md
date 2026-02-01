# PHSI-msc-project
MATLAB code and results for Master's project on Polarization Hyperspectral Imaging (PHSI). 

## Repository Structure

- 'hsi/'  Contains code for **hyperspectral imaging analysis**, based on https://drive.google.com/drive/folders/1OmyS4ESn0ryNYrc_qjWoNNmOyTwHn_yj?usp=drive_link 
- 'pol/' - Contains code for **polarization imaging analysis**, based on https://github.com/benesprk/bc_project/tree/main
- 'pol_peanuts/' - Unrelated to the other codes, studies polarization states through Fourier analysis obtained with rotating QWP and polarizers, in a study similar to https://physlab.org/experiment/polarization-peanuts-with-fourier-analysis/
- 'pol_tests' - Contains saved polarization figures for individual standard polarization imaging analysis.
%- 'pol results+figures/' - Contains saved polarization figures and processed data .mat files.
%- 'hsi results+figures/' - Contains saved hsi figures and processed data .mat files.
- 'MAIN_example_hs_load_adapted2phsi_november.m' - Main script to run **HSI and polarization analysis** on data from the hyperspectral imaging system.
  - 'MAIN_example_hs_load_adapted2phsi_load_figures_november.m' - Loads data previously saved in .mat files and plots the respective DoLP and AoLP figures (avoiding the >1h long computation time).
  - 'MAIN_example_hs_load_adapted2phsi_compare_methods_november.m' - Plots polarization parameters in grayscale for each method used and compares DoLP and AoLP differences between those methods.
- 'MAIN_pol_analysis.m' - Main script to run **polarization-only analysis** on data from the polarization camera.
- 'MAIN_phsi_december.m' - Final script to **fuse HSI and polarization results**, after calibrating the HSI and polarization camera by an homographic transformation matrix.
  - 'MAIN_phsi_useMatrix_woSpectral_december.m' - Fuses HSI and polarization results, using an homographic transformation matrix previously computed for calibration and skips spectral analysis for specific points or regions available in the main og version.

