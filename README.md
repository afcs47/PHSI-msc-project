# PHSI-msc-project
MATLAB code and results for Master's project on Polarization Hyperspectral Imaging (PHSI). 

## Repository Structure

- 'hsi/'  Contains code for **hyperspectral imaging analysis**, based on https://drive.google.com/drive/folders/1OmyS4ESn0ryNYrc_qjWoNNmOyTwHn_yj?usp=drive_link 
- 'pol/' - Contains code for **polarization imaging analysis**, based on https://github.com/benesprk/bc_project/tree/main
- 'pol_peanuts/' - Unrelated to the other codes, studies polarization states through Fourier analysis obtained with rotating QWP and polarizers, in a study similar to https://physlab.org/experiment/polarization-peanuts-with-fourier-analysis/
- 'results/' - Contains saved polarization figures and processed data.
- 'figures/' - Contains saved hsi figures.
- 'MAIN_example_hs_load_adapted2phsi.m' - Main script to run **HSI and polarization analysis** on data from the hyperspectral imaging system.
- 'MAIN_pol_analysis.m' - Main script to run **polarization-only analysis** on data from the polarization camera.
- 'MAIN_phsi.m' - Final script to **fuse HSI and polarization results**, after calibrating the HSI and polarization camera by an homographic transformation matrix.
