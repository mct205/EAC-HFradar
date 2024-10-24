EAC HF radar scripts

This repository contains the MATLAB code for processing the HF radar data in the Eastern Australia. 

The Data-processing directory is used for processing the radar data using 2dVar approach. 
The Plots directory contains the code for plotting the data. 

For the code to run, some external MATLAB toolboxes and data are needed: 
- m_map toolbox (v1.4o): https://www-old.eoas.ubc.ca/~rich/map.html
- UTide toolbox (v1.0.0): https://au.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions

The codes contained here included the following: 

Data management: 
-	plot_radar_coverage_and_mooring_and_drifters.m: List all the data used in the study
-	plot_RadTemporal_stats.m: List all the radar data available from 4 radar sites in the eastern Australia

Plot and analysis scripts: 
-	plot_avg_current.m: Computes the multi-year and seasonal averaged surface current velocities in both regions and plot: NEWC and COF
-	plot_cross_section.m: Plot the cross-section and assess the interannual variability of the surface current velocities from both regions 
-	 get_core_vel_eac.m: Identifying the Eastern Australian current (EAC) by the jet following method based on the algorithm of Archer et al. (2017)
-	compute_spectral_density_both_radars.m: Calculate the kinetic spectral density of the radar-derived velocities.
-	get_tidal_analysis_Utide.m: Compute tidal harmonic of the surface current velocities, Utide analysis package is required (https://www.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions)
-	plot_tidal_current.m: For plotting tidal analysis results from get_tidal_analysis_Utide.m
-	compare_data_cutoff_COF.m and compare_data_cutoff_NEWC.m: For testing the performance of 2dVar approach based on several synthetic gap scenarios. The computation takes a long time so the results are kept at folder processed_data for analysis

Some functions required for the code to run: 
-	spectrelisse2.m: For computing the spectrum
-	func_get_core_vel_eac.m: algorithm for identifying the EAC using COF radar data 

Notes: The HF radar data contains the reprocessed radar dataset is available here: 10.5281/zenodo.13984639
