# Range-shift
Scripts for range-shift analysis using NMRS data

1_NMRS_data_cleaning: taking out unused columns, and scaling lat/lon and easting/northing to 10km scale at the bottom left hand corner of hectad and remove Isle of Man records. Returns: NMRS_hectad_cleaned

2_hectad_recording_levels: using NMRS_hectad_cleaned created in script 1, class each hectad as recorded, well-recorded or heavily recorded depending on regional species richness. Returns list of hectads with corresponding recording level for each time period (1975-1991 and 2012-2016): Hectad_recording_levels_1975_1991_2012_2016

3_elevation: Extract elevation data, change to 10km resolution and merge with NMRS_hectad_cleaned. Returns NMRS_hectad_elevation.

4_climate: Extract annual mean temperature and total precipitation data for each time period from the MetOffice CP18 dataset. Change to 10km resolution and match with hectads from NMRS_hectad_elevation. Returns NMRS_hectad_elevation_climate. 

5_cool_adapted_moths: calculate each species' 75th temperature percetile between 1975-1991 to determine which species are 'cool-adapted'. Using species selected as 'cool-adapted', filter NMRS_hectad_elevation_climate by those 76 species. Returns NMRS_cool_moths_final

6_Range_margin_shift_NMRS: Using NMRS_cool_moths_final, calculate trailing edge shift in elevation and latitude between two time periods for well recorded and heavily recorded hectads

7_Multidimensional_shift_NMRS: Using NMRS_cool_moths_final, calculate direction and distance of centroid shift of each species between two time periods and species shift based on historic climate envelope

8_Local_extinction_NMRS: Using NMRS_cool_moths_final, run GLMMs with local extinction as response and climate variables as explanatory variables. Plot results as heatmap. 
