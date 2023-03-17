
# Precipitation buffers temperature-driven local extinctions of moths at warm range margins

### L.A. Hordley, R. Fox, A.J. Suggitt, N.A.D. Bourn

To reference the data or methods here, please cite the published paper [here](https://doi.org/10.1111/ele.14195). Contact Lisbeth Hordley with any questions at lhordley@butterfly-conservation.org. 

Data repository DOI

[![DOI](https://zenodo.org/badge/424910381.svg)](https://zenodo.org/badge/latestdoi/424910381)

## Overview

This repository contains code and data to:
* Categorise species as being a cool-adapted moth
* Calculate trailing edge shifts in elevation and latitude
* Estimate multidirectional shifts in species' range centroids
* Determine the drivers of local extinction of cool-adapted moths

## Instructions

Scripts should be run in the following order:

1. ```1_NMRS_data_cleaning```: taking out unused columns, and scaling lat/lon and easting/northing to 10km scale at the bottom left hand corner of hectad and remove Isle of Man records. Returns: NMRS_hectad_cleaned

2. ```2_hectad_recording_levels```: using NMRS_hectad_cleaned created in script 1, class each hectad as recorded, well-recorded or heavily recorded depending on regional species richness. Returns list of hectads with corresponding recording level for each time period (1975-1991 and 2012-2016): Hectad_recording_levels_1975_1991_2012_2016

3. ```3_elevation```: Extract elevation data, change to 10km resolution and merge with NMRS_hectad_cleaned. Returns NMRS_hectad_elevation.

4. ```4_climate```: Extract annual mean temperature and total precipitation data for each time period from the MetOffice CP18 dataset. Change to 10km resolution and match with hectads from NMRS_hectad_elevation. Returns NMRS_hectad_elevation_climate. 

5. ```5_cool_adapted_moths```: calculate each species' 75th temperature percetile between 1975-1991 to determine which species are 'cool-adapted'. Using species selected as 'cool-adapted', filter NMRS_hectad_elevation_climate by those 76 species. Returns NMRS_cool_moths_final

6. ```6_Range_margin_shift_NMRS```: Using NMRS_cool_moths_final, calculate trailing edge shift in elevation and latitude between two time periods for well recorded and heavily recorded hectads

7. ```7_Multidimensional_shift_NMRS```: Using NMRS_cool_moths_final, calculate direction and distance of centroid shift of each species between two time periods and species shift based on historic climate envelope

8. ```8_Local_extinction_NMRS```: Using NMRS_cool_moths_final, run GLMMs with local extinction as response and climate variables as explanatory variables. Plot results as heatmap. 
