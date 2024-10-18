# Overview

This repository contains files for “Forecasting Range Shifts of Dioecious Plant Species Under Climate Change.”

# Files

Note: Due to the large number of files in the repository, not all files are described below.

## data 

Note: *Raster files are too large to be provided in a public repository and are stored on a local machine. [Here](https://www.dropbox.com/scl/fo/em8fok5gqwyhsum1hmier/ALoREAEAcgsADRWDyCqR4FA?rlkey=d92vrqi4ue5osdd857qcjlb9r&dl=1) is  is the link to download the rasters.*

**Poa_pra** *(.csv)* - Precipitation data.

**Poa_tas** *(.csv)* - Temperature data.

**poar_ms_quantities** *(.rds)* - Data for field and experiment data.

**poar_season** *(.rds)* - Vital models output.

## Scripts

**01_Map**  *(.R)* - Generate all the maps used in the manuscript.

**02_Vital_rates_models** *(.R)* - Vital rate models (survival, growth, flowering, fertility) for subsequent use in MPMs. This is computationally intensive and may take a while.

**03_Population viability_vs_niche** *(.R)* - Matrix Projection Models from vital (survival, growth, flowering,fertility) to estimate species niche . This is computationally intensive and may take a while.

**04_Population viability_vs_geography**  *(.R)* - Projection of the Probabibility of viability into geographic space. This is computationally intensive and may take a while.

**05_climate_vs_Sex_ratio** *(.R)* - Sex ratio models. This is computationally intensive and may take a while.

**06_LTRE** *(.R)* - Life Table Response Experiment. This is computationally intensive and may take a while.

**twosexMPMLTRE** *(.R)* - Function used in Matrix Projection Models.

**poar_season** *(.stan)* - code stan used for the vital rates models.

## Others

**Manuscript** *(folder)* - A folder with code and various other objects used to generate the manuscript, including figures.
