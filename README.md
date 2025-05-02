LUC impacts Brazil
-------------------------------------------------- 
Created on 02/05/2025
@author: Thomas GERARD
Version: 1

Python Environement set-up:
- Python Version: 3.11.5
- pandas: 2.0.3
- numpy: 1.24.3
- spyder: 5.4.3
- GDAL: 3.6.2
- matplotlib: 3.7.2
- netCDF4: 1.6.2

R Environement set-up:
- CoordinateCleaner version: 3.0.1 
- readxl version: 1.4.3 
- writexl version: 1.5.0 
- rotl version: 3.1.0 
- spThin version: 0.2.0
- caret version: 6.0.94 
- pROC version: 1.18.5 
- flexsdm version: 1.3.4 
- sf version: 1.0.16 
- spdep version: 1.3.4 
- ROCR version: 1.0.11 
- MASS version: 7.3.60.2 

Impacts: 
 - Impact 1: Change in Land Carbon Stock (ton/ha)
 - Impact 2: Change in Mammal Species Richness (Nbr of species/ grid cell)
 - Impact 3: Change in Agricultural Revenues (USD/ha/year)

Land use scenarios designed by Silva Bezerra et al., (2022):
 - SSP1-1.9
 - SSP2-4.5
 - SSP3-7.0
(https://doi.org/10.1371/journal.pone.0256052)


Script directory
-------------------------------------------------- 

The 'Script'  directory contain all the following script:
1) "Main.py, in which the following elements are calculated across the LU projection:
    -  Carbon stock
    -  Mammal Species Richness
    -  Agricultural revenue
    
2) "Tool_functions.py", contains some functions called in other .py file

3) "Create_input.py", python script to:
   - Generate input data to calculate the Carbon Stock
   - Generate input data to calculate the agricultural revenue
   - Extract the land use map


Input directory
-------------------------------------------------- 
Contain all the necceassary input data

With data to calculate the:
- Biodiveristy objective stored in the Biodiveristy directory
- Climate objective stored in the Climate directory
- Economic objective stored in the Economic directory

These input data were created using "Create_input.py", 
which itslef use input data stored in the Premilinary_Input directory 

The input of the Biodiveristy model were generated using: 
- InputBiodiveristy1.R             # Taxonomic and geographic cleaning of biodiveristy database + thinning (1 records per 10km)
- InputBiodiveristy2.py            # Create Pseudo-absence data
- InputBiodiveristy3.R             # Calibrate and validate the SDMs
These file can be found in the Premilinary_Input directory 





