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


Impacts: 
 - Impact 1: Change in Land Carbon Stock (ton/ha)
 - Impact 2: Change in Mammal Species Richness (Nbr of species/ grid cell)
 - Impact 3: Change in Agricultural Revenues (USD/ha/year)

Land use projections designed by Silva Bezerra et al., (2022)(https://doi.org/10.1371/journal.pone.0256052):
 - SSP1-1.9
 - SSP2-4.5
 - SSP3-7.0



Script
-------------------------------------------------- 

1) "Main.py, in which the following elements are calculated across the LU projection:
    -  Change in Carbon stock across SSP1-1.9, SSP2-4.5, SSP3-7.0
    -  Change in Mammal Species Richness across SSP1-1.9, SSP2-4.5, SSP3-7.0
    -  Change in Agricultural revenue across SSP1-1.9, SSP2-4.5, SSP3-7.0
    The output data are stored in the Output directory
    
2) "Tool_functions.py", contains some functions called in other .py file

3) "Create_input.py", python script to:
   - Generate input data to calculate the Carbon Stock
   - Generate input data to calculate the agricultural revenue
   - Extract the land use map


Input directory
------------------------------------------------------------------ 
Contain all the necceassary input data to run the script in "Main.py", including: 

- Input data to calculate changes in carbon stock stored in Input/Carbon (See Generation of Input: carbon Stock, Agricultural Revenue, Land-use)
- Input data to calculate changes in agricultural revenue stored in Input/Economy (See Generation of Input: carbon Stock, Agricultural Revenue, Land-use)
- Land use projections under SSP1-1.9, SSP2-4.5 and SSP3-7.0 stored in Input_Land use (See Generation of Input: carbon Stock, Agricultural Revenue, Land-use)
- Input data to calculate changes in Species Richness stored in Input/Biodiveristy (See Generation of Input: Species Richness)


Generation of Input: carbon Stock, Agricultural Revenue, Land-use
------------------------------------------------------------------ 
"Create_input.py" use Secondary data to create some input data.

The secondary data are stored in Input/Preliminary_Input, as follow:
- Secondary data calculate changes in carbon stock stored in Preliminary_Input/Carbon
- Secondary data to calculate changes in agricultural revenue stored in Preliminary_Input/Economy
- Land use projections designed by Silva Bezerra et al., (2022) (also available in: https://doi.org/10.1371/journal.pone.0256052)


Generation of Input: Species Richness
------------------------------------------------------------------ 

Performance.xlsx and Species_coef.xlsx,i.e., input of the Species Richness calculation, were generated using the following script:

	- InputBiodiveristy1.R # Taxonomic and geographic cleaning of biodiveristy database + thinning (1 records per 10km)
	
	- InputBiodiveristy2.py # Create Pseudo-absence data

	- InputBiodiveristy3.R # Calibrate and validate the SDMs These files can be found in Premilinary_Input/Biodiveristy.


The input of these script can be downloaded at : 10.5281/zenodo.15497375









