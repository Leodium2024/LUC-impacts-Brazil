# LUC impacts Brazil
1. Description

The main file is MAIN.py, in which the following elements are calculated across the LU projection:
- Carbon model: Total carbon stock, SOC stock, Biomass carbon stock
- Biodiversity model: SR, WE, CWE and Average species area per order
- Agriculture model: Total agricultural revenue, Pasture revenue, Cropland revenue

2. Inputs

The input of the carbon model were generated using:
- InputCarbon1.py         # Calculate reference SOC ref based on LU, soil type, climate zone and Mapbiomass SOC map (input of SOC stock calculation)
- InputCarbon2.py         # Calculate the distribution of biome and phytosiome in each cell of the LU projection (input of SOC stock calculation)
- InputCarbon3.py         # Calculate the distribution of Soil type and Climate in each cell of the LU projection (input of Biomass C stock calculation)

The input of the Biodiveristy model were generated using: 
- InputBiodiveristy1.R             # Taxonomic and geographic cleaning of biodiveristy database + thinning (1 records per 10km)
- InputBiodiveristy2.py            # Create Pseudo-absence data
- InputBiodiveristy3.R             # Calibrate and validate the SDMs

The input  of the Agricultural revenue model are from IMAGE 3.2

5. Output






