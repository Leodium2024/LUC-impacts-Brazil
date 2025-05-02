"""
Created on 20/03/2025

@author: Thomas GERARD

Version: 1

# --------------------------------------------------
#  Main 
# --------------------------------------------------  

"""

from Tool_function import arrayToMap, rasterToArray, plot_map_with_mask, get_distribution_map, get_ottid
import os
import pandas as pd
import numpy as np


Main_directory = "E:\Brazil\Scenarios\Git_hub"
Input_directory = os.path.join(Main_directory,'Input')
Output_directory = os.path.join(Main_directory,'Output') 

#Define scale colors for mapping
colors1 = [(1, 1, 1), (0, 1, 0)]  # White to Green
colors2 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0)]  # White to Green to Dark Green 
colors3 = [(1, 0, 0), (1, 1, 1), (0, 1, 0)] # red to white to Green
colors4 = [(1, 0, 0), (1, 1, 1), (0, 0, 1)] # red to white to Blue
colors5 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors5 = [(1, 1, 1),(0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors6 = [(1, 1, 1), (0.5, 1, 0.5), (0, 1, 0), (0, 0.5, 0), (0, 0.25, 0), (0, 0, 0)]

#Select which scenario SSP1, SSP2 or SSP3
#scenarios = ['SSP2','SSP1','SSP2','SSP3', 'Observed'] 
scenarios = ['Observed','SSP1','SSP2','SSP3']
years_id=[0,1,2,3,4,5,6,7] 
years = ['2015', '2020', '2025', '2030', '2035', '2040','2045', '2050']

Scale_names = ['National', 'Amazon', 'Cerrado', 'Caatinga', 'Mata', 'Pantanal', 'Pampa']

land_use_types = [
    'Forest', 'Grassland', 'Pasture', 'Cropland', 'Forestry', 
    'Other', 'Annual', 'Perennial', 'Semiperennial'
]

env_variables = {
    'Elevation': 'Elevation.tif',
    'Slope': 'Slope.tif',
    'Bio_1': 'Bio_1.tif',
    'Bio_12': 'Bio_12.tif',
    'Bio_4': 'Bio_4.tif',
    'Bio_15': 'Bio_15.tif'
}


# --------------------------------------------------
#  Calculation of the objectives
# --------------------------------------------------  

Result_Total_carbon = {}
Result_Carbon_map = {}
Result_Total_revenue = {}
Result_revenue_map = {}
Result_Area = {}
Result_SR_map = {}

def Calculate_carbon(land_use, Crefs):
    
    Carbon_land_use = {
                        land_type: land_use[land_type] * Crefs[land_type]
                        for land_type in land_use_types if land_type not in {'Cropland', 'Other'}
                       }
    
    Carbon_map = sum(Carbon_land_use.values(), np.zeros_like(next(iter(Carbon_land_use.values()))))
    
    Total_Carbon = np.sum(Carbon_map*10000) #Multiply by 10000 to convert the ton/ha in ton per cell (10000 ha)

    Total_Carbon = Total_Carbon/1000000000 # to convert in Gt   
    
    return Total_Carbon, Carbon_map
        
         

def Calculate_Economy(land_use,Erefs):

    Revenue_land_use = {
                        land_type: land_use[land_type] * Erefs[land_type]
                        for land_type in ['Pasture','Annual', 'Perennial', 'Semiperennial']
                       }
    
    
    Revenue_map = sum(Revenue_land_use.values(), np.zeros_like(next(iter(Revenue_land_use.values()))))
    
    Total_Revenue = np.sum(Revenue_map*10000) #Multiply by 10000 to convert the USD/ha in ton per cell (10000 ha)
    
    Total_Revenue = Total_Revenue/1000000000 # to convert in billion USD
    
    return Total_Revenue, Revenue_map



def Calculate_Biodiversity (land_use, Coeficient, Performance, env, Ott_id):
    
    SR_map = pd.DataFrame(0, index=range(464), columns=range(544))
    Areas={}
    
    for Species_ID in Ott_id:
        # Get species distribution map
        Species_distribution = get_distribution_map(Input_directory, Species_ID, Coeficient, Performance, land_use, env)
    
        DistributionArea = Species_distribution.sum().sum()
        Areas[Species_ID] = DistributionArea
        
        # Update species richness (SR)
        SR_map += Species_distribution
   
    return Areas, SR_map



#Loop throught the scenarios
for Sc in scenarios:
    # Limit the years if the scenario is 'Observed'
    years_to_process = years[:2] if Sc == 'Observed' else years
    years_id_to_process = years_id[:2] if Sc == 'Observed' else years_id
    
    # Loop through the years
    for year, year_id in zip(years_to_process, years_id_to_process):
        
        #Import Land use
        land_use = {
                     land_type: rasterToArray(os.path.join(Input_directory, 'Land_Use', f'{Sc}', f'{year}_{land_type}.tiff'))
                     for land_type in land_use_types
                    }
        
        #Import Carbon reference
        Crefs = {
                  land_type: rasterToArray(os.path.join(Input_directory, 'Carbon', f'Cref_{land_type}.tiff'))
                  for land_type in land_use_types if land_type not in {'Cropland', 'Other'}
                 }
        
        #Import Economic reference
        Erefs = {
                  land_type: rasterToArray(os.path.join(Input_directory, 'Economy', f'{Sc}', f'{year}_Eref_{land_type}.tiff'))
                  for land_type in ['Pasture','Annual', 'Perennial', 'Semiperennial']
                 }
        
        #Import SDM parameter
        Coeficient = pd.read_excel(os.path.join(Input_directory,'Biodiversity', 'Species_coef.xlsx'))
        Order = pd.read_excel(os.path.join(Input_directory,'Biodiversity', 'Range_Ott_id.xlsx'))
        Performance = pd.read_excel(os.path.join(Input_directory,'Biodiversity', 'Performance.xlsx'))
        Ott_id = get_ottid(Performance)
        env = {
            key: rasterToArray(os.path.join(Input_directory, 'Biodiversity', filename))
            for key, filename in env_variables.items()
            }

        
        #Calculation of the objectives
        Total_carbon,carbon_map = Calculate_carbon(land_use,Crefs)
        Total_revenue, revenue_map = Calculate_Economy(land_use,Erefs)
        Areas, SR_map = Calculate_Biodiversity (land_use, Coeficient, Performance, env, Ott_id)
        
        #Store the result in dictionary per scenarios and year
        Result_Total_carbon.setdefault(Sc, {})[year] = Total_carbon
        Result_Carbon_map.setdefault(Sc, {})[year] = carbon_map
        Result_Total_revenue.setdefault(Sc, {})[year] = Total_revenue
        Result_revenue_map.setdefault(Sc, {})[year] = revenue_map
        Result_Area.setdefault(Sc, {})[year] = Areas
        Result_SR_map.setdefault(Sc, {})[year] = SR_map
        
        #Plot the maps
        plot_map_with_mask(carbon_map,0,200,colors5,"C Stock (t/ha)", f'Land carbon stock in {year} under {Sc}')
        plot_map_with_mask(revenue_map,0,1000,colors5,"Revenue (USD/ha)", f'Agricultural revenue in {year} under {Sc}')
        plot_map_with_mask(SR_map,0,100,colors5,"Species Richness (nbr_species)", f'Species Richness in {year} under {Sc}')

        
# ---------------------------------------------------------------
#  Spatial change
# ---------------------------------------------------------------  

mask=rasterToArray(os.path.join(Input_directory,'Mask_Albers.tif'))

for Sc in ['SSP1','SSP2','SSP3']:
    
    #Spatial changes
    Carbon_change_map = Result_Carbon_map[Sc]['2050'] - Result_Carbon_map[Sc]['2015']
    Revenue_change_map = Result_revenue_map[Sc]['2050'] - Result_revenue_map[Sc]['2015']
    SR_change_map = Result_SR_map[Sc]['2050'] - Result_SR_map[Sc]['2015']
    
    arrayToMap(Carbon_change_map, os.path.join(Output_directory, f'{Sc}_Change_carbon'),mask,Input_directory)
    arrayToMap(Revenue_change_map, os.path.join(Output_directory, f'{Sc}_Change_Revenue'),mask,Input_directory)
    arrayToMap(SR_change_map, os.path.join(Output_directory, f'{Sc}_Change_SR'),mask,Input_directory)
    arrayToMap(Result_Carbon_map[Sc]['2015'], os.path.join(Output_directory, f'{Sc}_2015_carbon'),mask,Input_directory)
    arrayToMap(Result_revenue_map[Sc]['2015'], os.path.join(Output_directory, f'{Sc}_2015_Revenue'),mask,Input_directory)
    arrayToMap(Result_SR_map[Sc]['2015'], os.path.join(Output_directory, f'{Sc}_2015_SR'),mask,Input_directory)

    plot_map_with_mask(Carbon_change_map,-50,50,colors4,"C Stock (t/ha)", f'Change in Land carbon stock under {Sc}')
    plot_map_with_mask(Revenue_change_map,-500,500,colors4,"Revenue (USD/ha)", f'Change in revnue under {Sc}')
    plot_map_with_mask(SR_change_map,-20,20,colors4,"Species Richness (Number of species)", f'Change in Species Richness under {Sc}')
    plot_map_with_mask(Result_Carbon_map[Sc]['2015'],0,200,colors5,"C Stock (t/ha)", f'Land carbon stock in 2015 under {Sc}')
    plot_map_with_mask(Result_revenue_map[Sc]['2015'],0,1000,colors5,"Revenue (USD/ha)", f'Agricultural revenue in 2015 under {Sc}')
    plot_map_with_mask(Result_SR_map[Sc]['2015'],0,100,colors5,"Species Richness (nbr_species)", f'Species Richness in 2015 under {Sc}')


# ---------------------------------------------------------------
#  Change at national scale
# ---------------------------------------------------------------  

   
#Carbon
Carbon = [list(ssp_data.values()) for ssp_data in Result_Total_carbon.values()]
Carbon = pd.DataFrame(Carbon, columns=years, index=Result_Total_carbon.keys())
Carbon_change_absolute = Carbon.sub(Carbon['2015'], axis=0)
    
with pd.ExcelWriter(os.path.join(Output_directory, 'Carbon.xlsx')) as writer:
    Carbon.to_excel(writer, sheet_name='Total')
    Carbon_change_absolute.to_excel(writer, sheet_name='Change')

#Economy      
Revenue = [list(ssp_data.values()) for ssp_data in Result_Total_revenue.values()]
Revenue = pd.DataFrame(Revenue, columns=years, index=Result_Total_revenue.keys())
Revenue_change_absolute = Revenue.sub(Revenue['2015'], axis=0)

with pd.ExcelWriter(os.path.join(Output_directory, 'Revenue.xlsx')) as writer:
     Revenue.to_excel(writer, sheet_name='Total')
     Revenue_change_absolute.to_excel(writer, sheet_name='Change')


#Biodiveristy 
Area_change = pd.DataFrame(columns=years)
for Sc in scenarios:
    SSP_Area = pd.DataFrame(Result_Area[Sc])  # Convert to DataFrame
    SSP_Area_change = SSP_Area.sub(SSP_Area['2015'], axis=0)  # Absolute change
    SSP_Area_change_relative = SSP_Area_change.div(SSP_Area['2015'], axis=0) * 100  # Relative change
    Average_Change_per_year = SSP_Area_change_relative.mean(axis=0)  # Mean per column (year)
    # Store the result in the final DataFrame (each row is a scenario)
    Area_change.loc[Sc] = Average_Change_per_year


with pd.ExcelWriter(os.path.join(Output_directory, 'Biodiversity.xlsx')) as writer:
     Area_change.to_excel(writer, sheet_name='Relative change')
     
     
     
     
     