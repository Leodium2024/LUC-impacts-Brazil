"""
Created on 20/03/2025

@author: Thomas GERARD

Version: 1

# --------------------------------------------------
#  Create Input for objectives calculations
# --------------------------------------------------  

Land use map:(!! input file have to be unzipped!!)
    => Extract all the land use map
    => Split the mosaic category
    => Split the cropland category


Carbon Sequestration objectives:
    - C Stock ref map Forest
    - C Stock ref map Grassland
    - C Stock ref map Pasture
    - C Stock ref map Perennial cropland
    - C Stock ref map Annual Cropland
    - C Stock ref map SemiPerennial Cropland
    - C Stock ref map Forest Plantation
    
Agro-economic development Objectives:
    - Revenue ref map Pasture
    - Revenue ref map Perennial cropland
    - Revenue ref map SemiPerennial Cropland
    - Revenue ref map Annual Cropland

"""

import numpy as np
import pandas as pd
import os
import netCDF4 as nc
from Tool_function import arrayToMap, rasterToArray, reclass_dataframe, get_raster_combi, plot_map_with_mask


Main_directory = "E:\Brazil\Scenarios\Git_hub"
Input_directory = os.path.join(Main_directory,'Input\Preliminary_Input')
Output_directory = os.path.join(Main_directory,'Input') 

#Define scale colors for mapping
colors1 = [(1, 1, 1), (0, 1, 0)]  # White to Green
colors2 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0)]  # White to Green to Dark Green 
colors3 = [(1, 0, 0), (1, 1, 1), (0, 1, 0)] # red to white to Green
colors4 = [(1, 0, 0), (1, 1, 1), (0, 0, 1)] # red to white to Blue
colors5 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors5 = [(1, 1, 1),(0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors6 = [(1, 1, 1), (0.5, 1, 0.5), (0, 1, 0), (0, 0.5, 0), (0, 0.25, 0), (0, 0, 0)]

#Select which scenario SSP1, SSP2 or SSP3
scenarios = ['SSP1','SSP2','SSP3'] 
years_id=[0,1,2,3,4,5,6,7] 
years =['2015', '2020', '2025', '2030', '2035', '2040','2045', '2050']

# --------------------------------------------------
# Extraction of land use maps
# --------------------------------------------------

def Get_Year_NetCDF (Scenario,lu,year,mask):
    """ Upload of the LU layer for a specific LU and year
    Scenario is the path leading to the .nc file 
    lu is the name of the land use type to be extracted (string)
    Year to be extracted (int)"""
    data = np.zeros((464, 544))
    
    if lu == 'pastp': 
        # In the nc provided by Bezera et al. (2022), the layer pasture is the same than the layer Grassland
        # Therefore here the pature layer is deducted by doing 1 - all the others LU
        nc_file = nc.Dataset(Scenario, 'r+')
    
        # Define the layers to replace pastp with
        layers_to_replace_with = ['veg', 'gveg', 'agric', 'mosc', 'fores', 'others']
                
        # Iterate over each layer to replace pastp with
        for layer in layers_to_replace_with:
            # Read the data from the respective layer and add to pastp
            data += nc_file.variables[layer][year, :, :]
                
        # Invert the data to match the logic "1 - 'veg' - 'gveg' -  'agric'- 'mosc'- 'fores' -'others'"
        data = 1 - data
        data = np.nan_to_num(data, nan=0)
        data = data * mask
        # Replace the values in 'pastp' layer with the calculated data
        nc_file.variables['pastp'][year, :, :] = data
            
            
    else:
        #If not pasture, read the netcdf as usual
        NetCDF = nc.Dataset(Scenario, 'r')
        Var = NetCDF.variables[lu]
        data = Var[year,:,:]
        data =  0 + data
        data = np.nan_to_num(data, nan=0)
    
    return data

def get_agri_weight():
    State = rasterToArray( os.path.join(Input_directory,'State_Albers.tif')) 
    Weight = pd.read_excel(os.path.join(Input_directory,'Land_Use','Percentage_agri_mos.xlsx'),sheet_name='Agriculture')
    Wpe = reclass_dataframe(Weight, State, 0 , 2)
    Ws = reclass_dataframe(Weight, State, 0 , 3)
    Wt = reclass_dataframe(Weight, State, 0 , 1)
    return Wpe,Ws,Wt



def Split_Mosaic ():
    ''' This function will return a array containing the proportion in each cell of 
        Perenial, Temporary, semiperenial agriculture and 
        Mosaic Agri-Forest, mosaic  Agri-Grassland
        All those proportion are state-specific'''
    Mesoregion = rasterToArray( os.path.join(Input_directory,'Land_Use','Mesoregion_Albers.tif'))
    Weight = pd.read_excel(os.path.join(Input_directory,'Land_Use','Percentage_agri_mos.xlsx'),sheet_name='Mosaic')
    Wf = reclass_dataframe(Weight, Mesoregion,0,  1)
    Wg = reclass_dataframe(Weight, Mesoregion,0,  2)
    Wp = reclass_dataframe(Weight, Mesoregion,0,  3)
    Wa = reclass_dataframe(Weight, Mesoregion,0 , 4)
    Wfo = reclass_dataframe(Weight, Mesoregion,0 , 5)
    Wot = reclass_dataframe(Weight, Mesoregion,0 , 6)

    return  Wf, Wg, Wp, Wa, Wfo, Wot


def get_lu():
    
    mask=rasterToArray(os.path.join(Input_directory,'Mask_Albers.tif'))
    Area = {}
    for sc in scenarios:
        Area [sc] = []
        Lu_path = os.path.join(Input_directory,'Land_Use',f'{sc}.nc')  
        for year, year_id in zip(years, years_id):
            if sc == 'SSP1' and year_id == 3: #Missing data
            #if  year_id == 3: #Missing data
                Cropland = (Get_Year_NetCDF(Lu_path,'agric',2, mask)+Get_Year_NetCDF(Lu_path,'agric',4, mask))/2
                Forest = (Get_Year_NetCDF(Lu_path,'veg',2, mask)+Get_Year_NetCDF(Lu_path,'veg',4, mask))/2
                Forestry = (Get_Year_NetCDF(Lu_path,'fores',2, mask)+Get_Year_NetCDF(Lu_path,'fores',4, mask))/2
                Other = (Get_Year_NetCDF(Lu_path,'others',2, mask)+Get_Year_NetCDF(Lu_path,'others',4, mask))/2
                Mosaic = (Get_Year_NetCDF(Lu_path,'mosc',2, mask)+Get_Year_NetCDF(Lu_path,'mosc',4, mask))/2
                Pasture = (Get_Year_NetCDF(Lu_path,'pastp',2, mask)+Get_Year_NetCDF(Lu_path,'pastp',4, mask))/2
                Grassland = (Get_Year_NetCDF(Lu_path,'gveg',2, mask)+Get_Year_NetCDF(Lu_path,'gveg',4, mask))/2
                 
            else:
                Cropland = Get_Year_NetCDF(Lu_path,'agric',year_id, mask)
                Forest = Get_Year_NetCDF(Lu_path,'veg',year_id, mask)
                Forestry = Get_Year_NetCDF(Lu_path,'fores',year_id, mask)
                Other = Get_Year_NetCDF(Lu_path,'others',year_id, mask)
                Mosaic = Get_Year_NetCDF(Lu_path,'mosc',year_id, mask)
                Pasture = Get_Year_NetCDF(Lu_path,'pastp',year_id, mask)
                Grassland = Get_Year_NetCDF(Lu_path,'gveg',year_id, mask)
                
            #Split Mosaic 
            Wf, Wg, Wp, Wa, Wfo, Wot = Split_Mosaic() 
            Forest = Forest + (Mosaic*Wf)
            Pasture = Pasture + (Mosaic*Wp)
            Forestry = Forestry + (Mosaic*Wfo)
            Cropland = Cropland + (Mosaic*Wa)
            Grassland = Grassland + (Mosaic*Wg)
            Other = Other + (Mosaic*Wot)
                
            #Split Cropland
            Wpe,Ws,Wt = get_agri_weight()
            Temporary=Wt*Cropland
            Perennial= Wpe*Cropland
            Semiperennial= Ws*Cropland
                

            arrayToMap(Forest, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Forest'),mask,Input_directory)
            arrayToMap(Grassland, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Grassland'),mask,Input_directory)
            arrayToMap(Pasture, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Pasture'),mask,Input_directory)
            arrayToMap(Cropland, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Cropland'),mask,Input_directory)
            arrayToMap(Forestry, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Forestry'),mask,Input_directory)
            arrayToMap(Other, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Other'),mask,Input_directory)
            arrayToMap(Temporary, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Annual'),mask,Input_directory)
            arrayToMap(Perennial, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Perennial'),mask,Input_directory)
            arrayToMap(Semiperennial, os.path.join(Output_directory,'Land_Use',f'{sc}', f'{year}_Semiperennial'),mask,Input_directory)
            
            # Calculate sums per category
            category_sums = {
                'Year': year,
                'Forest': 100*np.sum(Forest),
                'Grassland': 100*np.sum(Grassland),
                'Pasture': 100*np.sum(Pasture),
                'Cropland': 100*np.sum(Cropland),
                'Forestry': 100*np.sum(Forestry),
                'Other': 100*np.sum(Other),
                'Annual': 100*np.sum(Temporary),
                'Perennial': 100*np.sum(Perennial),
                'Semiperennial': 100* np.sum(Semiperennial),
            }
            Area[sc].append(category_sums)
            
    # Export to Excel with a sheet per scenario
    output_file = os.path.join(Output_directory,'Land_Use', 'Land_Use_Area.xlsx')
    with pd.ExcelWriter(output_file) as writer:
        for sc, data in Area.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=sc, index=False)

get_lu()



# --------------------------------------------------
# Creation of the input for Carcbon Objective
# --------------------------------------------------


def get_Fc (LU, StateBiome_Map):
     
     ''' This function will return a array containing the Fc in each cell
         The Fc value are State and biome specific'''
         
     Fc_LU = pd.read_excel(os.path.join(Input_directory,'Carbon','Fc.xlsx'),sheet_name=LU)
     F = reclass_dataframe(Fc_LU, StateBiome_Map, 0 , 1)
     return F
  
def get_SOC_ref (LU,ID_raster,SoilClim_Per):
    
    ''' This function will return an array containing the SOCref for the different LU.
        For each cell the SOCref have been weighted by the proportion of each type of SoilClim.
        SoilClim is a dataframe contaning the proportion of each Soil climate type in every cells
        SOC is a dataframe containing the SOC reference value for each soil type combination snf Land use type
        ID_raster is a np.array contating the Id of each cell
        1)SoilClim and SOC: The two DataFrames being merged in one dataframe that contained the proportion of each soilclim and the corresponding SOC reference value
          how='left': it's a left join, meaning that all the rows from the SoilClim DataFrame will be retained, and matching rows from the SOC DataFrame will be appended where possible.
          left_on='ClimSoil': This parameter specifies the column in the left DataFrame (SoilClim) to join on. 
          right_on='SoilClim': This parameter specifies the column in the right DataFrame (SOC) to join on.
        2)Creation of a new column that contained the SOC ref weighted by the proportion of each SoilClim in each cell
        3)Sum of the weighted SOC ref to obtain the Final SOC ref per cell
        4)Reclass of the ID dataframe with the weighted SOC ref 
    
    '''
    
    # reference value of SOC per climSoil type for this specific LU
    SOCref = pd.read_excel(os.path.join(Input_directory,'Carbon','SOCref_Subordem.xlsx'),sheet_name=LU) 
    
    # Merge the two DataFrames based on the condit
    merged_df = pd.merge(SoilClim_Per, SOCref, how='left', left_on='ClimSoil', right_on='ClimSoil')
    
    #Create a new column that is equal to the weight SOCref (weight = proportion of the area covered by this ClimSoil type)
    merged_df['SOC'] = merged_df['%'] * merged_df['Average']
    
    # Give a table that is the soil ref per cell (sum of all the weight SOCref)
    summed_table = merged_df.groupby('ID')['SOC'].sum().reset_index()
    
    output = reclass_dataframe(summed_table,ID_raster,0,1)
    

    return output        


def get_biomass_ref(LU, ID_raster, State_ID, BiomePhyto_per,Clim_raster):

    if LU in {'Forest','Grassland'}:
        
            merged = pd.merge(BiomePhyto_per, State_ID , how='left', left_on='ID', right_on='ID_cell')
            merged.drop(columns=['ID_cell'], inplace=True)
            
            
            # Reference value for the biome of Amazone, Caatinga, Pampa, Mata and Pantanal
            filtered_OtherBiome = merged[(merged['BiomePhyto'] < 200) | (merged['BiomePhyto'] > 300)] #As MCTI has state specific value for Cerrado, we filter to keep all the other biome
            ref= pd.read_excel(os.path.join(Input_directory,'Carbon', 'MCTI_ref.xlsx'), sheet_name='OtherBiome')
            merged_df = pd.merge(filtered_OtherBiome, ref, how='left', left_on='BiomePhyto', right_on='BiomePhyto_ID') #We associate the cell of the other biome to their reference value
            merged_df.drop(columns=['BiomePhyto_ID'], inplace=True)
            
            
            # Reference calue for the biome Cerrado (depend on the state)
            for i, state in enumerate([5, 7, 9, 10, 11, 12, 13, 17, 18, 26, 27]): # State in the cerrado
                ref = pd.read_excel(os.path.join(Input_directory,'Carbon', 'MCTI_ref.xlsx'), sheet_name=f'Cerrado{state}')
                filtered_State = merged[merged['State'].isin([state])]
                filtered_State = filtered_State[(filtered_State['BiomePhyto'] > 200) & (filtered_State['BiomePhyto'] < 300)]
                merged_Cerrado = pd.merge(filtered_State, ref, how='left', left_on='BiomePhyto', right_on='BiomePhyto_ID')
                merged_Cerrado.drop(columns=['BiomePhyto_ID'], inplace=True)
                merged_df = pd.concat([merged_df, merged_Cerrado], ignore_index=True)
        
        
            merged_df.drop(columns=['State'], inplace=True)
            
            #Calculation of the reference value in each cell
            #Create a new column that is equal to the weight SOCref (weight = proportion of the area covered by this ClimSoil type)
            merged_df['AGBref'] = merged_df['%'] * merged_df['AGB']
            summed_AGB = merged_df.groupby('ID')['AGBref'].sum().reset_index()
            AGBmap = reclass_dataframe(summed_AGB,ID_raster,0,1)
        
            merged_df['BGBref'] = merged_df['%'] * merged_df['BGB']
            summed_BGB = merged_df.groupby('ID')['BGBref'].sum().reset_index()
            BGBmap = reclass_dataframe(summed_BGB,ID_raster,0,1)
        
   
    else:
        
        ref= pd.read_excel(os.path.join(Input_directory,'Carbon','IPCC_ref.xlsx'), sheet_name=LU)
        AGBmap= reclass_dataframe(ref,Clim_raster,0,1) 
        BGBmap= reclass_dataframe(ref,Clim_raster,0,2) 
                
    
    Biomass = AGBmap+BGBmap
    
    return Biomass 



def get_C() :
    
    #Input data for mapping
    ID = rasterToArray(os.path.join(Input_directory,'ID_cell_Albers.tif')) 
    mask=rasterToArray(os.path.join(Input_directory,'Mask_Albers.tif'))
    
    #Input data SOC
    ClimSoil = pd.read_excel(os.path.join(Input_directory,'Carbon','ID_ClimSoil_Subordem.xlsx')) # Proportion of each ClimSoil type in each Cell
    Biome = 100*(rasterToArray(os.path.join(Input_directory,'Biome_Albers.tif')))
    State = rasterToArray( os.path.join(Input_directory,'State_Albers.tif')) 
    StateBiome = get_raster_combi(State,Biome)
    
    #Input data biomass
    State_ID = pd.read_excel(os.path.join(Input_directory,'Carbon', 'ID_State.xlsx'))
    BiomePhyto_per = pd.read_excel(os.path.join(Input_directory,'Carbon','Id_BiomePhyto.xlsx'))
    Clim_raster = rasterToArray(os.path.join(Input_directory,'Carbon','Climate_IPCC.tif'))
    
    # Calculation of referecence stock maps
    Ff= get_Fc('Forest', StateBiome)
    SOC_Forest =Ff*get_SOC_ref('Forest',ID,ClimSoil)
    Bio_Forest = get_biomass_ref('Forest', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Forest = SOC_Forest +  Bio_Forest
    plot_map_with_mask(C_Forest,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Forest')
    arrayToMap(C_Forest, os.path.join(Output_directory,'Carbon','Cref_Forest'),mask,Input_directory)
    
    Fg= get_Fc('Grassland', StateBiome)
    SOC_Grassland = Fg*get_SOC_ref('Grassland',ID,ClimSoil)
    Bio_Grassland = get_biomass_ref('Grassland', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Grassland = SOC_Grassland+Bio_Grassland
    plot_map_with_mask(C_Grassland,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Grassland')
    arrayToMap(C_Grassland, os.path.join(Output_directory,'Carbon','Cref_Grassland'),mask,Input_directory)
    
    Fp= get_Fc('Pasture', StateBiome)
    SOC_Pasture = Fp*get_SOC_ref('Pasture',ID,ClimSoil)
    Bio_Pasture = get_biomass_ref('Pasture', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Pasture = SOC_Pasture + Bio_Pasture
    plot_map_with_mask(C_Pasture,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Pasture')
    arrayToMap(C_Pasture, os.path.join(Output_directory,'Carbon','Cref_Pasture'),mask,Input_directory)
    
    Fpe= get_Fc('Perennial', StateBiome)
    SOC_Perennial = Fpe* get_SOC_ref('Perennial',ID,ClimSoil)
    Bio_Perennial = get_biomass_ref('Perennial', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Perennial= SOC_Perennial+Bio_Perennial
    plot_map_with_mask(C_Perennial,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Perennial')
    arrayToMap(C_Perennial, os.path.join(Output_directory,'Carbon','Cref_Perennial'),mask,Input_directory)
   
    Fs= get_Fc('Semiperennial', StateBiome)
    SOC_Semiperennial = Fs* get_SOC_ref('Semiperennial',ID,ClimSoil)
    Bio_Semiperennial = get_biomass_ref('Semiperennial', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Semiperennial = SOC_Semiperennial + Bio_Semiperennial
    plot_map_with_mask(C_Semiperennial,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Semiperennial')
    arrayToMap(C_Semiperennial, os.path.join(Output_directory,'Carbon','Cref_Semiperennial'),mask,Input_directory)
    
    Ft= get_Fc('Temporary',StateBiome)
    SOC_Annual = Ft * get_SOC_ref('Temporary',ID,ClimSoil)
    Bio_Annual = get_biomass_ref('Temporary', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Annual = SOC_Annual + Bio_Annual
    plot_map_with_mask(C_Annual,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Annual')
    arrayToMap(C_Annual, os.path.join(Output_directory,'Carbon','Cref_Annual'),mask,Input_directory)
    
    Ffo= get_Fc('Forestry', StateBiome)
    SOC_Forestry = Ffo*get_SOC_ref('Forestry',ID,ClimSoil)
    Bio_Forestry = get_biomass_ref('Forestry', ID, State_ID, BiomePhyto_per,Clim_raster)
    C_Forestry = SOC_Forestry + Bio_Forestry
    plot_map_with_mask(C_Forestry,0,200,colors5,"C Stock (t/ha)", 'Reference C Stock Forestry')
    arrayToMap(C_Forestry, os.path.join(Output_directory,'Carbon','Cref_Forestry'),mask,Input_directory)


    return C_Forest, C_Grassland, C_Pasture, C_Perennial, C_Annual, C_Semiperennial, C_Forestry # return le SOC of that year and sc


get_C()


# --------------------------------------------------
# Revenue
# --------------------------------------------------


def get_Revenue_Ref():
    
    
    #Input data
    ID = rasterToArray(os.path.join(Input_directory,'ID_cell_Albers.tif'))
    mask=rasterToArray(os.path.join(Input_directory,'Mask_Albers.tif'))
    State = rasterToArray( os.path.join(Input_directory,'State_Albers.tif')) 
    
    
    for Sc in scenarios:
        
        for year, year_id in zip(years, years_id):
            # --------------------------------------------------
            # Yield
            # --------------------------------------------------
            
            #Rainfed Yield in 2050 (IMAGE)
            Managment_Factor = pd.read_excel(os.path.join(Input_directory,'Economy', 'Managment_Factor.xlsx'), sheet_name=f'{Sc}')
            Feed_Efficiency = pd.read_excel(os.path.join(Input_directory,'Economy', 'Feed_Efficiency.xlsx'), sheet_name=f'{Sc}')
            Intensity_Factor = pd.read_excel(os.path.join(Input_directory,'Economy', 'Intensity_Factor.xlsx'), sheet_name=f'{Sc}')
            
            
            Grass_yield = (rasterToArray(os.path.join(Input_directory,'Economy', f'{Sc}', f'{year}_Grass.tif'))*Intensity_Factor.iloc[year_id, 1])/100
            Bovine_Yield = Grass_yield/Feed_Efficiency.iloc[year_id,1] #ton of non dairy cattle per ha/year
            Soybean_yield = (rasterToArray(os.path.join(Input_directory,'Economy',f'{Sc}',  f'{year}_Soybean.tif'))*Managment_Factor.iloc[year_id, 1]*Intensity_Factor.iloc[year_id, 2])/100
            Sugar_yield = (rasterToArray(os.path.join(Input_directory,'Economy', f'{Sc}', f'{year}_Sugarcane.tif'))*Managment_Factor.iloc[year_id, 2]*Intensity_Factor.iloc[year_id, 3])/100
            Coffee_yield = (rasterToArray(os.path.join(Input_directory,'Economy', f'{Sc}', f'{year}_OtherNonFood.tif'))*Managment_Factor.iloc[year_id, 3]*Intensity_Factor.iloc[year_id, 4])/100
            
        
            #plot_map_with_mask(Bovine_Yield, 0, 1, colors5, f"ton/ha/year", f'Yield livestock')
            #plot_map_with_mask(Soybean_yield, 0, 6, colors5, f"ton/ha/year", f'Yield Soybeans')
            #plot_map_with_mask(Coffee_yield, 0, 5, colors5, f"ton/ha/year", f'Yield Coffee')
            #plot_map_with_mask( Sugar_yield, 0, 140, colors5, f"ton/ha/year", f'Yield Sugarcane')
            
            
            # --------------------------------------------------
            # Price
            # --------------------------------------------------
            
            #Price in 2020 (CONAB)
            Price = pd.read_excel(os.path.join(Input_directory,'Economy', 'PRICE.xlsx'), sheet_name=f'PRICE_CONAB')
            
            Soybean_Price_map = reclass_dataframe(Price, State, 1, 2)
            Sugar_Price_map = reclass_dataframe(Price, State, 1, 3)
            Coffee_Price_map = reclass_dataframe(Price, State, 1, 4)
            Bovine_Price_map = reclass_dataframe(Price, State, 1, 5)
            
            #plot_map_with_mask(Bovine_Price_map, 0, 5000, colors5, f"USD/ton", f'Price Bovine')
            #plot_map_with_mask(Soybean_Price_map, 0, 400, colors5, f"USD/ton", f'Price Soybeans')
            #plot_map_with_mask(Coffee_Price_map, 0, 2500, colors5, f"USD/ton", f'Price Coffee')
            #plot_map_with_mask(Sugar_Price_map, 0, 30, colors5, f"USD/ton", f'Price Sugarcane')
            
            
            #Price index (IMAGE)
            Price_index = pd.read_excel(os.path.join(Input_directory,'Economy', 'Price_Index.xlsx'), sheet_name=f'{Sc}')
            
            Soybean_Price_index = Price_index.iloc[year_id, 1]
            Sugar_Price_index = Price_index.iloc[year_id, 2]
            Coffee_Price_index = Price_index.iloc[year_id, 3]
            Livestock_Price_index = Price_index.iloc[year_id, 4]
            
            # --------------------------------------------------
            # Revenue
            # --------------------------------------------------
            
            #Revenue reference map
            Rev_temporary = Soybean_yield * Soybean_Price_map * Soybean_Price_index
            Rev_Perennial = Coffee_yield * Coffee_Price_map * Coffee_Price_index 
            Rev_SemiPerennial = Sugar_yield * Sugar_Price_map * Sugar_Price_index
            Rev_Pasture = Bovine_Yield * Bovine_Price_map * Livestock_Price_index
    
            
            plot_map_with_mask(Rev_temporary, 0, 4000, colors5, "USD/ha", 'Ref Revenue Temporary cropland')
            plot_map_with_mask(Rev_Perennial, 0, 4000, colors5, "USD/ha", 'Ref Revenue Perennial cropland')
            plot_map_with_mask(Rev_SemiPerennial, 0, 4000, colors5, "USD/ha", 'Ref Revenue Semiperennial cropland')
            plot_map_with_mask( Rev_Pasture, 0, 500, colors5, "USD/ha", 'Ref Revenue Planted pasture')
            
            arrayToMap(Rev_temporary, os.path.join(Output_directory, 'Economy',f'{Sc}', f'{year}_Eref_Annual'),mask,Input_directory)
            arrayToMap(Rev_Perennial, os.path.join(Output_directory, 'Economy',f'{Sc}', f'{year}_Eref_Perennial'),mask,Input_directory)
            arrayToMap(Rev_SemiPerennial, os.path.join(Output_directory, 'Economy',f'{Sc}', f'{year}_Eref_Semiperennial'),mask,Input_directory)
            arrayToMap(Rev_Pasture, os.path.join(Output_directory, 'Economy',f'{Sc}', f'{year}_Eref_Pasture'),mask,Input_directory)
     
get_Revenue_Ref()      




# --------------------------------------------------
# Revenue for 2020
# --------------------------------------------------



def get_Revenue_Observed_Ref():
    Sc='SSP2'
    years_id=[0,1] 
    years =['2015', '2020']
    
    for year, year_id in zip(years, years_id):
        #Input data
        ID = rasterToArray(os.path.join(Input_directory,'ID_cell_Albers.tif'))
        mask=rasterToArray(os.path.join(Input_directory,'Mask_Albers.tif'))
        State = rasterToArray( os.path.join(Input_directory,'State_Albers.tif')) 
        
        
        #Rainfed Yield in 2050 (IMAGE)
        Managment_Factor = pd.read_excel(os.path.join(Input_directory,'Economy', 'Managment_Factor.xlsx'), sheet_name=f'{Sc}')
        Feed_Efficiency = pd.read_excel(os.path.join(Input_directory,'Economy', 'Feed_Efficiency.xlsx'), sheet_name=f'{Sc}')
        Intensity_Factor = pd.read_excel(os.path.join(Input_directory,'Economy', 'Intensity_Factor.xlsx'), sheet_name=f'{Sc}')
        
        
        Grass_yield = (rasterToArray(os.path.join(Input_directory,'Economy', f'{Sc}', f'{year}_Grass.tif'))*Intensity_Factor.iloc[year_id, 1])/100
        Bovine_Yield=Grass_yield/Feed_Efficiency.iloc[year_id,1] #ton of non dairy cattle per ha/year
        Soybean_yield = (rasterToArray(os.path.join(Input_directory,'Economy',f'{Sc}',  f'{year}_Soybean.tif'))*Managment_Factor.iloc[year_id, 1]*Intensity_Factor.iloc[year_id, 2])/100
        Sugar_yield = (rasterToArray(os.path.join(Input_directory,'Economy', f'{Sc}', f'{year}_Sugarcane.tif'))*Managment_Factor.iloc[year_id, 2]*Intensity_Factor.iloc[year_id, 3])/100
        Coffee_yield = (rasterToArray(os.path.join(Input_directory,'Economy', f'{Sc}', f'{year}_OtherNonFood.tif'))*Managment_Factor.iloc[year_id, 3]*Intensity_Factor.iloc[year_id, 4])/100
        
        
        #plot_map_with_mask(Bovine_Yield, 0, 1, colors5, f"ton/ha/year", f'Yield livestock')
        #plot_map_with_mask(Soybean_yield, 0, 6, colors5, f"ton/ha/year", f'Yield Soybeans')
        #plot_map_with_mask(Coffee_yield, 0, 5, colors5, f"ton/ha/year", f'Yield Coffee')
        #plot_map_with_mask( Sugar_yield, 0, 140, colors5, f"ton/ha/year", f'Yield Sugarcane')
        
        
        # --------------------------------------------------
        # Price
        # --------------------------------------------------
        
        if year_id==0:
            Price = pd.read_excel(os.path.join(Input_directory,'Economy', 'PRICE.xlsx'), sheet_name=f'PRICE_CONAB')
        else: 
            Price = pd.read_excel(os.path.join(Input_directory,'Economy', 'PRICE_2020.xlsx'), sheet_name=f'PRICE_CONAB')
        
        Soybean_Price_map = reclass_dataframe(Price, State, 1, 2)
        Sugar_Price_map = reclass_dataframe(Price, State, 1, 3)
        Coffee_Price_map = reclass_dataframe(Price, State, 1, 4)
        Bovine_Price_map = reclass_dataframe(Price, State, 1, 5)
        
        #plot_map_with_mask(Bovine_Price_map, 0, 5000, colors5, f"USD/ton", f'Price Bovine')
        #plot_map_with_mask(Soybean_Price_map, 0, 400, colors5, f"USD/ton", f'Price Soybeans')
        #plot_map_with_mask(Coffee_Price_map, 0, 2500, colors5, f"USD/ton", f'Price Coffee')
        #plot_map_with_mask(Sugar_Price_map, 0, 30, colors5, f"USD/ton", f'Price Sugarcane')
        
        
        #Revenue reference map
        Rev_temporary = Soybean_yield * Soybean_Price_map
        Rev_Perennial = Coffee_yield * Coffee_Price_map 
        Rev_SemiPerennial = Sugar_yield * Sugar_Price_map 
        Rev_Pasture = Bovine_Yield * Bovine_Price_map 
        
        
        plot_map_with_mask(Rev_temporary, 0, 4000, colors5, "USD/ha", 'Ref Revenue Temporary cropland')
        plot_map_with_mask(Rev_Perennial, 0, 4000, colors5, "USD/ha", 'Ref Revenue Perennial cropland')
        plot_map_with_mask(Rev_SemiPerennial, 0, 4000, colors5, "USD/ha", 'Ref Revenue Semiperennial cropland')
        plot_map_with_mask( Rev_Pasture, 0, 500, colors5, "USD/ha", 'Ref Revenue Planted pasture')
        
        arrayToMap(Rev_temporary, os.path.join(Output_directory, 'Economy','Observed', f'{year}_Eref_Annual'),mask,Input_directory)
        arrayToMap(Rev_Perennial, os.path.join(Output_directory, 'Economy','Observed', f'{year}_Eref_Perennial'),mask,Input_directory)
        arrayToMap(Rev_SemiPerennial, os.path.join(Output_directory, 'Economy','Observed', f'{year}_Eref_Semiperennial'),mask,Input_directory)
        arrayToMap(Rev_Pasture, os.path.join(Output_directory, 'Economy','Observed', f'{year}_Eref_Pasture'),mask,Input_directory)

get_Revenue_Observed_Ref()
