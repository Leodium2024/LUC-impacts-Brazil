import pandas as pd
from osgeo import gdal
from osgeo import osr
import glob
import os
import math
import numpy as np
import pandas as pd
from scipy.ndimage import zoom
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from os.path import join as opj
import warnings


########################################################################################
########################################################################################
#############################       Import/Export      #################################
########################################################################################
########################################################################################


def rasterToArray(inpRaster):
    """Return a numpy array based on a raster file. The 'No data' of the raster 
    are converted to zero"""
    inpRaster = gdal.Open(inpRaster)
    unqBand = inpRaster.GetRasterBand(1)
    noData = unqBand.GetNoDataValue()
    npArr = inpRaster.ReadAsArray()
    npArr[npArr == noData] = 0  # Converting noData values to zero
    if npArr[npArr < 0].any() == True:  # Making sure no neg value is in array
        npArr[npArr < 0] = 0
    return npArr

def setInputPath():
    # Define the directory where the results will be stored
    inputDir = "E:/Brazil/Model/Input"
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)
    return inputDir

def setOutputPath():
    """Define the directory where the results will be stored"""
    resultsDir = "E:/Brazil/Model/Output"
    if not os.path.exists(resultsDir):
        os.makedirs(resultsDir)
    return resultsDir

def referenceRaster():
    """Set the reference raster in which its geoparameters will be used to 
    create/save new rasters/numpy arrays."""
    refRaster = os.path.join(inputDir,'Mask_Albers.tif')
    return refRaster

def Get_Year_NetCDF (Scenario,lu,year):
    """ Upload of the LU layer for a specific LU and year
    Scenario is the path leading to the .nc file 
    lu is the name of the land use type to be extracted (string)
    Year to be extracted (int)
    The land use "pastp" need to be caluclated based on the other land use
    """
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
        data =data*mask
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


def arrayToMap(inpArray, outRaster):
    """Convert an array to a raster map and save to disk (tiff format). 
    Arguments: inpArray = array to be converted; outRaster = string with the 
    filename representing the raster map"""
    # Getting geoparameters of the original raster
    refRaster = gdal.Open(referenceRaster())
    geotransform = refRaster.GetGeoTransform()
    originX = geotransform[0]  # top left-x
    originY = geotransform[3]  # top left-y
    pixelWidth = geotransform[1]  # w-e pixel resolution
    pixelHeight = geotransform[5]  # n-s pixel resolution
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    inpArray=inpArray*mask
    # Creating new raster, then adding geoparameters
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(
        opj(resultsDir, outRaster + ".tiff"), cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(inpArray)
    outband.SetNoDataValue(0)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(refRaster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    

########################################################################################
########################################################################################
#############################       Tool Functions      ################################
########################################################################################
########################################################################################
    
def reclass_dataframe(dataframe, raster, col_with_ID, Col_with_Value):
    """Reclass an array based on the values in the DataFrame."""
    # Extract code-value pairs from the DataFrame
    code = dataframe.iloc[:, col_with_ID].values
    value = 10000*dataframe.iloc[:, Col_with_Value].values  # Multiply by 10000 to make sure the value are not rounded while int32 tranformation
    
    # Lookup table with all the codes of npRaster. '+1' encompasses all codes
    max_code = np.max(code)
    lut = np.zeros(max(max_code + 1, raster.max() + 1), dtype=np.float32)
    
    # Assign a value for the lookup table based on the associated code
    lut[code] = value
    
    # Reclassify npRaster filled w/ codes according to lookup table's values
    raster_reclassif = (lut[raster]).astype(np.int32)
    raster_reclassif =  raster_reclassif/10000 # Divided by 10000 to go back to the same order of magnitude
    
    return raster_reclassif

def get_raster_combi(raster1, raster2):
 '''Function that combine two categorical array into one
    The two input arrays have to be reclassify  beforehand'''

 #getting climate + soil combination
 combi = raster1 + raster2
 # Creating mask where soil == 0 (meaning that the cell represents water)
 #mask = np.ma.masked_where(raster2 == 0,raster2)
 # Replacing the output cells with the mask 
 #combi[mask.mask] = 0 
 return combi

def get_agri_weight(State_raster):
    Weight = pd.read_excel(os.path.join(inputDir,'Percentage_agri_mos.xlsx'),sheet_name='Agriculture')
    Wpe = reclass_dataframe(Weight, State_raster, 0 , 2)
    Ws = reclass_dataframe(Weight, State_raster, 0 , 3)
    Wt = reclass_dataframe(Weight, State_raster, 0 , 1)
    return Wpe,Ws,Wt

def get_lu(sc, year):
    
    State = rasterToArray( os.path.join(inputDir,'State_Albers.tif'))
    Agriculture = Get_Year_NetCDF(os.path.join(inputDir,sc),'agric',year)
    Forest = Get_Year_NetCDF(os.path.join(inputDir,sc),'veg',year)
    Forestry = Get_Year_NetCDF(os.path.join(inputDir,sc),'fores',year)
    Other = Get_Year_NetCDF(os.path.join(inputDir,sc),'others',year)
    Mosaic = Get_Year_NetCDF(os.path.join(inputDir,sc),'mosc',year)
    Pasture = Get_Year_NetCDF(os.path.join(inputDir,sc),'pastp',year)
    Grassland = Get_Year_NetCDF(os.path.join(inputDir,sc),'gveg',year)
    
    #Split Mosaic 
    Wf, Wg, Wp, Wa, Wfo, Wot = Split_Mosaic (State)
    
    Forest = Forest + (Mosaic*Wf)
    Pasture = Pasture + (Mosaic*Wp)
    Forestry = Forestry + (Mosaic*Wfo)
    Agriculture = Agriculture + (Mosaic*Wa)
    Grassland = Grassland + (Mosaic*Wg)
    Other = Other + (Mosaic*Wot)
    
    return Forest, Pasture, Forestry, Agriculture, Grassland, Other    

def Split_Mosaic (State_raster):
    ''' This function will return a array containing the proportion in each cell of 
        Perenial, Temporary, semiperenial agriculture and 
        Mosaic Agri-Forest, mosaic  Agri-Grassland
        All those proportion are state-specific'''
    Weight = pd.read_excel(os.path.join(inputDir,'Percentage_agri_mos.xlsx'),sheet_name='Mosaic')
    Wf = reclass_dataframe(Weight, State_raster,0,  1)
    Wg = reclass_dataframe(Weight, State_raster,0,  2)
    Wp = reclass_dataframe(Weight, State_raster,0,  3)
    Wa = reclass_dataframe(Weight, State_raster,0 , 4)
    Wfo = reclass_dataframe(Weight, State_raster,0 , 5)
    Wot = reclass_dataframe(Weight, State_raster,0 , 6)

    return  Wf, Wg, Wp, Wa, Wfo, Wot
    
    

########################################################################################
########################################################################################
#############################       Carbon Functions      ##############################
########################################################################################
########################################################################################


def get_Fc (LU, StateBiome_Map):
     
     ''' This function will return a array containing the Fc in each cell
         The Fc value are State and biome specific'''
         
     Fc_LU = pd.read_excel(os.path.join(inputDir,'Carbon','Fc.xlsx'),sheet_name=LU)
     F = reclass_dataframe(Fc_LU, StateBiome_Map, 0 , 1)
     return F
  
def get_SOC_ref (LU,ID_raster,SoilClim_Per):
    
    ''' This function will return an array containing the SOCref for a specific LU.
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
    SOCref = pd.read_excel(os.path.join(inputDir,'Carbon','SOCref_Subordem.xlsx'),sheet_name=LU) 
    
    # Merge the two DataFrames based on the condit
    merged_df = pd.merge(SoilClim_Per, SOCref, how='left', left_on='ClimSoil', right_on='ClimSoil')
    
    #Create a new column that is equal to the weight SOCref (weight = proportion of the area covered by this ClimSoil type)
    merged_df['SOC'] = merged_df['%'] * merged_df['Average']
    
    # Give a table that is the soil ref per cell (sum of all the weight SOCref)
    summed_table = merged_df.groupby('ID')['SOC'].sum().reset_index()
    
    output = reclass_dataframe(summed_table,ID_raster,0,1)

    return output        

def get_SOC(sc, year, ClimSoil, StateBiome, ID, State) :

    ''' This function will return a array containing the SOC for a specific year
        and scenario'''
    
    Forest, Pasture, Forestry, Agriculture, Grassland, Other = get_lu(sc, year)
    
    Fpe= get_Fc('Perennial', StateBiome)
    Fs= get_Fc('Semiperennial', StateBiome)
    Ft= get_Fc('Temporary',StateBiome)
    Ff= get_Fc('Forest', StateBiome)
    Fg= get_Fc('Grassland', StateBiome)
    Fp= get_Fc('Pasture', StateBiome)
    Ffo= get_Fc('Forestry', StateBiome)
    
    Wpe, Ws, Wt = get_agri_weight (State)
    SOC_Perenial = Fpe* Wpe * get_SOC_ref('Perennial',ID,ClimSoil)
    SOC_Semiperenial = Fs* Ws * get_SOC_ref('Semiperennial',ID,ClimSoil)
    SOC_Temporary = Ft * Wt * get_SOC_ref('Temporary',ID,ClimSoil)
    SOC_Agri=(SOC_Perenial+SOC_Semiperenial+SOC_Temporary)*Agriculture
    
    SOC_Pasture= Fp*get_SOC_ref('Pasture',ID,ClimSoil)*Pasture
    SOC_Forest=Ff*get_SOC_ref('Forest',ID,ClimSoil)*Forest
    SOC_Grassland= Fg*get_SOC_ref('Grassland',ID,ClimSoil)*Grassland
    SOC_Forestry= Ffo*get_SOC_ref('Forestry',ID,ClimSoil)*Forestry
    
    SOC_Stock=SOC_Agri+SOC_Forest+SOC_Forestry+SOC_Grassland+SOC_Pasture
    
    return SOC_Stock 

def get_biomass_ref(LU, ID_raster, State_ID, BiomePhyto_per, Clim_raster):
    ''' This function will return an array containing the biomass ref for a specific LU.
        For Grassland and forest, we used the MCTI ref values varying with the biome, phytophisiome and sometimes the the federal state.
        For other land use, we took IPPC reference values
    '''
    if LU in {'Forest','Grassland'}:
        
            merged = pd.merge(BiomePhyto_per, State_ID , how='left', left_on='ID', right_on='ID_cell')
            merged.drop(columns=['ID_cell'], inplace=True)
            
            
            # Reference value for the biome of Amazone, Caatinga, Pampa, Mata and Pantanal
            filtered_OtherBiome = merged[(merged['BiomePhyto'] < 200) | (merged['BiomePhyto'] > 300)] #As MCTI has state specific value for Cerrado, we filter to keep all the other biome
            ref= pd.read_excel(os.path.join(inputDir,'Carbon', 'MCTI_ref.xlsx'), sheet_name='OtherBiome')
            merged_df = pd.merge(filtered_OtherBiome, ref, how='left', left_on='BiomePhyto', right_on='BiomePhyto_ID') #We associate the cell of the other biome to their reference value
            merged_df.drop(columns=['BiomePhyto_ID'], inplace=True)
            
            
            # Reference calue for the biome Cerrado (depend on the state)
            for i, state in enumerate([5, 7, 9, 10, 11, 12, 13, 17, 18, 26, 27]): # State in the cerrado
                ref = pd.read_excel(os.path.join(inputDir,'Carbon', 'MCTI_ref.xlsx'), sheet_name=f'Cerrado{state}')
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
        
        ref= pd.read_excel(os.path.join(inputDir,'Carbon','IPCC_ref.xlsx'), sheet_name=LU)
        AGBmap= reclass_dataframe(ref,Clim_raster,0,1) 
        BGBmap= reclass_dataframe(ref,Clim_raster,0,2) 
        
    return AGBmap,BGBmap 

def get_biomass(sc, year, Clim_raster, BiomePhyto_per, State_ID, ID, State):
        
    ''' This function will return a array containing the Biomass carbon for a specific year
        and scenario'''
    
        Forest, Pasture, Forestry, Agriculture, Grassland, Other = get_lu(sc, year)
        
        #Agriculture
        Wpe, Ws, Wt = get_agri_weight (State)
        AGB_pere,BGB_pere =  Wpe * get_biomass_ref('Perennial', ID, State_ID, BiomePhyto_per,Clim_raster)
        AGB_semi,BGB_semi = Ws * get_biomass_ref('Semiperennial', ID, State_ID, BiomePhyto_per,Clim_raster)
        AGB_temp,BGB_temp =  Wt * get_biomass_ref('Temporary', ID, State_ID, BiomePhyto_per,Clim_raster)
        Bio_Agri=(AGB_pere+BGB_pere+AGB_semi+BGB_semi+AGB_temp+BGB_temp)*Agriculture
        
        #Other land use
        AGB_Forest,BGB_Forest= get_biomass_ref('Forest', ID, State_ID, BiomePhyto_per,Clim_raster)
        Bio_Forest = (AGB_Forest+BGB_Forest)*Forest
        AGB_Grassland,BGB_Grassland= get_biomass_ref('Grassland', ID, State_ID, BiomePhyto_per,Clim_raster)
        Bio_Grassland= (AGB_Grassland+BGB_Grassland)*Grassland
        AGB_Pasture,BGB_Pasture= get_biomass_ref('Pasture', ID, State_ID, BiomePhyto_per,Clim_raster)
        Bio_Pasture = (AGB_Pasture+BGB_Pasture)*Pasture
        AGB_Forestry,BGB_Forestry= get_biomass_ref('Forestry', ID, State_ID, BiomePhyto_per,Clim_raster)
        Bio_Forestry = (AGB_Forestry+BGB_Forestry)*Forestry
       
        Bio_stock = Bio_Agri + Bio_Forest + Bio_Forestry + Bio_Grassland + Bio_Pasture
       
        return Bio_stock    

def Calculate_Carbon(sc,Years):

    ''' This function will return dictionaries containing arrays of the total Carbon stock, Biomass carbon stock
        SOC stock across the land use projections. 
        It also a dataframe with the total Carbon stock, biomass carbon and SOC for the study area '''
    
    ID = rasterToArray(os.path.join(inputDir,'Carbon','ID_cell_Albers.tif')) #Table with the ID of each cell
    State = rasterToArray( os.path.join(inputDir,'State_Albers.tif')) 
    #Input data SOC
    ClimSoil = pd.read_excel(os.path.join(inputDir,'Carbon','ID_ClimSoil_Subordem.xlsx')) # Proportion of each ClimSoil type in each Cell
    Biome = 100*(rasterToArray(os.path.join(inputDir,'Carbon','Biome_Albers.tif')))
    StateBiome = get_raster_combi(State,Biome)
    #Input data biomass
    State_ID = pd.read_excel(os.path.join(inputDir,'Carbon', 'ID_State.xlsx'))
    BiomePhyto_per = pd.read_excel(os.path.join(inputDir,'Carbon','Id_BiomePhyto.xlsx'))
    Clim_raster = rasterToArray(os.path.join(inputDir,'Carbon','Climate_IPCC.tif'))
    
    #Output data
    TotalStockSOC = {}
    TotalStockBio = {}
    TotalStockCarbon={}
    SOC_total = []
    Bio_total = []
    Carbon_total = []
    
    for year in Years:
        TotalStockSOC[year] = get_SOC(sc, year, ClimSoil, StateBiome, ID, State)
        TotalStockBio[year] = get_biomass(sc, year, Clim_raster, BiomePhyto_per, State_ID, ID, State)
        TotalStockCarbon[year] = TotalStockSOC[year] + TotalStockBio[year]
        SOC = np.sum(10000 * TotalStockSOC[year]) / 1000000
        SOC_total.append(SOC)
        Bio = np.sum(10000 * TotalStockBio[year]) / 1000000
        Bio_total.append(Bio)
        Carbon=SOC+Bio
        Carbon_total.append(Carbon)
        
    Evolution_Carbon= pd.DataFrame([SOC_total, Bio_total, Carbon_total],
                                   index=['SOC', 'Bio', 'Total'])    
    
    return TotalStockSOC,TotalStockBio,TotalStockCarbon,Evolution_Carbon
    
    
        
########################################################################################
########################################################################################
#############################       Biodiversity Functions      ########################
########################################################################################
########################################################################################

def get_perfomant_ottid(Performance):
    ''' This function select the SDMs to be included and return the corresponding ott_id '''
    
    Performance=Performance[(Performance['N_Pre'] >=10)]  #More than 10 records
    Performance=Performance[(Performance['NNI'] >=0.5)]   #Not clustered records
    Performance=Performance[(Performance['AUC'] >=0.8)]   #AUC >= 0.8
    Performance=Performance[(Performance['TSS'] >=0.6)]   #TSS >= 0.6
    Ott_id =Performance['Ott_id'].unique()
    return Ott_id

def get_distribution_map(Species_ID, Coeficient, Performance, land_use, env):
    
    ''' This function dtermine the distribution area of a species based on the
        the logistic regression coeficient, land use and envrionmental variables'''
    
    #Upload of the range of the species
    Species_range_directory = os.path.join(inputDir,'Biodiversity','Species_Range')
    Species_range = rasterToArray(os.path.join(Species_range_directory, f'ott_id_{Species_ID}.tif'))
    
    
    #Upload of the regression coeficient of the species
    Coeficient_ID = Coeficient[Coeficient['Ott_id'] == Species_ID]
    
    #Upload of the threshold of the species
    Performance_ID = Performance[Performance['Ott_id'] == Species_ID]
    Threshold = Performance_ID.iloc[0,12]
    
    ######Logistic Model######
    #1) Dataframe to store the result of the logistic model
    Total_fitted = pd.DataFrame(0, index=range(464), columns=range(544))
    #2) Calculation of the polynome of the logistic model
    Forest_poly= Coeficient_ID.iloc[0, 1]*land_use["Forest"]
    Pasture_poly= Coeficient_ID.iloc[0, 2]*land_use["Pasture"]
    Grassland_poly= Coeficient_ID.iloc[0, 3]*land_use["Grassland"]
    Agriculture_poly= Coeficient_ID.iloc[0, 4]*land_use["Agriculture"]
    Forestry_poly= Coeficient_ID.iloc[0, 5]*land_use["Forestry"]
    Bio_1_poly= Coeficient_ID.iloc[0, 6]*env["Bio_1"]
    Bio_4_poly= Coeficient_ID.iloc[0, 7]*env["Bio_4"]
    Bio_12_poly= Coeficient_ID.iloc[0, 8]*env["Bio_12"]
    Bio_15_poly= Coeficient_ID.iloc[0, 9]*env["Bio_15"]
    Elevation_poly= Coeficient_ID.iloc[0, 10]*env["Elevation"]
    Slope_poly= Coeficient_ID.iloc[0, 11]*env["Slope"]
    Intercept= Coeficient_ID.iloc[0, 12]

    #3) Calculate the logistic model
    warnings.filterwarnings('ignore')
    exponent = Forest_poly + Pasture_poly + Grassland_poly + Agriculture_poly + Forestry_poly + Bio_1_poly + Bio_4_poly + Bio_12_poly + Bio_15_poly + Elevation_poly + Slope_poly + Intercept
    z = np.exp(-exponent)
    Model_fitted = 1/(1+z)
    Total_fitted=Total_fitted+ Model_fitted
      
    
    ##### Distribution Map #####
    Species_distribution = Total_fitted.applymap(lambda x: 1 if x > Threshold else 0)
    Species_distribution = Species_distribution*Species_range
    
    
    return Species_distribution

def get_N_SR_WE_CWE (sc, year, Coeficient, Performance, env, Ott_id, Order):
    
   ''' This function calculate  N,SR,WE and CWE for a specific year and scenario'''
    
    #Upload of land use variable for a specific year into a dictionary
    Forest, Pasture, Forestry, Agriculture, Grassland, Other = get_lu(sc,year)
    land_use = {'Agriculture': Agriculture,'Forest': Forest,'Forestry': Forestry,'Other': Other,'Pasture': Pasture,'Grassland': Grassland}
    
    
    # Output data
    Order_name = ['Artiodactyla', 'Carnivora', 'Chiroptera', 'Cingulata','Didelphimorphia', 'Lagomorpha', 'Primates', 'Rodentia']
    
    N = {order: 0 for order in Order_name}  # Create a dictionary of 0
    SR = {order: pd.DataFrame(0, index=range(464), columns=range(544)) for order in Order_name}
    WE = {order: pd.DataFrame(0, index=range(464), columns=range(544)) for order in Order_name}
    CWE = {order: pd.DataFrame(0, index=range(464), columns=range(544)) for order in Order_name}
    distribution_areas = {}
    
    N['All'] = 0
    SR['All'] = pd.DataFrame(0, index=range(464), columns=range(544))
    WE['All'] = pd.DataFrame(0, index=range(464), columns=range(544))
    CWE['All'] = pd.DataFrame(0, index=range(464), columns=range(544))
    
    for Species_ID in Ott_id:
        # Get species distribution map
        Species_distribution = get_distribution_map(Species_ID, Coeficient, Performance, land_use, env)
        # Check if the species has a distribution area and move to the next species if not
        DistributionArea = Species_distribution.sum().sum()
        distribution_areas[Species_ID] = DistributionArea
        
        if DistributionArea == 0:
            print(f"No distribution area for species {Species_ID}")
            continue
        
        # Update Number of species (N)
        N['All'] += 1
        
        # Update species richness (SR)
        SR['All'] += Species_distribution * mask
        
        # Update weighted endemism (WE)
        Inverse_DistributionArea = 1 / DistributionArea
        WE_species = Species_distribution * Inverse_DistributionArea * mask
        WE['All'] += WE_species
        

        # Get the order of the species
        Order_species = Order.loc[Order['Ott_id'] == Species_ID, 'Order'].values[0]

        if Order_species in N:
            SR[Order_species] += Species_distribution * mask
            WE[Order_species] += WE_species
            N[Order_species] += 1
    
    # Calculate CWE (WE/SR) 
    for order in SR:
        # Avoid division by zero
        CWE[order] = WE[order] / SR[order].replace(0, np.nan)
        CWE[order] = CWE[order].fillna(0)  # Replace NaN values with 0
    
    return N, SR, WE, CWE, distribution_areas


def get_change_distribution_Area (TotalArea,Order,Years):
    '''This function calculate the average distribution area of each mammals order
    '''
    
    #Dataframe to stock the change of distribution area for each species
    # Convert TotalArea to a DataFrame
    dataframe_area = pd.DataFrame(TotalArea)
    
    # Reset the index and update column names
    dataframe_area.reset_index(inplace=True)
    dataframe_area.columns = ['Ott_id'] + list(dataframe_area.columns[1:])
    
    # Merge with the Order DataFrame on 'Ott_id'
    dataframe_area = pd.merge(dataframe_area, Order, on='Ott_id', how='left')
    
    # Ensure the columns corresponding to 'Years' are selected
    year_columns = [dataframe_area.columns[i+1] for i in Years]
    
    # Calculate the average of the columns corresponding to the years for each order
    average = dataframe_area.groupby('Order')[year_columns].mean().reset_index()
    
    return average
                        
        
           
def Calculate_Biodiversity (sc, Years):
     ''' This function will return dictionaries containing arrays of SR, WE, CWE across time for each mammal order. 
        It also return a dictionaries with the number of species per order and time 
        The average distribution areas of each order is also caluclated for each timestep and returned in a dataframe'''
    
    #Upload of environemental variable into a dictionary
    Elevation=rasterToArray(os.path.join(inputDir,'Biodiversity','Elevation_Mean.tif'))
    Slope =rasterToArray(os.path.join(inputDir,'Biodiversity','Slope_Mean.tif'))
    biovar_1 = rasterToArray(os.path.join(inputDir,'Biodiversity','Bio_1.tif'))
    biovar_12 = rasterToArray(os.path.join(inputDir,'Biodiversity','Bio_12.tif'))
    biovar_4 = rasterToArray(os.path.join(inputDir,'Biodiversity','Bio_4.tif'))
    biovar_15 = rasterToArray(os.path.join(inputDir,'Biodiversity','Bio_15.tif'))
    env = {'Bio_1': biovar_1,'Bio_12': biovar_12,'Bio_4': biovar_4,'Bio_15': biovar_15,'Elevation': Elevation,'Slope': Slope}
    
    #Other 
    Coeficient = pd.read_excel(os.path.join(inputDir,'Biodiversity', 'Species_coef_LR_AIC_1000.xlsx'))
    Order = pd.read_excel(os.path.join(inputDir,'Biodiversity', 'Ott_order.xlsx'))
    
    #Selection of performant model
    Performance = pd.read_excel(os.path.join(inputDir,'Biodiversity', 'Performance_LR_AIC_1000.xlsx'))
    Ott_id=get_perfomant_ottid(Performance)
    
    #Output Data
    TotalN = {}
    TotalSR = {}
    TotalWE = {}
    TotalCWE = {}
    TotalArea = {}
    
    for year in Years:
        N, SR, WE, CWE, distribution_areas = get_N_SR_WE_CWE(sc, year, Coeficient, Performance, env, Ott_id, Order)
        
        TotalN[year] = N
        TotalSR[year] = SR
        TotalWE[year] = WE
        TotalCWE[year] = CWE
        TotalArea[year] = distribution_areas
    
    Evolution_Area= get_change_distribution_Area (TotalArea,Order, Years)
    
    return TotalN, TotalSR, TotalWE, TotalCWE, Evolution_Area
  
########################################################################################
########################################################################################
#############################       Revenue Functions      ########################
########################################################################################
########################################################################################

def Calculate_Revenue (sc_path,sc,Years):

    Agri_Revenue = Revenue = pd.read_excel(os.path.join(inputDir,'Agriculture', 'Revenue.xlsx'),sheet_name='Agriculture')
    Beef_Revenue = Revenue = pd.read_excel(os.path.join(inputDir,'Agriculture', 'Revenue.xlsx'),sheet_name='Pasture')
    
    #Output Data
    TotalRevenue = {}
    Revenue_From_Agri_total = []
    Revenue_From_Pasture_total = []
    Revenue_Total = []
    
    for year in Years:
        Forest, Pasture, Forestry, Agriculture, Grassland, Other = get_lu(sc_path,year)
        Revenue_From_Agri = Agriculture*Agri_Revenue.iloc[sc-1, year+1]*100
        Revenue_From_Pasture = Pasture*Beef_Revenue.iloc[sc-1, year+1]*100
        Revenue =  Revenue_From_Agri + Revenue_From_Pasture
        TotalRevenue[year] = Revenue
        
        Rev_From_Agri = np.sum(Revenue_From_Agri)
        Revenue_From_Agri_total.append(Rev_From_Agri)
        
        Rev_From_Pasture = np.sum(Revenue_From_Pasture)
        Revenue_From_Pasture_total.append(Rev_From_Pasture)
        
        R = np.sum(Revenue)
        Revenue_Total.append(R)
        
        
    Evolution_Revenue= pd.DataFrame([Revenue_From_Agri_total, Revenue_From_Pasture_total, Revenue_Total],
                                   index=['Cropland', 'Pasture', 'Total'])    
        
    return TotalRevenue,Evolution_Revenue


########################################################################################
########################################################################################
###################################       Plot      ####################################
########################################################################################
########################################################################################

def plot_map_with_mask(data_array, vmin, vmax, colors, scale_title, title):
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors=colors, N=256)
    
    plt.imshow(data_array, cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    
    # Remove the frame around the map
    plt.gca().set_frame_on(False)
    
    # Remove axis ticks and labels
    plt.xticks([])
    plt.yticks([])
    
    # Add colorbar with scale title
    cbar = plt.colorbar()
    cbar.set_label(scale_title)
    
    # Add title
    plt.title(title)
    
    plt.show()
   
def replace_with_nan(data_array, mask_array):
    data_array[mask_array == 0] = np.nan
    return data_array

    
########################################################################################
########################################################################################
###################################       Command      #################################
########################################################################################
########################################################################################

#Define scale colors for mapping
colors1 = [(1, 1, 1), (0, 1, 0)]  # White to Green
colors2 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0)]  # White to Green to Dark Green 
colors3 = [(1, 0, 0), (1, 1, 1), (0, 1, 0)] # red to white to Green
colors4 = [(1, 0, 0), (1, 1, 1), (0, 0, 1)] # red to white to Blue
colors5 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors5 = [(1, 1, 1),(0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors6 = [(1, 1, 1), (0.5, 1, 0.5), (0, 1, 0), (0, 0.5, 0), (0, 0.25, 0), (0, 0, 0)]

#Define directory
resultsDir = setOutputPath()
inputDir = setInputPath()
refRaster = referenceRaster()
scenariosDict = {}

#Configuration of the model (you must turn on/turn off, depending of which part you want to run) 

#1) Which type of impacts should be estimated? (0=no, 1=Yes)
CarbonComponent = 0
BiodiveristyComponent = 1 
AgricultureComponent = 0

#2) Which scenarios? (0=no, 1=Yes)
sc1 = 1
sc2 = 1
sc3 = 1

if sc1 == 1:
    scenariosDict[1] = os.path.join(inputDir,'SSP1_RCP19_Albers.nc')
if sc2 == 1:
    scenariosDict[2] = os.path.join(inputDir,'SSP2_RCP45_Albers.nc')
if sc3 == 1:
    scenariosDict[3] = os.path.join(inputDir,'SSP3_RCP70_Albers.nc')
    
#3) Exporting/Ploting the result? (0=no, 1=Yes)
Export_result = 1 
Plot_results = 1

#4) For all year of the projection ('All_year'), for only the first and the last ('First_and_Last_year')
Type ='All_year'

if Type == 'First_and_Last_year':
    Years=[0,7] #2015 and 2050
    Columns=['2015', '2050']
    
if Type == 'All_year':
    Years=[0,1,2,3,4,5,6,7] #2015, 3020, 2025, 2030, 2035, 2040, 2045 and 2050
    Columns=['2015', '2020', '2025', '2030', '2035', '2040','2045', '2050']
    
########################################################################################
########################################################################################
###########################      Process the code      #################################
########################################################################################
########################################################################################


# Input rasters
mask=rasterToArray(refRaster)

for sc, sc_path in scenariosDict.items():
    print(f'SC{sc}')
    if CarbonComponent == 1:
        TotalStockSOC,TotalStockBio,TotalStockCarbon,Evolution_Carbon = Calculate_Carbon(sc_path,Years)
        Diff_SOC = TotalStockSOC[7] - TotalStockSOC[0]
        Diff_Bio = TotalStockBio[7] - TotalStockBio[0]
        Diff_Carbon = TotalStockCarbon[7] - TotalStockCarbon[0]
        Evolution_Carbon.columns=Columns
        
        if Export_result == 1:
            
            arrayToMap(TotalStockSOC[0], os.path.join('Carbon', f'SC{sc}_SOC2015'))
            arrayToMap(TotalStockSOC[7], os.path.join('Carbon', f'SC{sc}_SOC2050'))
            arrayToMap(Diff_SOC, os.path.join('Carbon', f'SC{sc}_SOC_change'))
            
            arrayToMap(TotalStockBio[0], os.path.join('Carbon', f'SC{sc}_Bio2015'))
            arrayToMap(TotalStockBio[7], os.path.join('Carbon', f'SC{sc}_Bio2050'))
            arrayToMap(Diff_Bio, os.path.join('Carbon', f'SC{sc}_Bio_change'))
            
            arrayToMap(TotalStockCarbon[0], os.path.join('Carbon', f'SC{sc}_Carbon2015'))
            arrayToMap(TotalStockCarbon[7], os.path.join('Carbon', f'SC{sc}_Carbon2050'))
            arrayToMap(Diff_Carbon, os.path.join('Carbon', f'SC{sc}_Carbon_change'))
            
            Evolution_Carbon.to_excel(os.path.join(resultsDir, 'Carbon', f'SC{sc}_Evolution_Carbon.xlsx'))
            
        if Plot_results == 1:
        
            plot_map_with_mask(Diff_SOC,-5,5,colors3,"T/ha", f'Change of SOC under scenario SSP{sc}')
            plot_map_with_mask(Diff_Bio,-50,50,colors3,"T/ha", f'Change of carbon biomass under scenario SSP{sc}')
            plot_map_with_mask(Diff_Carbon,-50,50,colors3, "Ton/ha", f'Change of carbon stock under scenario SSP{sc}')
            

    if BiodiveristyComponent == 1:
        
        N,SR,WE,CWE,Evolution_Area=Calculate_Biodiversity(sc_path, Years)
        
        mammals_orders = ['All','Primates', 'Rodentia', 'Chiroptera']
        Diff_SR = {order: SR[7][order] - SR[0][order]for order in mammals_orders}
        Diff_WE = {order: WE[7][order] - WE[0][order]for order in mammals_orders}
        Diff_CWE = {order: CWE[7][order] - CWE[0][order] for order in mammals_orders}
        
        if Export_result == 1:
            for order in mammals_orders:
                # Export SR maps
                arrayToMap(SR[0][order].to_numpy(), os.path.join('Biodiversity', f'SC{sc}_SR2015_{order}'))
                arrayToMap(SR[7][order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_SR2050_{order}'))
                arrayToMap(Diff_SR[order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_SR_change_{order}'))
                
                # Export WE maps
                arrayToMap(WE[0][order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_WE2015_{order}'))
                arrayToMap(WE[7][order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_WE2050_{order}'))
                arrayToMap(Diff_WE[order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_WE_change_{order}'))
                
                # Export CWE maps
                arrayToMap(CWE[0][order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_CWE2015_{order}'))
                arrayToMap(CWE[7][order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_CWE2050_{order}'))
                arrayToMap(Diff_CWE[order].to_numpy() , os.path.join('Biodiversity', f'SC{sc}_CWE_change_{order}'))
                
                Evolution_Area.to_excel(os.path.join(resultsDir, 'Biodiversity', f'SC{sc}_Evolution_Area.xlsx'))
            
            
        if Plot_results == 1:
            
            #All
            plot_map_with_mask(Diff_SR['All'],-10,10,colors3,"Number of species", f'Change of SR under scenario SSP{sc}')
            plot_map_with_mask(Diff_WE['All'],-0.001,0.001,colors3,"Number of species", f'Change of WE under scenario SSP{sc}')
            plot_map_with_mask(Diff_CWE['All'],-0.00001,0.00001,colors3,"Number of species", f'Change of CWE under scenario SSP{sc}')
            
            # Primates
            plot_map_with_mask(Diff_SR['Primates'], -5, 5, colors3, "Number of species", f'Change of SR for Primates under scenario SSP{sc}')
            plot_map_with_mask(Diff_WE['Primates'], -0.001, 0.001, colors3, "Number of species", f'Change of WE for Primates under scenario SSP{sc}')
            plot_map_with_mask(Diff_CWE['Primates'], -0.0001, 0.0001, colors3, "Number of species", f'Change of CWE for Primates under scenario SSP{sc}')
            
            # Rodentia
            plot_map_with_mask(Diff_SR['Rodentia'], -5, 5, colors3, "Number of species", f'Change of SR for Rodentia under scenario SSP{sc}')
            plot_map_with_mask(Diff_WE['Rodentia'], -0.001, 0.001, colors3, "Number of species", f'Change of WE for Rodentia under scenario SSP{sc}')
            plot_map_with_mask(Diff_CWE['Rodentia'], -0.0001, 0.0001, colors3, "Number of species", f'Change of CWE for Rodentia under scenario SSP{sc}')
            
            # Chiroptera
            plot_map_with_mask(Diff_SR['Chiroptera'], -10, 10, colors3, "Number of species", f'Change of SR for Chiroptera under scenario SSP{sc}')
            plot_map_with_mask(Diff_WE['Chiroptera'], -0.001, 0.001, colors3, "Number of species", f'Change of WE for Chiroptera under scenario SSP{sc}')
            plot_map_with_mask(Diff_CWE['Chiroptera'], -0.00001, 0.00001, colors3, "Number of species", f'Change of CWE for Chiroptera under scenario SSP{sc}')

    if AgricultureComponent == 1:
        
        Total_revenue,Evolution_revenue=Calculate_Revenue(sc_path,sc,Years)
        Diff_revenue = Total_revenue[7] - Total_revenue[0]
        Evolution_revenue.columns=Columns
        
        
        if Export_result == 1:
           arrayToMap(Total_revenue[0], os.path.join('Agriculture', f'SC{sc}_Revenue2015'))
           arrayToMap(Total_revenue[7], os.path.join('Agriculture', f'SC{sc}_Revenue2050'))
           arrayToMap(Diff_revenue, os.path.join('Agriculture', f'SC{sc}_Revenue_change')) 
           Evolution_revenue.to_excel(os.path.join(resultsDir, 'Agriculture', f'SC{sc}_Evolution_Revenue.xlsx'))
            
        if Plot_results == 1:

            plot_map_with_mask(Diff_revenue,-10000,10000,colors3,"R$2015", f'Change of revenue under scenario SSP{sc}')




