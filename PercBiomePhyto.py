from osgeo import gdal
from osgeo import osr
import glob
import os
import math
import numpy as np
import pandas as pd
from scipy.ndimage import zoom



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




def setOutputPath():
    # Define the directory where the results will be stored
    resultsDir = "D:/Data/Results"
    if not os.path.exists(resultsDir):
        os.makedirs(resultsDir)
    return resultsDir


def arrayToMap(inpArray, outRaster,refRaster):
    """Convert an array to a raster map and save to disk (tiff format). 
    Arguments: inpArray = array to be converted; outRaster = string with the 
    filename representing the raster map"""
    # Getting geoparameters of the original raster
    
    geotransform = refRaster.GetGeoTransform()
    originX = geotransform[0]  # top left-x
    originY = geotransform[3]  # top left-y
    pixelWidth = geotransform[1]  # w-e pixel resolution
    pixelHeight = geotransform[5]  # n-s pixel resolution
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    # Generating a masked array representing cells with null value 
    # (further used when converting array to raster)
    # reference-array for masking null value
    nvArr = (refRaster.ReadAsArray()).astype(np.int32)
    # Converts to -9999 any value that doesn't represent a LU type code
    nvArr[(nvArr < 0) | (nvArr > 20)] = -9999
    nvMsk = np.ma.masked_where(nvArr != -9999, nvArr)  # null value mask
    # Updating inpArray with the null value of the MskArr
    inpArray[~nvMsk.mask] = -9999
    # Creating new raster, then adding geoparameters
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(
        os.path.join(resultsDir, outRaster + ".tif"), cols, rows, 1, gdal.GDT_Int32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(inpArray)
    outband.SetNoDataValue(-9999)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(refRaster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def reclass (txtFile,raster):
    """Reclass an array based on the value in the txtFile"""
    #This function is used in preparation of get_raster_combi
    #Soil type reclass between 1 and 15 
    #Climate type reclassed between 100 and 900
    #LU reclassed between 3000 and 62 000
    
    
    # 'Code-value' reference matrix: convert the .txt into array
    cvrm = np.loadtxt(txtFile, dtype=np.int32, skiprows=1)  
    if cvrm.ndim == 1:
        cvrm = cvrm.reshape((cvrm.size // 2, 2))
    code, value = cvrm.T  # Transposed array
    # Lookup table with all the codes of npRaster. '+1' encompasses all codes
    #lut = np.arange(raster.max() +1)
    max_code = max(code)
    lut = np.zeros(max(max_code + 1, raster.max() + 1), dtype=np.int32)
    # Assign a value for the lookup table based on the associated code
    lut[code] = value
    # Reclassify npRaster filled w/ codes according to lookup table's values
    raster_reclassif = (lut[raster]).astype(np.int32)
    return raster_reclassif

def get_raster_combi(raster1, raster2):
    """Combine two array. The output array is the sum of the two input array"""
    
    combi = raster1 + raster2
    
    # Creating mask where soil == 0 (meaning that the cell represents water) When enough RAM
    #mask = np.ma.masked_where(raster2 == 0,raster2 )
    # Replacing the output cells with the mask 
    #combi[mask.mask] = 0 
    return combi




# Output and input directory
resultsDir = r'E:/Brazil/Model/Carbon/Output/Result_priliminaire'
Biome_directory = r'E:/Brazil/Data/Zone/Biome/Biome_reprojected'
Phyto_directory = r'E:/Brazil/Data/Zone/Phytophisiome/Phyto_reprojected'
ID_directory=r'E:/Brazil/Data/Land use/Bezera/ID_reprojected'

Biome_reclass_path = r'E:/Brazil/Model/Carbon/Input/Biome_reclass.txt'





dfs = []



for i in range(0 , 35): # Brasil divided into 36 array
    if i not in [0, 1, 4, 5, 6, 7, 11, 12, 17, 34, 35]: #blank array
       print (i)
       
       # Climate type and preparation before combining
       Biome_file_path = os.path.join(Biome_directory, f'Biome{i}.tif')
       Biome = rasterToArray(Biome_file_path)
       Biome_reclass=reclass(Biome_reclass_path, Biome)
       Biome = None
       
       # Climate type and preparation before combining
       Phyto_file_path = os.path.join(Phyto_directory, f'Phyto{i}.tif')
       Phyto = rasterToArray(Phyto_file_path)
       
       
       # Creation of an array with a distinct ID (CSS) per soil (SS) and climate type (C)
       BiomePhyto=get_raster_combi(Biome_reclass, Phyto) 
       Biome_reclass = None
       Phyto = None
       
       ID_file_path = os.path.join(ID_directory,  f'ID{i}.tif')
       ID = rasterToArray(ID_file_path)
       
       # Extraction of all the ID in this part of Brasil
       ID_All, count_ID= np.unique(ID, return_counts=True)
       
       # Initialize lists to store results for this iteration
       iteration_data = []
       
       for idx, zone_ID in enumerate(ID_All):
           
           # Extract values from the raster for the current zone
           ID_in_zone = BiomePhyto[ID == zone_ID].astype(np.int64)
           # Count occurrences of each unique value in values_in_zone
           BiomePhyto_ID, counts_BiomePhyto_ID = np.unique(ID_in_zone, return_counts=True)
           
           # Add rows for each unique value within the zone
           for j, BiomePhyto_value in enumerate(BiomePhyto_ID):
               iteration_data.append([i, zone_ID, BiomePhyto_value, counts_BiomePhyto_ID[j], count_ID[idx], (counts_BiomePhyto_ID[j]/count_ID[idx])])
               print(f"Part = {i}, ID = {zone_ID}, ClimSoil = {BiomePhyto_value}, Count_BiomePhyto = {counts_BiomePhyto_ID[j]},Count_ID = {count_ID[idx]} % = {(counts_BiomePhyto_ID[j]/count_ID[idx])}")
           
           # Create DataFrame for this iteration's data
           iteration_df = pd.DataFrame(iteration_data, columns=['Part', 'ID', 'BiomePhyto', 'CountBiomePhyto', 'CountID','%']) 
           
       # one dataframe per part (list)
       dfs.append(iteration_df)
           
       BiomePhyto = None
       ID = None
       
# To merge the dataframe of each part       
all_results_dfs = pd.concat(dfs, ignore_index=True)


# Group by the values in the second and third columns and sum the values in the fourth and fifth columns
summed_df = all_results_dfs.groupby(['ID', 'BiomePhyto'])[['CountBiomePhyto', 'CountID']].sum().reset_index()

# Filter out rows where 'ClimSoil' is different from 0 and bigger than 10
filtered_df = summed_df[summed_df['BiomePhyto'] != 0]
filtered_df = filtered_df[filtered_df['BiomePhyto'] > 100]                                                                                                                                                                                    

# Filter out rows where 'ClimSoil' is different 10,20,30,40,50,50,60,70,80,90
filtered_df = filtered_df[~filtered_df['BiomePhyto'].isin([100, 200, 300, 400, 500, 600])]

# Filter out rows where 'ID' is different 0
filtered_df = filtered_df[filtered_df['ID'] != 0]

# Group by 'ID' and sum the 'CountSoilClim' column
summed_count_df = filtered_df.groupby('ID')['CountBiomePhyto'].sum().reset_index()

# Rename the 'CountSoilClim' column to 'CountID'
summed_count_df.rename(columns={'CountBiomePhyto': 'CountID'}, inplace=True)

# Merge the summed 'CountID' DataFrame with the filtered DataFrame on 'ID'
filtered_df = pd.merge(filtered_df, summed_count_df, on='ID')

# Now 'CountID' will be equal to the sum of 'CountSoilClim' for each 'ID'

filtered_df.drop(columns=['CountID_x'], inplace=True)

filtered_df['%'] = filtered_df['CountBiomePhyto'] / filtered_df['CountID_y']



excel_file_path = os.path.join(resultsDir, 'Id_stats_BiomePhyto.xlsx')
filtered_df.to_excel(excel_file_path, index=False, engine='openpyxl')
    