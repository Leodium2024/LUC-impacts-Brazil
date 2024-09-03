from osgeo import gdal
#import rasterio
import os
import numpy as np
import pandas as pd
import statistics


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

'''
def rasterToArray(inpRaster):
    """Return a numpy array based on a raster file. The 'No data' of the raster 
    are converted to zero"""
    with rasterio.open(inpRaster) as src:
        npArr = src.read(1)  # Assuming single band raster
        noData = src.nodata
        npArr[npArr == noData] = 0  # Converting noData values to zero
        npArr = np.where(npArr < 0, 0, npArr)  # Setting negative values to zero
    return npArr

'''
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


def get_raster_combi(raster1, raster2):
    """Combine two array. The output array is the sum of the two input array"""
    
    combi = raster1 + raster2
    
    # Creating mask where soil == 0 (meaning that the cell represents water) When enough RAM
    #mask = np.ma.masked_where(raster2 == 0,raster2 )
    # Replacing the output cells with the mask 
    #combi[mask.mask] = 0 
    return combi



def setOutputPath():
    """Create necessary folders to run the script and return the folder in which
    the results will be stored"""
    results = os.path.join('Output')
    paths = [results]
    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)
    return results





# Output and input directory
resultsDir = setOutputPath()
Clim_directory = r'C:\Users\Gerar007\Documents\Brazil\Data\Zone\Clim\KoppenAlvares\Clim_Albers'
Soil_directory = r'C:/Users/Gerar007/Documents/Brazil/Data/Zone/Soil/IBGE_albers'
SOC_directory = r'C:/Users/Gerar007/Documents/Brazil/Data/Carbon/SOC/2015Albers'
lu_directory = r'C:\Users\Gerar007\Documents\Brazil\Data\Land use\Mapbiomass\LU15_albers'

Clim_reclass_path = r'C:/Users/Gerar007/Documents/Brazil/Model/Carbon/Input/climate_reclass.txt'
Soil_reclass_path = r'C:/Users/Gerar007/Documents/Brazil/Model/Carbon/Input/soil_reclass.txt'
LU_reclass_path = r'C:/Users/Gerar007/Documents/Brazil/Model/Carbon/Input/LU_reclass.txt'

# Initialize an empty dictionary
zone_dict = {i: [] for i in range(101,801)}  # IDs from 0 to 800
zone_dict = {key: value for key, value in zone_dict.items() if key % 10 >= 1 and key % 10 not in [0, 8, 9]}


for i in range(16, 36): # Brasil divided into 36 array
    if i not in [0, 1, 4, 5, 6, 7, 11, 12, 17, 34, 35]: #blank array
       print (i)
    # Climate type and preparation before combining
       Clim_file_path = os.path.join(Clim_directory, f'clim{i}.tif')
       Clim = rasterToArray(Clim_file_path)
       Clim_reclass=reclass(Clim_reclass_path, Clim)
       Clim = None 
    
    # Soil type and preparation before combining
       Soil_file_path = os.path.join(Soil_directory,  f'Soil{i}.tif')
       Soil = rasterToArray(Soil_file_path)
       Soil_reclass=reclass(Soil_reclass_path, Soil)
       Soil = None
    
    # Creation of an array with a distinct ID (CS) per soil (S) and climate type (C)
       ClimSoil=get_raster_combi(Clim_reclass, Soil_reclass) 
       Soil_reclass = None
       Clim_reclass= None
       print ('ClimSoil done')
    
    # LU type and preparation before combining
       lu_file_path = os.path.join(lu_directory,  f'LU{i}.tif')
       LU = rasterToArray(lu_file_path)
       LU_reclass= reclass(LU_reclass_path, LU)
       LU = None
           
    # Creation of an array with a distinct ID (XXCSS) per LU (XX), soil (SS) and climate type(C)
       LUClimSoil=get_raster_combi(LU_reclass, ClimSoil)
       LU_reclass = None
       ClimSoil = None
       print ('LuClimSoil done')
       
    #Import of Mapbiomas SOC
       SOC_file_path = os.path.join(SOC_directory,  f'SOC{i}.tif')
       SOC = rasterToArray(SOC_file_path)
       print ('SOC done')
       
    

    # Iterate through the keys of the zone_dict
       for zone_value in zone_dict.keys():
          if zone_value != 0:
           # Extract SOC values where LUclimSoil equals the zone ID
               soc_values = SOC[LUClimSoil == zone_value]
           # Store SOC values in the corresponding list in the dictionary
               zone_dict[zone_value].extend(soc_values.tolist())
       
       SOC=None
       LUClimSoil=None
       
    
    #To empty the memory before the next iteration
       variables_to_delete = [
           'Clim_file_path', 'Clim', 'Clim_reclass_file_path', 'Clim_reclass',
           'Soil_file_path', 'Soil', 'Soil_reclass_file_path',
           'ClimVege', 'lu_file_path', 'LU', 'LU_reclass',
           'SOC_file_path', 'SOC', 'unique_zones',
           'result_array', 'idx', 'zone_value', 'values_in_zone',
           'column_names', 'df', 'excel_file_path']
    # Delete variables using a loop
       for variable_name in variables_to_delete:
           try:
               del globals()[variable_name]
           except KeyError:
               pass  # Variable may not exist
  
            
# Define sheet names
sheet_names = {
    1: 'Forest',
    2: 'Grassland',
    3: 'Pasture',
    4: 'Semiperennial',
    5: 'Temporary',
    6: 'Perenial',
    7: 'Forestry'
}

# Create a dictionary to hold DataFrames for each first digit
dfs_per_first_digit = {}

# Calculate statistics for each list and store them in the dictionary of DataFrames
for key, lst in zone_dict.items():
    first_digit = int(str(key)[0])  # Get the first digit of the key
    rest_of_key = int(str(key)[1:])  # Get the rest of the key after the first digit as integer
    average = np.mean(lst)
    median = np.median(lst)
    std_dev = np.std(lst)
    
    # Create a DataFrame for this key if it doesn't exist yet
    if first_digit not in dfs_per_first_digit:
        dfs_per_first_digit[first_digit] = pd.DataFrame(columns=['ID', 'Average', 'Median', 'Standard Deviation'])
    
    # Append the statistics to the corresponding DataFrame
    new_row = {'ID': rest_of_key, 'Average': average, 'Median': median, 'Standard Deviation': std_dev}
    dfs_per_first_digit[first_digit] = pd.concat([dfs_per_first_digit[first_digit], pd.DataFrame([new_row])], ignore_index=True)

# Export each DataFrame to a separate sheet in an Excel file
excel_file_path = os.path.join(resultsDir, 'SOC_Stat.xlsx')
with pd.ExcelWriter(excel_file_path, engine='openpyxl') as writer:
    for first_digit, df in dfs_per_first_digit.items():
        sheet_name = sheet_names.get(first_digit, f'Sheet_{first_digit}')
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print(f"All iterations data have been exported to {excel_file_path}")
  