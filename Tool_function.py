"""
Created on 20/03/2025

@author: Thomas GERARD

Version: 1

"""

import numpy as np
import pandas as pd
from osgeo import gdal
import os
from os.path import join as opj
import osgeo.osr as osr
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings


# --------------------------------------------------
# Import & Export Functions
# --------------------------------------------------

def referenceRaster(Input_directory):
    """Set the reference raster in which its geoparameters will be used to 
    create/save new rasters/numpy arrays."""
    refRaster = os.path.join(Input_directory,'Mask_Albers.tif')
    return refRaster

def rasterToArray(inpRaster):
    """Return a numpy array based on a raster file. The 'No data' of the raster 
    are converted to zero"""
    inpRaster = gdal.Open(inpRaster)
    unqBand = inpRaster.GetRasterBand(1)
    noData = unqBand.GetNoDataValue()
    npArr = inpRaster.ReadAsArray()
    npArr[npArr == noData] = 0  # Converting noData values to zero
    #if npArr[npArr < 0].any() == True:  # Making sure no neg value is in array
        #npArr[npArr < 0] = 0
    return npArr

def arrayToMap(inpArray, outRaster, mask, Input_directory):
    """Convert an array to a raster map and save to disk (tiff format). 
    Arguments: inpArray = array to be converted; outRaster = string with the 
    filename representing the raster map"""
    # Getting geoparameters of the original raster
    refRaster = gdal.Open(referenceRaster(Input_directory))
    geotransform = refRaster.GetGeoTransform()
    originX = geotransform[0]  # top left-x
    originY = geotransform[3]  # top left-y
    pixelWidth = geotransform[1]  # w-e pixel resolution
    pixelHeight = geotransform[5]  # n-s pixel resolution
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    # Apply the mask to the input array and set no-data value outside the mask
    noDataValue = -9999  # Define a suitable no-data value
    inpArray = np.where(mask, inpArray, noDataValue)
    
    # Creating new raster, then adding geoparameters
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(
        opj(outRaster + ".tiff"), cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(inpArray)
    outband.SetNoDataValue(noDataValue)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(refRaster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


# --------------------------------------------------
# Plot functions
# --------------------------------------------------

#Define scale colors for mapping
colors1 = [(1, 1, 1), (0, 1, 0)]  # White to Green
colors2 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0)]  # White to Green to Dark Green 
colors3 = [(1, 0, 0), (1, 1, 1), (0, 1, 0)] # red to white to Green
colors4 = [(1, 0, 0), (1, 1, 1), (0, 0, 1)] # red to white to Blue
colors5 = [(1, 1, 1), (0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors5 = [(1, 1, 1),(0, 1, 0), (0, 0.5, 0), (0, 0, 0)]  # White to Green to Dark Green to black
colors6 = [(1, 1, 1), (0.5, 1, 0.5), (0, 1, 0), (0, 0.5, 0), (0, 0.25, 0), (0, 0, 0)]


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


# --------------------------------------------------
# Tool functions
# --------------------------------------------------

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


# --------------------------------------------------
# Biodiversity function
# --------------------------------------------------

def get_ottid(Performance):
    
    Performance=Performance[(Performance['N_Pre'] >=10)]  #More than 10 records
    Performance=Performance[(Performance['N_Abs'] >=10)] 
    Performance=Performance[(Performance['NNI'] >=0.5)]   #Not clustered records
    Performance=Performance[(Performance['AUC'] >=0.5)]   #AUC >= 0.8
    Performance=Performance[(Performance['TSS'] >=0)]   #TSS >= 0.6
    Ott_id =Performance['Ott_id'].unique()
    return Ott_id


def get_distribution_map(Input_directory,Species_ID, Coeficient, Performance, land_use, env):
    
    #Upload of the range of the species
    Species_expert_range = rasterToArray(os.path.join(Input_directory,'Biodiversity','Species_Range', f'{Species_ID}.tiff'))
    
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
    Grassland_poly= Coeficient_ID.iloc[0, 2]*(land_use["Grassland"]+land_use["Pasture"])
    Cropland_poly= Coeficient_ID.iloc[0, 3]*land_use["Cropland"]
    Forestry_poly= Coeficient_ID.iloc[0, 4]*land_use["Forestry"]
    Bio_1_poly= Coeficient_ID.iloc[0, 5]*env["Bio_1"]
    Bio_4_poly= Coeficient_ID.iloc[0, 6]*env["Bio_4"]
    Bio_12_poly= Coeficient_ID.iloc[0, 7]*env["Bio_12"]
    Bio_15_poly= Coeficient_ID.iloc[0, 8]*env["Bio_15"]
    Elevation_poly= Coeficient_ID.iloc[0, 9]*env["Elevation"]
    Slope_poly= Coeficient_ID.iloc[0, 10]*env["Slope"]
    Intercept= Coeficient_ID.iloc[0, 11]

    #3) Calculate the logistic model
    warnings.filterwarnings('ignore')
    exponent = Forest_poly +  Grassland_poly + Cropland_poly + Forestry_poly+  Bio_1_poly + Bio_4_poly + Bio_12_poly + Bio_15_poly + Elevation_poly + Slope_poly + Intercept
    z = np.exp(-exponent)
    Model_fitted = 1/(1+z)
    Total_fitted=Total_fitted+ Model_fitted
      
    
    ##### Distribution Map #####
    Species_distribution = Total_fitted.applymap(lambda x: 1 if x > Threshold else 0)
    Species_distribution = Species_distribution*Species_expert_range
    
    
    return Species_distribution
