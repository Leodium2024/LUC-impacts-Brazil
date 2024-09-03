import geopandas as gpd
from shapely.geometry import Point
import random

#Pseudo-absence data resulting from a random sampling in areas with habitat unsuitable for the species. 
#This method has been shown to be the most reliable for the logistic model, especially when using a large number of pseudo-absence data with equal weighting for presence and absence (51). Those recommendations were applied in this research and 1000 pseudo-absence data were randomly sampled in unsuitable areas. 
#Expert species range maps were used to identify unsuitable areas

########################################################################################
##################     1) Generation of Pseudo-Absence records     #####################
########################################################################################

def generate_random_points_excluding_areas(study_area_shp, exclusion_areas_gpkg, output_shp, num_points):
    # Load the study area from the shapefile
    study_area_gdf = gpd.read_file(study_area_shp)
    
    # Ensure the study area contains only one polygon
    if len(study_area_gdf) != 1:
        raise ValueError("The study area shapefile should contain only one polygon.")
    
    # Get the study area polygon
    study_area = study_area_gdf.geometry.iloc[0]
    
    # Get the bounding box of the study area
    minx, miny, maxx, maxy = study_area.bounds
    
    # Load the exclusion areas from the GPKG file
    exclusion_areas_gdf = gpd.read_file(exclusion_areas_gpkg)
    
    # List to store all generated points
    all_points = []
    
    for idx, exclusion_area in exclusion_areas_gdf.iterrows():
        exclusion_polygon = exclusion_area.geometry
        points = []
        
        while len(points) < num_points:
            # Generate a random point within the bounding box
            random_point = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
            
            # Check if the point is within the study area and not within the exclusion area
            if study_area.contains(random_point) and not exclusion_polygon.contains(random_point):
                points.append({
                    'geometry': random_point,
                    'exclusion_id': idx,  # Add an ID to differentiate based on the exclusion area
                    **exclusion_area.drop('geometry').to_dict()  # Add the attributes from the exclusion area
                })
        
        # Append the points to the list of all points
        all_points.extend(points)
    
    # Create a GeoDataFrame with all the points
    points_gdf = gpd.GeoDataFrame(all_points, crs=study_area_gdf.crs)
    
    # Save the GeoDataFrame to a shapefile
    points_gdf.to_file(output_shp, driver="ESRI Shapefile")

# Example usage:
study_area_shp = "Brazil.shp"
exclusion_areas_gpkg = "HMW_Mammal.gpkg"
output_shp = "Absence.shp" 
num_points = 1000
generate_random_points_excluding_areas(study_area_shp, exclusion_areas_gpkg, output_shp, num_points)

########################################################################################
##################     2) Extraction of environemental variables     ###################
########################################################################################

# Environmental Variable Extraction in ArcGIS: Using the generated SHP file in ArcGIS to extract environmental variable values for each species absence record.


