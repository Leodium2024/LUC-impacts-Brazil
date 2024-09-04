import geopandas as gpd
from shapely.geometry import Point
import random

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
study_area_shp = "C:/Users/Gerar007/Documents/Brazil/Data/Zone/State/Federal/Brazil.shp"
exclusion_areas_gpkg = "C:/Users/Gerar007/Documents/Brazil/Data/Biodiversity/HSW_Mammalia/Brasil/HMW_Artiodactyla.gpkg"
output_shp = "C:/Users/Gerar007/Documents/Brazil/Data/Biodiversity/HSW_Mammalia/RandomPoint/RandomPoint.shp" 
num_points = 100
generate_random_points_excluding_areas(study_area_shp, exclusion_areas_gpkg, output_shp, num_points)



