import geopandas as gpd
from shapely.geometry import Point

def distance_meters(place_gdf, proj_crs = "EPSG:32619"):
    bounds = place_gdf.to_crs(proj_crs).total_bounds # Projected crs - meters
    
    x1 = bounds[0] + (bounds[2] - bounds[0])/3
    x2 = bounds[2] - (bounds[2] - bounds[0])/3
    y = bounds[1] + (bounds[3] - bounds[1])/2

    points = gpd.GeoSeries([Point(x1, y), Point(x2, y)], crs=proj_crs)
    return points[0].distance(points[1])
