import geopandas as gpd
from shapely.geometry import Point

def distance_meters(place_gdf, proj_crs = "EPSG:32619"):
    bounds = place_gdf.to_crs(proj_crs).total_bounds # Projected crs - meters
    
    x1 = bounds[0] + (bounds[2] - bounds[0])/3
    x2 = bounds[2] - (bounds[2] - bounds[0])/3
    y = bounds[1] + (bounds[3] - bounds[1])/2

    points = gpd.GeoSeries([Point(x1, y), Point(x2, y)], crs=proj_crs)
    return points[0].distance(points[1])

def result_rounding(results):
    # reduces results size for export
    # score columns such as there.there_index, employment_subtotal etc -> 3 decimal places
    # distance columns such as employment1 -> 0 decimal places (nearest metre)
    # avoid doing anything to connect_id, x or y
    rounding_dict = {**{k:3 for k in results.columns if "Index" in k or "." in k
                        and 'connect_id' not in k},
                     **{k:0 for k in results.columns if "Index" not in k and "." not in k
                        and k != 'x'
                        and k != 'y'}}
    return results.round(rounding_dict)

def export(folder, result, area, resultname, crs=proj_crs):
    result_gdf = gpd.GeoDataFrame(result_rounding(result), 
                                    geometry = gpd.GeoSeries.from_xy(result.x, result.y, crs=crs))
    result_gdf.to_file((folder + area + resultname + ".gpkg"), driver="GPKG")
    result_rounding(result.filter(regex='x|y$|_')).to_csv(folder + area + resultname + ".csv")