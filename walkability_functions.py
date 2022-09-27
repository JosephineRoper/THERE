# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 12:15:25 2021

@author: z3258367
All POI and network data should be in a projected CRS

"""

import numpy as np
import geopandas as gpd
import pandas as pd
import math
import osmnx as ox
from shapely.geometry import MultiPoint

#%% functions

def nearest_to_centroids(polygons, points, **kwargs):
    orig_geometry = polygons.geometry
    # convert polygons to centroids
    polygons.geometry = polygons.centroid
    polygons = polygons.sjoin_nearest(points, **kwargs)
    # convert back to polygons
    polygons.geometry = orig_geometry
    return polygons
    
def poi_downloader(place, poi_dictionary, proj_crs, timeout=None):
    # Download WalkTHERE points of interest from OSM using OSMnx
    # place should either be a string suitable for finding a place by name from the Nominatim API
    # or a geodataframe containing a single polygon of the place's boundaries

    # changing timeout here changes setting for any subsequent use of OSMNx
    # this seems unavoidable
    if isinstance(timeout, int):
        ox.settings.timeout = timeout
    
    tags = {}
    for category, values in poi_dictionary.items():
        tags = {x: values.get(x,[]) + tags.get(x,[]) for x in set(values).union(tags)}

    if type(place) == str:
        gdf = ox.geometries_from_place(place, tags).to_crs(proj_crs)
    elif type(place) == gpd.GeoDataFrame:
        place = place.to_crs("EPSG:4326")
        bbox = place.bounds
        bbox_pois = ox.geometries.geometries_from_bbox(bbox['maxy'][0], bbox['miny'][0], bbox['maxx'][0], bbox['minx'][0], tags)
        gdf = gpd.clip(bbox_pois, place, keep_geom_type=False).to_crs(proj_crs)
        #? does this work with multiple polygons
    else:
        print("'place' should be a string for querying OSM, or a geodataframe containing a polygon of the place boundaries.")
        return

    # OSM POIs include domestic swimming pools in some areas. This line removes swimming pools less than 100m2.
    # Same for domestic tennis courts appearing as 'pitches'. Remove pitches below 450m2.
    gdf = gdf[~((gdf['leisure']=='swimming_pool') & (gdf.area < 100))]
    gdf = gdf[~((gdf['leisure']=='pitch') & (gdf.area < 450))]
        
    gdf['orig_geometry'] = gdf.geometry
    # convert all to centroids
    gdf.geometry = gdf.centroid
    return gdf

def single_points(data_gdf, area_gdf=None):
    # similar to poly_vertices but can be used where polygon centroids
    # are wanted rather than vertices
    orig_crs = data_gdf.crs
    
    if area_gdf is not None:
        area_gdf.to_crs(data_gdf.crs, inplace = True)
        data_gdf = gpd.clip(data_gdf, area_gdf)
        
    exploded = data_gdf.explode(index_parts=True)
    exploded.geometry = (exploded.geometry
                         .to_crs('+proj=wintri')
                         .centroid
                         .to_crs(orig_crs))
    return exploded

def poly_vertices(data_gdf):
    # Add area column, convert polygons to vertices as Multipoints, explode to Points.
    # Thus original polygon area will be retained with the vertices,
    # which is useful for some analyses

    exploded = data_gdf.copy()
    exploded['poly_area'] = exploded.to_crs('+proj=wintri').area
    exploded.geometry = exploded.geometry.apply(lambda x: MultiPoint(list(x.exterior.coords)))
    exploded = exploded.explode(index_parts=True)
    
    return exploded

def remove_duplicate_pois(datasets, buffer=10):
    # input is a list of maybe one or several geodataframes
    # buffer default is 10m = assumes projected CRS
    
    pois = pd.concat(datasets).reset_index(drop=True)
    pois_buffer = pois.copy()
    pois_buffer.geometry = pois.geometry.buffer(buffer)
    
    joined_pois = gpd.sjoin(
        pois, pois_buffer, how="left", predicate='intersects', 
        lsuffix='1', rsuffix='2')
    
    joined_pois = joined_pois[joined_pois['category_1'] == joined_pois['category_2']].copy()
    joined_pois['unique_id'] = joined_pois.apply(
        lambda row: np.nanmin([row.name, row['index_2']]), axis=1)
    joined_pois = joined_pois.drop_duplicates(subset='unique_id', keep=False)
    
    print("Removed " 
          + "{0:.2f}".format((1-len(joined_pois)/len(pois))*100) 
          + "% duplicate points from dataframes")  
    
    return pois[pois.index.isin(joined_pois.index)]

def access_weight(x, distance, beta=0.001):
    if x == distance:
        return 0
    else:
        return math.exp(-beta*x)

def there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5):
    # return_no = 5 is the default setting, to return individual distance results
    # for maximum 5 closest points in each category. Enables some debugging &
    # visualisation options, but not returning an overly large results matrix
    # (as poi_nums searched may be 100s per category)
    
    results = distance_network.nodes_df.copy()

    total_weight = sum(poi_weights)
    
    for category in poi_weights.index:
        results = results.copy()  # to remove fragmentation warning
        
        dim_const = poi_lambdas[category]
        dist_const = poi_gammas[category]
        num_pois = poi_nums[category]
        weight = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(weight)))

        if category not in poi_dictionary:
            print("Category", category, "is not in the POI dictionary")

        else:
            relevant_pois = gpd.GeoDataFrame()
            for key in poi_dictionary[category]:
                if key in pois:
                    relevant_pois = pd.concat([relevant_pois, pois.loc[(pois[key].isin(poi_dictionary[category][key]))]])
            
            if len(relevant_pois) == 0:
                print("No pois in category: "+ category)
                results[str(cat_name)] = 0

            elif poi_variables[category] == 'count':
                
                x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)
                distance_network.set_pois(category, distance, num_pois, x, y)
                
                access = distance_network.nearest_pois(
                    distance=distance, category=category, num_pois=num_pois)

                discounted = access.applymap(access_weight,distance=distance,beta=dist_const).sum(axis=1)

                results[cat_name] = weight*(1-np.exp(-dim_const*discounted))
                
                # this provides columns with the distance of the return_no closest destinations in the category
                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                    
            else:
                # makes eg. a job count column the index column, so that 'include_poi_ids' returns job counts
                relevant_pois = relevant_pois.set_index(poi_variables[category])
                
                x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)

                distance_network.set_pois(category, distance, num_pois, x, y)

                access = distance_network.nearest_pois(
                    distance=distance, category=category, num_pois=num_pois, include_poi_ids=True)

                results[poi_variables[category]] = ((access.iloc[:,0:num_pois].applymap(access_weight, distance=distance,beta=dist_const))*
                                    access.iloc[:,num_pois:2*num_pois].values
                                    ).sum(axis=1)

                results[cat_name] = weight*(1-np.exp(-dim_const*results[poi_variables[category]]))

                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                    
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary]   
    
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
    
    return results
