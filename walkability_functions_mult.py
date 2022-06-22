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
from shapely.geometry import MultiPoint

#%% functions

def single_points(data_gdf, area_gdf=None):
    # mostly should be replaced by poly_vertices but can be used where polygon centroids
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
    # Add area column, convert polygons to vertices Multipoints, explode to Points.
    # Thus original polygon area will be retained with the vertices,
    # which is useful for some analyses

    exploded = data_gdf.copy()
    exploded['poly_area'] = exploded.to_crs('+proj=wintri').area
    exploded.geometry = exploded.geometry.apply(lambda x: MultiPoint(list(x.exterior.coords)))
    exploded = exploded.explode(index_parts=True)
    
    return exploded

def categorise_pois(pois, categories, old_column, new_column):
    # pois as df or gdf, categories as dictionary, old_column, new_column as string
    
    categories_dict = {
        newkey: key 
        for key, values in categories.items() 
        for newkey in values
        }
    
    remain_pois = pois[~pois[old_column].isin(categories_dict)].copy()
    print("Tags present in the dataset but not categorised:")
    print(pd.unique(remain_pois[old_column]))
    
    pois = pois[pois[old_column].isin(categories_dict)].copy()
    pois[new_column] = pois[old_column]
    pois.replace({new_column: categories_dict},inplace=True)
    
    return pois

def remove_duplicate_pois(datasets, buffer=10):
    # input is a list of maybe one or several datasets
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

def access_weight(x, distance):
    beta = -0.001
    if x == distance:
        return 0
    else:
        return math.exp(beta*x)

def walk_index(distance_network, pois, poi_weights, poi_numbers, distance):
    ########### need to change this to keep names of pois but not much else? 
    # otherwise need to drop columns above
    # need to intialise results probably?
    
    results = distance_network.nodes_df.copy()
    
    for category in poi_weights:
        results = results.copy()  # to remove fragmentation warning
        
        num_pois = poi_numbers[category]
        x, y = (pois[pois.category == category]['geometry'].x, 
                pois[pois.category == category]['geometry'].y)
        cat_total = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(cat_total)))
        
        if len(x) == 0:
            print("Category "+ category + " is empty")
            results[str(cat_name)] = 0
            
        else:
            distance_network.set_pois(category, distance, num_pois, x, y)
            
            access = distance_network.nearest_pois(
                distance=distance, category=category, num_pois=num_pois)
            
            for i in range(num_pois):
                col_name = ''.join((str(category),str(i+1)))
                results[col_name] = access[i+1]
            
            results[cat_name] = ((access.applymap(access_weight,distance=distance))*
                                 poi_weights[category]/poi_numbers[category]).sum(axis=1) 
            
            print("Finished category: " + category)
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_weights]   
    weights = sum(poi_weights.values())
                           
    results['Walk_Index'] = 100/weights*(results[col_list].sum(axis=1))
    
    return results

def primal(distance_network, pois, distance, max_items, categories):
    results = distance_network.nodes_df.copy()
    
    for cat in categories:
        x, y = (pois[pois.category == cat]['geometry'].x, 
                    pois[pois.category == cat]['geometry'].y)
        
        if len(x) == 0:
            print("Category "+ cat + " is empty")
            results[cat] = 0
            
        else:
            distance_network.set_pois(category=cat, maxdist=distance, maxitems=max_items, x_col=x, y_col=y)

            access = distance_network.nearest_pois(
                distance=distance, category=cat, num_pois=max_items, include_poi_ids=False)

            number_accessible = max_items - (access == 2400).sum(axis=1)

            results[cat] = number_accessible
            print("Finished category: " + cat)

    return results

def distance_weight(x):
    if x < 400:
        return 1
    if x >= 2399:
        return 0
    if x > 2000:
        return 0.05
    else:
        return 1 - ((1-0.05)/(2000-400))*(x-400)


