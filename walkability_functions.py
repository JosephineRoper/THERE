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

#%% functions

# clip to area, convert to single part, get centroids
def single_points(data_gdf, area_gdf=None):
    
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

def categorise_pois(pois, categories, category_column):
    # pois as df or gdf, categories as dictionary, category_column as string
    
    categories_dict = {
        newkey: key 
        for key, values in categories.items() 
        for newkey in values
        }
    
    remain_pois = pois[~pois[category_column].isin(categories_dict)].copy()
    print("Tags present in the dataset but not categorised:")
    print(pd.unique(remain_pois[category_column]))
    
    pois = pois[pois[category_column].isin(categories_dict)].copy()
    pois.replace({category_column: categories_dict},inplace=True)
    
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
    
    joined_pois = joined_pois[joined_pois['fclass_1'] == joined_pois['fclass_2']].copy()
    joined_pois['unique_id'] = joined_pois.apply(
        lambda row: np.nanmin([row.name, row['index_2']]), axis=1)
    joined_pois = joined_pois.drop_duplicates(subset='unique_id', keep=False)
    
    print("Removed " 
          + "{0:.2f}".format((1-len(joined_pois)/len(pois))*100) 
          + "% duplicate points from dataframes")  
    
    return pois[pois.index.isin(joined_pois.index)]


def walk_index(distance_network, pois, poi_weights, distance):
    ########### need to change this to keep names of pois but not much else? 
    # otherwise need to drop columns above
    # need to intialise results probably?
    
    def access_weight(x):
        beta = -0.001
        if x == distance:
            return 0
        else:
            return math.exp(beta*x)
    
    results = distance_network.nodes_df.copy()
    
    for category in poi_weights:  
        num_pois = len(poi_weights[category])
        x, y = (pois[pois.fclass == category]['geometry'].x, 
                pois[pois.fclass == category]['geometry'].y)
        cat_total = int(sum(poi_weights[category]))
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
            
            results[cat_name] = ((access.applymap(access_weight))*
                                 poi_weights[category]).sum(axis=1) 
            
            print("Finished category: " + category)
            
    col_list = [''.join((str(category),"_",str(int(sum(poi_weights[category])))))
                for category in poi_weights]   
    weightslist = sum(list(poi_weights.values()),[])
                           
    results['Walk_Index'] = 100/sum(weightslist)*(results[col_list].sum(axis=1))
    
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


