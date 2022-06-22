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

# Download WalkTHERE points of interest from OSM using OSMnx
def poi_downloader(place, poi_dictionary, proj_crs):
    # place should either be a string suitable for finding a place by name from the Nominatim API
    # or a geodataframe
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
    else:
        print("'place' should be a string for querying OSM, or a geodataframe containing a polygon of the place.")
        return
        
    gdf['orig_geometry'] = gdf.geometry
    # convert all to centroids
    gdf.geometry = gdf.centroid
    return gdf

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

def categorise_pois_new(pois, categories, old_column, new_column='category'):
    # pois as df or gdf, categories as dictionary, old_column, new_column as string
    
    categories_dict = {
        newkey: key 
        for key, values in categories.items() 
        for newkey in values
        }
    
    remain_pois = pois[~pois[old_column].isin(categories_dict)].copy()
    print("Some tags are present in the dataset but not in the category dictionary. " +
        "POIs with these tags have been removed:")
    print(pd.unique(remain_pois[old_column]))
    
    pois = pois[pois[old_column].isin(categories_dict)].copy()
    pois['category'] = pois[old_column]
    
    # taken this out because it no longer makes sense if some pois
    # can count in multiple categories
    #pois.replace({new_column: categories_dict},inplace=True)
    
    return pois

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

def remove_duplicate_pois_new(datasets, buffer=10):
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

def access_weight(x, distance, beta=0.001):
    if x == distance:
        return 0
    else:
        return math.exp(-beta*x)

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
    
def walk_index_new(distance_network, pois, poi_dictionary, poi_weights, poi_numbers, distance):
    ########### need to change this to keep names of pois but not much else? 
    # otherwise need to drop columns above
    
    results = distance_network.nodes_df.copy()
    
    for category in poi_weights:
        results = results.copy()  # to remove fragmentation warning
        
        num_pois = poi_numbers[category]
        x, y = (pois[pois['category'].isin(poi_dictionary[category])]['geometry'].x, 
                pois[pois['category'].isin(poi_dictionary[category])]['geometry'].y)
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

def walk_index_m(distance_network, pois, poi_weights, distance):
    ########### need to change this to keep names of pois but not much else? 
    # otherwise need to drop columns above
    # need to intialise results probably?
    
    results = distance_network.nodes_df.copy()
    
    for cat in poi_weights:  
        num_pois = len(poi_weights[cat])
        x, y = (pois[pois.category == cat]['geometry'].x, 
                pois[pois.category == cat]['geometry'].y)
        cat_total = int(sum(poi_weights[cat]))
        cat_name = ''.join((str(cat),"_",str(cat_total)))
        
        if len(x) == 0:
            print("Category "+ cat + " is empty")
            results[str(cat_name)] = 0
            
        else:
            distance_network.set_pois(cat, distance, num_pois, x, y)
            
            access = distance_network.nearest_pois(
                distance=distance, category=cat, num_pois=num_pois)
            
            for i in range(num_pois):
                col_name = ''.join((str(cat),str(i+1)))
                results[col_name] = access[i+1]
            
            results[cat_name] = ((access.applymap(access_weight, distance=distance))*
                                 poi_weights[cat]).sum(axis=1) 
            
            print("Finished category: " + cat)
            
    col_list = [''.join((str(category),"_",str(int(sum(poi_weights[category])))))
                for category in poi_weights]   
    weightslist = sum(list(poi_weights.values()),[])
                           
    results['Walk_Index'] = 100/sum(weightslist)*(results[col_list].sum(axis=1))
    
    return results

def walk_index_exp(distance_network, pois, poi_dictionary, poi_weights, poi_lambdas, distance):
    
    results = distance_network.nodes_df.copy()

    num_pois = 300
    
    for category in poi_weights:
        results = results.copy()  # to remove fragmentation warning
        
        constant = poi_lambdas[category]
        x, y = (pois[pois['category'].isin(poi_dictionary[category])]['geometry'].x, 
                pois[pois['category'].isin(poi_dictionary[category])]['geometry'].y)
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

            discounted = access.applymap(access_weight,distance=distance)
                     
            discounted = discounted*[constant*math.exp(-constant*(x-0.5)) for x in range(1,num_pois+1)]
            
            results[cat_name] = (discounted*poi_weights[category]).sum(axis=1) 
            
            print("Finished category: " + category)
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_weights]   
    weights = sum(poi_weights.values())
                           
    results['Walk_Index'] = 100/weights*(results[col_list].sum(axis=1))
    
    return results

def walk_index_full(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=10):
    # return_no = 10 is a default setting to return distance results for maximum
    # 10 closest points in each category. Enables some debugging & visualisation options,
    # but not returning an overly large results matrix
    
    results = distance_network.nodes_df.copy()

    total_weight = sum(poi_weights)
    
    for category in poi_weights.index:
        results = results.copy()  # to remove fragmentation warning
        
        dim_const = poi_lambdas[category]
        dist_const = poi_gammas[category]
        num_pois = poi_nums[category]
        weight = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(weight)))
        
        if len(pois[pois['category'].isin(poi_dictionary[category])]) == 0:
            print("Category "+ category + " is empty")
            results[str(cat_name)] = 0

        elif poi_variables[category] == 'count':
            
            x, y = (pois[pois['category'].isin(poi_dictionary[category])]['geometry'].x, 
                pois[pois['category'].isin(poi_dictionary[category])]['geometry'].y)
            distance_network.set_pois(category, distance, num_pois, x, y)
            
            access = distance_network.nearest_pois(
                distance=distance, category=category, num_pois=num_pois)
          
            for i in range(return_no):
                col_name = ''.join((str(category),str(i+1)))
                results[col_name] = access[i+1]

            discounted = access.applymap(access_weight,distance=distance,beta=dist_const)
                     
            discounted = discounted*[dim_const*math.exp(-dim_const*(x-0.5)) for x in range(1,num_pois+1)]
            
            results[cat_name] = (discounted*weight).sum(axis=1) 
            
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])))

        else:
            relevant_pois = pois[pois['category'].isin(poi_dictionary[category])]

            relevant_pois = relevant_pois.set_index(poi_variables[category])

            x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)

            distance_network.set_pois(category, distance, num_pois, x, y)

            access = distance_network.nearest_pois(
                distance=distance, category=category, num_pois=num_pois, include_poi_ids=True)

            results[poi_variables[category]] = ((access.iloc[:,0:num_pois].applymap(access_weight, distance=distance,beta=dist_const))*
                                access.iloc[:,num_pois:2*num_pois].values
                                ).sum(axis=1)

            for i in range(return_no):
                col_name = ''.join((str(category),str(i+1)))
                results[col_name] = access[i+1]

            discounted = np.array([1 - math.exp(-dim_const*x) for x in results[poi_variables[category]]])
            
            results[cat_name] = (discounted*weight)

            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])))
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_weights.index]   
    
    results['Walk_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
    
    return results

def there_index_osmnx(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5):
    # return_no = 5 is the default setting to return individual distance results for maximum
    # 5 closest points in each category. Enables some debugging & visualisation options,
    # but not returning an overly large results matrix (as poi_nums searched may be 100+ per category)
    
    results = distance_network.nodes_df.copy()

    total_weight = sum(poi_weights)
    
    for category in poi_weights.index:
        results = results.copy()  # to remove fragmentation warning
        
        dim_const = poi_lambdas[category]
        dist_const = poi_gammas[category]
        num_pois = poi_nums[category]
        weight = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(weight)))

        if category in poi_dictionary:
            relevant_pois = gpd.GeoDataFrame()
            for key in poi_dictionary[category]:
                relevant_pois = pd.concat([relevant_pois, pois.loc[(pois[key].isin(poi_dictionary[category][key]))]])
            
            if len(relevant_pois) == 0:
                print("No pois in category: "+ category)
                results[str(cat_name)] = 0

            elif poi_variables[category] == 'count':
                
                x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)
                distance_network.set_pois(category, distance, num_pois, x, y)
                
                access = distance_network.nearest_pois(
                    distance=distance, category=category, num_pois=num_pois)
              
                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]

                discounted = access.applymap(access_weight,distance=distance,beta=dist_const)
                         
                discounted = discounted*[dim_const*math.exp(-dim_const*(x-0.5)) for x in range(1,num_pois+1)]
                
                results[cat_name] = (discounted*weight).sum(axis=1) 

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

                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]

                discounted = np.array([1 - math.exp(-dim_const*x) for x in results[poi_variables[category]]])
                
                results[cat_name] = (discounted*weight)

            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
            
        else: print("Category", category, "is not in the POI dictionary")
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary]   
    
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
    
    return results

def itermerge(access_results_df, variable_df):
    i=0
    for column in access_results_df:
        access_results_df = access_results_df.merge(variable_df, how='left', left_on = column, right_index = True, suffixes = [None, i])
        i = i + 1
    return access_results_df

def walk_index_aggregates(distance_network, pois, poi_dictionary, poi_weights, poi_numbers, distance, variable):

    for cat in category:

        x, y = (pois[pois['category'].isin(poi_dictionary[category])]['geometry'].x, 
                pois[pois['category'].isin(poi_dictionary[category])]['geometry'].y)

        num_pois = poi_numbers[category]

        distance_network.set_pois(category=cat, maxdist=distance, maxitems=num_pois, x_col=x, y_col=y)

        access = distance_network.nearest_pois(
            distance=distance, category=cat, num_pois=num_pois, include_poi_ids=True)

        counts = itermerge(access.iloc[:,num_pois:num_pois*2], access['Jobs_count'])

        results_3['jobs'] = ((employment_access.iloc[:,0:100].applymap(access_weight, distance=maximum_dist))*
                                    jobcounts.iloc[:,100:200].values
                                    ).sum(axis=1)
        weight = 100*poi_weights['employment']/sum(poi_weights.values())

        results_3['employment'] = weight*results_3['jobs']/max(results_3['jobs'])

        results_3['Walk_Index'] = results_3['Walk_Index'] + results_3['employment']

        if results['Walk Index']:
            results['Walk_Index'] = results['Walk_Index'] + aggregate
        else:
            results['Walk_Index'] = aggregate

def primal(distance_network, pois, distance, max_items, categories, poi_dictionary):
    # what we are doing here could be replaced by the Pandana aggregate function
    # except 
    results = distance_network.nodes_df.copy()
    
    for cat in categories:
        x, y = (pois[pois['category'].isin(poi_dictionary[cat])]['geometry'].x, 
                    pois[pois['category'].isin(poi_dictionary[cat])]['geometry'].y)
        
        if len(x) == 0:
            print("Category "+ cat + " is empty")
            results[cat] = 0
            
        else:
            distance_network.set_pois(category=cat, maxdist=distance, maxitems=max_items, x_col=x, y_col=y)

            access = distance_network.nearest_pois(
                distance=distance, category=cat, num_pois=max_items, include_poi_ids=False)

            number_accessible = max_items - (access == distance).sum(axis=1)

            results[cat] = number_accessible
            print("Finished category: " + cat)

    return results

def primal_2(distance_network, pois, distance, categories, poi_dictionary):
    # what we are doing here could be replaced by the Pandana aggregate function
    # except 
    results = distance_network.nodes_df.copy()
    
    for cat in categories:
        x, y = (pois[pois['category'].isin(poi_dictionary[cat])]['geometry'].x, 
                    pois[pois['category'].isin(poi_dictionary[cat])]['geometry'].y)
        
        if len(x) == 0:
            print("Category "+ cat + " is empty")
            results[cat] = 0
            
        else:
            nodes = distance_network.get_node_ids(x_col=x, y_col=y, mapping_distance=100)
            distance_network.set(nodes, name=cat)

            access = distance_network.aggregate(
                distance=distance, type='count', decay='flat',
                name=cat)

            #number_accessible = max_items - (access == distance).sum(axis=1)

            results[cat] = access
            print("Finished category: " + cat)
            print(access.describe())

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


