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
        # I'm not sure why I'm doing this instead of using geometries_from_polygon
        bbox = place.bounds
        bbox_pois = ox.geometries.geometries_from_bbox(bbox['maxy'][0], bbox['miny'][0], bbox['maxx'][0], bbox['minx'][0], tags)
        gdf = gpd.clip(bbox_pois, place, keep_geom_type=False).to_crs(proj_crs)
        #this doesn't work with multiple polygons, it dissolves them
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

def cluster_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=1/3, loop=2):
    # Run the index once to give POIs an initial accessibility
    # if going to recurse further, need to not reset the index each time, it causes an error
    # think it's fixed by setting the index back to original again
    # why am I even doing this? to attribute THERE to the pois correctly
    print(chain_frac)
    # this will run a number of times equal to loop, updating the POI THERE Index each time
    # so far, I haven't tested what difference multiple rounds makes.
    for i in range(loop):
        poi_results = there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no, chain_frac)
        pois['node_id'] = distance_network.get_node_ids(pois['geometry'].x, pois['geometry'].y)
        pois = pois.reset_index().set_index('node_id').copy()
        pois['THERE'] = poi_results['THERE_Index']
        pois = pois.reset_index().set_index('index').copy()
        print(pois['THERE'])

    # having this out of the loop means can skip the final assignment of pois['THERE'], but it seems cleaner to have it all inside
    #results = there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, #poi_variables, distance, return_no, chain_frac)
    return poi_results
                   
def there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=1/3):
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
                    # some duplicates arise with this method of constructing relevant_pois,
                    # because the same POI may be tagged "shopping:supermarket" and "building:supermarket"
                    # in an OSM dataset for example. We are relying on the poi index to be an unique identifier.
                    relevant_pois = relevant_pois[~relevant_pois.index.duplicated()]

            if len(relevant_pois) == 0:
                print("No pois in category: "+ category)
                results[str(cat_name)] = 0

            else:
                if poi_variables[category] == 'count':
                        relevant_pois['attract'] = [1 for i in range(len(relevant_pois))]
                else: 
                    # makes eg. a job count column the index column, so that 'include_poi_ids' returns job counts
                    relevant_pois['attract'] = relevant_pois[poi_variables[category]]
                        
                # check if relevant_pois has a 'THERE' column
                if 'THERE' in relevant_pois.columns:
                #try: # if 'THERE' exists, 
                    relevant_pois['chain_weight'] = (1-chain_frac) + chain_frac*relevant_pois['THERE']/100
                #except:
                else:
                    relevant_pois['chain_weight'] = 1
              
                x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)
                distance_network.set_pois(category, distance, num_pois, x, y)
                
                access = distance_network.nearest_pois(
                    distance=distance, category=category, num_pois=num_pois, include_poi_ids=True)

                impedance = np.exp(-dist_const*access.iloc[:,0:num_pois])
                impedance[access.iloc[:,0:num_pois] == distance] = 0

                attract_dict = pd.Series(relevant_pois.attract.values,index=relevant_pois.index).to_dict()
                there_dict = pd.Series(relevant_pois.chain_weight.values,index=relevant_pois.index).to_dict()

                attract = access.iloc[:,num_pois:2*num_pois].copy().applymap(lambda y: attract_dict.get(y,y))
                there_weight = access.iloc[:,num_pois:2*num_pois].copy().applymap(lambda y: there_dict.get(y,y))

                # nan to 0 is necessary because where eg. for some origin point, only 4 POIs of a certain category
                # are found (within distance), the returned poi_ids will be NaN for the remaining num_pois columns
                attractiveness = np.nan_to_num(attract.values, nan=0.0)
                there = np.nan_to_num(there_weight.values, nan=0.0)
                attractiveness_sum = attractiveness.cumsum(axis=1).astype(np.int64)

                #max_opps = np.max(attractiveness_sum)
                #dim = np.array([(1-np.exp(-dim_const*(x+1))) - (1-np.exp(-dim_const*x)) for x in range(max_opps+1)])
                dim = np.array([(1-np.exp(-dim_const*(x+1))) - (1-np.exp(-dim_const*x)) for x in attractiveness_sum])

                # this is no longer useful/meaningful as an output column. could output attract_sum instead
                results[poi_variables[category]] = (dim*attractiveness*impedance*there).sum(axis=1)
                
                results[cat_name] = weight*results[poi_variables[category]]

                #results[cat_name] = weight*(impedance*dim_opp_weight).sum(axis=1)
                
                # this provides columns with the distance of the return_no closest destinations in the category
                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                    
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary]   
    
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
    
    return results

def there_index_transit(distance_network, walk_node_ids,
                        pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5):
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

        dim_opp_weight = [(1-np.exp(-dim_const*x)) - (1-np.exp(-dim_const*(x-1))) for x in range(1,num_pois+1)]

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

                impedance = np.exp(-dist_const*access.iloc[:,0:num_pois])
                impedance[access.iloc[:,0:num_pois] == distance] = 0

                results[cat_name] = weight*(impedance*dim_opp_weight).sum(axis=1)
                
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

                impedance = np.exp(-dist_const*access.iloc[:,0:num_pois])
                impedance[access.iloc[:,0:num_pois] == distance] = 0
                
                attractiveness = np.nan_to_num(access.iloc[:,num_pois:2*num_pois].values, nan=0.0)
                attractiveness_sum = attractiveness.cumsum(axis=1).astype(np.int64)

                max_opps = np.max(attractiveness_sum)
                dim = np.array([(1-np.exp(-dim_const*(x+1))) - (1-np.exp(-dim_const*x)) for x in range(max_opps+1)])

                # this is no longer useful/meaningful as an output column. could output attract_sum instead
                results[poi_variables[category]] = (dim[attractiveness_sum]*attractiveness*impedance).sum(axis=1)
                
                results[cat_name] = weight*results[poi_variables[category]]

                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                    
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary]   
    
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
    
    return results

