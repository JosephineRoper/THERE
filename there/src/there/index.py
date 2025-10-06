# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 12:15:25 2021

@author: z3258367
All POI and network data should be in a projected CRS, except for transit which is always in WGS84.

"""

import numpy as np
import geopandas as gpd
import pandas as pd
import pandana as pdna
import math
import osmnx as ox

def cluster_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, 
                  return_no=5, chain_frac=1/3, loop=2, prev_results=None, trace=False):
    # Run the index once to give POIs an initial accessibility
    # if going to recurse further, need to not reset the index each time, it causes an error
    # think it's fixed by setting the index back to original again. Why not leave it as the node_id? I don't remember.
    # why am I even doing this? to attribute THERE to the pois correctly
    # this will run a number of times equal to loop, updating the POI THERE Index each time
    # prev_results is for convenience when testing clustering effects, reuse non-cluster results to reduce runtime
    if prev_results is not None:
        pois['node_id'] = distance_network.get_node_ids(pois['geometry'].x, pois['geometry'].y)
        pois = pois.reset_index().set_index('node_id').copy()
        pois['THERE'] = prev_results['THERE_Index']
        pois = pois.reset_index().set_index('index').copy()
        loop -= 1

    for i in range(loop):
        results = there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no, chain_frac)
        print("loop:", i)
        if i == 0:
            result_cache = results
        pois['node_id'] = distance_network.get_node_ids(pois['geometry'].x, pois['geometry'].y)
        pois = pois.reset_index().set_index('node_id').copy()
        pois['THERE'] = results['THERE_Index']
        if trace:
            pois['THERE_'+str(i)] = results['THERE_Index']
            result_cache['THERE_'+str(i)] = results['THERE_Index']
        pois = pois.reset_index().set_index('index').copy()

    return result_cache, pois

def transit_cluster_index(walk_network, integrated_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=1/3, loop=2, prev_results=None, trace=False):

    if pois.crs is not "EPSG:4326":
        pois = pois.to_crs("EPSG:4326")  # transit networks are always in WGS84 

    if prev_results is not None:
        pois['node_id'] = walk_network.get_node_ids(pois['geometry'].x, pois['geometry'].y)
        pois = pois.reset_index().set_index('node_id').copy()
        pois['THERE'] = prev_results['THERE_Index']
        pois = pois.reset_index().set_index('index').copy()
        loop -= 1

    for i in range(loop):
        results = transit_index(walk_network, integrated_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no, chain_frac)
        print("loop:", i)
        if i == 0:
            result_cache = results
        pois['node_id'] = walk_network.get_node_ids(pois['geometry'].x, pois['geometry'].y)
        pois = pois.reset_index().set_index('node_id').copy()
        pois['THERE'] = results['THERE_Index']
        if trace:
            pois['THERE_'+str(i)] = results['THERE_Index']
            result_cache['THERE_'+str(i)] = results['THERE_Index']
        pois = pois.reset_index().set_index('index').copy()

    return result_cache, pois

def extract_pois(poi_dictionary, pois, category):
    relevant_pois = gpd.GeoDataFrame()
    for key in poi_dictionary[category]:
        if key in pois:
            relevant_pois = pd.concat([relevant_pois, pois.loc[(pois[key].isin(poi_dictionary[category][key]))]])
            # some duplicates arise with this method of constructing relevant_pois,
            # because the same POI may be tagged "shopping:supermarket" and "building:supermarket"
            # in an OSM dataset for example. We are relying on the poi index to be an unique identifier.
            relevant_pois = relevant_pois[~relevant_pois.index.duplicated()]
    return relevant_pois

def _prepare_poi_attractiveness(relevant_pois, poi_variables, category, chain_frac):
    """Prepare POI attractiveness and chain weight values."""
    if poi_variables[category] == 'count':
        relevant_pois['attract'] = [1 for i in range(len(relevant_pois))]
    else: 
        # makes eg. a job count column the index column, so that 'include_poi_ids' returns job counts
        relevant_pois['attract'] = relevant_pois[poi_variables[category]]
            
    # check if pois already have a 'THERE' column
    if 'THERE' in relevant_pois.columns: 
        relevant_pois['chain_weight'] = (1-chain_frac) + chain_frac*relevant_pois['THERE']/100
    else:
        relevant_pois['chain_weight'] = 1 if chain_frac != 0 else (1-chain_frac)
    
    return relevant_pois

def _calculate_accessibility_score(access, relevant_pois, dist_const, dim_const, num_pois, distance):
    """Calculate the accessibility score using impedance, attractiveness, and diminishing returns."""
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

    dim = np.array([(1-np.exp(-dim_const*(x+1))) - (1-np.exp(-dim_const*x)) for x in attractiveness_sum])

    # Calculate the final accessibility score
    accessibility_score = (dim*attractiveness*impedance*there).sum(axis=1)
    
    return accessibility_score

def _process_category_standard(distance_network, pois, poi_dictionary, category, poi_variables, 
                              poi_gammas, poi_nums, poi_lambdas, distance, chain_frac):
    """Process a single category for standard (non-transit) networks."""
    relevant_pois = extract_pois(poi_dictionary, pois, category)
    
    if len(relevant_pois) == 0:
        return None, None
    
    relevant_pois = _prepare_poi_attractiveness(relevant_pois, poi_variables, category, chain_frac)
    
    # Set POIs on network
    x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)
    distance_network.set_pois(category, distance, poi_nums[category], x, y)
    
    # Get accessibility
    access = distance_network.nearest_pois(
        distance=distance, category=category, num_pois=poi_nums[category], include_poi_ids=True)
    
    # Calculate score
    accessibility_score = _calculate_accessibility_score(
        access, relevant_pois, poi_gammas[category], poi_lambdas[category], poi_nums[category], distance)
    
    return accessibility_score, access

def _process_category_transit(walk_network, integrated_network, pois, poi_dictionary, category, 
                            poi_variables, poi_gammas, poi_nums, poi_lambdas, distance, chain_frac):
    """Process a single category for transit networks."""
    relevant_pois = extract_pois(poi_dictionary, pois, category)
    
    if len(relevant_pois) == 0:
        return None, None
    
    relevant_pois = _prepare_poi_attractiveness(relevant_pois, poi_variables, category, chain_frac)
    
    # Map POIs to walk network nodes
    relevant_pois['node_id'] = walk_network.get_node_ids(relevant_pois['geometry'].x, relevant_pois['geometry'].y)
    relevant_pois.set_index('node_id', inplace=True, drop=False)
    
    # Set POIs on integrated network using walk network node locations
    integrated_network.set_pois(category=category, 
            maxdist=distance, 
            maxitems=poi_nums[category], 
            x_col=walk_network.nodes_df.loc[relevant_pois['node_id']].x, 
            y_col=walk_network.nodes_df.loc[relevant_pois['node_id']].y)
    
    # Find nearest POIs in category
    access = integrated_network.nearest_pois(
        distance=distance, category=category, num_pois=poi_nums[category], include_poi_ids=True)
    
    # Calculate score for category
    accessibility_score = _calculate_accessibility_score(
        access, relevant_pois, poi_gammas[category], poi_lambdas[category], poi_nums[category], distance)
    
    return accessibility_score, access
                   
def there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=0, walk_network=None):
    # return_no = 5 is the default setting, to return individual distance results
    # for maximum 5 closest points in each category. Enables some debugging &
    # visualisation options, but not returning an overly large results matrix
    # (as poi_nums searched may be 100s per category)
    
    results = distance_network.nodes_df.copy()

    total_weight = sum(poi_weights)
    
    for category in poi_weights.index:
        results = results.copy()  # to remove fragmentation warning
        
        weight = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(weight)))

        if category not in poi_dictionary:
            print("Category", category, "is not in the POI dictionary")
        else:
            accessibility_score, access = _process_category_standard(
                distance_network, pois, poi_dictionary, category, poi_variables, 
                poi_gammas, poi_nums, poi_lambdas, distance, chain_frac)

            if accessibility_score is None:
                print("No pois in category: "+ category)
                results[str(cat_name)] = 0
            else:
                # Store the variable score and weighted category score
                results[poi_variables[category]] = accessibility_score
                results[cat_name] = weight * accessibility_score
                
                # Store distance to closest destinations for debugging/visualization
                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                    
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
            
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary]   
    
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
    
    return results

def transit_index(walk_network, integrated_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=1/3):
    # differences compared to there_index: uses two networks in order to ascribe POIs only to walk network node, and converts to EPSG:4326 because UrbanAccess creates integrated networks direct from GTFS data in WGS84.

    results = integrated_network.nodes_df.copy()

    if pois.crs is not "EPSG:4326":
        pois = pois.to_crs("EPSG:4326")  # transit networks are always in WGS84 

    total_weight = sum(poi_weights)

    for category in poi_weights.index:
        results = results.copy()  # to remove fragmentation warning

        weight = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(weight)))
        
        if category not in poi_dictionary:
            print("Category", category, "is not in the POI dictionary")
        else:
            accessibility_score, access = _process_category_transit(
                walk_network, integrated_network, pois, poi_dictionary, category, 
                poi_variables, poi_gammas, poi_nums, poi_lambdas, distance, chain_frac)

            if accessibility_score is None:
                print("No pois in category: "+ category)
                results[str(cat_name)] = 0
            else:
                # Store the variable score and weighted category score
                results[poi_variables[category]] = accessibility_score
                results[cat_name] = weight * accessibility_score

                # Store distance to closest destinations for debugging/visualization
                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                                                
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
                
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary]   
        
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))

    return results

def test_transit_index(walk_network, integrated_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=1/3):

    results = integrated_network.nodes_df.copy()

    total_weight = sum(poi_weights)

    for category in poi_weights.index:
        results = results.copy()  # to remove fragmentation warning

        weight = poi_weights[category]
        cat_name = ''.join((str(category),"_",str(weight)))
        
        if category not in poi_dictionary:
            print("Category", category, "is not in the POI dictionary")
        else:
            # Use the same logic as transit_index but also extract intermediate values for debugging
            relevant_pois = extract_pois(poi_dictionary, pois, category)
            
            if len(relevant_pois) == 0:
                print("No pois in category: "+ category)
                results[str(cat_name)] = 0
                # Set debug variables to None
                access, impedance, attractiveness, prev_there, attractiveness_sum, dim = [None] * 6
            else:
                relevant_pois = _prepare_poi_attractiveness(relevant_pois, poi_variables, category, chain_frac)
                
                # Map POIs to walk network nodes
                relevant_pois['node_id'] = walk_network.get_node_ids(relevant_pois['geometry'].x, relevant_pois['geometry'].y)
                relevant_pois.set_index('node_id', inplace=True, drop=False)
                
                # Set POIs on integrated network using walk network node locations
                integrated_network.set_pois(category=category, 
                        maxdist=distance, 
                        maxitems=poi_nums[category], 
                        x_col=walk_network.nodes_df.loc[relevant_pois['node_id']].x, 
                        y_col=walk_network.nodes_df.loc[relevant_pois['node_id']].y)
                
                # Get accessibility
                access = integrated_network.nearest_pois(
                    distance=distance, category=category, num_pois=poi_nums[category], include_poi_ids=True)
                
                # Calculate intermediate values for debugging (copied from _calculate_accessibility_score)
                impedance = np.exp(-poi_gammas[category]*access.iloc[:,0:poi_nums[category]])
                impedance[access.iloc[:,0:poi_nums[category]] == distance] = 0

                attract_dict = pd.Series(relevant_pois.attract.values,index=relevant_pois.index).to_dict()
                there_dict = pd.Series(relevant_pois.chain_weight.values,index=relevant_pois.index).to_dict()

                attract = access.iloc[:,poi_nums[category]:2*poi_nums[category]].copy().applymap(lambda y: attract_dict.get(y,y))
                there_weight = access.iloc[:,poi_nums[category]:2*poi_nums[category]].copy().applymap(lambda y: there_dict.get(y,y))
                
                attractiveness = np.nan_to_num(attract.values, nan=0.0)
                prev_there = np.nan_to_num(there_weight.values, nan=0.0)
                attractiveness_sum = attractiveness.cumsum(axis=1).astype(np.int64)

                dim = np.array([(1-np.exp(-poi_lambdas[category]*(x+1))) - (1-np.exp(-poi_lambdas[category]*x)) for x in attractiveness_sum])

                # Calculate final score
                accessibility_score = (dim*attractiveness*impedance*prev_there).sum(axis=1)
                
                results[poi_variables[category]] = accessibility_score
                results[cat_name] = weight * accessibility_score

                # Store distance to closest destinations for debugging/visualization
                for i in range(return_no):
                    col_name = ''.join((str(category),str(i+1)))
                    results[col_name] = access[i+1]
                                                
            print("Finished category: " + category)
            print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
                
    col_list = [''.join((str(category),"_",str(poi_weights[category])))
                for category in poi_dictionary if category in poi_weights.index]   
        
    results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))

    return results, access, impedance, attractiveness, prev_there, attractiveness_sum, dim, relevant_pois

def point_index(x,y, distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no=5):
    # this is a version of the index that can be used for a single point
    location = gpd.points_from_xy(x,y)
    
    
    results = there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, return_no)
    
    return results
