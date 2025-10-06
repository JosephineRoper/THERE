# -*- coding: utf-8 -*-
"""
Object-oriented refactoring of the WalkTHERE index calculation

This is Option 3: Class-based approach with composition
"""

import numpy as np
import geopandas as gpd
import pandas as pd
import pandana as pdna
import math
import osmnx as ox


class POIExtractor:
    """Handles POI extraction and preparation."""
    
    @staticmethod
    def extract_pois(poi_dictionary, pois, category):
        """Extract relevant POIs for a given category."""
        relevant_pois = gpd.GeoDataFrame()
        for key in poi_dictionary[category]:
            if key in pois:
                relevant_pois = pd.concat([relevant_pois, pois.loc[(pois[key].isin(poi_dictionary[category][key]))]])
                # Remove duplicates that may arise from multiple OSM tags
                relevant_pois = relevant_pois[~relevant_pois.index.duplicated()]
        return relevant_pois
    
    @staticmethod
    def prepare_attractiveness(relevant_pois, poi_variables, category, chain_frac):
        """Prepare POI attractiveness and chain weight values."""
        if poi_variables[category] == 'count':
            relevant_pois['attract'] = [1 for i in range(len(relevant_pois))]
        else: 
            relevant_pois['attract'] = relevant_pois[poi_variables[category]]
                
        # Calculate chain weights
        if 'THERE' in relevant_pois.columns: 
            relevant_pois['chain_weight'] = (1-chain_frac) + chain_frac*relevant_pois['THERE']/100
        else:
            relevant_pois['chain_weight'] = 1 if chain_frac != 0 else (1-chain_frac)
        
        return relevant_pois


class AccessibilityCalculator:
    """Handles the core accessibility calculations."""
    
    @staticmethod
    def calculate_impedance(access, dist_const, num_pois, distance):
        """Calculate distance decay impedance."""
        impedance = np.exp(-dist_const * access.iloc[:, 0:num_pois])
        impedance[access.iloc[:, 0:num_pois] == distance] = 0
        return impedance
    
    @staticmethod
    def prepare_attraction_matrices(access, relevant_pois, num_pois):
        """Prepare attractiveness and chain weight matrices."""
        attract_dict = pd.Series(relevant_pois.attract.values, index=relevant_pois.index).to_dict()
        there_dict = pd.Series(relevant_pois.chain_weight.values, index=relevant_pois.index).to_dict()

        attract = access.iloc[:, num_pois:2*num_pois].copy().applymap(lambda y: attract_dict.get(y, y))
        there_weight = access.iloc[:, num_pois:2*num_pois].copy().applymap(lambda y: there_dict.get(y, y))
        
        attractiveness = np.nan_to_num(attract.values, nan=0.0)
        there = np.nan_to_num(there_weight.values, nan=0.0)
        
        return attractiveness, there
    
    @staticmethod
    def calculate_diminishing_returns(attractiveness, dim_const):
        """Calculate diminishing returns based on cumulative attractiveness."""
        attractiveness_sum = attractiveness.cumsum(axis=1).astype(np.int64)
        dim = np.array([(1-np.exp(-dim_const*(x+1))) - (1-np.exp(-dim_const*x)) for x in attractiveness_sum])
        return dim
    
    @classmethod
    def calculate_score(cls, access, relevant_pois, dist_const, dim_const, num_pois, distance):
        """Calculate the complete accessibility score."""
        impedance = cls.calculate_impedance(access, dist_const, num_pois, distance)
        attractiveness, there = cls.prepare_attraction_matrices(access, relevant_pois, num_pois)
        dim = cls.calculate_diminishing_returns(attractiveness, dim_const)
        
        accessibility_score = (dim * attractiveness * impedance * there).sum(axis=1)
        return accessibility_score


class NetworkStrategy:
    """Base class for different network processing strategies."""
    
    def set_pois_on_network(self, network, category, relevant_pois, poi_nums, distance):
        """Set POIs on the network. To be implemented by subclasses."""
        raise NotImplementedError
    
    def get_nearest_pois(self, network, category, poi_nums, distance):
        """Get nearest POIs from the network. To be implemented by subclasses."""
        raise NotImplementedError


class StandardNetworkStrategy(NetworkStrategy):
    """Strategy for standard (non-transit) networks."""
    
    def set_pois_on_network(self, network, category, relevant_pois, poi_nums, distance):
        """Set POIs directly on the network using their coordinates."""
        x, y = (relevant_pois['geometry'].x, relevant_pois['geometry'].y)
        network.set_pois(category, distance, poi_nums, x, y)
    
    def get_nearest_pois(self, network, category, poi_nums, distance):
        """Get nearest POIs from the network."""
        return network.nearest_pois(
            distance=distance, category=category, num_pois=poi_nums, include_poi_ids=True)


class TransitNetworkStrategy(NetworkStrategy):
    """Strategy for transit networks."""
    
    def __init__(self, walk_network):
        self.walk_network = walk_network
    
    def set_pois_on_network(self, integrated_network, category, relevant_pois, poi_nums, distance):
        """Set POIs on integrated network using walk network node locations."""
        # Map POIs to walk network nodes
        relevant_pois['node_id'] = self.walk_network.get_node_ids(
            relevant_pois['geometry'].x, relevant_pois['geometry'].y)
        relevant_pois.set_index('node_id', inplace=True, drop=False)
        
        # Set POIs on integrated network using walk network node locations
        integrated_network.set_pois(
            category=category, 
            maxdist=distance, 
            maxitems=poi_nums, 
            x_col=self.walk_network.nodes_df.loc[relevant_pois['node_id']].x, 
            y_col=self.walk_network.nodes_df.loc[relevant_pois['node_id']].y)
    
    def get_nearest_pois(self, integrated_network, category, poi_nums, distance):
        """Get nearest POIs from the integrated network."""
        return integrated_network.nearest_pois(
            distance=distance, category=category, num_pois=poi_nums, include_poi_ids=True)


class THEREIndexCalculator:
    """Main class for calculating THERE indices."""
    
    def __init__(self, network_strategy):
        self.network_strategy = network_strategy
        self.poi_extractor = POIExtractor()
        self.accessibility_calculator = AccessibilityCalculator()
    
    def process_category(self, network, pois, poi_dictionary, category, poi_variables, 
                        poi_gammas, poi_nums, poi_lambdas, distance, chain_frac):
        """Process a single POI category."""
        # Extract relevant POIs
        relevant_pois = self.poi_extractor.extract_pois(poi_dictionary, pois, category)
        
        if len(relevant_pois) == 0:
            return None, None
        
        # Prepare attractiveness values
        relevant_pois = self.poi_extractor.prepare_attractiveness(
            relevant_pois, poi_variables, category, chain_frac)
        
        # Set POIs on network
        self.network_strategy.set_pois_on_network(
            network, category, relevant_pois, poi_nums[category], distance)
        
        # Get accessibility
        access = self.network_strategy.get_nearest_pois(
            network, category, poi_nums[category], distance)
        
        # Calculate accessibility score
        accessibility_score = self.accessibility_calculator.calculate_score(
            access, relevant_pois, poi_gammas[category], poi_lambdas[category], 
            poi_nums[category], distance)
        
        return accessibility_score, access
    
    def calculate_index(self, network, pois, poi_dictionary, poi_weights, poi_gammas, 
                       poi_nums, poi_lambdas, poi_variables, distance, return_no=5, 
                       chain_frac=0, pois_crs_conversion=None):
        """Calculate the complete THERE index."""
        
        # Handle CRS conversion if needed (for transit networks)
        if pois_crs_conversion is not None:
            if pois.crs is not pois_crs_conversion:
                pois = pois.to_crs(pois_crs_conversion)
        
        results = network.nodes_df.copy()
        total_weight = sum(poi_weights)
        
        for category in poi_weights.index:
            results = results.copy()  # to remove fragmentation warning
            
            weight = poi_weights[category]
            cat_name = ''.join((str(category), "_", str(weight)))

            if category not in poi_dictionary:
                print("Category", category, "is not in the POI dictionary")
                results[str(cat_name)] = 0
            else:
                accessibility_score, access = self.process_category(
                    network, pois, poi_dictionary, category, poi_variables, 
                    poi_gammas, poi_nums, poi_lambdas, distance, chain_frac)

                if accessibility_score is None:
                    print("No pois in category: " + category)
                    results[str(cat_name)] = 0
                else:
                    # Store results
                    results[poi_variables[category]] = accessibility_score
                    results[cat_name] = weight * accessibility_score
                    
                    # Store distance to closest destinations for debugging/visualization
                    for i in range(return_no):
                        col_name = ''.join((str(category), str(i+1)))
                        results[col_name] = access[i+1]
                        
                print("Finished category: " + category)
                print("Maximum score: " + str(max(results[cat_name])) + " out of " + str(weight))
                
        # Calculate final index
        col_list = [''.join((str(category), "_", str(poi_weights[category])))
                    for category in poi_dictionary]   
        
        results['THERE_Index'] = 100/total_weight*(results[col_list].sum(axis=1))
        
        return results


class ClusterIndexCalculator:
    """Handles clustering calculations that iterate the index multiple times."""
    
    def __init__(self, index_calculator, network_for_poi_mapping):
        self.index_calculator = index_calculator
        self.network_for_poi_mapping = network_for_poi_mapping
    
    def calculate(self, network, pois, poi_dictionary, poi_weights, poi_gammas, 
                 poi_nums, poi_lambdas, poi_variables, distance, return_no=5, 
                 chain_frac=1/3, loop=2, prev_results=None, trace=False, 
                 pois_crs_conversion=None):
        """Calculate clustered index with multiple iterations."""
        
        # Handle previous results
        if prev_results is not None:
            pois['node_id'] = self.network_for_poi_mapping.get_node_ids(pois['geometry'].x, pois['geometry'].y)
            pois = pois.reset_index().set_index('node_id').copy()
            pois['THERE'] = prev_results['THERE_Index']
            pois = pois.reset_index().set_index('index').copy()
            loop -= 1

        for i in range(loop):
            results = self.index_calculator.calculate_index(
                network, pois, poi_dictionary, poi_weights, poi_gammas, 
                poi_nums, poi_lambdas, poi_variables, distance, return_no, 
                chain_frac, pois_crs_conversion)
            
            print("loop:", i)
            if i == 0:
                result_cache = results
                
            # Update POI THERE values for next iteration
            pois['node_id'] = self.network_for_poi_mapping.get_node_ids(pois['geometry'].x, pois['geometry'].y)
            pois = pois.reset_index().set_index('node_id').copy()
            pois['THERE'] = results['THERE_Index']
            
            if trace:
                pois['THERE_'+str(i)] = results['THERE_Index']
                result_cache['THERE_'+str(i)] = results['THERE_Index']
                
            pois = pois.reset_index().set_index('index').copy()

        return result_cache, pois


# Factory functions to create the different index calculators
def create_standard_index_calculator():
    """Create a calculator for standard (non-transit) networks."""
    strategy = StandardNetworkStrategy()
    return THEREIndexCalculator(strategy)


def create_transit_index_calculator(walk_network):
    """Create a calculator for transit networks."""
    strategy = TransitNetworkStrategy(walk_network)
    return THEREIndexCalculator(strategy)


# Wrapper functions that maintain the original API
def there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, 
               poi_nums, poi_lambdas, poi_variables, distance, return_no=5, chain_frac=0):
    """Calculate THERE index using standard network."""
    calculator = create_standard_index_calculator()
    return calculator.calculate_index(
        distance_network, pois, poi_dictionary, poi_weights, poi_gammas, 
        poi_nums, poi_lambdas, poi_variables, distance, return_no, chain_frac)


def transit_index(walk_network, integrated_network, pois, poi_dictionary, poi_weights, 
                 poi_gammas, poi_nums, poi_lambdas, poi_variables, distance, 
                 return_no=5, chain_frac=1/3):
    """Calculate THERE index using transit network."""
    calculator = create_transit_index_calculator(walk_network)
    return calculator.calculate_index(
        integrated_network, pois, poi_dictionary, poi_weights, poi_gammas, 
        poi_nums, poi_lambdas, poi_variables, distance, return_no, chain_frac, 
        pois_crs_conversion="EPSG:4326")


def cluster_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, 
                 poi_nums, poi_lambdas, poi_variables, distance, return_no=5, 
                 chain_frac=1/3, loop=2, prev_results=None, trace=False):
    """Calculate clustered THERE index."""
    index_calculator = create_standard_index_calculator()
    cluster_calculator = ClusterIndexCalculator(index_calculator, distance_network)
    return cluster_calculator.calculate(
        distance_network, pois, poi_dictionary, poi_weights, poi_gammas, 
        poi_nums, poi_lambdas, poi_variables, distance, return_no, 
        chain_frac, loop, prev_results, trace)


def transit_cluster_index(walk_network, integrated_network, pois, poi_dictionary, 
                         poi_weights, poi_gammas, poi_nums, poi_lambdas, poi_variables, 
                         distance, return_no=5, chain_frac=1/3, loop=2, prev_results=None, 
                         trace=False):
    """Calculate clustered THERE index for transit networks."""
    index_calculator = create_transit_index_calculator(walk_network)
    cluster_calculator = ClusterIndexCalculator(index_calculator, walk_network)
    return cluster_calculator.calculate(
        integrated_network, pois, poi_dictionary, poi_weights, poi_gammas, 
        poi_nums, poi_lambdas, poi_variables, distance, return_no, 
        chain_frac, loop, prev_results, trace, pois_crs_conversion="EPSG:4326")


def point_index(x, y, distance_network, pois, poi_dictionary, poi_weights, poi_gammas, 
               poi_nums, poi_lambdas, poi_variables, distance, return_no=5):
    """Calculate THERE index for a single point."""
    location = gpd.points_from_xy(x, y)
    return there_index(distance_network, pois, poi_dictionary, poi_weights, poi_gammas, 
                      poi_nums, poi_lambdas, poi_variables, distance, return_no)
