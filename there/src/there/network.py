import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
import pandana as pdna
import urbanaccess as ua
from urbanaccess.config import settings
from urbanaccess.gtfsfeeds import feeds
from urbanaccess import gtfsfeeds
from urbanaccess.gtfs.gtfsfeeds_dataframe import gtfsfeeds_dfs
from urbanaccess.network import ua_network, load_network

def transit_cost(edges):
    # adjusting for differential value of time for in-vehicle time vs waiting vs walking. 
    # Current values based on Wardman 2004
    # note UA default headway is half headway, thus using 1.2 instead of 0.6
    type_conds = [edges['net_type'] == 'walk',
                  edges['net_type'] == 'transit',
                 edges['net_type'] == 'transit to osm',
                 edges['net_type'] == 'osm to transit']
    type_choices = [2, 1, 2, 1.2]
    edges['adjustment'] = np.select(type_conds, type_choices, default = 1)
    return edges

def cycle_cost(edges):    
    def int_50(n):
    # some values of OSM tag maxspeed are unwieldy, eg text containing kmh, lists (produced by OSMNx simplification) etc. 
    # This tries to strip text, but for lists, NaN etc gives up and returns the default speed in NSW - 50km/hr.
    # Note that OSMNx simplification already collapses lists when values are the same.
    # The most common case of lists is speed limit changes on arterial roads, eg the list is [60,80],  
    # which is all equivalent to >50 for the purposes of this bicycle network anyway.
        try:
            return int(n)
        except:
            try:
                return int(filter(lambda x: x in '0123456789.', n))
            except:
                return 50
        
    edges['inferred_speed'] = edges['maxspeed'].apply(int_50)

    type_conds = [edges['bicycle'] == 'no',
                 edges['cycleway'] == 'track',
                 edges['cycleway'] == 'lane',
                (edges['highway'] == 'cycleway') & (edges['bicycle'] == 'designated') & (edges['foot'] == 'designated'),
                 edges['cycleway'].isin(['shared_path', 'crossing']),
                 edges['cycleway'] == 'shared_busway',
                 (edges['cycleway'].isin(['shared_lane', 'doorzone'])) | (edges['bicycle']=='yes')]
    type_choices = ['prohibited',
                   'track',
                   'lane',
                   'shared_path', 'shared_path',
                   'shared_busway',
                   'sharrow']
    edges['cycleway_type'] = np.select(type_conds, type_choices, default = 'no_cycleway')

    impedance_conds = [edges['cycleway_type'] == 'prohibited',
                edges['cycleway_type'] == 'track',
                edges['cycleway_type'] == 'shared_path',
                edges['cycleway_type'] == 'shared_busway',
                edges['inferred_speed'] <= 30,
                (edges['highway'].isin(['residential','unclassified'])) & (edges['cycleway_type'] == 'lane'),
                (edges['highway'].isin(['residential','unclassified'])) & (edges['cycleway_type'] == 'sharrow'),
                (edges['highway'].isin(['residential','unclassified'])) & (edges['cycleway_type'] == 'no_cycleway'),
                (edges['highway'] == 'tertiary') & (edges['cycleway'] == 'lane'),
                (edges['highway'] == 'tertiary') & (edges['cycleway'] == 'sharrow'),
                (edges['highway'] == 'tertiary') & (edges['cycleway'] == 'no_cycleway'),
                (edges['highway'] == 'secondary') & (edges['cycleway'] == 'lane'),
                (edges['highway'] == 'secondary') & (edges['cycleway'] == 'sharrow'),
                (edges['highway'] == 'secondary') & (edges['cycleway'] == 'no_cycleway'),
                (edges['highway'] == 'primary') & (edges['cycleway'] == 'lane'),
                (edges['highway'] == 'primary') & (edges['cycleway'] == 'sharrow'),
                (edges['highway'] == 'primary') & (edges['cycleway'] == 'no_cycleway'),]
    impedance_choice = [5, 1, 1.5, 1.5, 1.1, 
                    1.1, 1.2, 1.3, # residential and unclassified
                    2, 2.5, 3, # tertiary
                    3, 3.5, 4, # secondary
                    4, 4.5, 5] # primary
    edges['cycle_impedance'] = np.select(impedance_conds, impedance_choice, default = 5)

    return edges['length'] * edges['cycle_impedance']

def two_way_edges(edges):
    # 2017 David Wasserman via Pandana github forum
    # copy all two way edges and duplicate with reversed to and from nodes
    edges = edges.reset_index()
    twoway_edges = edges.loc[(edges['oneway'] == False)|(edges['oneway'] == 'False')].copy()
    if len(twoway_edges) < 1:
        print("No two-way edges found")
    twoway_u = twoway_edges["u"].copy()
    twoway_edges["u"] = twoway_edges["v"]
    twoway_edges["v"] = twoway_u
    return pd.concat([edges, twoway_edges])

def prepare_cycle_net(place_gdf, proj_crs):
    # To create a 'cycling friendly' network it is useful to have some additional tags
    ox.settings.useful_tags_way = list(set(ox.settings.useful_tags_way + 
                                           ['bicycle'] + ['foot'] + ['cycleway'] + ['cycleway:left'] + ['cycleway:right']))
    
    G = ox.graph.graph_from_polygon(place_gdf.to_crs('EPSG:4326').geometry[0], network_type='bike', simplify=False)
    # because we use attributes of the edges for the cycling impedance, it is better to avoid combining incident
    # edges that come from different OSM ways. Thus simplify separately with strict=False. 
    G = ox.simplify_graph(G, strict=False)

    # Get nodes and edges as geodataframes (gdfs) from OSMNX network
    graph_df = ox.graph_to_gdfs(G)
    nodes_gdfs = graph_df[0].to_crs(proj_crs)
    edges_gdfs = graph_df[1].to_crs(proj_crs)

    #For walking it's ok to assume all edges are two-way, but for cycling and driving it's not. So we need to create Pandana networks with twoway=False, but first, duplicate & flip the OSM edges that are two-way.
    edges_gdfs = two_way_edges(edges_gdfs)

    edges_gdfs['cycleway'] = edges_gdfs['cycleway'].fillna(edges_gdfs['cycleway:left'])
    edges_gdfs['cycleway'] = edges_gdfs['cycleway'].fillna(edges_gdfs['cycleway:right'])

    # Setting indices of Edges gdfs to match expected dataframe for Pandana
    edges_gdfs['from_idx'] = edges_gdfs['u']
    edges_gdfs['to_idx'] = edges_gdfs['v']
    edges_gdfs= edges_gdfs.set_index(['from_idx', 'to_idx'])
    edges_gdfs.index.names= ['','']

    # Setting indices of Nodes gdfs to match expected dataframe for Pandana
    nodes_gdfs.index.name = 'id'

    edges_gdfs['cycle_cost'] = cycle_cost(edges_gdfs)

    cycle_network_imp = pdna.Network(nodes_gdfs.geometry.x, nodes_gdfs.geometry.y, 
                                edges_gdfs['u'], edges_gdfs['v'], 
                                edge_weights=edges_gdfs[['cycle_cost']],                           
                                twoway=False)
    return cycle_network_imp

def prepare_walk_net(place_gdf, proj_crs):
    # downloads a network from OSM within the polygon of the place_gdf
    # and prepares a pandana network for walking
    G = ox.graph.graph_from_polygon(place_gdf.to_crs('EPSG:4326').geometry[0], network_type='walk')

    # Get nodes and edges as geodataframes (gdfs) from OSMNX network
    graph_df = ox.graph_to_gdfs(G)
    nodes_gdfs = graph_df[0]
    edges_gdfs = graph_df[1]

    edges_gdfs = edges_gdfs.to_crs(proj_crs)
    nodes_gdfs = nodes_gdfs.to_crs(proj_crs)

    # with new OSMnx graph from polygon seems to be different
    edges_gdfs = edges_gdfs.reset_index()
    # Setting indices of Edges gdfs to match expected dataframe for Pandana
    edges_gdfs['from_idx'] = edges_gdfs['u']
    edges_gdfs['to_idx'] = edges_gdfs['v']
    # Pandana expects edges to have a two item index based on the same IDs as the node index. 
    # (with thanks to https://github.com/shriv/accessibility-series/blob/master/Accounting%20for%20hills%20in%20accessibility%20analyses.ipynb)
    edges_gdfs= edges_gdfs.set_index(['from_idx', 'to_idx'])
    edges_gdfs.index.names= ['','']

    # Setting indices of Nodes gdfs to match expected dataframe for Pandana
    nodes_gdfs.index.name = 'id'
    # Create a pandana network with data extracted from an OSMNX graph
    distance_network = pdna.Network(nodes_gdfs.geometry.x, nodes_gdfs.geometry.y,
                                    edges_gdfs['u'], edges_gdfs['v'], 
                                    edges_gdfs[['length']])
    return distance_network

def prepare_drive_net(place_gdf, proj_crs):
    G = ox.graph.graph_from_polygon(place_gdf.to_crs('EPSG:4326').geometry[0], network_type='drive')

    # Get nodes and edges as geodataframes (gdfs) from OSMNX network
    graph_df = ox.graph_to_gdfs(G)
    nodes_gdfs = graph_df[0]
    edges_gdfs = graph_df[1]

    edges_gdfs = edges_gdfs.to_crs(proj_crs)
    nodes_gdfs = nodes_gdfs.to_crs(proj_crs)

    # with new OSMnx graph from polygon seems to be different
    edges_gdfs = edges_gdfs.reset_index()
    # Setting indices of Edges gdfs to match expected dataframe for Pandana
    edges_gdfs['from_idx'] = edges_gdfs['u']
    edges_gdfs['to_idx'] = edges_gdfs['v']
    edges_gdfs= edges_gdfs.set_index(['from_idx', 'to_idx'])
    edges_gdfs.index.names= ['','']

    # Setting indices of Nodes gdfs to match expected dataframe for Pandana
    nodes_gdfs.index.name = 'id'
    # Create a pandana network with data extracted from an OSMNX graph
    distance_network = pdna.Network(nodes_gdfs.geometry.x, nodes_gdfs.geometry.y,
                                    edges_gdfs['u'], edges_gdfs['v'], 
                                    edges_gdfs[['length']])
    return distance_network

def prepare_transit_net(place_gdf, proj_crs, gtfs_folder, walk_network=None, day='monday', timerange=['07:00:00', '12:00:00']):
    # need to return ua_net as much more flexible.
    bbox = place_gdf.to_crs("EPSG:4326").bounds
    bbox_points = (bbox['minx'][0], bbox['miny'][0], bbox['maxx'][0], bbox['maxy'][0])

    loaded_feeds = ua.gtfs.load.gtfsfeed_to_df(gtfsfeed_path=(gtfs_folder),
                                                validation=True, 
                                                verbose=True, 
                                                bbox=bbox_points, 
                                                remove_stops_outsidebbox=True, 
                                                append_definitions=True)
    
    ua.gtfs.network.create_transit_net(gtfsfeeds_dfs=loaded_feeds,
                                       day=day,
                                       timerange=timerange,
                                       calendar_dates_lookup=None)
    
    urbanaccess_net = ua.network.ua_network

    ua.gtfs.headways.headways(gtfsfeeds_df=loaded_feeds,
                              headway_timerange=timerange)
    
    if walk_network:
        # possibility to supply a previously created Pandana network based on osmnx (or any other loading method)
        # which retains a better street network than default OSMnet loader, in my opinion
        # note that edges should have already been flipped and appended in the previous Pandana walk_network, by building the
        # network with twoway=True (default). The OSMnet loader also does this.
        # But, I use walk network in projected CRS and UA networks need to be lat-lon, so have to convert.
        nodes = gpd.GeoDataFrame(walk_network.nodes_df, geometry=gpd.GeoSeries.from_xy(walk_network.nodes_df.x, walk_network.nodes_df.y, crs=proj_crs)).to_crs("EPSG:4326")
        nodes['x'] = nodes.geometry.x
        nodes['y'] = nodes.geometry.y
        edges = walk_network.edges_df

        edges['distance'] = edges['length']
        ua.osm.network.create_osm_net(osm_edges=edges, 
                                      osm_nodes=nodes, 
                                      travel_speed_mph=1.875)

        # for some reason this prevents a memoryerror in the network integration step
        urbanaccess_net.osm_nodes['id'] = urbanaccess_net.osm_nodes.index  

    else:
        nodes, edges = ua.osm.load.ua_network_from_bbox(bbox=bbox, 
                                                        remove_lcn=True)
        ua.osm.network.create_osm_net(osm_edges=edges, 
                                      osm_nodes=nodes, 
                                      travel_speed_mph=1.875)
    
    ua.network.integrate_network(urbanaccess_network=urbanaccess_net,
                                 headways=True,
                                 urbanaccess_gtfsfeeds_df=loaded_feeds,
                                 headway_statistic='mean')
    return urbanaccess_net

def pandana_transit_net(urbanaccess_net, place_gdf):
    #place = place_gdf.to_crs("EPSG:4326").geometry[0]
    nodes_in_place = urbanaccess_net.net_nodes #gpd.sjoin(urbanaccess_net.net_nodes, place, how='inner', op='intersects')    
    transit_ped_net = pdna.Network(nodes_in_place["x"],
                                    nodes_in_place["y"],
                                    urbanaccess_net.net_edges["from_int"],
                                    urbanaccess_net.net_edges["to_int"],
                                    urbanaccess_net.net_edges[["weight"]], 
                                    twoway=False)
    
    osm_node_ids = urbanaccess_net.net_nodes.loc[urbanaccess_net.net_nodes.net_type == 'walk'].index.values
    osm_nodes = urbanaccess_net.net_nodes.reindex(osm_node_ids)  # filters for rows with index in list of values
    osm_edges = urbanaccess_net.net_edges.loc[urbanaccess_net.net_edges.from_int.isin(osm_node_ids) &
                                    urbanaccess_net.net_edges.to_int.isin(osm_node_ids)]

    # outstanding question: can this be replaced by the original walk network? Not very simply because of the CRS I think
    osm_net = pdna.Network(osm_nodes["x"], 
                            osm_nodes["y"],
                            osm_edges["from_int"],
                            osm_edges["to_int"],
                            osm_edges[["weight"]], 
                            twoway=False)
    
    return transit_ped_net, osm_net