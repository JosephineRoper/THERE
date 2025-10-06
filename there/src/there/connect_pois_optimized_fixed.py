# Fixed version of connect_pois function
# Key performance improvements over original implementation + bug fixes

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiPoint
from shapely.ops import snap, split
from shapely import wkt
import rtree
import itertools

pd.options.mode.chained_assignment = None

def connect_poi_optimized_fixed(pois, nodes, edges, key_col=None, path=None, threshold=200, knn=5, meter_epsg=3857, prefix=9990000000):
    """
    Fixed version of connect_poi_optimized that properly returns all edges (original + new).
    
    FIXES:
    1. Proper variable scoping for update_edges function
    2. Ensures all original edges are preserved with proper connect_type marking
    3. Better handling of edge concatenation
    """
    
    # Helper functions (optimized versions)
    def find_kne(point, lines):
        """Find the nearest edge to a point"""
        # Handle both Series and list inputs
        if hasattr(lines, 'values'):
            # It's a pandas Series, convert to list of geometries
            line_geoms = lines.values
        else:
            line_geoms = lines
            
        dists = np.array([geom.distance(point) for geom in line_geoms])
        kne_pos = dists.argmin()
        
        if hasattr(lines, 'iloc'):
            # pandas Series - use iloc for positional indexing
            kne = lines.iloc[[kne_pos]]
            kne_idx = lines.index[kne_pos]
        else:
            # regular list/array
            kne_idx = kne_pos
            kne = line_geoms[kne_pos]
            
        return kne_idx, line_geoms[kne_pos]

    def get_pp_vectorized(points, lines):
        """Vectorized projection point calculation"""
        return [line.interpolate(line.project(point)) for point, line in zip(points, lines)]

    def split_line(line, pps):
        """Split 'line' by all intersecting 'pps'."""
        line = snap(line, pps, 1e-8)
        try:
            new_lines = list(split(line, pps).geoms)
            return new_lines
        except TypeError as e:
            print('Error when splitting line: {}\n{}\n{}\n'.format(e, line, pps))
            return []

    def update_poi_nodes(nodes, new_points, meter_epsg=3857):
        """Update nodes with POI points"""
        new_nodes = new_points.copy()
        new_nodes['connect_type'] = node_type_poi
        new_nodes['connect_id'] = new_nodes[key_col].astype(int)
        new_nodes['new'] = 1
        new_nodes.set_crs(epsg=meter_epsg, inplace=True)

        gdfs = [nodes, new_nodes]
        nodes = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True, sort=False), crs=gdfs[0].crs)
        nodes.reset_index(drop=True, inplace=True)
        nodes['connect_id'] = nodes.index
        return nodes, new_nodes

    def update_pap_nodes(nodes, new_points, meter_epsg=3857):
        """Update nodes with projected access points"""
        new_nodes = gpd.GeoDataFrame(new_points, columns=['geometry'], crs=f'epsg:{meter_epsg}')
        n = len(new_nodes)
        new_nodes['connect_type'] = node_type_pp
        new_nodes['connect_id'] = [int(connect_id_prefix + i) for i in range(n)]
        new_nodes['new'] = 1
        new_nodes.set_crs(epsg=meter_epsg, inplace=True)

        gdfs = [nodes, new_nodes]
        nodes = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True, sort=False), crs=gdfs[0].crs)
        nodes.reset_index(drop=True, inplace=True)
        nodes['connect_id'] = nodes.index
        return nodes, new_nodes

    def extract_coords_efficient(geometries):
        """Extract coordinates more efficiently without WKT serialization"""
        start_coords = []
        end_coords = []
        
        for geom in geometries:
            coords = list(geom.coords)
            start_coords.append(coords[0])
            end_coords.append(coords[-1])
        
        return start_coords, end_coords

    def create_nodes_lookup_efficient(nodes_gdf, precision=2):
        """Create efficient coordinate lookup without WKT serialization"""
        factor = 10 ** precision
        coords_dict = {}
        
        for idx, geom in enumerate(nodes_gdf['geometry']):
            x, y = geom.x, geom.y
            x_round = round(x * factor) / factor
            y_round = round(y * factor) / factor
            coord_key = (x_round, y_round)
            coords_dict[coord_key] = nodes_gdf.iloc[idx]['connect_id']
        
        return coords_dict

    def update_edges(edges, new_lines, replace, line_pps_dict, nodes_id_dict, pois, key_col, edge_type, threshold, meter_epsg):
        """Update edge info by adding new_lines or replacing existing ones"""
        if replace:
            kne_idxs = list(line_pps_dict.keys())
            lens = [len(item) for item in new_lines]
            new_lines_gdf = gpd.GeoDataFrame({
                'kne_idx': np.repeat(kne_idxs, lens),
                'geometry': list(itertools.chain.from_iterable(new_lines))
            })
            
            cols = list(edges.columns)
            cols.remove('geometry')
            new_edges = new_lines_gdf.merge(edges[cols], how='left', left_on='kne_idx', right_index=True)
            new_edges.drop('kne_idx', axis=1, inplace=True)
            new_lines = new_edges['geometry']
        else:
            new_edges = gpd.GeoDataFrame(pois[[key_col]], geometry=new_lines, columns=[key_col, 'geometry'])
            new_edges['oneway'] = False
            new_edges['connect_type'] = edge_type

        # Optimized coordinate extraction and lookup
        new_edges['length'] = [l.length for l in new_lines]
        
        start_coords, end_coords = extract_coords_efficient(new_edges['geometry'])
        
        # Round coordinates for matching
        precision = 1
        factor = 10 ** precision
        start_coords_rounded = [(round(x * factor) / factor, round(y * factor) / factor) for x, y in start_coords]
        end_coords_rounded = [(round(x * factor) / factor, round(y * factor) / factor) for x, y in end_coords]
        
        new_edges['from'] = [nodes_id_dict.get(coord, None) for coord in start_coords_rounded]
        new_edges['to'] = [nodes_id_dict.get(coord, None) for coord in end_coords_rounded]
        new_edges['connect_id'] = ['_'.join(list(map(str, s))) for s in zip(new_edges['from'], new_edges['to'])]
        
        print("Missing 'to' nodes:", len(new_edges[new_edges['to'].isna()]))
        
        start = edges.index[-1] + 1
        stop = start + len(new_edges)
        new_edges.index = range(start, stop)

        if replace:
            edges = edges.drop(kne_idxs, axis=0)
        else:
            valid_pos = np.where(new_edges['length'] <= threshold)[0]
            n = len(new_edges)
            n_fault = n - len(valid_pos)
            f_pct = n_fault / n * 100
            print("Remove faulty projections: {}/{} ({:.2f}%)".format(n_fault, n, f_pct))
            new_edges = new_edges.iloc[valid_pos]

        new_edges.set_crs(epsg=meter_epsg, inplace=True)
        dfs = [edges, new_edges]
        edges = gpd.GeoDataFrame(pd.concat(dfs, ignore_index=False, sort=False), crs=dfs[0].crs)
        return edges, new_edges

    # Configuration
    orig_crs = pois.crs
    node_type_pp = 'projected_pap'
    node_type_poi = 'poi'
    edge_type = 'projected_footway'
    connect_id_prefix = prefix

    # Convert CRS
    pois_meter = pois.to_crs(epsg=meter_epsg)
    nodes_meter = nodes.to_crs(epsg=meter_epsg)
    edges_meter = edges.to_crs(epsg=meter_epsg)
    
    # IMPORTANT FIX: Mark original edges with connect_type if they don't have it
    if 'connect_type' not in edges_meter.columns:
        edges_meter['connect_type'] = 'original_network'

    # OPTIMIZED: Use BallTree instead of individual R-tree queries
    print("Building spatial index...")
    try:
        from sklearn.neighbors import BallTree
        use_balltree = True
    except ImportError:
        print("scikit-learn not available, falling back to R-tree, this will be much slower.")
        use_balltree = False
        Rtree = rtree.index.Index()
        [Rtree.insert(fid, geom.bounds) for fid, geom in edges_meter['geometry'].items()]

    # STAGE 1: interpolation
    print("Updating external nodes...")
    nodes_meter, _ = update_poi_nodes(nodes_meter, pois_meter, meter_epsg=meter_epsg)
    
    print("Projecting POIs to the network...")
    
    if use_balltree:
        # OPTIMIZED: Batch spatial query using BallTree
        print("Using optimized BallTree spatial indexing...")
        
        # Extract edge midpoints for spatial indexing
        edge_midpoints = np.array([[geom.centroid.x, geom.centroid.y] for geom in edges_meter['geometry']])
        poi_coords = np.array([[pt.x, pt.y] for pt in pois_meter['geometry']])
        
        # Build BallTree for fast spatial queries
        tree = BallTree(edge_midpoints, metric='euclidean')
        distances, near_indices = tree.query(poi_coords, k=knn)
        
        # Convert indices to edge IDs and get corresponding geometries
        edge_ids = edges_meter.index.values
        pois_meter['near_idx'] = [[edge_ids[idx] for idx in indices] for indices in near_indices]
        pois_meter['near_lines'] = [edges_meter['geometry'].loc[edge_ids[indices]] for indices in near_indices]
        
        # Optimized find_kne for better performance
        kne_indices = []
        knes = []
        
        for i, (point, near_lines) in enumerate(zip(pois_meter['geometry'], pois_meter['near_lines'])):
            # Calculate distances to candidate edges only
            distances = np.array([line.distance(point) for line in near_lines.values])
            min_pos = distances.argmin()
            
            kne_idx = near_lines.index[min_pos]
            kne_indices.append(kne_idx)
            knes.append(near_lines.iloc[min_pos])
        
        pois_meter['kne_idx'] = kne_indices
        knes = knes
        
    else:
        # Original R-tree method (slower)
        print("Using R-tree spatial indexing...")
        pois_meter['near_idx'] = [list(Rtree.nearest(point.bounds, knn)) for point in pois_meter['geometry']]
        pois_meter['near_lines'] = [edges_meter['geometry'][near_idx] for near_idx in pois_meter['near_idx']]
        pois_meter['kne_idx'], knes = zip(*[find_kne(point, near_lines) for point, near_lines in 
                                          zip(pois_meter['geometry'], pois_meter['near_lines'])])
    
    # OPTIMIZED: Vectorized projection point calculation
    pois_meter['pp'] = get_pp_vectorized(pois_meter['geometry'], knes)

    print("Updating internal nodes...")
    nodes_meter, _ = update_pap_nodes(nodes_meter, list(pois_meter['pp']), meter_epsg=meter_epsg)
    
    # OPTIMIZED: Create coordinate lookup more efficiently
    print("Creating coordinate lookup...")
    nodes_id_dict = create_nodes_lookup_efficient(nodes_meter, precision=1)

    print("Updating internal edges...")
    line_pps_dict = {k: MultiPoint(list(v)) for k, v in pois_meter.groupby(['kne_idx'])['pp']}
    
    # Fix: Handle case where idx might return Series instead of single geometry
    new_lines = []
    for idx, pps in line_pps_dict.items():
        # Get the geometry, handling both single values and Series
        geom = edges_meter['geometry'][idx]
        if hasattr(geom, 'iloc'):
            # It's a Series, take the first element
            geom = geom.iloc[0]
        #same for pps - convert multiPoint to single Point and take first
        if isinstance(pps, MultiPoint):
            points = [p for p in pps.geoms]
            pps = points[0]
        split_result = split_line(geom, pps)
        new_lines.append(split_result)
    
    print(len(edges_meter), len(new_lines))
    edges_meter, _ = update_edges(edges_meter, new_lines, replace=True, 
                                line_pps_dict=line_pps_dict, nodes_id_dict=nodes_id_dict,
                                pois=pois_meter, key_col=key_col, edge_type=edge_type,
                                threshold=threshold, meter_epsg=meter_epsg)
    print(len(edges_meter))

    # STAGE 2: connection
    print("Updating external links...")
    pps_gdf = nodes_meter[(nodes_meter['connect_type'] == node_type_pp) & (nodes_meter['new'] == 1)]
    new_lines = [LineString([p1, p2]) for p1, p2 in zip(pois_meter['geometry'], pps_gdf['geometry'])]
    edges_meter, _ = update_edges(edges_meter, new_lines, replace=False,
                                line_pps_dict=line_pps_dict, nodes_id_dict=nodes_id_dict,
                                pois=pois_meter, key_col=key_col, edge_type=edge_type,
                                threshold=threshold, meter_epsg=meter_epsg)

    # STAGE 3: output
    edges_meter['length'] = edges_meter.length
    nodes = nodes_meter.to_crs(orig_crs).drop('new', axis=1)
    edges = edges_meter.to_crs(orig_crs)

    # Preprocess for pandana
    nodes.index = nodes['connect_id']
    nodes['x'] = [p.x for p in nodes['geometry']]
    nodes['y'] = [p.y for p in nodes['geometry']]
    edges['length'] = edges['length'].astype(float)

    # Report issues
    if len(nodes_meter) != len(nodes_id_dict):
        print("NOTE: duplication in node coordinates keys")
        print("Nodes count:", len(nodes_meter))
        print("Node coordinates key count:", len(nodes_id_dict))

    print("Missing 'from' nodes:", len(edges[edges['from'].isna()]))
    print("Missing 'to' nodes:", len(edges[edges['to'].isna()]))
    
    # VERIFICATION: Print edge type breakdown
    print("\nEdge breakdown:")
    print(edges['connect_type'].value_counts())

    if path:
        nodes.to_file(path+'/nodes.shp')
        edges.to_file(path+'/edges.shp')

    return nodes, edges
