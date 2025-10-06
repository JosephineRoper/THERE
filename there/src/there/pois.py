import osmnx as ox
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPoint

def get_poi_dict():
    poi_dictionary = {
        'employment':{
            'category':['employment']
        },
        'shopping':{
            'shop':['bakery', 'clothes', 'supermarket', 'mall', 'greengrocer',
                    'seafood', 'wine', 'butcher','convenience',
                    'beverages', 'alcohol', 'bicycle_shop', 'department_store', 
                    'doityourself', 'beauty_shop', 'outdoor_shop', 
                    'stationery', 'bookshop', 'gift_shop', 'newsagent', 
                    'car_dealership', 'furniture_shop', 'sports_shop',
                    'garden_centre', 'computer_shop', 'shoe_shop', 'florist', 
                    'video_shop', 'toy_shop', 'mobile_phone_shop', 'jeweller'],
            # possibly we could pick up all shop=True excluding a few. but not sure how
            # and many options to exclude
            'amenity':['marketplace'],
            'building':['kiosk', 'supermarket',],
        },
        'errands':{
            'amenity':['atm','bank','courthouse','post_box', 'post_office',
                    'clinic', 'dentist', 'doctors', 'hospital',
                    'pharmacy', 'veterinary', 'travel_agent',
                    'place_of_worship'],
            'shop':['optician', 'hairdresser', 'laundry',],
            'healthcare':['physiotherapist'],
            'office':['government'], #### further refine ?
        },
        'recreation':{
            'leisure':['dog_park', 'ice_rink', 'park', 'pitch', 'playground',
                    "fitness_centre","sports_centre", 'stadium', 'swimming_pool',
                    'swimming_area', 'track', 'water_park','golf_course',],
            'club':['social'],
            'amenity':['bar', 'biergarten', 'cafe', 'fast_food', 'food_court',
                    'ice_cream', 'pub', 'restaurant', 'nightclub',
                    'library', 'arts_centre', 'cinema', 'community_centre',
                    'social_centre', 'theatre',],
            'building':['stadium', 'castle', 'ruins',],
            'tourism':['aquarium', 'artwork', 'attraction', 'gallery',
                    'museum', 'picnic_site', 'theme_park', 'viewpoint',
                    'zoo'],
            'natural':['beach'],
        },
        'education':{
            'amenity':['college', 'kindergarten', 'music_school',
                    'school', 'university', 'childcare'],
        }
    }
    return poi_dictionary

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
        try:
            gdf = ox.geometries_from_place(place, tags).to_crs(proj_crs)
        except Exception as e:
            print(f"Error downloading POIs for {place}: {e}")
            raise ValueError("Check your internet connection and the place name. Try using osmnx.geocode_to_gdf(place) first and using the result with poi_downloader.")
    elif type(place) == gpd.GeoDataFrame:
        place = place.to_crs("EPSG:4326")
        # I'm not sure why I'm doing this instead of using geometries_from_polygon. I think because it errors a lot.
        bbox = place.bounds
        bbox_pois = ox.geometries.geometries_from_bbox(bbox['maxy'][0], bbox['miny'][0], bbox['maxx'][0], bbox['minx'][0], tags)
        gdf = gpd.clip(bbox_pois, place, keep_geom_type=False).to_crs(proj_crs)
        #this doesn't work with multiple polygons, it dissolves them
    else:
        raise TypeError("'place' should be a string for querying OSM, or a geodataframe containing a polygon of the place boundaries.")

    # OSM POIs include domestic swimming pools in some areas. This line removes swimming pools less than 100m2.
    # Same for domestic tennis courts appearing as 'pitches'. Remove pitches below 450m2.
    gdf = gdf[~((gdf['leisure']=='swimming_pool') & (gdf.area < 100))]
    gdf = gdf[~((gdf['leisure']=='pitch') & (gdf.area < 450))]
        
    gdf['orig_geometry'] = gdf.geometry
    # convert all to centroids
    gdf.geometry = gdf.centroid

    # sometimes there is a multiindex with element_type returned. need to work out if this always happens.
    gdf.index = gdf.index.droplevel('element_type')
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

def nearest_to_centroids(polygons, points, **kwargs):
    orig_geometry = polygons.geometry
    # convert polygons to centroids
    polygons.geometry = polygons.centroid
    polygons = polygons.sjoin_nearest(points, **kwargs)
    # convert back to polygons
    polygons.geometry = orig_geometry
    return polygons

def poly_vertices(data_gdf):
    # Add area column, convert polygons to vertices as Multipoints, explode to Points.
    # Thus original polygon area will be retained with the vertices,
    # which is useful for some analyses

    exploded = data_gdf.copy()
    exploded['poly_area'] = exploded.to_crs('+proj=wintri').area
    exploded.geometry = exploded.geometry.apply(lambda x: MultiPoint(list(x.exterior.coords)))
    exploded = exploded.explode(index_parts=True)
    
    return exploded

def remove_duplicate_pois(datasets, buffer=10, variable=None):
    # input is a list of maybe one or several geodataframes
    # buffer default is 10m = assumes projected CRS
    # datasets must have a 'category' column defining which duplicates will be removed - only pois that have the same category and are within the buffer distance of each other will be removed
    # variable is used where the pois have attributes that need to be summed eg. job counts
    
    pois = pd.concat(datasets).reset_index(drop=True)
    pois_buffer = pois.copy()
    pois_buffer.geometry = pois.geometry.buffer(buffer)
    
    joined_pois = gpd.sjoin(
        pois, pois_buffer, how="left", predicate='intersects', 
        lsuffix='1', rsuffix='2')

    joined_pois = joined_pois[joined_pois['category_1'] == joined_pois['category_2']].copy()
    joined_pois['unique_id'] = joined_pois.apply(lambda row: np.nanmin([row.name, row['index_2']]), axis=1)

    final_pois = joined_pois.drop_duplicates(subset='unique_id', keep=False)

    # sum variable using unique_id
    if variable is not None:
        joined_pois['index_1'] = joined_pois.index
        final_pois[variable + '_1'] = joined_pois.groupby(['index_1'])[variable + '_2'].sum()

    print("Removed " 
          + "{0:.2f}".format((1-len(final_pois)/len(pois))*100) 
          + "% duplicate points from dataframes - "
          + str((len(pois)-len(final_pois))) + " points")
    
    # might need to drop duplicate columns with 2 suffixes
    final_pois = final_pois.loc[:, ~final_pois.columns.str.endswith('_2')]
    # then remove _1 suffixes on remaining columns
    final_pois.columns = final_pois.columns.str.rstrip('_1')

    return final_pois  #pois[pois.index.isin(final_pois.index)]

def employ_points(dzns, meshblocks, meshblock_shape, proj_crs):

    employ_mbs = meshblocks[meshblocks['MB_CATEGORY_NAME_2016'].isin(
        ['Commercial','Primary Production','Hospital/Medical','Education','Other','Industrial'])]
    employ_areas = pd.DataFrame(employ_mbs.groupby('DZN_CODE_2016')['AREA_ALBERS_SQKM'].sum())
    DZN_areas = dzns.join(employ_areas, on='DZN (POW)', how='left')

    # sometimes 'Count', sometimes 'Number' or 'Jobs'
    excluded = DZN_areas[(DZN_areas['Jobs']>0) & (DZN_areas['AREA_ALBERS_SQKM'].isna())]['Jobs'].sum()/DZN_areas['Jobs'].sum()
    print(['Proportion of jobs in excluded meshblocks:'], excluded)

    DZN_areas['JobDensity'] = DZN_areas['Jobs']/DZN_areas['AREA_ALBERS_SQKM']

    employ_mbs = employ_mbs.join(DZN_areas.set_index('DZN (POW)'), on='DZN_CODE_2016', how='inner', rsuffix='_DZN')

    employ_mbs['Jobs'] = employ_mbs['JobDensity']*employ_mbs['AREA_ALBERS_SQKM']

    meshblock_shape['MB_CODE16'] = meshblock_shape['MB_CODE16'].astype('int64')
    employ_mbs['MB_CODE_2016'] = employ_mbs['MB_CODE_2016'].astype('int64')

    employ_shapes = meshblock_shape.join(employ_mbs.set_index('MB_CODE_2016')[['DZN_CODE_2016','Jobs']], how='right', on='MB_CODE16')

    centroids = employ_shapes.copy()
    centroids.geometry = (centroids.geometry
                            .to_crs(proj_crs)
                            .centroid)
    return centroids

def default_poi_params():
    #this returns the default parameters for the default POI categories based on Sydney household travel survey data and walking mode for the distance constant.
    data = {'category': ['employment', 'education', 'shopping', 'errands', 'recreation'],
        'weight': [31.9, 14.3, 21.6, 8.9, 23.3],
        'distance_constant': [0.001, 0.001, 0.001, 0.001, 0.001],
        'diminishing_returns_constant': [0.000017, 0.3, 0.3, 0.3, 0.3],
        'variable': ['Jobs', 'count', 'count', 'count', 'count'],
        'num_pois': [100, 25, 25, 25, 25]}
    df = pd.DataFrame(data)
    df.set_index('category', inplace=True)
    return df