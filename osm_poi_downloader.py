# Download WalkTHERE points of interest from OSM using OSMnx

def poi_downloader_(place, tags):
    tags = {

    # shopping
    'shop':['bakery', 'clothes', 'supermarket', 'mall', 'greengrocer',
            'seafood', 'wine', 'butcher','convenience',
            'beverages', 'alcohol', 'bicycle_shop', 'department_store', 
            'doityourself', 'beauty_shop', 'outdoor_shop', 
            'stationery', 'bookshop', 'gift_shop', 'newsagent', 
            'car_dealership', 'furniture_shop', 'sports_shop',
            'garden_centre', 'computer_shop', 'shoe_shop', 'florist', 
            'video_shop', 'toy_shop', 'mobile_phone_shop', 'jeweller'],
    # possibly we could pick up all shop=True excluding a few. but so many options.
    'amenity':['marketplace'],
    'building':['kiosk', 'supermarket',],

    # errands
    'amenity':['atm','bank','courthouse','post_box', 'post_office',
               'clinic', 'dentist', 'doctors', 'hospital',
               'pharmacy', 'veterinary', 'travel_agent',
               'place_of_worship'],
    'shop':['optician', 'hairdresser', 'laundry',],
    'healthcare':['physiotherapist'],
    'office':['government'], #### further refine ?

    # recreation
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

    # education
    'amenity':['college', 'kindergarten', 'music_school',
               'school', 'university', 'childcare'],

    }


    # convert all to centroids?
    
    gdf = ox.geometries_from_place(place, tags)

    return gdf
