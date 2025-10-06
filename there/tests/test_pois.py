import pytest
import geopandas as gpd
import there.pois as tp

valid_name = "Sydney, Australia"
invalid_name = "Invalid Place, Nowhere"
list_name = ["New York, USA"]

def test_get_poi_dict():
    poi_dict = tp.get_poi_dict(valid_name)
    assert isinstance(poi_dict, dict)

@pytest.fixture
def poi_dict():
    return tp.get_poi_dict(valid_name)

@pytest.fixture
def proj_crs():
    return "EPSG:4326"

def test_downloader_valid(valid_name, poi_dict, proj_crs):
    # Test the poi_downloader function with a valid name
    pois = tp.poi_downloader(valid_name, poi_dict, proj_crs)
    assert isinstance(pois, gpd.GeoDataFrame)
    assert not pois.empty
    assert 'geometry' in pois.columns
    assert len(pois) > 100

def test_downloader_invalid(invalid_name, poi_dict, proj_crs):
    # Test the poi_downloader function with an invalid name
    pois = tp.poi_downloader(invalid_name, poi_dict, proj_crs)
    assert isinstance(pois, gpd.GeoDataFrame)
    assert pois.empty
    with pytest.raises(ValueError):
        tp.poi_downloader(invalid_name, poi_dict, proj_crs)

def test_downloader_list(list_name, poi_dict, proj_crs):
    # Test the poi_downloader function with a list instead of string or gdf
    with pytest.raises(TypeError):
        tp.poi_downloader(list_name, poi_dict, proj_crs)
    