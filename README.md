# THERE (Travel from Home to Everywhere Required Everyday)

THERE is a collection of scripts and functions using the Python open-source geospatial ecosystem to calculate accessibility across multiple modes, for any area in the world. The model was developed during my PhD research and the theoretical basis for its design will be published in an upcoming journal article.

The principle of the model is to include all categories of destinations that are regular sources of travel, weighted according to trip frequency in travel surveys, preferably surveys specific to the area in question.

The model therefore includes employment as a major travel destination. Employment data at a scale appropriate to the mode under analysis (but preferably as fine as possible, eg. jobs per block) is needed.
Alternately the model can be run without employment as a category, by adjusting the weightings of other categories, but this should be very clearly specified anywhere you use the results "THERE minus employment" or similar.

If employment or job numbers are available at a coarse scale (eg Destination Zones in Australia), and you also know which blocks or even parcels are possible employment land and which are not, the Employment Points notebook can be used to produce approximate finer-scale employment data by distributing the DZN jobs across the relevant blocks.

Data for other types of destinations is downloaded from OpenStreetMap, in the example notebooks provided, or can be uploaded from any other data sources you have available (or a combination used).
