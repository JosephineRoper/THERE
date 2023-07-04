# THERE (Travel from Home to Everywhere Required Everyday)

THERE is a collection of scripts and functions using the Python open-source geospatial ecosystem to calculate accessibility across multiple modes, for any area in the world. The model is based on my PhD research, the following preprint explains the theoretical basis for its design. To be cited if using this code:
Roper, J., Ng, M., & Pettit, C. (2022). WalkTHERE: Travel from Home to Everywhere Required Everyday - An open-source walkability index based on multi-activity accessibility. Unpublished (under review @ JTLU). https://doi.org/10.13140/RG.2.2.34461.59365

### Index principles
The index aims to answer the question "what percentage of people's total needs can they access by a given mode of transport from each location?". The principle of the model is to include all categories of destinations that are regular sources of travel, weighted according to trip frequency in travel surveys, preferably surveys specific to the area in question. The model therefore includes employment as a major travel destination.

The model incorporates potentially infinite numbers of destinations, at any distance from each origin point, using a theory of diminishing returns to opportunity and increasing generalised cost of travel, rather than using a hard limit for number of destinations per category, or a distance threshold. (In these scripts 300-1500 points per category are used, in order to keep running time to a few minutes, but this can be adjusted depending on the accuracy required).

The generalised cost of travel in these scripts is based on network distance with an exponential distance decay for walking, cycling and driving, and time for public transport. All these can be adjusted if other data is available - for example if congested travel time information is available for driving, or detailed walking environment data (gradient, shade, attractiveness) can be used to weight network links.

### Data required to use the index
Employment data at a scale appropriate to the mode under analysis (but preferably as fine as possible, eg. jobs per block) is needed. Alternately the model can be run without employment as a category, by adjusting the weightings of other categories, but this should be very clearly specified anywhere you use the results "THERE minus employment" or similar.

If employment or job numbers are available at a coarse scale (eg Destination Zones in Australia), and you also know which blocks or even parcels are possible employment land and which are not, the Employment Points notebook can be used to produce approximate finer-scale employment data by distributing the DZN jobs across the relevant blocks.

Data for other types of destinations is downloaded from OpenStreetMap, in the example notebooks provided, or can be uploaded from any other data sources you have available (or a combination used).

### Results
The results of the index for walking for some Australian cities can currently be viewed at building level on the Colouring platform: for [Sydney](https://sydney.colouringaustralia.org/view/context/493004), [Adelaide](https://adelaide.colouringaustralia.org/view/context), and [Melbourne](https://melbourne.colouringaustralia.org/view/community/).

An example of comparative results for four modes, in Edinburgh:
![Edinburgh 4 modes](https://user-images.githubusercontent.com/61174184/232359132-4ee3f094-df7e-4ff3-8f1f-c94d0fc7666a.jpeg)
