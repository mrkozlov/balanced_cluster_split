# Balanced Spatial Cluster Thinning

An R implementation for spatially balanced data splitting using k-means clustering, designed for species distribution modeling. This tool splits presence and pseudo-absence data into training and testing sets while ensuring spatial balance using environmental covariates and Surface Range Envelope (SRE) masking.

## Features
- **Cluster-based spatial splitting**: Uses k-means clustering to group spatial points and ensure balanced representation in training and testing sets.
- **Automated environmental masking**: Applies Surface Range Envelope (SRE) to generate pseudo-absence points within environmentally suitable areas.
- **Visualization**: Produces diagnostic plots showing the spatial distribution of training and testing sets with nearest-neighbor connections.
- **Customizable parameters**: Allows tuning of clustering and splitting parameters for optimal spatial balance.

## Requirements
- R (version 4.0 or higher recommended)
- Required R packages:
  ```R
  install.packages(c("sf", "terra", "FNN", "ggplot2", "svglite"))
Input files:
Species presence points (Shapefile, e.g., new_species_presence.shp)
Study area polygon (Shapefile, e.g., 15_11Layers.shp)
Environmental layers (GeoTIFFs in a folder, e.g., t30/)
Installation
Install R and the required packages listed above.
Download or clone this repository to your local machine.
Ensure all input files (Shapefiles and GeoTIFFs) are accessible and paths are correctly specified in the script.
Usage
Prepare input files:
A Shapefile with species presence points (e.g., new_species_presence.shp).
A Shapefile defining the study area polygon (e.g., 15_11Layers.shp).
A folder containing environmental raster layers in GeoTIFF format (e.g., t30/ with files like wc2.1_30s_bio_11.tif).
Modify the script: Update file paths in balanced_cluster_split.R:
R

Копировать
train_balanced <- st_read("path/to/new_species_presence.shp")
polygon <- st_read("path/to/15_11Layers.shp")
path_30s <- "path/to/t30/"
Adjust parameters: Customize the splitting parameters in the params list:
R

Копировать
params <- list(
  test_size = 0.2,              # Proportion of data for testing (e.g., 20%)
  n_clusters = 12,              # Number of clusters for k-means
  max_points_per_cluster = 2,   # Maximum points per cluster in test set
  max_attempts = 90000,         # Maximum attempts for optimal split
  w_cross_dist = 7.0,           # Weight for cross-group distance
  w_coverage = 9.0,             # Weight for cluster coverage
  w_balance = 60.0,             # Weight for cluster balance
  w_dist_diff = 0.1             # Weight for distance difference penalty
)
Run the script: Execute the script in R:
R

Копировать
source("balanced_cluster_split.R")
Outputs
Training and testing sets:
train_data.txt: Training set with coordinates, environmental covariates, and presence/absence labels.
test_data.txt: Testing set with the same structure.
Diagnostic plots:
balanced_cluster_split_plot.png: PNG visualization of the split.
balanced_cluster_split_plot.tiff: TIFF visualization.
balanced_cluster_split_plot.pdf: PDF visualization.
Splitting report:
splitting_report.txt: Text file with metrics (e.g., dataset sizes, presence ratios, spatial metrics, cluster coverage).
Example Output
Plot: A scatter plot showing presence (grey50) and pseudo-absence (grey90) points, with training (blue) and testing (red) sets overlaid. Dashed lines connect test points to their nearest training points.
Report: Example metrics from splitting_report.txt:
text

Копировать
=== SPLITTING RESULTS ===
Parameters:
- Test size: 0.2
- Max attempts: 90000
- Number of clusters: 12
- Cross-distance weight: 7.0
- Coverage weight: 9.0
- Balance weight: 60.0
- Distance difference weight: 0.1
- Max points per cluster: 2
Dataset sizes:
- Training set: 800 points (80.0%)
- Testing set: 200 points (20.0%)
Presence/Absence distribution:
- Training set presence: 50.2%
- Testing set presence: 49.8%
Spatial metrics:
- Mean distance in train: 0.123
- Mean distance in test: 0.119
- Cross-group distance: 0.115
- Cluster coverage: 100.0%
Notes
The script creates temporary files in D:/temp_terra1. Ensure this directory exists or modify temp_dir in the script.
The script assumes all rasters are in the WGS84 coordinate system (+proj=longlat +datum=WGS84). Verify your input data matches this CRS.
Adjust the selected_layers variable to choose specific environmental layers for SRE masking.