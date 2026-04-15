# User Manual — FSLAM Calibration Tool

## Purpose

This tool calibrates input parameters for FSLAM, a physically based model used to estimate landslide occurrence on slopes. The calibration compares model predictions against points where landslide occurrence is known, then searches for parameter adjustments that improve model performance.

The main file is `nsga2_FSLAM_v1.R`. It contains the code required to load data, run FSLAM, calibrate parameters, review metrics, and optionally visualize results. This manual explains how to prepare and use the tool without repeating the code already included in the script.

## What The User Gets

After a calibration run, the tool provides:

- A calibrated parameter set for adjusting soil, land-cover, rainfall, and infiltration properties.
- Performance indicators for pre-event and post-event conditions.
- A Pareto curve when using multi-objective calibration.
- A table with the best parameter combinations according to the selected criteria.
- An optional interactive map with stability/instability predictions and derived hydrological variables.

## Requirements

You need R installed. RStudio is recommended because it makes it easier to open the project, inspect tables, run the script by sections, and view interactive maps.

The script uses R packages for spatial data handling, data manipulation, calibration, classification metrics, interactive plots, and maps. If R reports that a package is missing, install that package from the R console before continuing. The main packages are `sf`, `dplyr`, `caret`, `nsga2R`, `plotly`, and, for maps, `mapview`.

The single-objective calibration option uses a function named `dds`. If you want to use that option, the function must already be available in your R environment. If it is not available, use the multi-objective calibration workflow, which is the main option in the script.

By default, the script keeps the optional single-objective calibration and interactive map steps disabled. Enable those options inside the R file only when the required dependencies are installed and you want to run those stages.

## Project Structure

The project contains the main script and two sample data folders:

- `nsga2_FSLAM_v1.R`: main calibration and visualization script.
- `MORLE_input_data`: sample data for areas with observed landslides.
- `noMORLE_input_data`: sample data for areas without observed landslides.

Each data folder includes three files:

- `points_sample.csv`: terrain points or cells.
- `soil.csv`: soil-type properties.
- `hmtu.csv`: land-cover/use properties and curve numbers.

## Input Data

### Terrain Points

The points file represents the terrain cells or locations where the model will be evaluated. Each row must include coordinates, observed class, slope, contributing area, soil type, land-cover type, rainfall, antecedent recharge, and curve number.

The observation column must indicate whether a landslide occurred. In the sample file, `des = 1` means landslide and `des = 0` means stable condition.

Slope must be provided in radians. Coordinates must match the coordinate reference system that will later be used for mapping.

### Soil Properties

The soil file contains, for each soil type, hydraulic conductivity, minimum and maximum cohesion, minimum and maximum internal friction angle, soil depth, density, porosity, and hydrologic soil group.

The soil identifier must match the `soil` column in the points file.

### Land-Cover Properties

The `hmtu.csv` file contains root cohesion and curve numbers associated with hydrologic soil groups A, B, C, and D.

The land-cover identifier must match the `lucl` column in the points file.

### Units Row

The `soil.csv` and `hmtu.csv` files include a units row immediately after the header. The script is prepared to ignore that row during calculations, so do not remove it if you are working with the same format as the sample data.

## Preparation Before Running

Open the project from the root folder, meaning the folder that contains `nsga2_FSLAM_v1.R`, `MORLE_input_data`, and `noMORLE_input_data`.

At the beginning of the script, check that the input paths point to the dataset you want to use. By default, the script is configured to use `MORLE_input_data`.

If you want to calibrate with your own data, replace the sample files or update the paths in the script to point to your data folder. Keep the same column names and expected units.

Also review the fixed model parameters:

- `g`: gravity.
- `dw`: water density.
- `b`: raster or grid cell size.
- `ds`: soil density.

The parameter that usually needs the most attention is `b`, because it must match the spatial resolution of your data.

## Parameters Calibrated By The Model

The calibration adjusts ten parameters:

| Parameter | Meaning | Adjustment type |
|-----------|---------|-----------------|
| `Cs` | Soil cohesion | Multiplier |
| `Cr` | Root cohesion | Multiplier |
| `phi_min` | Minimum internal friction angle | Offset in degrees |
| `phi_max` | Maximum internal friction angle | Offset in degrees |
| `Pa` | Antecedent recharge | Multiplier |
| `z` | Soil depth | Offset in meters |
| `Pe` | Rainfall event | Multiplier |
| `k` | Hydraulic conductivity | Multiplier |
| `por` | Porosity | Offset |
| `CN` | Curve number | Multiplier |

The search ranges are defined in the script through lower and upper limits. If a parameter has the same lower and upper limit, it remains fixed during calibration.

## Recommended Workflow

1. Open the project in RStudio.
2. Check that the input data are in the correct folder.
3. Confirm that cell size and units match your study area.
4. Run a short calibration first to verify that the data load correctly.
5. Review the pre-event and post-event metrics.
6. If the result is reasonable, increase the number of generations or population size for a more complete search.
7. Select a parameter combination with a good balance between pre-event and post-event performance.
8. Visualize the result on a map and check whether the predictions make spatial sense.
9. Export the result only after you are satisfied with the calibration.

## Multi-Objective Calibration

This is the recommended option. It searches for two objectives at the same time:

- Good performance before the rainfall event, when cells should ideally remain stable.
- Good performance after the event, when the model should identify cells with observed landslides.

The main output is a Pareto front. Each point on that front represents a different parameter combination. There is not always a single best solution; in most cases, you choose a combination that provides an acceptable balance between both metrics.

In the script, the filtered table keeps only the combinations that exceed minimum performance thresholds. If the table is empty, temporarily reduce the thresholds or increase the algorithm search effort.

## Single-Objective Calibration

This option focuses only on post-event performance. It can be useful when the main goal is to detect observed landslides, but it should be used carefully because it may accept parameter sets that predict too much instability before the event.

Use this option only if the `dds` function is available and you understand the difference between optimizing one metric and optimizing the balance between pre-event and post-event conditions.

## Interpreting Results

The main metrics are:

- `acc_ini`: accuracy before the rainfall event.
- `acc_fin`: accuracy after the rainfall event.

A good parameter set should maintain high stability before the event while also detecting unstable zones after the event.

Do not automatically select the row with the highest post-event accuracy. Review the pre-event accuracy and the resulting map as well. A calibration that marks too many cells as unstable may look good for landslide detection, but it may be physically unrealistic.

## Map Visualization

The map function converts the final table into a spatial object and displays interactive layers. The main layers are:

- Post-event probability-of-failure prediction.
- Pre-event probability-of-failure prediction.
- Ratios between water level, antecedent recharge, event contribution, and soil depth.

Before mapping, confirm the coordinate reference system. The CRS must match the `x` and `y` coordinates in the points file. If the CRS does not match, the map may appear shifted or in the wrong location.

Map visualization runs only when the map option is enabled inside the script. This prevents a normal calibration run from failing on systems where `mapview` is not installed yet.

## Export

Spatial export is optional. Enable it only after selecting a calibrated solution. The script exports using modern `sf` functions, so it does not depend on `rgdal`.

Use a clear output path and a dedicated results folder so generated files are not mixed with input data.

## Common Issues

### The Script Cannot Find The Files

Check that you are working from the project root folder and that the paths at the beginning of the script match the actual location of the CSV files.

### An R Package Is Missing

Install the package named in the error message and run the script again. Map visualization requires `mapview`; NSGA-II calibration requires `nsga2R`.

### The Filtered Table Is Empty

The accuracy thresholds may be too strict for a short run. Try a longer search or check whether the observed labels and input units are correct.

### The Map Does Not Appear

Check that `mapview` is installed, that a calibrated solution exists, and that the CRS matches the coordinates of your points.

### Results Look Physically Inconsistent

Review units, slope in radians, cell size, soil depth, porosity, and curve numbers. Small unit errors can strongly affect the results.

## Recommendations

Start with the sample data to confirm that your R environment works. Then move to your own data.

Use short runs for technical checks and longer runs for final results.

Document the data folder used, parameter ranges, selection thresholds, and run date. This makes the calibration easier to reproduce and compare across scenarios.

Do not interpret the calibration using only one numerical metric. Combine metrics, maps, and knowledge of the study area.
