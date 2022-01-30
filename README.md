# paper-d18o-phosphate-riparian

Data and code repository for the paper titled *Phosphate Oxygen Isotope Ratios in Vegetated Riparian Buffer Strip Soils* in Vadose Zone Journal.

This README.md file was generated on 2022-01-30 by Sheila Saia.

[![DOI](https://zenodo.org/badge/452524774.svg)](https://zenodo.org/badge/latestdoi/452524774)

## General Information

**Title of Dataset**<br>

"paper-d18o-phosphate-riparian"

**Brief Dataset Description**

This GitHub repository was created to provide access to collected data, analysis code, and other information associated with the paper by Bauke et al. titled, *CPhosphate Oxygen Isotope Ratios in Vegetated Riparian Buffer Strip Soils* in Vadose Zone Journal (in Press).

**Dataset Contact Information**<br>
Name: Sheila Saia<br>
Institution: State Climate Office of North Carolina, North Carolina State University<br>
Address: 1005 Capability Drive, Raleigh, NC 27606, USA<br>
Email: ssaia at ncsu dot edu<br>
Role: author, maintainer<br>

**Date of Data Collection**<br>

We collected all soil samples for this study at the University of Bonn Research Farm (Frankenforst) in October 2016. We brough soil samples back to the University of Bonn and analyzed them for d18O-phosphate values and Hedley phosphorus (P) extraction concentrations as described in the associated paper by Bauke et al. (2022).<br>

**Geographic location of data collection**<br>

We collected all soil samples for this study at the University of Bonn Research Farm (Frankenforst) in Germany (7° 12' 22'' E, 50° 42' 49'' N).

**Information about funding sources that supported the collection of the data**<br>

We have no specific funding sources and support to declare here.

## Sharing & Access Information ##

**Licenses/restrictions placed on the data**<br>

Please cite this work (see recommended citation below and in the [`CITATION.md` file](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/CITATION.md) in this repository) and use/distribute according to the CC-BY v4.0 license.

For a full view of the license visit the [`LICENSE.md` file](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/LICENSE.md) file in this repository. For a human readable version of this license visit https://creativecommons.org/licenses/by/4.0/.

**Links to publications that cite or use the data**<br>

As of 2022-01-30 there are no other publications that cite or use these data.

**Links to other publicly accessible locations of the data**<br>

This dataset and associated R code are available at https://github.com/sheilasaia/paper-d18o-phosphate-riparian and Zenodo (https://zenodo.org/record/5920851#.Yfb04mBOlYs). The associated publication is available via Vadose Zone Journal (LINK_HERE).

**Links/relationships to ancillary data**<br>

There are no ancillary data associated with this paper.

**Data derived from another source**<br>

This study does not rely on data derived from another source.

**Additional related data collected that was not included in the current data package**<br>

There are no other related data collected besides what is included here.

**Are there multiple versions of the dataset?**<br>

There are no other versions of the data associated with this paper.

**Recommended citation for the data**<br>

See the [`CITATION.md` file](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/CITATION.md) for the recommended citation for these data and code.

**Paper Availability**<br>

The associated publication is available via Vadose Zone Journal (LINK_HERE).

## Methodological Information ##

**Description of methods used for collection/generation of data:**<br>

Briefly, we collected soil samples and analyzed them for d18O-phosphate values and various Hedley phosphorus (P) extraction concentrations. We provide a detailed description of the methods used to collect and analyze these data in the associated code and publication. The publication is available via Vadose Zone Journal (LINK_HERE).

**Methods for processing the data:**<br>

We included the raw data (as analyzed by Dr. Christian von Sperber and Carina Popp) in the [`raw_data` directory](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/tree/main/data/raw_data) of this repository and all post-processed data are included in the [`data` directory](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/tree/main/data) of this repository.

We post-processed raw data by manually reorganizing data into a csv file and quality controlled these data by checking against the raw data several times. We also used table sorts and basic plots in R to visually inspect these data and ensure there were no typos/transpositions of data values.

**Instrument- or software-specific information needed to interpret the data:**<br>

All data collection, processing, and analysis was done in R (v4.2.1) using RStudio desktop (v2021.09.1+372) on a Macbook Air laptop (macOS Monterey v12.1, Apple M1 2020 processor, 16 GB RAM, 500 GB SSD hard drive).

We also used the following R packages: `tidyverse` (v1.3.1), `forecats` (v0.5.1), `here` (v1.0.1), `car` (v3.0.12), `lme4` (v1.1.27.1), `sjPlot` (v2.8.10), and `effects` (v4.2.0), and `renv` (v0.15.2).

*NOTE:* Users of these scripts can use the [`renv` R package](https://rstudio.github.io/renv/articles/renv.html) to activate the R packages used in this R project repository. Run the code below to activate the R environment that we used for all analysis in your RStudio session. **Make sure you are running code as an *RStudio project***. To do this, download the whole GitHub repository and double click on the "paper-d18o-phosphate-riparian.Rproj" file. This will open RStudio and anchor your RStudio session to this project. Once you have done this, follow the steps below.

```
# Steps to Setup Your RStudio Environment

# 1. install the renv R package
install.packages("renv")

# 2. load the library
libary(renv)

# 3. initialize the environment in your RStudio session
renv::activate()

# 4. Run the scripts as you normally would.
```

**Standards and calibration information, if appropriate:**<br>

There are no applicable standards and calibration information.

**Environmental/experimental conditions:**<br>

There are no applicable environmental/experimental conditions.

**Describe any quality-assurance procedures performed on the data:**<br>

We post-processed raw data by manually reorganizing data into a csv file and quality controlled these data by checking against the raw data several times. We also used table sorts and basic plots in R to visually inspect these data and ensure there were no typos/transpositions of data values.

**People involved with sample collection, processing, analysis and/or submission:**<br>

Sheila Saia wrangled, analyzed, and submitted the data using the code associated with this repository. If you find any errors, please submit [an issue](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/issues) or contact Sheila Saia at ssaia at ncsu dot edu.

## Data & File Overview

### 1. main directory ###

**File List**

The main directory contains the following documentation files:

* `CITATION.md` - This file contains information on how to cite the data, code, and paper associated with this repository.
* `CONTRIBUTING.md` - This file contains information on how to notify the maintainer of issues with the code.
* `LICENSE.md` - This file contains information regarding how this work is licensed and how to give attribution.
* `README.md` - This file (here) provides overall information on about the content and organization of this repository.

The main directory contains the following R-related files:

* `paper-d18o-phosphate-riparian.Rproj` - This is the RStudio project file.
* `analysis_script.R` - The purpose of this R script is to process, analysis, and plot results of data collected in this study.

The main directory contains the following sub-directories, which are explained in further detail below.

* `data` - This directory contains the raw and post-processed data associated with the paper.
* `figures_tables` - This is the output directory for figures and tables generated by the `analysis_script.R script`.

**Relationship Between Files**<br>

Outputs (e.g., figures) from `analysis_script.R` are found in the `figures_tables` directory. Raw data files used to generate `d18o_data.csv` and `hedley_data.csv` in the `data` directory can be found in the `raw_data` directory, included with the `data` directory. See the `README.md` file in the `raw_data` directory for more information on these raw data.

### 2. data directory ###

#### 2.1 data directory ####

**File List & Relationship Between Files**

This directory contains the post-processed (i.e., tidied) data used in this study. This directory also includes the `raw_data` sub-directory, which is described in Section 2.2 below. Both csv files in this directory are inputs to `analysis_script.R`.

* `d18o_data.csv` - This dataset represents the d18O-phosphate values from samples taken along the study transect at different soil profile depths.

* `hedley_data.csv` - This dataset represents the Hedley phosphorus (P) extraction concentrations (including inorganic P and organic P) from samples taken along the study transect at different soil profile depths.

## Data-Specific Information For: `d18o_data.csv` ##

**Number of columns/variables**

Number of columns: 7

**Number of rows**

Number of rows: 48

**Variable list**<br>

* row_num - unique dataset row number
* bonn_id - University of Bonn sample id number
* transect_dist_id - transect distance id number (ranges from 1 to 6, where 1 is closest to the stream and 6 is furthest away from the stream)
* transect_dist_m - transect distance in meters (ranges from 1 to 55 meters away from the stream)
* replicate_num - sample replicate id number
* depth_cm - soil sample profile depth (cm)
* d18o_value - d18O-phosphate value (per mil)

**Missing data codes**<br>

There is no missing data so it is not necessary to explain missing data codes.

**Specialized formats or other abbreviations used:**

All specialized formats and abbreviations are described here.


## Data-Specific Information For: `hedley_data.csv` ##

**Number of columns/variables**

Number of columns: 9

**Number of rows**

Number of rows: 768

**Variable list**<br>

* row_num - unique dataset row number
* transect_dist_id - transect distance id number (ranges from 1 to 6, where 1 is closest to the stream and 6 is furthest away from the stream)
* depth_cm - soil sample profile depth (cm)
* extract - Hedley phosphorus (P) extraction pool description
* fraction - description of fraction, either Po for organic P or Pi for inorganic P
* replicate - sample replicate id number
* soil_moist_type - description of soil moisture conditions of sample analysis, either dried for dried soil sample preparation or field_fresh for field fresh preparation
* p_conc_mgperkg - phosphorus (P) concentration of Hedley extraction (mg of P per kg of soil)
* notes - any specific notes about the sample

**Missing data codes**<br>

There is no missing data so it is not necessary to explain missing data codes.

**Specialized formats or other abbreviations used:**

All specialized formats and abbreviations are described here.


#### 2.2 raw_data sub-directory ####

**File List & Relationship Between Files**

We included the raw data (as analyzed by Dr. Christian von Sperber and Carina Popp) in this sub-directory for reference. Further description of these files is contained in the [`README.md` file](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/data/raw_data/README.md) with this sub-directory as well as the README tabs of each file. 

This sub-directory contains the following three files:

* `Frankenforst d18O values.xlxs` - Raw d18O-phosphate data. See [d18o_data.csv](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/data/d18o_data.csv) for the final tidied data version.

* `Hedley field fresh vs dry inorganic P pools.xlsx` - Raw Hedley inorganic P pools data. See [hedley_data.csv](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/data/hedley_data.csv) for the final tided data version.

* `Hedley field fresh vs dry organic P pools.xlsx` - Raw Hedley organic P pools data. See [hedley_data.csv](https://github.com/sheilasaia/paper-d18o-phosphate-riparian/blob/main/data/hedley_data.csv) for the final tided data version.


#### 2.3 figure_table directory ####

We generated the files included in this directory via `analysis_script.R`. Please see that script for more information on these image files. The exception to this is the `table.xlxs` file which was generated manally by Sheila Saia based on model outputs in `analysis_script.R`.

Sheila Saia did some minor manual figure editing in Inkscape (v1.1) to generate .svg and .png files and publication-ready, high resolution (i.e., 300 dpi) pdfs.