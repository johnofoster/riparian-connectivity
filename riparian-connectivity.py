#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
riparian_connectivity.py.

Description
-----------
This script determines the location and extent of riparian connectivity in a watershed
then outputs several riparian connectivity statistics and an interactive map of the 
watershed.

Parameters (Inputs)
------
a projected watershed raster file, 
a watercourse line vector file, 
a waterbodies polygon vector file and 
a watershed polygon vector file
a user defined threshold

Process
-------
The script takes the above inputs, creates a riparian buffer around the waterbodies and 
watercourses. The Normalized Difference Vegetation Index is calculated on the buffer. 
sA threshold is used to separate vegetated areas from non-vegetated areas and statistics
computed on the result. The statistics and vegetated to non-vegetated buffers are
integrated into an interactive map, which aids in visualization of non-vegetated areas/
areas of low riparian connectivity.

Returns (Outputs)
-------
Riparian Statistics
Interactive Map

About
-----
This script was created by a group of student from Carleton University, in Ottawa, 
Canada and was originally written for the Ottawa Riverkeeper to aid in their assessment
of riparian connectivity, one of their fourteen (14) indicators of watershed health

Authors: John Foster (Lead), Benjamin Colbourne, Taji Hamilton, Haley Nicholson
"""

# %% Import modules

# built-ins
from datetime import datetime
import os
from pathlib import Path
import sys

# data manipulation
import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rxr

# plotting
import matplotlib.pyplot as plt
import folium
from folium import plugins

# raster to vector function
from ast import literal_eval
from joblib import Parallel, delayed
import multiprocessing
from rasterio import features
from shapely.geometry import shape

# Suppress FutureWarning
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


# %% User input function - John


def user_input(input_message):
    """Get user input or q to quit."""
    raw_input = input(input_message)
    if raw_input.lower() == "q" or raw_input.lower() == "quit":
        print("\nGoodbye!")
        sys.exit()
    elif raw_input.lower() == "":
        return None
    else:
        return raw_input


# %% 1) Function to read the vector and raster data - John


def load_data_ui():
    """
    Provide a textual user interface for the loading of the required data.

    Additionally, prevent errors later in the script by controlling for the following:

        Watershed boundary:
            - Geometry type must be 'Polygon'
            - Maximum of 1 feature

        Water bodies
            - Geometry type must be 'Polygon'

        Water courses
            - Geometry type must be 'LineString'

        Sentinel-2 imagery
            - Must have a corodinate reference system (CRS)
            - CRS must be projected
            - CRS must use linear units of 'metres'
            - Must have band 4 (Red)
            - Must have band 8 (NIR)

    Returns
    -------
    data_dict : dictionary
        Returns a dictionary with each of the data objects as the value of an
        appropriately named key.
    """
    print("\nRiparian Connectivity")
    print("-------------------------------------------------")
    print("Please enter required information or 'q' to quit.")
    print("-------------------------------------------------")

    # Get the watershed name
    loaded = False
    while loaded == False:
        watershed_name = user_input("Watershed name: ")
        if watershed_name:
            loaded = True

    # Read the watershed boundary vector data
    loaded = False
    while loaded == False:
        input_path = user_input(
            "\nPath to the watershed boundary polygon shapefile/geopackage: "
        )
        if input_path:
            try:
                watershed_path = os.path.join(*Path(input_path).parts)
                print("\nReading file...", end="")
                watershed_gdf = gpd.read_file(watershed_path)
                # Check for errors
                if watershed_gdf["geometry"].any().geom_type not in [
                    "Polygon",
                    "MultiPolygon",
                ]:
                    raise Exception(
                        "Error: Wrong geometry type. Expected 'Polygon' or "
                        "'MultiPolygon' but received"
                        f"{watershed_gdf['geometry'].all().geom_type}."
                    )
                elif len(watershed_gdf) > 1:
                    raise Exception(
                        "Error: Too many features. Expected 1 feature but"
                        f" received {len(watershed_gdf)} features."
                    )
                else:
                    print("done")
                    loaded = True
            except Exception as e:
                print("\n")
                print(e, "\n")
                print(
                    "Watershed boundary polygon file not loaded successfully."
                    " Please enter a path to a valid file or 'q' to quit."
                )

    # Read the water bodies vector data
    loaded = False
    while loaded == False:
        input_path = user_input(
            "\nPath to the water bodies polygon shapefile/geopackage: "
        )
        if input_path:
            try:
                waterbodies_path = os.path.join(*Path(input_path).parts)
                print("\nReading file...", end="")
                waterbodies_gdf = gpd.read_file(waterbodies_path)
                # Check for errors
                if waterbodies_gdf["geometry"].any().geom_type not in [
                    "Polygon",
                    "MultiPolygon",
                ]:
                    raise Exception(
                        "Error: Wrong geometry type. Expected 'Polygon' but"
                        f" received '{waterbodies_gdf['geometry'].all().geom_type}'."
                    )
                else:
                    print("done")
                    loaded = True
            except Exception as e:
                print("\n")
                print(e, "\n")
                print(
                    "Water bodies boundary polygon file not loaded successfully."
                    " Please enter a path to a valid file or 'q' to quit."
                )

    # Read the water courses vector data
    loaded = False
    while loaded == False:
        input_path = user_input(
            "\nPath to the water courses line shapefile/geopackage: "
        )
        if input_path:
            try:
                watercourses_path = os.path.join(*Path(input_path).parts)
                print("\nReading file...", end="")
                watercourses_gdf = gpd.read_file(watercourses_path)
                # Check for errors
                if watercourses_gdf["geometry"].any().geom_type not in [
                    "LineString",
                    "MultiLineString",
                ]:
                    raise Exception(
                        "Error: Wrong geometry type. Expected 'LineString' but"
                        f" received '{watercourses_gdf['geometry'].all().geom_type}'."
                    )
                else:
                    print("done")
                    loaded = True
            except Exception as e:
                print("\n")
                print(e, "\n")
                print(
                    "Water courses boundary line file not loaded successfully."
                    " Please enter a path to a valid file or 'q' to quit."
                )

    # Read the Sentinel-2 multispectral GeoTiff imagery file
    loaded = False
    while loaded == False:
        input_path = user_input(
            "\nPath to the Sentinel-2 multispectral GeoTiff imagery file: "
        )
        if input_path:
            try:
                imagery_path = os.path.join(*Path(input_path).parts)
                print("\nReading file...", end="")
                imagery_da = rxr.open_rasterio(filename=imagery_path)
                imagery_crs = imagery_da.rio.crs
                # Check for errors
                if imagery_crs == None:
                    raise Exception(
                        "Error: Missing coordinate reference system. Expected imagery"
                        f" with a valid CRS but received '{imagery_crs}'."
                    )
                elif imagery_crs.is_projected == False:
                    raise Exception(
                        "Error: Coordinate reference system is not projected. Expected"
                        f" imagery with a projected CRS but received '{imagery_crs}'."
                    )
                elif imagery_crs.linear_units != "metre":
                    raise Exception(
                        "Error: Coordinate reference system linear units are not in"
                        " meters. Expected imagery with linear units of 'metre' but"
                        f" received '{imagery_crs.linear_units}'."
                    )

                elif 4 not in imagery_da["band"]:
                    raise Exception(
                        "Error: Missing Band 4 (Red). Band 4 required for the"
                        " calculation of the Normalized Vegetation Difference Index."
                    )
                elif 8 not in imagery_da["band"]:
                    raise Exception(
                        "Error: Missing Band 8 (NIR). Band 8 required for the"
                        " calculation of the Normalized Vegetation Difference Index."
                    )
                else:
                    print("done")
                    loaded = True
            except Exception as e:
                print("\n")
                print(e, "\n")
                print(
                    "Sentinel-2 multispectral imagery file not loaded successfully."
                    " Please enter a path to a valid file or 'q' to quit."
                )

    # Input the buffer width
    loaded = False
    while loaded == False:
        buffer_width = user_input("\nRiparian buffer width in meters: ")
        if buffer_width:
            try:
                buffer_width = float(buffer_width)
                if buffer_width > 0:
                    loaded = True
                else:
                    raise Exception
            except:
                print("\nError: Buffer width must be a valid integer or float > 0\n")

    print("\nAll input data loaded successfully.\n")

    # Create a directory to hold the results
    try:
        os.mkdir("results")
    except FileExistsError:
        pass
    os.chdir("results")

    # Create a unique directory to hold the results from the current analysis
    # Modified from https://stackoverflow.com/a/56680778
    results_dir_name = watershed_name.replace(" ", "_")
    try:
        os.mkdir(f"{results_dir_name}")
    except FileExistsError:
        counter = 2
        while os.path.isdir(f"{results_dir_name}_{counter}"):
            counter += 1
        os.mkdir(f"{results_dir_name}_{counter}")
        results_dir_name = f"{results_dir_name}_{counter}"

    print("\nA log file and outputs will be found here:")
    print(os.path.abspath(results_dir_name) + "\n")

    user_input("Press 'Enter' or 'Return' to proceed or 'q' to quit: ")

    # Change the working directory to the results directory for the current session
    os.chdir(results_dir_name)

    # Write initial entries to the log file
    log_filename = watershed_name.replace(" ", "_") + "-log.txt"
    log_filepath = os.path.abspath(log_filename)
    with open(log_filepath, "a") as file:
        file.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S\n\n"))
        file.write(f"Riparian connectivity log:\n{watershed_name}\n\n")
        file.write(f"Results directory:\nf{os.path.abspath(results_dir_name)}\n\n")
        file.write(f"Watershed boundary input:\n{os.path.abspath(watershed_path)}\n\n")
        file.write(f"Water bodies input:\n{os.path.abspath(waterbodies_path)}\n\n")
        file.write(f"Water courses input:\n{os.path.abspath(watercourses_path)}\n\n")
        file.write(f"Imagery input:\n{os.path.abspath(imagery_path)}\n\n")
        file.write(f"Imagery coordinate reference system: {imagery_crs}\n\n")
        file.write(f"Riparian buffer width (m) input: {buffer_width}\n\n")

    # Add the inputs to a dictionary
    data_dict = {
        "watershed_name": watershed_name,
        "watershed": watershed_gdf,
        "waterbodies": waterbodies_gdf,
        "watercourses": watercourses_gdf,
        "imagery": imagery_da,
        "imagery_crs": imagery_crs,
        "buffer_width": buffer_width,
        "log_filepath": log_filepath,
    }

    return data_dict


# %% 2) Function to perform vector operations to create riparian zone - Ben


def vector_operations(
    watershed_gdf,
    waterbodies_gdf,
    watercourses_gdf,
    imagery_crs,
    buffer_width,
    log_filepath,
):
    """
    Perform vector operations to create a riparian buffer.

    Create a riparian buffer from the watershed's waterbodies and watercourses.

    Parameters
    ----------
    watershed_gdf : GeoPandas GeoDataFrame
        The watershed boundary. Geometry type must be Polygon.

    waterbodies_gdf : GeoPandas GeoDataFrame
        The water bodies within the watershed. Geometry type must be Polygon.

    watercourses_gdf : GeoPandas GeoDataFrame
        The water courses within the watershed. Geometry type must be LineString.

    imagery_crs : CRS
        The coordinate reference system of the Sentinel-2 imagery. The GeoDataFrames
        will be reprojected to this CRS.

    buffer_width : float
        The width in meters of the riparian buffer. The buffer will be applied to the
        water bodies and water courses and will be clipped to the watershed boundary.

    log_filepath : str
        Path to the log file.

    Returns
    -------
    riparian_buff_geom : GeoPandas GeoSeries
        The geometry of the watershed's riparian buffer.

    """
    print("\nCreating riparian buffer...", end="")

    with open(log_filepath, "a") as file:
        file.write(
            "Computations started @" f"{datetime.now().strftime('%H:%M:%S')}\n\n"
        )

    # Reproject the GeoDataFrames
    watershed_gdf = watershed_gdf.to_crs(imagery_crs)
    waterbodies_gdf = waterbodies_gdf.to_crs(imagery_crs)
    watercourses_gdf = watercourses_gdf.to_crs(imagery_crs)

    # Dissolve the GeoDataFrames
    watershed_gdf = watershed_gdf[["geometry"]].dissolve()
    waterbodies_gdf = waterbodies_gdf[["geometry"]].dissolve()
    watercourses_gdf = watercourses_gdf[["geometry"]].dissolve()

    # Access the dissolved geometries
    watershed_geom = watershed_gdf["geometry"]
    waterbodies_geom = waterbodies_gdf["geometry"]
    watercourses_geom = watercourses_gdf["geometry"]

    # Create water bodies and water courses buffer geometry
    waterbodies_buff_geom = waterbodies_geom.buffer(buffer_width)
    watercourses_buff_geom = watercourses_geom.buffer(buffer_width)

    # Perform union between the water bodies and water courses buffer geometry
    water_buff_geom = waterbodies_buff_geom.union(other=watercourses_buff_geom)

    # Clip to watershed boundary
    water_buff_geom = water_buff_geom.clip(mask=watershed_geom)

    # Get difference of unioned buffers and water bodies to get the outer buffer
    riparian_buff_geom = water_buff_geom.difference(waterbodies_geom)

    print("done\n")

    # Write intermediate data to file
    print("Exporting riparian buffer...", end="")

    # Create new GeoDataFrame containing the difference of the geometries
    # NOTE: This GeoDataFrame is just for saving this intermediate result
    riparian_buff_gdf = gpd.GeoDataFrame(geometry=riparian_buff_geom)

    riparian_buff_gdf.to_file(filename="1_riparian_buffer.gpkg", driver="GPKG")
    print("done\n")

    print("Riparian buffer can be found here:")
    print(os.path.abspath("1-riparian_buffer.gpkg") + "\n")

    # Write a log entry
    with open(log_filepath, "a") as file:
        file.write(
            "Riparian buffer completed @ " f" {datetime.now().strftime('%H:%M:%S')}\n\n"
        )

    return riparian_buff_geom


# %% 3) Function to perform NDVI image processing - Haley


def create_ndvi(imagery_da, riparian_buff_geom, log_filepath):
    """
    Create a Normalized Difference Vegetation Index (NDVI) image of the riparian buffer.

    Parameters
    ----------
    imagery_da : RioXarray DataArray
        Sentinel-2 DataArray where band 4 is the Red band and band 8 is the NIR band.

    riparian_buff_geom : GeoPandas GeoSeries
        Geometry of the riparian buffer.

    log_filepath : str
        Path to the log file.

    Returns
    -------
    ndvi_da : RioXarray DataArray
        Normalized Difference Vegetation Index (NDVI) image of the riparian buffer.

        NDVI = (NIR band - Red band) / (NIR band + Red band)
    """
    print("\nCreating NDVI image of the riparian buffer...")
    # Set the nodata value to be nan
    imagery_da.attrs["_FillValue"] = np.nan

    # Clip the DataArray
    riparian_buff_da = imagery_da.rio.clip(geometries=riparian_buff_geom)

    # Select the Red band and keep only its long name attribute
    red_da = riparian_buff_da.sel(band=4)
    red_da.attrs["long_name"] = red_da.attrs["long_name"][3]

    # Select the NIR band and keep only its long name attribute
    nir_da = riparian_buff_da.sel(band=8)
    nir_da.attrs["long_name"] = nir_da.attrs["long_name"][7]

    # Calculate NDVI
    ndvi_da = (nir_da - red_da) / (nir_da + red_da)
    ndvi_da.attrs["long_name"] = "NDVI (Normalized Difference Vegetation Index)"

    print("done\n")

    # Write intermediate data to file
    print("Exporting NDVI image of riparian buffer...", end="")
    ndvi_da.rio.to_raster("2-riparian_buffer-NDVI.tiff")
    print("done\n")

    print("Riparian buffer NDVI image can be found here:")
    print(os.path.abspath("2-riparian_buffer-NDVI.tiff") + "\n")

    # Write a log entry
    with open(log_filepath, "a") as file:
        file.write(
            "Riparian buffer NDVI image completed @"
            f" {datetime.now().strftime('%H:%M:%S')}\n\n"
        )
    # ---------------------------------------------------------------------

    return ndvi_da


# %% 4) Function to create Riparian Vegetation DataArray - Taji


def create_binary_riparian_da(ndvi_da, log_filepath):
    """
    Create a binary DataArray of riparian vegetation and not-vegetation.

    Use the user input NDVI threshold as the break point. Pixel values of the NDVI
    DataArray that are equal to or over the threshold will be considered vegetation.
    Pixel values of the NDVI DataArray that are under the threshold will be considered
    not-vegetation.

    Parameters
    ----------
    ndvi_da : RioXarray DataArray
        Normalized Difference Vegetation Index (NDVI) image of the riparian buffer.

    log_filepath : str
        Path to the log file.

    Returns
    -------
    riparian_da : RioXarray DataArray

        Riparian buffer pixels classified into vegetation (1) and non-vegetation (2).
        Nodata is 0.

    """
    # Plot riparian NDVI histogram
    ndvi_da.plot.hist(bins=50, figsize=(10, 10))
    plt.title("Histogram of Riparian Zone NDVI Values")
    plt.savefig("3-riparian_buffer-NDVI_histogram.png")

    # Provide threshold suggestions to the user
    print("A Histogram of the riparian buffer NDVI values has been saved here:")
    print(os.path.abspath("4_NDVI_histogram.png") + "\n")

    # Input the NDVI threshold
    loaded = False
    while loaded == False:
        threshold = user_input("Please enter the NDVI threshold: ")
        try:
            threshold = float(threshold)
            if threshold > 0 and threshold < 1:
                loaded = True
            else:
                raise Exception
        except:
            print("\nError: NDVI threshold must be a decimal number between 0 and 1")

    # Write a log entry
    with open(log_filepath, "a") as file:
        file.write(f"NDVI threshold input: {threshold}\n\n")

    # Specify a threshold and create a binary riparian vegetation DataArray
    # ------------------------------------------------------------------------

    print("Creating image of riparian buffer vegetation...", end="")

    # Boolean vegetation DataArray
    veg_da = ndvi_da >= threshold

    # Reassign vegetation pixels from True to 1
    veg_da.values = np.where(veg_da.values == True, 1, 0)
    veg_da = veg_da.astype("uint8")

    # Boolean not-vegetation DataArray
    not_veg_da = ndvi_da < threshold

    # Reassign not-vegetation pixels from True to 2
    not_veg_da.values = np.where(not_veg_da.values == True, 2, 0)
    not_veg_da = not_veg_da.astype("uint8")

    # Create a riparian vegetation/not vegeation DataArray
    riparian_da = veg_da + not_veg_da

    print("done\n")

    # Write intermediate data to file
    print("Exporting image of riparian buffer vegetation...", end="")
    riparian_da.rio.to_raster("3-riparian_buffer-vegetation.tiff")
    print("done\n")

    print("Riparian buffer vegetation image can be found here:")
    print(os.path.abspath("3-riparian_buffer-vegetation.tiff") + "\n")

    # Write a log entry
    with open(log_filepath, "a") as file:
        file.write(
            "Riparian buffer vegetation image completed @"
            f" {datetime.now().strftime('%H:%M:%S')}\n\n"
        )

    return riparian_da, threshold


# %% 5) Function to convert Riparian Vegetation DataArray to GeoDataFrame - John


def extract_raster_features(da, log_filepath, n_jobs=-1):
    """
    Convert a RioXarray DataArray to a GeoPandas GeoDataFrame.

    Contiguous pixels of the same value are turned into polygons with a column named
    'value' which represents the values of the source pixels. The boundary of the
    polygons matches the boundary of the pixels.

    Based on: https://pysal.org/tobler/generated/tobler.dasymetric.extract_raster_features.html

    Parameters
    ----------
    da : DataArray
        The RioXarray DataArray to be converted to a GeoPandas GeoDataFrame. Must
        have a valid Rasterio CRS (e.g. da.rio.crs).

    n_jobs : int, optional
        The default is -1. This will then count the number of available processes for
        multiprocessing.

    log_filepath : str
        Path to the log file.

    Returns
    -------
    gdf : GeoDataFrame
        A GeoDataFrame containing Shapely polygons and pixel values.
    """
    # Print initial log message
    print("Converting the riparian vegetation DataArray to a GeoDataFrame...", end="")

    def _chunk_dfs(geoms_to_chunk, n_jobs):
        chunk_size = geoms_to_chunk.shape[0] // n_jobs + 1
        for i in range(n_jobs):
            start = i * chunk_size
            yield geoms_to_chunk.iloc[start : start + chunk_size]

    def _apply_parser(df):
        return df.apply(_parse_geom)

    def _parse_geom(geom_str):
        return shape(literal_eval(geom_str))

    # CRS of the DataArray
    raster_crs = da.rio.crs

    # Convert regions to polygons
    shapes = list(features.shapes(da.values, transform=da.rio.transform()))

    # Get number of processes for multiprocessing
    if n_jobs == -1:
        n_jobs = multiprocessing.cpu_count()

    res = list(zip(*shapes))
    geoms = pd.Series(res[0], name="geometry").astype(str)
    pieces = _chunk_dfs(geoms, n_jobs)
    geoms = pd.concat(
        Parallel(n_jobs=n_jobs)(delayed(_apply_parser)(i) for i in pieces)
    )

    geoms = gpd.GeoSeries(geoms).buffer(0)  # we sometimes get self-intersecting rings
    vals = pd.Series(res[1], name="value")
    riparian_gdf = gpd.GeoDataFrame(vals, geometry=geoms, crs=raster_crs)

    # Write a log entry
    with open(log_filepath, "a") as file:
        file.write(
            "Riparian buffer vegetation image converted to vector @"
            f" {datetime.now().strftime('%H:%M:%S')}\n\n"
        )

    # Subset by vegetation and not-vegetation features
    vegetation_gdf = riparian_gdf.loc[riparian_gdf["value"] == 1]
    not_vegetation_gdf = riparian_gdf.loc[riparian_gdf["value"] == 2]

    # Dissolve all riparian features into unique contiguous features
    riparian_buffer_gdf = (
        riparian_gdf.loc[(riparian_gdf["value"] == 1) | (riparian_gdf["value"] == 2)]
        .dissolve()
        .explode(ignore_index=True)
    )

    # Print completed log message
    print("done\n")

    # Save the results to a dictionary

    riparian_dict = {
        "riparian_buffer_gdf": riparian_buffer_gdf,
        "vegetation_gdf": vegetation_gdf,
        "not_vegetation_gdf": not_vegetation_gdf,
    }

    # Write intermediate data to file
    print("Exporting riparian buffer vegetation geopackages...", end="")
    vegetation_gdf.to_file(filename="4-riparian_buffer-vegetation.gpkg", driver="GPKG")
    not_vegetation_gdf.to_file(
        filename="4-riparian_buffer-not_vegetation.gpkg", driver="GPKG"
    )
    riparian_buffer_gdf.to_file(filename="4-riparian_buffer.gpkg", driver="GPKG")

    print("done\n")

    return riparian_dict


# %% 6) Function to calculate Riparian Connectivity Statistics - John / team


def riparian_stats(
    watershed_name,
    buffer_width,
    ndvi_threshold,
    watershed_gdf,
    riparian_buffer_gdf,
    vegetation_gdf,
    not_vegetation_gdf,
    log_filepath,
):
    """
    Generate statistics related to the riparian buffer.

    Parameters
    ----------
    watershed_name : str
        Name of the watershed being evaluated.

    buffer_width : float
        The width in meters of the riparian buffer. Included in the table of statistics
        for reference purposes.

    ndvi_threshold : float
        The NDVI threshold used to classify pixels as vegetation and not-vegetation.
        Included in the table of statistics for reference purposes.

    watershed_gdf : GeoPandas GeoDataFrame
        The watershed being evaluated. Needed for the calculation of the watershed's
        total area.

    riparian_buffer_gdf : GeoPandas GeoDataFrame
        The riparian buffer features. Needed for the calculation of the number of
        riparian buffer features, area, and perimeter.

    vegetation_gdf : GeoPandas GeoDataFrame
        The vegetation features that fall within the riparian buffer. Needed for the
        calculation of the number of riparian vegetation features, area, and perimeter.

    not_vegetation_gdf : GeoPandas GeoDataFrame
        The not-vegetation features that fall within the riparian buffer. Needed for the
        calculation of the number of riparian not-vegetation features, area, and
        perimeter.

    log_filepath : str
        Path to the log file.


    Returns
    -------
    stats_df : Pandas DataFrame
        A DataFrame with the following column names and associated values:

        'Watershed name'
            The name of the watershed being evaluated.

        'Watershed area (km2)'
            The area of the watershed (km2).

        'Riparian buffer area (km2)'
            The area of the riparian buffer (km2).

        'Vegetation area (km2)'
            The area of the riparian buffer classified as vegetation (km2).

        'Not-vegetation area (km2)'
            The area of the riparian buffer classified as not-vegetation (km2).

        'Vegetation coverage (%)'
            The percentage of the riparian buffer's area that is covered by vegetation.
            (Vegetation area / Riparian area) * 100

        'Not-vegetation coverage (%)'
            The percentage of the riparian buffer's area that is not covered by
            vegetation.
            (Not-vegetation area / Riparian area) * 100

        'Mean area of not-vegetation patches (km2)'
            The mean area (km2) of the not-vegetation features in the riparian buffer.

        'Number of riparian buffer features'
            The number of separate riparian buffer features.

        'Number of vegetation features'
            The number of separate vegetation features in the riparian buffer.

        'Number of not-vegetation features'
            The number of separate not-vegetation features in the riparian buffer.

        'Perimeter of riparian buffer (km)'
            The total perimeter (km) of the riparian buffer features.

        'Perimeter of vegetation features (km)'
            The total perimeter (km) of the vegetation features.

        'Perimeter of not-vegetation features (km)'
            The total perimeter (km) of the not-vegetation features.

        'Vegetation Connectivity'
            A rough metric to quantify the connectivity of the riparian features.
            Number of riparian buffer features / Number of vegetation features

        'Vegetation Compactness'
            An attempt to quantify the the compactness of the riparian vegetation
            features. The long name of this could be 'normalized isoperimetric ratio'.
            See script README.md for more info.

            # Reference compactness i.e. riparian buffer compactness
            riparian_compactness =  riparian_perimeter**2 / riparian_area

            # Actual compactness i.e. riparian vegetation compactness
            vegetation_compactness = veg_perimeter**2 / veg_area

            # Normalized compactness
            compactness =  riparian_compactness / actual_compactness

    """
    # Print initial log message
    print("Calculating riparian statistics...", end="")

    # Area of features
    veg_area = vegetation_gdf["geometry"].area.sum() / 1_000_000
    not_veg_area = not_vegetation_gdf["geometry"].area.sum() / 1_000_000
    riparian_area = riparian_buffer_gdf["geometry"].area.sum() / 1_000_000
    watershed_area = watershed_gdf["geometry"].area.sum() / 1_000_000

    # Mean patch size of not-vegetation features
    not_veg_mean_size = not_vegetation_gdf["geometry"].area.mean() / 1_000_000

    # Perimeter of features
    veg_perimeter = vegetation_gdf["geometry"].length.sum() / 1000
    not_veg_perimeter = not_vegetation_gdf["geometry"].length.sum() / 1000
    riparian_perimeter = riparian_buffer_gdf["geometry"].length.sum() / 1000

    # Number of features
    veg_n_feature = len(vegetation_gdf)
    not_veg_n_feature = len(not_vegetation_gdf)
    riparian_n_feature = len(riparian_buffer_gdf)  # Not used yet

    # Coverage of total riparian area
    veg_coverage = veg_area / riparian_area * 100
    not_veg_coverage = not_veg_area / riparian_area * 100

    # Connectivity
    connectivity = riparian_n_feature / veg_n_feature

    # Compactness aka normalized isoperimetric ratio
    # Adapted from p21 of: http://www.cyto.purdue.edu/cdroms/micro2/content/education/wirth10.pdf
    riparian_compactness = riparian_perimeter ** 2 / riparian_area
    vegetation_compactness = veg_perimeter ** 2 / veg_area
    compactness = riparian_compactness / vegetation_compactness

    data = {
        "Watershed name": [watershed_name],
        "Buffer width:": [buffer_width],
        "NDVI threshold": [ndvi_threshold],
        "Watershed area (km2)": [watershed_area],
        "Riparian buffer area (km2)": [riparian_area],
        "Vegetation area (km2)": [veg_area],
        "Not-vegetation area (km2)": [not_veg_area],
        "Vegetation coverage (%)": [veg_coverage],
        "Not-vegetation coverage (%)": [not_veg_coverage],
        "Mean area of not-vegetation patches (km2)": [not_veg_mean_size],
        "Number of riparian buffer features": [riparian_n_feature],
        "Number of vegetation features": [veg_n_feature],
        "Number of not-vegetation features": [not_veg_n_feature],
        "Perimeter of riparian buffer (km)": [riparian_perimeter],
        "Perimeter of vegetation features (km)": [veg_perimeter],
        "Perimeter of not-vegetation features (km)": [not_veg_perimeter],
        "Vegetation Connectivity": [connectivity],
        "Vegetation Compactness": [compactness],
    }

    stats_df = pd.DataFrame.from_dict(data=data)

    print("done")

    # Write a log entry
    with open(log_filepath, "a") as file:
        file.write(
            "Riparian statistics completed @"
            f" {datetime.now().strftime('%H:%M:%S')}\n\n"
        )
        file.write("Results:\n")
        file.write(stats_df.transpose().to_string())

    # FOR DEBUGGING ------------------------------------------------------
    # Write intermediate data to file for debugging
    print("\nExporting results...", end="")

    stats_df.to_csv("4-statistics_table.csv")

    print("done\n")

    return stats_df


# %% 7) Function to produce a report - John & Haley


def report(stats_df, vegetation_gdf, not_vegetation_gdf, log_filepath):
    """
    Write an HTML report to the results directory.

    The report contains the statistical results and an interactive folium map.

    Parameters
    ----------
    stats_df : Pandas DataFrame
        The DataFrame containing the riparian connectivity statistics. Displayed as an
        HTML table in the report.

    vegetation_gdf : GeoPandas GeoDataFrame
        The GeoDataFrame containing the riparian buffer vegetation features. Displayed
        in an interactive Folium map.

    not_vegetation_gdf : GeoPandas GeoDataFrame
        The GeoDataFrame containing the riparian buffer not-vegetation features.
        Displayed in an interactive Folium map.

    log_filepath : str
        Path to the log file.

    Returns
    -------
    None.

    """
    vegetation_gdf["value"] = "vegetation"
    not_vegetation_gdf["value"] = "not-vegetation"

    m = vegetation_gdf.explore(
        width="85%",
        color="green",
        tooltip="value",
        name="Vegetation",
        legend_kwds={"caption": "Riparian vegetation"},
        style_kwds={"stroke": False, "opacity": 0.3},
    )

    not_vegetation_gdf.explore(
        m=m,
        color="black",
        tooltip="value",
        name="Not-vegetation",
        style_kwds={
            "stroke": True,
            "opacity": 1,
            "fillColor": "black",
            "weight": 1,
        },
    )

    minimap = plugins.MiniMap()
    m.add_child(minimap)
    folium.LayerControl().add_to(m)

    # 1. Set up multiple variables to store the titles, text within the report
    page_title_text = "Riparian Connectivity"
    title_text = "Riparian Connectivity in the Ottawa River Watershed"

    stats_text = "Riparian Connectivity Summary Statistics"
    stats_description = "Statistic Name Descriptions"

    # 2. Combine them together using a long f-string

    html = f"""
    <html>
        <head>
            <title>{page_title_text}</title>
        </head>
        <body>
        <body style="background-color:WhiteSmoke;">
            <h1><b><u>{title_text}</u></b></h1>
            
            <h4>{stats_text}</h4>
           {stats_df.transpose().to_html()}
          
            <h4><u>{stats_description}</u></h4>
         
             <ul>
      <li>Watershed area (km2) - The area of the watershed (km2).</li>
      <li>Riparian buffer area (km2) - The area of the riparian buffer.</li>
      <li>Vegetation area (km2) - The area of the riparian buffer classified as vegetation.</li>
      <li>Not-vegetation area (km2) - The area of the riparian buffer classified as not-vegetation</li>
      <li>Vegetation coverage (%) - The percentage of the riparian buffer's area that is covered by vegetation</li>
      <li>Not-vegetation coverage (%) - The percentage of the riparian buffer's area that is not covered by vegetation</li>
      <li>Mean area of not-vegetation patches (km2) - The mean area (km2) of the not-vegetation features in the riparian buffer</li>
      <li>Number of riparian buffer features - The number of separate riparian buffer features</li>
      <li>Number of vegetation features - The number of separate vegetation features in the riparian buffer</li>
      <li>Number of not-vegetation features - The number of separate not-vegetation features in the riparian buffer</li>
      <li>Perimeter of riparian buffer (km) - The total perimeter (km) of the riparian buffer features.</li>
      <li>Perimeter of riparian buffer (km) - The total perimeter (km) of the riparian buffer features.</li>
      <li>Perimeter of not-vegetation features (km) - The total perimeter of not-vegetation features.</li>
      <li>Vegetation Connectivity - The quotient of the number of riparian buffer features divided by the number of vegetation features.</li>
      <li>Vegetation Compactness - The normalized isoperimetric ratio (see README.md for more info)</li>
        </ul>
        </head>

    <body>
    
    <h4><u>Map Legend</u></h4>
    
    <table>
    <table style="width:15%"
      <tr>
        <th>Class</th>
        <th>Colour</th>
      </tr>
      
      
      <tr>
        <td>Vegetation</td>
        <td>Green</td>
      </tr>
      <tr>
        <td>Non-Vegetation</td>
        <td>Black</td>
      </tr>
    </table>
    
    </body>
            
            {m.get_root().render()}
            </body>
      
        </html>
    """

    # 3. Write the html string as an HTML file
    with open("5-report.html", "w") as f:
        f.write(html)


# %% def main function


def main():
    """
    Is the main function of the script.

    Calls all other functions and executes the program.

    Parameters (inputs)
    -------------------
    None.

    Returns (outputs)
    -----------------
    None.
    """
    # 1) Function to read the vector and raster data - John
    data_dict = load_data_ui()

    # 2) Function to perform vector operations to create riparian zone - Ben
    riparian_buff_geom = vector_operations(
        watershed_gdf=data_dict["watershed"],
        waterbodies_gdf=data_dict["waterbodies"],
        watercourses_gdf=data_dict["watercourses"],
        imagery_crs=data_dict["imagery_crs"],
        buffer_width=data_dict["buffer_width"],
        log_filepath=data_dict["log_filepath"],
    )

    # 3) Function to perform NDVI image processing - Haley
    ndvi_da = create_ndvi(
        imagery_da=data_dict["imagery"],
        riparian_buff_geom=riparian_buff_geom,
        log_filepath=data_dict["log_filepath"],
    )

    # 4) b. Function to create Riparian Vegetation DataArray - Taji
    riparian_da, ndvi_threshold = create_binary_riparian_da(
        ndvi_da=ndvi_da,
        log_filepath=data_dict["log_filepath"],
    )

    # 5) Function to convert Riparian Vegetation DataArray to GeoDataFrame - John
    riparian_dict = extract_raster_features(
        da=riparian_da, log_filepath=data_dict["log_filepath"]
    )

    # 6) Function to calculate Riparian Connectivity Statistics - Taji / John
    stats_df = riparian_stats(
        watershed_name=data_dict["watershed_name"],
        buffer_width=data_dict["buffer_width"],
        ndvi_threshold=ndvi_threshold,
        watershed_gdf=data_dict["watershed"],
        riparian_buffer_gdf=riparian_dict["riparian_buffer_gdf"],
        vegetation_gdf=riparian_dict["vegetation_gdf"],
        not_vegetation_gdf=riparian_dict["not_vegetation_gdf"],
        log_filepath=data_dict["log_filepath"],
    )
    print(stats_df.transpose())

    # 7) Function to produce a report - John & Haley
    report(
        stats_df=stats_df,
        vegetation_gdf=riparian_dict["vegetation_gdf"],
        not_vegetation_gdf=riparian_dict["not_vegetation_gdf"],
        log_filepath=data_dict["log_filepath"],
    )


# %% run main function
if __name__ == "__main__":
    main()
