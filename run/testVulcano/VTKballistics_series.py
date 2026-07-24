#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate VTK files and a ParaView .series file from OpenFOAM cloudInfo CSV data.

This script reads particle data from:
    postProcessing/cloudInfo1/<time>/output.csv

and writes:
    postProcessing/cloudInfo1/particles_<time>.vtk
    postProcessing/cloudInfo1/particles.vtk.series

The .series file is written with a configurable time scaling factor so that
ParaView uses the correct physical time values.

Authors
-------
S. Giansante
M. de' Michieli Vitturi

Affiliation
-----------
INGV Pisa
"""

import json
import os

import pandas as pd


TIME_SCALE = 1.0
CLOUDINFO_SUBPATH = os.path.join("postProcessing", "cloudInfo1")
SERIES_FILENAME = "particles.vtk.series"
REQUIRED_COLUMNS = {"x", "y", "z", "d", "Ux", "Uy", "Uz"}


def is_number(string):
    """
    Return True if the input string can be converted to float, else False.
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def write_vtk(points, diameters, velocities, output_file):
    """
    Write a legacy ASCII VTK POLYDATA file.

    Parameters
    ----------
    points : array-like
        Particle coordinates, shape (n, 3).
    diameters : array-like
        Particle diameters, shape (n,).
    velocities : array-like
        Particle velocity vectors, shape (n, 3).
    output_file : str
        Output VTK file path.
    """
    n_points = len(points)

    with open(output_file, "w", encoding="utf-8") as vtk_file:
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write("Ballistics data\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET POLYDATA\n")
        vtk_file.write(f"POINTS {n_points} float\n")

        for point in points:
            vtk_file.write(f"{point[0]} {point[1]} {point[2]}\n")

        vtk_file.write(f"\nPOINT_DATA {n_points}\n")
        vtk_file.write("SCALARS diameter float 1\n")
        vtk_file.write("LOOKUP_TABLE default\n")

        for diameter in diameters:
            vtk_file.write(f"{diameter}\n")

        vtk_file.write("VECTORS velocity float\n")

        for velocity in velocities:
            vtk_file.write(f"{velocity[0]} {velocity[1]} {velocity[2]}\n")


def process_simulation(sim_path, time_scale=TIME_SCALE):
    """
    Process one simulation directory.

    The function reads all valid time folders inside:
        <sim_path>/postProcessing/cloudInfo1/

    For each folder containing output.csv, it generates a corresponding VTK file
    and adds an entry to the ParaView .series file.

    Parameters
    ----------
    sim_path : str
        Base path of the simulation.
    time_scale : float, optional
        Scaling factor applied to folder times when writing the .series file.
        For example, if folders are 116, 116.1, ..., 120 but the physical time
        is 11.6, 11.61, ..., 12.0, then use time_scale=0.1.
    """
    cloudinfo_path = os.path.join(sim_path, CLOUDINFO_SUBPATH)
    print(f"Searching in: {cloudinfo_path}")

    if not os.path.exists(cloudinfo_path):
        print(f"Warning: cloudInfo1 directory not found in {sim_path}")
        return

    series_files = []

    folders = sorted(
        os.listdir(cloudinfo_path),
        key=lambda name: float(name) if is_number(name) else float("inf"),
    )

    for folder in folders:
        folder_path = os.path.join(cloudinfo_path, folder)

        if not os.path.isdir(folder_path):
            continue

        if not is_number(folder):
            continue

        csv_path = os.path.join(folder_path, "output.csv")
        if not os.path.exists(csv_path):
            print(f"Skipping {folder}: output.csv not found")
            continue

        try:
            dataframe = pd.read_csv(csv_path)
        except Exception as exc:
            print(f"Warning: could not read {csv_path}: {exc}")
            continue

        if not REQUIRED_COLUMNS.issubset(dataframe.columns):
            print(
                f"Warning: missing required columns in {csv_path}. "
                f"Required columns are: {sorted(REQUIRED_COLUMNS)}"
            )
            continue

        points = dataframe[["x", "y", "z"]].values
        diameters = dataframe["d"].values
        velocities = dataframe[["Ux", "Uy", "Uz"]].values

        vtk_filename = f"particles_{folder}.vtk"
        vtk_path = os.path.join(cloudinfo_path, vtk_filename)

        write_vtk(points, diameters, velocities, vtk_path)
        print(f"Written: {vtk_path}")

        raw_time = float(folder)
        physical_time = raw_time * time_scale

        series_files.append(
            {
                "name": vtk_filename,
                "time": physical_time,
            }
        )

    if not series_files:
        print("No VTK files were generated. The .series file was not created.")
        return

    series_path = os.path.join(cloudinfo_path, SERIES_FILENAME)
    series_data = {
        "file-series-version": "1.0",
        "files": series_files,
    }

    with open(series_path, "w", encoding="utf-8") as series_file:
        json.dump(series_data, series_file, indent=4)

    print(f"\nSeries file created: {series_path}")
    print(
        "You can now open 'particles.vtk.series' in ParaView "
        "to visualize the animation with the correct time values."
    )


if __name__ == "__main__":
    BASE_PATH = os.getcwd()
    process_simulation(BASE_PATH, time_scale=TIME_SCALE)
