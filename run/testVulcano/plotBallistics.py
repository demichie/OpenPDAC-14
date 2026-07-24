import os
import re
import glob
import csv
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from linecache import getline
import pandas as pd
from matplotlib.colors import LightSource
from matplotlib.colors import BoundaryNorm
import argparse  # Added for command-line arguments
import sys  # Added for sys.exit with argparse errors
from matplotlib.ticker import FuncFormatter

# from sklearn.preprocessing import StandardScaler
# from sklearn.cluster import DBSCAN

# Global default values that can be overridden by CLI args
DEFAULT_RASTER_RESOLUTION = 50.0  # Was 'step_dens' globally
PARTICLE_IMPACT_VELOCITY_THRESHOLD = 1.0  # Was 'toll' globally
DEFAULT_STOP_DISTANCE_TOLERANCE = 1.0e-6


def log_to_percent_formatter(x, pos):
    """
    Converte un valore log10(percentuale) in una stringa di etichetta
    percentuale.
    x: il valore del tick (log10 della percentuale).
    pos: la posizione del tick (non usata qui).
    """
    if np.isneginf(x):  # Gestisce log10(0)
        return "0%"

    percentage = 10**x

    if percentage == 0:  # Potrebbe accadere per underflow vicino a zero
        return "0%"
    # Formattazione per leggibilità
    elif percentage < 0.001:  # Percentuali molto piccole
        return f"{percentage:.1e}%"  # Es. 1.2e-4%
    elif percentage < 0.1:  # Percentuali piccole
        return f"{percentage:.3f}%"  # Es. 0.012%
    elif percentage < 1:  # Percentuali < 1%
        return f"{percentage:.2f}%"  # Es. 0.12%
    elif percentage < 10:  # Percentuali tra 1% e 10%
        return f"{percentage:.1f}%"  # Es. 5.3%
    elif percentage <= 100:  # Percentuali fino al 100%
        return f"{percentage:.0f}%"  # Es. 75%
    else:  # Non dovrebbe accadere se la percentuale massima è 100
        return f"{percentage:.0f}%"


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def read_csv_to_arrays(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        columns = {field: [] for field in reader.fieldnames}
        for row in reader:
            for key in reader.fieldnames:
                columns[key].append(float(row[key]))
    arrays = {key.strip(): np.array(values) for key, values in columns.items()}
    return arrays


def readASC(DEM_file_path):  # Renamed variable for clarity
    print('Reading DEM file: ' + DEM_file_path)
    hdr = [getline(DEM_file_path, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell_size_dem, nd_val_dem = values
    cols = int(cols)
    rows = int(rows)

    xs_DEM = lx + 0.5 * cell_size_dem + \
        np.linspace(0, (cols - 1) * cell_size_dem, cols)
    ys_DEM = ly + 0.5 * cell_size_dem + \
        np.linspace(0, (rows - 1) * cell_size_dem, rows)
    dem_extent = lx, lx + cols * cell_size_dem, ly, ly + rows * cell_size_dem

    DEM_data = pd.read_table(DEM_file_path,
                             sep=r'\s+',
                             header=None,
                             skiprows=6).astype(float).values

    DEM_data = np.flipud(DEM_data)
    # Assuming 0 is a safe fill value for NODATA in DEM for hillshade
    DEM_data[DEM_data == nd_val_dem] = 0.0

    x_coords_dem_centers = xs_DEM
    y_coords_dem_centers = ys_DEM
    xmin_dem = np.amin(x_coords_dem_centers)
    xmax_dem = np.amax(x_coords_dem_centers)
    print('xmin_dem_center, xmax_dem_center (DEM cell centers):', xmin_dem,
          xmax_dem)
    ymin_dem = np.amin(y_coords_dem_centers)
    ymax_dem = np.amax(y_coords_dem_centers)
    print('ymin_dem_center, ymax_dem_center (DEM cell centers):', ymin_dem,
          ymax_dem)

    X_dem_centers, Y_dem_centers = np.meshgrid(x_coords_dem_centers,
                                               y_coords_dem_centers)
    Z_dem_data = DEM_data
    return X_dem_centers, Y_dem_centers, Z_dem_data, cell_size_dem, dem_extent


def read_topoGridDict(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    xVent = float(re.search(r'xVent\s+([\d\.\-eE]+)', content).group(1))
    yVent = float(re.search(r'yVent\s+([\d\.\-eE]+)', content).group(1))
    raster_match = re.search(r'rasterFile\s+(\S+);', content)
    DEMfile_path_from_dict = raster_match.group(
        1) if raster_match else None  # Renamed variable
    return xVent, yVent, DEMfile_path_from_dict


# --- Command Line Argument Parsing Functions ---


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=("Process OpenFOAM Lagrangian particle data to create "
                     "ballistic impact maps."))
    parser.add_argument(
        "--lx",
        type=float,
        default=None,
        help="Lower-left x-coordinate of the output raster grid (optional).")
    parser.add_argument(
        "--ly",
        type=float,
        default=None,
        help="Lower-left y-coordinate of the output raster grid (optional).")
    parser.add_argument(
        "--res",
        type=float,
        default=None,
        help="Output raster grid resolution (cell size) (optional).")
    parser.add_argument(
        "--nrows",
        type=int,
        default=None,
        help="Number of rows in the output raster grid (optional).")
    parser.add_argument(
        "--ncols",
        type=int,
        default=None,
        help="Number of columns in the output raster grid (optional).")
    parser.add_argument(
        "--sdt", "--stop_distance_tolerance",
        type=float,
        default=DEFAULT_STOP_DISTANCE_TOLERANCE,
        help=("Distance tolerance [m] used to decide whether a particle has "
              "the same final position in consecutive outputs. Default: "
              f"{DEFAULT_STOP_DISTANCE_TOLERANCE:g}."))
    parser.add_argument(
        "--map-mode",
        choices=["deposited", "depositedEscaped"],
        default="depositedEscaped",
        help=("Particle statuses to include in raster maps. Use 'deposited' "
              "to map only deposited particles, or 'depositedEscaped' to map "
              "both deposited and escaped particles. Default: depositedEscaped."))
    # --np argument removed as per request
    return parser.parse_args()


def validate_args(args):
    provided = {
        "lx": args.lx is not None,
        "ly": args.ly is not None,
        "res": args.res is not None,
        "nrows": args.nrows is not None,
        "ncols": args.ncols is not None
    }
    total_provided = sum(provided.values())

    if total_provided == 0:
        return "default"
    if total_provided == 1 and provided["res"]:
        return "res_only"
    if all(provided.values()):
        return "full"

    # If none of the above, it's an invalid combination
    error_message = (
        "Invalid combination of grid arguments (--lx, --ly, --res, --nrows, "
        "--ncols).\n"
        "Accepted usage patterns for grid definition:\n"
        "1. No grid arguments: Grid extents derived from DEM, "
        "resolution uses script default.\n"
        "2. Only --res: Grid extents derived from DEM, resolution from "
        "--res.\n"
        "3. All of --lx, --ly, --res, --nrows, --ncols: "
        "Fully user-defined grid.\n")
    # Using sys.exit(error_message) or parser.error(error_message) might
    # be cleaner but for now, a raised ValueError is fine.
    # To use parser.error, parse_arguments would need the parser instance.
    # For now, let's print and raise.
    print("ERROR: " + error_message, file=sys.stderr)
    raise ValueError("Invalid command-line arguments for grid definition.")


# --- Main Function ---


def main():
    args = parse_arguments()
    grid_definition_mode = validate_args(args)
    print(f"Grid definition mode selected: {grid_definition_mode}")

    xVent, yVent, DEMfile_name = read_topoGridDict("system/topoGridDict")

    DEMfile_path = os.path.join(os.getcwd(), "constant", "DEM", DEMfile_name)

    print('DEMfile_path from topoGridDict:', DEMfile_path)

    if not os.path.exists(DEMfile_path):
        raise FileNotFoundError(f"DEM file not found: {DEMfile_path}")
    print("DEMfile_path resolved to:", DEMfile_path)
    
    # X_dem_centers, Y_dem_centers are meshes of DEM cell centers
    # Z_dem_data is the elevation data
    # dem_cell_size is the resolution of the input DEM
    # dem_extent is (xmin_edge, xmax_edge, ymin_edge, ymax_edge) of the DEM
    (
        X_dem_centers,
        Y_dem_centers,
        Z_dem_data,
        dem_cell_size,
        dem_extent,
    ) = readASC(DEMfile_path)

    ls = LightSource(azdeg=45, altdeg=45)
    current_dir = os.getcwd()

    output_raster_root = os.path.join(current_dir, "postProcessing",
                                      "raster_maps", "ballistics")
    os.makedirs(output_raster_root, exist_ok=True)
    postprocessing_dir = os.path.join(current_dir, "postProcessing")
    os.makedirs(postprocessing_dir, exist_ok=True)

    search_pattern = os.path.join(current_dir, "postProcessing", "cloudInfo1",
                                  "*", "output.csv")
    cloud_files = glob.glob(search_pattern)

    if not cloud_files:
        print(f"No particle CSV files found with pattern: {search_pattern}")
        return

    try:
        timesteps = [
            float(os.path.basename(os.path.dirname(file)))
            for file in cloud_files
        ]
    except ValueError as e:
        print(f"Error extracting timesteps from file paths. "
              f"Example file: {cloud_files[0]}")
        print(f"Ensure files are in directories like "
              f"'postProcessing/cloudInfo1/NUMERICAL_TIMESTEP/output.csv'. "
              f"Error: {e}")
        return

    timevals = np.array(timesteps)
    sorted_indices = np.argsort(timevals)
    sorted_files = [cloud_files[i] for i in sorted_indices]
    sorted_times = timevals[sorted_indices]
    n_times = len(sorted_files)

    if n_times == 0:
        print("No particle CSV files found. Exiting.")
        return
    print(f"Found {n_times} timesteps.")

    print("Starting scan for unique particles...")
    all_particle_ids_set = set()
    particle_properties_map = {}

    for filename in sorted_files:
        data_scan = read_csv_to_arrays(filename)
        current_origProc = data_scan['origProc'].astype(int)
        current_origId = data_scan['origId'].astype(int)
        current_d = data_scan['d']
        current_rho = data_scan['rho']

        for j in range(len(current_origId)):
            particle_key = (current_origProc[j], current_origId[j])
            all_particle_ids_set.add(particle_key)
            if particle_key not in particle_properties_map:
                particle_properties_map[particle_key] = {
                    'd': current_d[j],
                    'rho': current_rho[j]
                }

    unique_global_particle_identifiers = sorted(list(all_particle_ids_set))
    particle_to_global_idx_map = {
        pid: i
        for i, pid in enumerate(unique_global_particle_identifiers)
    }

    nballistics_global = len(unique_global_particle_identifiers)
    print(f"Total unique particles found (nballistics_global): "
          f"{nballistics_global}")

    if nballistics_global == 0:
        print("No particles found in any file. Exiting.")
        return

    d_global = np.array([
        particle_properties_map[pid]['d']
        for pid in unique_global_particle_identifiers
    ])
    rho_global = np.array([
        particle_properties_map[pid]['rho']
        for pid in unique_global_particle_identifiers
    ])

    plot_titles = []
    unique_diams_global = np.unique(d_global)
    num_unique_sizes = len(unique_diams_global)

    if num_unique_sizes < 5:
        diams_for_plot_categories = unique_diams_global
        for i in range(len(diams_for_plot_categories) + 1):
            if i < len(diams_for_plot_categories):
                plot_titles.append('Diameter = ' +
                                   str(diams_for_plot_categories[i]) + 'm')
            else:
                plot_titles.append('All sizes')
    else:
        dMin_global = np.amin(d_global)
        dMax_global = np.amax(d_global)
        print('Global dMin:', dMin_global, 'Global dMax:', dMax_global)
        diams_for_plot_categories = np.linspace(dMin_global, dMax_global, 4)

        plot_titles.append('Diameter <= ' +
                           f"{diams_for_plot_categories[0]:.2e}" + 'm')
        for i in range(len(diams_for_plot_categories) - 1):
            plot_titles.append(f"{diams_for_plot_categories[i]:.2e}" +
                               'm < Diameter <= ' +
                               f"{diams_for_plot_categories[i+1]:.2e}" + 'm')
        plot_titles.append('Diameter > ' +
                           f"{diams_for_plot_categories[-1]:.2e}" + 'm')
        plot_titles.append('All sizes')

    print('Diameters for plot categorization (diams_for_plot_categories):',
          diams_for_plot_categories)

    A_velocity = np.full((nballistics_global, 3, n_times), np.nan)
    B_position = np.full((nballistics_global, 3, n_times), np.nan)
    matr_data = np.full((n_times, 8, nballistics_global), np.nan)

    print("Starting population of arrays A_velocity, B_position, matr_data...")
    for i_time, filename in enumerate(sorted_files):
        data = read_csv_to_arrays(filename)
        current_origProc = data['origProc'].astype(int)
        current_origId = data['origId'].astype(int)

        x_curr, y_curr, z_curr = data['x'], data['y'], data['z']
        Ux_curr, Uy_curr, Uz_curr = data['Ux'], data['Uy'], data['Uz']

        for local_particle_idx in range(len(current_origId)):
            proc, oid = current_origProc[local_particle_idx], current_origId[
                local_particle_idx]
            global_idx = particle_to_global_idx_map[(proc, oid)]

            A_velocity[global_idx, :, i_time] = [
                Ux_curr[local_particle_idx], Uy_curr[local_particle_idx],
                Uz_curr[local_particle_idx]
            ]
            B_position[global_idx, :, i_time] = [
                x_curr[local_particle_idx] + xVent,
                y_curr[local_particle_idx] + yVent, z_curr[local_particle_idx]
            ]

            matr_data[i_time, 1:4, global_idx] = B_position[global_idx, :,
                                                            i_time]
            matr_data[i_time, 4:7, global_idx] = A_velocity[global_idx, :,
                                                            i_time]
            if not np.isnan(A_velocity[global_idx, 0, i_time]):
                matr_data[i_time, 7,
                          global_idx] = LA.norm(A_velocity[global_idx, :,
                                                           i_time])

    print("Calculating last-seen status and final-position plateau...")
    stop_distance_tolerance = args.sdt
    print(f"Stop distance tolerance: {stop_distance_tolerance:g} m")

    # For every particle we store its last available record in the CSV series.
    # The status column then distinguishes particles that are still present at
    # the final output from particles that disappeared earlier, which are
    # interpreted as escaped/removed from the cloud.
    last_seen_idx = np.full(nballistics_global, -1, dtype=int)
    first_valid_idx = np.full(nballistics_global, -1, dtype=int)
    final_valid_idx = np.full(nballistics_global, -1, dtype=int)
    first_stationary_idx = np.full(nballistics_global, -1, dtype=int)
    last_moving_idx = np.full(nballistics_global, -1, dtype=int)
    n_valid_outputs = np.zeros(nballistics_global, dtype=int)
    last_seen_time = np.full(nballistics_global, np.nan, dtype=float)
    first_stationary_time = np.full(nballistics_global, np.nan, dtype=float)
    final_plateau_length = np.zeros(nballistics_global, dtype=int)
    last_moving_to_final_distance = np.full(nballistics_global, np.nan)
    status = np.full(nballistics_global, "noValid", dtype=object)

    final_output_idx = n_times - 1

    for s_global_idx in range(nballistics_global):
        pos_series = B_position[s_global_idx, :, :]
        valid_position_mask = np.all(~np.isnan(pos_series), axis=0)
        valid_indices = np.where(valid_position_mask)[0]

        if valid_indices.size == 0:
            continue

        first_idx = valid_indices[0]
        final_idx = valid_indices[-1]
        first_valid_idx[s_global_idx] = first_idx
        final_valid_idx[s_global_idx] = final_idx
        last_seen_idx[s_global_idx] = final_idx
        n_valid_outputs[s_global_idx] = valid_indices.size
        last_seen_time[s_global_idx] = sorted_times[final_idx]

        final_position = pos_series[:, final_idx]

        # Walk backward from the last seen particle position and find the
        # beginning of the final constant-position plateau. For particles still
        # present at the final output this is the deposition time. For particles
        # that disappear earlier it is only the beginning of the last observed
        # stationary segment, if any.
        plateau_start_idx = final_idx
        for idx in valid_indices[-2::-1]:
            dist_from_final = LA.norm(pos_series[:, idx] - final_position)
            if dist_from_final <= stop_distance_tolerance:
                plateau_start_idx = idx
            else:
                last_moving_idx[s_global_idx] = idx
                last_moving_to_final_distance[s_global_idx] = dist_from_final
                break

        first_stationary_idx[s_global_idx] = plateau_start_idx
        first_stationary_time[s_global_idx] = sorted_times[plateau_start_idx]
        final_plateau_length[s_global_idx] = int(
            np.sum(valid_indices >= plateau_start_idx))

        disappeared_before_final_output = final_idx < final_output_idx

        if disappeared_before_final_output:
            # A parcel that is not present in the last CSV output has been
            # removed from the cloud before the final output. In the current
            # wall-interaction setup this is interpreted as escaped/removed,
            # independently of whether it had a short stationary segment before
            # disappearing.
            status[s_global_idx] = "escaped"
        else:
            # The parcel is still present in the final CSV output. If the last
            # observed position differs from the previous available position,
            # the final position is not part of a stationary plateau, so the
            # parcel is still moving at the final time. Otherwise it belongs to
            # the final stationary plateau and is treated as deposited.
            if valid_indices.size == 1:
                status[s_global_idx] = "stillMoving"
            elif final_plateau_length[s_global_idx] <= 1:
                status[s_global_idx] = "stillMoving"
            else:
                status[s_global_idx] = "deposited"

    n_with_last_seen = int(np.sum(last_seen_idx >= 0))
    n_with_motion_before_last_seen = int(np.sum(last_moving_idx >= 0))
    print(f"Particles with a last-seen position: {n_with_last_seen} / {nballistics_global}")
    print("Particles with at least one earlier different position: "
          f"{n_with_motion_before_last_seen} / {nballistics_global}")

    print("Calculating min, mean and max velocities before last-seen/deposit...")
    velocities_summary = np.full((nballistics_global, 5), np.nan)
    for k_global_idx in range(nballistics_global):
        velocities_summary[k_global_idx, 0] = k_global_idx
        velocities_summary[k_global_idx, 1] = d_global[k_global_idx]

        # Use velocities up to the last moving index when there is a final
        # stationary plateau. If the particle disappears while still moving,
        # use all valid velocity samples up to the last seen output.
        if last_moving_idx[k_global_idx] >= 0:
            velocity_end_idx = last_moving_idx[k_global_idx]
        else:
            velocity_end_idx = last_seen_idx[k_global_idx]

        if velocity_end_idx >= 0:
            particle_vel_norms_trajectory = matr_data[:velocity_end_idx + 1, 7,
                                                      k_global_idx]
            valid_vels = particle_vel_norms_trajectory[
                ~np.isnan(particle_vel_norms_trajectory)]
            if len(valid_vels) > 0:
                velocities_summary[k_global_idx, 2] = np.amin(valid_vels)
                velocities_summary[k_global_idx, 3] = np.mean(valid_vels)
                velocities_summary[k_global_idx, 4] = np.amax(valid_vels)

    C_vel_headers = [
        'global_index', 'diameter [m]',
        'min_velocity_before_last_seen_or_deposit [m/s]',
        'mean_velocity_before_last_seen_or_deposit [m/s]',
        'max_velocity_before_last_seen_or_deposit [m/s]'
    ]
    df_velocities = pd.DataFrame(velocities_summary)
    df_velocities[0] = df_velocities[0].astype(int)
    velocities_csv_path = os.path.join(postprocessing_dir, "velocities.csv")
    df_velocities.to_csv(velocities_csv_path,
                         header=C_vel_headers,
                         index=False,
                         na_rep='NaN')
    print(f"Velocities summary saved to: {velocities_csv_path}")

    print("Calculating last-seen particle properties...")
    r_global = d_global / 2.0
    V_global = (4.0 / 3.0) * np.pi * (r_global**3)
    m_global = rho_global * V_global

    impact_rows = []
    for s_global_idx, particle_id in enumerate(unique_global_particle_identifiers):
        orig_proc, orig_id = particle_id
        row = {
            'global_index': int(s_global_idx),
            'origProc': int(orig_proc),
            'origId': int(orig_id),
            'status': status[s_global_idx],
            'diameter [m]': d_global[s_global_idx],
            'density [kg/m3]': rho_global[s_global_idx],
            'first_timestep_index': int(first_valid_idx[s_global_idx]),
            'first_time [s]': np.nan,
            'last_seen_timestep_index': int(last_seen_idx[s_global_idx]),
            'last_seen_time [s]': last_seen_time[s_global_idx],
            'x_last_seen [m]': np.nan,
            'y_last_seen [m]': np.nan,
            'z_last_seen [m]': np.nan,
            'last_seen_velocity_norm [m/s]': np.nan,
            'last_seen_kinetic_energy [J]': np.nan,
            'stationary_start_timestep_index': int(first_stationary_idx[s_global_idx]),
            'stationary_start_time [s]': first_stationary_time[s_global_idx],
            'last_moving_timestep_index': int(last_moving_idx[s_global_idx]),
            'last_moving_time [s]': np.nan,
            'last_moving_velocity_norm [m/s]': np.nan,
            'last_moving_kinetic_energy [J]': np.nan,
            'final_plateau_length_outputs': int(final_plateau_length[s_global_idx]),
            'last_moving_to_last_seen_distance [m]': last_moving_to_final_distance[s_global_idx],
            'n_csv_appearances': int(n_valid_outputs[s_global_idx]),
        }

        if first_valid_idx[s_global_idx] >= 0:
            row['first_time [s]'] = sorted_times[first_valid_idx[s_global_idx]]

        lsi = last_seen_idx[s_global_idx]
        if lsi >= 0:
            row['x_last_seen [m]'] = matr_data[lsi, 1, s_global_idx]
            row['y_last_seen [m]'] = matr_data[lsi, 2, s_global_idx]
            row['z_last_seen [m]'] = matr_data[lsi, 3, s_global_idx]
            last_seen_vel_norm = matr_data[lsi, 7, s_global_idx]
            row['last_seen_velocity_norm [m/s]'] = last_seen_vel_norm
            if not np.isnan(last_seen_vel_norm) and not np.isnan(m_global[s_global_idx]):
                row['last_seen_kinetic_energy [J]'] = 0.5 * m_global[s_global_idx] * (last_seen_vel_norm**2)

        lmi = last_moving_idx[s_global_idx]
        if lmi >= 0:
            row['last_moving_time [s]'] = sorted_times[lmi]
            last_moving_vel_norm = matr_data[lmi, 7, s_global_idx]
            row['last_moving_velocity_norm [m/s]'] = last_moving_vel_norm
            if not np.isnan(last_moving_vel_norm) and not np.isnan(m_global[s_global_idx]):
                row['last_moving_kinetic_energy [J]'] = 0.5 * m_global[s_global_idx] * (last_moving_vel_norm**2)

        impact_rows.append(row)

    df_impacts = pd.DataFrame(impact_rows)
    impacts_csv_path = os.path.join(postprocessing_dir, "impacts.csv")
    df_impacts.to_csv(impacts_csv_path,
                      index=False,
                      na_rep='NaN')
    print(f"Particle last-seen/status properties saved to: {impacts_csv_path}")
    print("Status counts:")
    for status_name, count in df_impacts['status'].value_counts().items():
        print(f"  {status_name}: {count}")

    if args.map_mode == "deposited":
        statuses_for_mapping = ["deposited"]
        mapping_label = "deposited particles"
    else:
        statuses_for_mapping = ["deposited", "escaped"]
        mapping_label = "deposited + escaped particles"

    print("Statuses included in raster maps: "
          + ", ".join(statuses_for_mapping))

    filtered_impact_data = df_impacts[
        df_impacts['status'].isin(statuses_for_mapping)
        & df_impacts['x_last_seen [m]'].notna()
        & df_impacts['y_last_seen [m]'].notna()
    ]

    if filtered_impact_data.shape[0] == 0:
        print(f"No {mapping_label} available for raster maps.")
    else:
        print(f"Number of {mapping_label} for mapping: "
              f"{filtered_impact_data.shape[0]}")
        x_impact_coords = filtered_impact_data['x_last_seen [m]'].to_numpy()
        y_impact_coords = filtered_impact_data['y_last_seen [m]'].to_numpy()
        diam_impacted = filtered_impact_data['diameter [m]'].to_numpy()

        if grid_definition_mode == "full":
            # ... (codice per la modalità "full") ...
            print("Using user-defined grid parameters for raster maps.")
            map_x_ll = args.lx
            map_y_ll = args.ly
            current_map_res = args.res
            map_ncols = args.ncols
            map_nrows = args.nrows
            x_density_grid_centers = map_x_ll + 0.5 * \
                current_map_res + np.arange(map_ncols) * current_map_res
            y_density_grid_centers = map_y_ll + 0.5 * \
                current_map_res + np.arange(map_nrows) * current_map_res
        else:
            data_extent_x_min, data_extent_x_max = dem_extent[0], dem_extent[1]
            data_extent_y_min, data_extent_y_max = dem_extent[2], dem_extent[3]
            if grid_definition_mode == "default":
                print(f"Using script default resolution "
                      f"({DEFAULT_RASTER_RESOLUTION}) and DEM-derived "
                      f"extents for raster maps.")
                current_map_res = DEFAULT_RASTER_RESOLUTION
            else:  # "res_only"
                print(f"Using user-defined resolution ({args.res}) "
                      f"and DEM-derived extents for raster maps.")
                current_map_res = args.res
            x_density_grid_centers = np.arange(
                data_extent_x_min + 0.5 * current_map_res, data_extent_x_max,
                current_map_res)
            y_density_grid_centers = np.arange(
                data_extent_y_min + 0.5 * current_map_res, data_extent_y_max,
                current_map_res)
            if len(x_density_grid_centers) == 0:
                x_density_grid_centers = np.array(
                    [data_extent_x_min + 0.5 * current_map_res])
            if len(y_density_grid_centers) == 0:
                y_density_grid_centers = np.array(
                    [data_extent_y_min + 0.5 * current_map_res])
            map_x_ll, map_y_ll = data_extent_x_min, data_extent_y_min

        xx_grid_centers, yy_grid_centers = np.meshgrid(x_density_grid_centers,
                                                       y_density_grid_centers)
        map_nrows, map_ncols = xx_grid_centers.shape
        map_x_ur = map_x_ll + map_ncols * current_map_res
        map_y_ur = map_y_ll + map_nrows * current_map_res
        density_map_plot_extent = (map_x_ll, map_x_ur, map_y_ll, map_y_ur)

        print('Density map grid definition:')
        print(f'  LL corner: ({map_x_ll:.2f}, {map_y_ll:.2f})')
        print(f'  Cell size: {current_map_res:.2f}')
        print(f'  NCols: {map_ncols}, NRows: {map_nrows}')
        print(f'  UR corner: ({map_x_ur:.2f}, {map_y_ur:.2f})')
        print(f'  Plot extent: {density_map_plot_extent}')

        num_plot_classes = len(plot_titles)
        count_ballistic_class = np.zeros(
            (num_plot_classes, map_nrows, map_ncols))

        for xi, yi, di in zip(x_impact_coords, y_impact_coords, diam_impacted):
            grid_ix_frac, grid_jy_frac = (xi - map_x_ll) / current_map_res, (
                yi - map_y_ll) / current_map_res
            i_col, j_row = int(
                np.clip(np.ceil(grid_ix_frac) - 1, 0, map_ncols - 1)), int(
                    np.clip(np.ceil(grid_jy_frac) - 1, 0, map_nrows - 1))

            count_ballistic_class[-1, j_row, i_col] += 1

            if num_unique_sizes < 5:
                for k_diam_idx in range(len(diams_for_plot_categories)):
                    if np.isclose(di, diams_for_plot_categories[k_diam_idx]):
                        count_ballistic_class[k_diam_idx, j_row, i_col] += 1
                        break
            else:
                assigned_to_specific_class = False
                if di <= diams_for_plot_categories[0]:
                    count_ballistic_class[0, j_row, i_col] += 1
                    assigned_to_specific_class = True
                else:
                    for k_bin_idx in range(len(diams_for_plot_categories) - 1):
                        if di <= diams_for_plot_categories[k_bin_idx + 1]:
                            count_ballistic_class[k_bin_idx + 1, j_row,
                                                  i_col] += 1
                            assigned_to_specific_class = True
                            break
                if not assigned_to_specific_class:
                    count_ballistic_class[len(diams_for_plot_categories),
                                          j_row, i_col] += 1

        # Loop per il plotting
        for i_class in range(num_plot_classes):
            fig, ax = plt.subplots()
            ax.imshow(ls.hillshade(np.flipud(Z_dem_data),
                                   vert_exag=1.0,
                                   dx=dem_cell_size,
                                   dy=dem_cell_size),
                      cmap='gray',
                      extent=dem_extent)
            ax.set_aspect('equal', 'box')

            zz_counts = np.squeeze(count_ballistic_class[i_class, :, :])
            sum_zz_counts = np.sum(zz_counts)

            if sum_zz_counts > 0:
                zz_percentage = zz_counts / sum_zz_counts * 100.0
                zz_log_percentage = np.log10(zz_percentage,
                                             out=np.full_like(
                                                 zz_percentage, -np.inf),
                                             where=(zz_percentage > 0))
            else:
                zz_log_percentage = np.full_like(zz_counts, -np.inf)

            finite_vals = zz_log_percentage[np.isfinite(zz_log_percentage)]
            if finite_vals.size > 0:
                min_val_plot = np.min(finite_vals)
                max_val_plot = np.max(finite_vals)
                if np.isclose(max_val_plot, min_val_plot):
                    max_val_plot = min_val_plot + 0.1
            else:
                min_val_plot = 0  # log10(1%)
                max_val_plot = 1  # log10(10%)

            # plot_levels sono ancora i valori logaritmici per la mappatura
            # dei colori
            plot_levels = np.linspace(min_val_plot, max_val_plot, 11)
            cmap_terrain = plt.get_cmap('terrain_r')
            norm_boundary = BoundaryNorm(plot_levels,
                                         ncolors=cmap_terrain.N,
                                         clip=True)

            im_overlay = ax.imshow(np.flipud(zz_log_percentage),
                                   cmap=cmap_terrain,
                                   norm=norm_boundary,
                                   interpolation='nearest',
                                   extent=density_map_plot_extent,
                                   alpha=0.65)

            ax.set_xlim(density_map_plot_extent[0], density_map_plot_extent[1])
            ax.set_ylim(density_map_plot_extent[2], density_map_plot_extent[3])

            # --- Modifica per la Colorbar ---
            # I ticks sono nelle posizioni dei livelli logaritmici,
            # ma le etichette mostrano la percentuale
            clb = plt.colorbar(im_overlay,
                               ticks=plot_levels,
                               format=FuncFormatter(log_to_percent_formatter))
            clb.set_label('% Ballistics (Color Scale is Logarithmic)',
                          labelpad=10)  # Etichetta aggiornata
            # --- Fine Modifica Colorbar ---

            ax.set_title(plot_titles[i_class])

            safe_title_part = plot_titles[i_class].replace(" ", "_").replace(
                "<=", "lte").replace(">", "gt").replace("=", "eq").replace(
                    ".", "p").replace(",", "")
            png_file = os.path.join(output_raster_root,
                                    f"map_{safe_title_part}_ballistic.png")
            plt.savefig(png_file, dpi=200)
            plt.close(fig)

            nodata_val_asc = -9999.0
            # Salviamo ancora i valori logaritmici nel file ASC
            asc_data_to_save = np.copy(zz_log_percentage)
            asc_data_to_save[np.isneginf(asc_data_to_save)] = nodata_val_asc

            asc_file = os.path.join(output_raster_root,
                                    f"map_{safe_title_part}_ballistic.asc")
            header_asc = f"ncols     {map_ncols}\n"
            header_asc += f"nrows    {map_nrows}\n"
            header_asc += f"xllcorner {map_x_ll:.6f}\n"
            header_asc += f"yllcorner {map_y_ll:.6f}\n"
            header_asc += f"cellsize {current_map_res:.6f}\n"
            header_asc += f"NODATA_value {nodata_val_asc}\n"
            np.savetxt(asc_file,
                       np.flipud(asc_data_to_save),
                       header=header_asc,
                       fmt='%.5f',
                       comments='')

    print("Processing completed.")


# Assicurati che if __name__ == '__main__': sia alla fine
if __name__ == '__main__':
    main()
