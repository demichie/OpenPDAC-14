import pandas as pd
from pathlib import Path
import sys
import re
from typing import List, Optional

# --- Configuration ---
# The pattern to match the ensemble directories.
# It looks for 'ensemble.' followed by any number of digits.
ENSEMBLE_DIR_PATTERN = r'ensemble\.(\d+)'

# The name of the output file expected inside each ensemble directory.
SOURCE_FILENAME = 'output.csv'

# The destination folder and filename for the aggregated results.
OUTPUT_FOLDER = 'CSV'
OUTPUT_FILENAME = 'results.csv'
# ---------------------


def extract_ensemble_id(directory_name: str) -> Optional[int]:
    """Extracts the numerical ID from an ensemble directory name.

    Uses a regular expression to find the digits following 'ensemble.'.

    Args:
        directory_name (str): The name of the directory (e.g., 'ensemble.001').

    Returns:
        Optional[int]: The extracted integer ID if the name matches the pattern,
                       otherwise None.
    """
    match = re.match(ENSEMBLE_DIR_PATTERN, directory_name)
    if match:
        # Return the first captured group (the digits) as an integer
        return int(match.group(1))
    return None


def read_ensemble_output(file_path: Path) -> Optional[pd.DataFrame]:
    """Reads a CSV file from an ensemble directory into a pandas DataFrame.

    This function attempts to read the specified file. If the file does not
    exist or is empty, it prints a warning and returns None.

    Args:
        file_path (Path): The full path to the CSV file to be read.

    Returns:
        Optional[pd.DataFrame]: The loaded data as a DataFrame, or None if
                                the file could not be read or was empty.
    """
    if not file_path.is_file():
        print(f"    Warning: File not found at '{file_path}'. Skipping.")
        return None

    try:
        # Read the CSV. Assuming standard format with a header.
        df = pd.read_csv(file_path)
        if df.empty:
            print(f"    Warning: File '{file_path}' is empty. Skipping.")
            return None
        return df
    except Exception as e:
        print(f"    Error reading '{file_path}': {e}")
        return None


def main():
    """
    Main script to aggregate results from multiple ensemble directories.

    Steps:
    1. Scans the current directory for folders matching 'ensemble.xxxxx'.
    2. For each folder, reads the 'output.csv' file.
    3. Adds an 'ensemble_id' column to the data based on the folder name.
    4. Concatenates all data into a single DataFrame.
    5. Saves the aggregated data to 'CSV/results.csv'.
    """
    base_dir = Path.cwd()
    print(f"Starting aggregation process in: {base_dir}\n")

    # 1. Find and sort all matching ensemble directories
    # We find all directories and filter them using our regex pattern.
    all_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
    ensemble_dirs = []
    for d in all_dirs:
        eid = extract_ensemble_id(d.name)
        if eid is not None:
            # Store tuple of (ID, Path) for easy sorting
            ensemble_dirs.append((eid, d))

    # Sort by the numerical ID to maintain logical order in the final file
    ensemble_dirs.sort(key=lambda x: x[0])

    if not ensemble_dirs:
        print(f"Error: No directories matching '{ENSEMBLE_DIR_PATTERN}' found.")
        sys.exit(1)

    print(f"Found {len(ensemble_dirs)} ensemble directories to process.")

    data_frames: List[pd.DataFrame] = []
    successful_reads = 0

    # 2. Process each directory
    print("\n--- Processing Directories ---")
    for ensemble_id, directory in ensemble_dirs:
        print(f"  Processing ID {ensemble_id} ({directory.name})...")
        file_path = directory / SOURCE_FILENAME

        # Read the data
        df = read_ensemble_output(file_path)

        if df is not None:
            # 3. Add the ensemble_id as the first column
            # We insert it at position 0 with the column name 'ensemble_id'
            df.insert(0, 'ensemble_id', ensemble_id)
            data_frames.append(df)
            successful_reads += 1

    # 4. Aggregate and Save
    print("\n--- Aggregation and Saving ---")
    if not data_frames:
        print("Error: No valid data could be read from any directory.")
        print(f"Could not create '{OUTPUT_FILENAME}'.")
        sys.exit(1)

    # Concatenate all individual DataFrames into one
    # ignore_index=True creates a fresh index for the master DataFrame
    aggregated_df = pd.concat(data_frames, ignore_index=True)

    # Ensure the output directory exists
    output_dir_path = base_dir / OUTPUT_FOLDER
    output_dir_path.mkdir(exist_ok=True)
    output_file_path = output_dir_path / OUTPUT_FILENAME

    # Save to CSV without the pandas index
    aggregated_df.to_csv(output_file_path, index=False)

    print(f"Successfully processed {successful_reads} out of {len(ensemble_dirs)} directories.")
    print(f"Aggregated data saved to: '{output_file_path}'")
    print(f"Total records: {len(aggregated_df)}")
    print("Script finished.")


if __name__ == '__main__':
    main()
