import pandas as pd
from pathlib import Path
import shutil
import sys


def replace_placeholders_in_directory(working_dir: Path, row_data: pd.Series):
    """Recursively finds files, replaces placeholders, and reports unused parameters.

    This function walks through the directory tree, replacing placeholders
    (e.g., 'ENSEMBLE_Volume') with values from the provided data row. It logs
    all successful replacements. After processing all files, it prints a
    warning if any parameters from the data row were not found in any file.

    Args:
        working_dir (Path):
            The path to the target directory whose files will be processed.
        row_data (pd.Series):
            A Pandas Series representing one row of the samples.csv. The index
            contains parameter names and values are the sample values.

    Returns:
        None. Files are modified in place and results are printed to the console.
    """
    print(f"  Processing directory: {working_dir}")

    # A set to keep track of all parameters that are successfully found and replaced.
    found_parameters = set()

    # First pass: iterate through files to perform replacements
    for path in working_dir.rglob('*'):
        if path.is_file():
            try:
                content = path.read_text(encoding='utf-8')
                replacements_made_in_file = []

                for param_name, value in row_data.items():
                    placeholder = f"ENSEMBLE_{param_name}"

                    if placeholder in content:
                        content = content.replace(placeholder, str(value))
                        replacements_made_in_file.append(
                            f"Replaced '{placeholder}' with '{str(value)}'")
                        # Add the used parameter to our tracking set
                        found_parameters.add(param_name)

                if replacements_made_in_file:
                    print(f"    -> In file '{path.relative_to(working_dir)}':")
                    for replacement_log in replacements_made_in_file:
                        print(f"         - {replacement_log}")
                    path.write_text(content, encoding='utf-8')

            except UnicodeDecodeError:
                # Silently skip binary files that cannot be decoded as text.
                continue
            except Exception as e:
                print(
                    f"    -> An unexpected error occurred with file {path}: {e}"
                )

    # Second pass (after all files are processed) to check for unused parameters.
    all_parameters = set(row_data.index)
    unused_parameters = all_parameters - found_parameters

    if unused_parameters:
        print(
            f"  Warning: The following parameters were defined in the CSV but not found in any template files:"
        )
        # Sort for consistent output order
        for param in sorted(list(unused_parameters)):
            print(f"    - {param}")


def main():
    """Main script to generate an ensemble of simulation directories.

    This script performs the following steps:
    1. Reads the `samples.csv` file.
    2. For each row (each set of parameters) in the CSV:
        a. Creates a new uniquely named directory (e.g., 'ensemble.00000').
        b. Checks for existence and stops to prevent accidental data loss.
        c. Copies the entire contents of a 'templatedir' into the new directory.
        d. Scans all files within the new directory (including subdirectories),
           replaces placeholder strings, and reports unused parameters.
    """
    base_dir = Path.cwd()
    samples_file = base_dir / 'CSV/samples.csv'
    template_dir = base_dir / 'templatedir'

    # --- Pre-run Checks ---
    if not samples_file.is_file():
        print(
            f"Error: Input file '{samples_file}' not found. Please run the sampling script first."
        )
        sys.exit(1)

    if not template_dir.is_dir():
        print(
            f"Error: Template directory '{template_dir}' not found. Please create it."
        )
        sys.exit(1)

    # --- Load Data ---
    df = pd.read_csv(samples_file)
    print(f"Loaded {len(df)} samples from '{samples_file}'.\n")

    # --- Check for 'sample_id' column ---
    if 'sample_id' not in df.columns:
        print(f"Error: The required 'sample_id' column was not found in '{samples_file}'.")
        print("Please re-run the sampling script to include this identifier column.")
        sys.exit(1)

    # --- Main Loop ---
    for index, row in df.iterrows():
        working_dir_name = f"ensemble.{index:05d}"
        working_dir = base_dir / working_dir_name

        print(f"Preparing ensemble member {index}: {working_dir_name}")

        # --- Safe Directory Handling ---
        # Check if the target directory already exists. If so, stop the script
        # with an error message to prevent overwriting existing data.
        if working_dir.exists():
            print(f"\nError: Directory '{working_dir_name}' already exists.")
            print(
                "Please remove or back up the existing ensemble directories before running again."
            )
            sys.exit(1)  # Exit the script with an error code

        # If the directory does not exist, proceed with copying.
        shutil.copytree(template_dir, working_dir, symlinks=True)

        params_to_replace = row.drop('sample_id')

        # --- Placeholder Replacement ---
        replace_placeholders_in_directory(working_dir, params_to_replace)
        print("-" * 20)

    print("\nScript finished successfully.")


if __name__ == '__main__':
    main()
