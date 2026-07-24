import pandas as pd
import numpy as np
import json
from pathlib import Path
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from SALib.analyze import sobol
from typing import Tuple, Dict, List
import os


def load_data(
        samples_file: Path, results_file: Path,
        config_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame, List[Dict]]:
    """Loads the input samples, output results, and parameter configuration.

    Args:
        samples_file (Path): Path to the input parameters CSV file (e.g., 'samples.csv').
        results_file (Path): Path to the simulation results CSV file (e.g., 'results.csv').
        config_file (Path): Path to the JSON file defining the parameter space.

    Returns:
        A tuple containing the input DataFrame (X), the output DataFrame (Y), and the problem configuration.
    """
    print("--- Loading Data ---")
    if not samples_file.is_file() or not results_file.is_file(
    ) or not config_file.is_file():
        print(f"Error: One or more required files not found.")
        print(f"  - Input samples expected at: {samples_file}")
        print(f"  - Simulation results expected at: {results_file}")
        print(f"  - Parameter configuration expected at: {config_file}")
        sys.exit(1)

    X = pd.read_csv(samples_file)    
    if 'sample_id' in X.columns:
        # Usiamo .drop() specificando il nome della colonna. 'axis=1' indica che stiamo
        # eliminando una colonna. 'inplace=True' modifica il DataFrame direttamente.
        X.drop('sample_id', axis=1, inplace=True)
        print("  - 'sample_id' column dropped from input samples.")
   
    
    Y = pd.read_csv(results_file)
    if 'ensemble_id' in Y.columns:
        Y.drop('ensemble_id', axis=1, inplace=True)
        print("  - 'ensemble_id' column dropped from results.")    
    
    with config_file.open('r') as f:
        param_config = json.load(f)

    if len(X) != len(Y):
        print(
            f"Error: Row count mismatch between '{samples_file}' ({len(X)}) and '{results_file}' ({len(Y)})."
        )
        sys.exit(1)

    print(f"Successfully loaded {len(X)} samples and results.")
    return X, Y, param_config


def plot_input_output_correlations(X: pd.DataFrame, Y: pd.DataFrame):
    """Generates scatter plots for each output vs. each input parameter.

    It assumes that any necessary log transformations have already been applied to the DataFrames
    and their column names have been updated accordingly.

    Args:
        X (pd.DataFrame): The (potentially transformed) input parameter samples.
        Y (pd.DataFrame): The (potentially transformed) simulation results.
    """
    print("\n--- Generating Input vs. Output Correlation Plots ---")
    numeric_inputs = X.select_dtypes(include=np.number).columns.tolist()

    for output_name in Y.columns:
        n_inputs = len(numeric_inputs)
        n_cols = 3
        n_rows = (n_inputs + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows,
                                 n_cols,
                                 figsize=(n_cols * 5, n_rows * 4),
                                 constrained_layout=True)
        axes = axes.flatten()

        for i, input_name in enumerate(numeric_inputs):
            ax = axes[i]
            # The regression is now calculated correctly on the (potentially log-transformed) data.
            sns.regplot(x=X[input_name],
                        y=Y[output_name],
                        ax=ax,
                        scatter_kws={'alpha': 0.6},
                        line_kws={'color': 'red'})
            ax.set_title(f"{output_name} vs. {input_name}")
            ax.set_xlabel(input_name)
            ax.set_ylabel(output_name)

        # Hide any unused subplots
        for i in range(n_inputs, len(axes)):
            axes[i].set_visible(False)

        fig.suptitle(f'Correlation Plots for Output: {output_name}',
                     fontsize=16)
        # Sanitize filename for robustness
        safe_output_name = "".join(c for c in output_name
                                   if c.isalnum() or c in (' ', '_')).rstrip()
        filename = f"ENSEMBLE_PLOTS/correlations_for_{safe_output_name.replace(' ', '_')}.png"
        plt.savefig(filename, dpi=150)
        filename = f"ENSEMBLE_PLOTS/correlations_for_{safe_output_name.replace(' ', '_')}.pdf"
        plt.savefig(filename)
        print(f"Saved correlation plot: {filename}")
        plt.close()


def perform_sobol_analysis(X: pd.DataFrame, Y: pd.DataFrame,
                           param_config: List[Dict],
                           log_scale_inputs: set) -> Tuple[Dict, Dict]:
    """Performs Sobol sensitivity analysis, transforming log-distributed variables.

    Args:
        X (pd.DataFrame): The raw, untransformed input parameter samples.
        Y (pd.DataFrame): The (potentially transformed) simulation results.
        param_config (List[Dict]): The parameter space configuration.
        log_scale_inputs (set): A set of input variable names that require log transformation.

    Returns:
        A tuple containing the analysis results and the problem definition for SALib.
    """
    print("\n--- Performing Sobol Sensitivity Analysis ---")

    problem_params = [
        p for p in param_config if p['type'] == 'continuous' and 'range' in p
    ]
    X_sobol_transformed = X.copy()
    problem_names, problem_bounds = [], []

    print("Transforming log-distributed inputs for Sobol analysis...")
    for p in problem_params:
        param_name = p['name']
        if param_name in log_scale_inputs:
            X_sobol_transformed[param_name] = np.log10(
                X_sobol_transformed[param_name])
            problem_bounds.append(np.log10(p['range']).tolist())
            problem_names.append(f"{param_name} (log10)")
            print(
                f"  - '{param_name}' was converted to log10 scale for Sobol.")
        else:
            problem_bounds.append(p['range'])
            problem_names.append(param_name)

    problem = {
        'num_vars': len(problem_params),
        'names': problem_names,
        'bounds': problem_bounds
    }

    original_names = [p['name'] for p in problem_params]
    X_sobol = X_sobol_transformed[original_names].to_numpy()

    analysis_results = {}

    for output_name in Y.columns:
        print(f"\nAnalysis for output: '{output_name}'")
        Y_sobol = Y[output_name].to_numpy()
        Si = sobol.analyze(problem,
                           Y_sobol,
                           print_to_console=False,
                           calc_second_order=True)
        analysis_results[output_name] = Si
        df_si = pd.DataFrame(
            {
                'S1': Si['S1'],
                'S1_conf': Si['S1_conf'],
                'ST': Si['ST'],
                'ST_conf': Si['ST_conf']
            },
            index=problem['names'])
        print(df_si.round(3))

    return analysis_results, problem


def plot_sobol_indices(analysis_results: Dict, problem: Dict):
    """Plots the calculated Sobol indices as bar charts with confidence intervals."""
    print("\n--- Generating Sobol Indices Plots ---")
    param_names = problem['names']
    for output_name, Si in analysis_results.items():
        df_sobol = pd.DataFrame({
            'S1': Si['S1'],
            'ST': Si['ST']
        },
                                index=param_names)
        errors = pd.DataFrame({
            'S1': Si['S1_conf'],
            'ST': Si['ST_conf']
        },
                              index=param_names)
        df_sobol.sort_values(by='ST', ascending=False, inplace=True)
        errors = errors.reindex(df_sobol.index)
        ax = df_sobol.plot(kind='bar',
                           y=['S1', 'ST'],
                           figsize=(12, 7),
                           yerr=errors,
                           capsize=4,
                           color=sns.color_palette("viridis", 2))
        ax.set_title(f'Sobol Indices for Output: {output_name}')
        ax.set_ylabel('Sensitivity Index')
        ax.set_xlabel('Parameter')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        ax.axhline(0, color='black', linewidth=0.8)
        plt.tight_layout()
        safe_output_name = "".join(c for c in output_name
                                   if c.isalnum() or c in (' ', '_')).rstrip()
        filename = f"ENSEMBLE_PLOTS/sobol_indices_for_{safe_output_name.replace(' ', '_')}.png"
        plt.savefig(filename, dpi=150)
        filename = f"ENSEMBLE_PLOTS/sobol_indices_for_{safe_output_name.replace(' ', '_')}.pdf"
        plt.savefig(filename)
        print(
            f"Saved Sobol indices plot with confidence intervals: {filename}")
        plt.close()


def main():
    """Main script for sensitivity analysis post-processing."""

    csv_folder = 'CSV'

    samples_file = Path(os.path.join(csv_folder, 'samples.csv'))
    results_file = Path(os.path.join(csv_folder, 'results.csv'))
    config_file = Path(os.path.join(csv_folder, 'parameters.json'))

    LOG_TRANSFORM_THRESHOLD = 1000

    # 1. Load raw data
    X_raw, Y_raw, param_config = load_data(samples_file, results_file,
                                           config_file)

    # 2. Create transformed copies of the data for both plotting and analysis
    X_plot_transformed = X_raw.copy()
    Y_plot_transformed = Y_raw.copy()

    # 3. Identify and transform input variables for plotting
    log_scale_inputs = {
        p['name']
        for p in param_config
        if p.get('distribution') in ['log', 'powerlaw', 'trunclognorm']
    }
    print("\n--- Transforming Input Variables for Plotting ---")
    for col_name in X_raw.columns:
        if col_name in log_scale_inputs:
            new_name = f"{col_name} (log10)"
            X_plot_transformed[new_name] = np.log10(
                X_plot_transformed[col_name])
            X_plot_transformed.drop(columns=[col_name], inplace=True)
            print(f"  - Input '{col_name}' transformed to '{new_name}'.")

    # 4. Identify and transform output variables for both plotting and analysis
    print("\n--- Transforming Output Variables for Plotting and Analysis ---")
    for col_name in Y_raw.columns:
        min_val, max_val = Y_raw[col_name].min(), Y_raw[col_name].max()
        if min_val > 0 and (max_val / min_val) > LOG_TRANSFORM_THRESHOLD:
            new_name = f"{col_name} (log10)"
            print(
                f"  - Output '{col_name}' spans > {LOG_TRANSFORM_THRESHOLD}x range. Transforming to '{new_name}'."
            )
            Y_plot_transformed[new_name] = np.log10(
                Y_plot_transformed[col_name])
            Y_plot_transformed.drop(columns=[col_name], inplace=True)
        else:
            if min_val <= 0:
                print(
                    f"  - Output '{col_name}' contains non-positive values. Keeping linear scale."
                )
            else:
                print(
                    f"  - Output '{col_name}' has a small dynamic range. Keeping linear scale."
                )

    # 5. Perform analyses
    # For plotting, use the fully transformed dataframes
    plot_input_output_correlations(X_plot_transformed, Y_plot_transformed)

    # For Sobol, use the raw X (it handles its own transformation) and the transformed Y
    sobol_results, problem_def = perform_sobol_analysis(
        X_raw, Y_plot_transformed, param_config, log_scale_inputs)
    plot_sobol_indices(sobol_results, problem_def)

    print("\nSensitivity analysis finished successfully.")


if __name__ == '__main__':
    main()
