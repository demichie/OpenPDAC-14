import pandas as pd
import numpy as np
from pathlib import Path
import json
from typing import List, Tuple, Dict, Any
from scipy.stats import truncnorm
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
from SALib.sample import saltelli

# Global or configuration constants
CSV_FOLDER = 'CSV'
OUTPUT_CSV_FILE = 'samples.csv'
OUTPUT_PLOT_FOLDER = 'ENSEMBLE_PLOTS'
OUTPUT_PLOT_BASE_FILENAME = 'samples_plot'
CONFIG_FILE = 'parameters.json'
GENERATE_PLOT = True

# N_SAMPLES is now the *base* number for the Saltelli sampler.
# The total number of samples will be N * (2*D + 2), where D is the number of continuous params.
# Powers of 2 are often recommended for N.
N_BASE_SAMPLES = 32
    

def read_asc_file(asc_file_path: Path) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Reads an ESRI ASCII grid file (.asc) and returns its data and metadata."""
    print(f"Reading .asc file: {asc_file_path}")
    if not asc_file_path.is_file():
        raise FileNotFoundError(f"File {asc_file_path} was not found.")
    header_lines = [getline(str(asc_file_path), i) for i in range(1, 7)]
    header = {
        key.lower(): float(value)
        for key, value in (line.split() for line in header_lines)
    }
    header['ncols'], header['nrows'] = int(header['ncols']), int(
        header['nrows'])
    data = pd.read_table(asc_file_path,
                         delim_whitespace=True,
                         header=None,
                         skiprows=6,
                         dtype=np.float64).values
    return np.flipud(data), header


def plot_pair_grid(df: pd.DataFrame, axis_vars: List, category_vars: List,
                   log_scale_vars: set, axis_limits: Dict, config: Dict):
    """Generates the main pair plot with unified histograms and a correct, manual legend.

    Args:
        df (pd.DataFrame): The DataFrame containing the sample data.
        axis_vars (List[str]): A list of column names to be used as plot axes.
        category_vars (List[str]): A list of column names to be used for styling (hue/style).
        log_scale_vars (set): A set of column names that require a logarithmic scale.
        config (Dict): A dictionary containing plot configuration, like 'base_filename'.

    Returns:
        None. The plot is saved to both PNG and PDF files.
    """
    plot_df = df.copy()
    hue_var, style_var, ordered_cat_vars = None, None, []
    if category_vars:
        cat_uniqueness = sorted([(var, plot_df[var].nunique())
                                 for var in category_vars],
                                key=lambda x: x[1],
                                reverse=True)
        hue_var = cat_uniqueness[0][0]
        ordered_cat_vars.append(hue_var)
        if len(cat_uniqueness) > 1:
            style_var = cat_uniqueness[1][0]
            ordered_cat_vars.append(style_var)

    print("\n--- Generating Main Pair Plot ---")
    print(f"Plotting axes: {axis_vars}")
    if hue_var:
        print(
            f"Variable for color (hue): '{hue_var}' ({plot_df[hue_var].nunique()} unique values)"
        )
    if style_var:
        print(
            f"Variable for marker (style): '{style_var}' ({plot_df[style_var].nunique()} unique values)"
        )

    hue_map, marker_map = None, None
    if hue_var:
        hue_categories = sorted(plot_df[hue_var].unique())
        palette = sns.color_palette("viridis", n_colors=len(hue_categories))
        hue_map = {cat: color for cat, color in zip(hue_categories, palette)}
    if style_var:
        style_categories = sorted(plot_df[style_var].unique())
        markers = ['o', 'X', 's', '^', 'P', 'D'][:len(style_categories)]
        marker_map = {
            cat: marker
            for cat, marker in zip(style_categories, markers)
        }

    n_vars = len(axis_vars)
    fig, axes = plt.subplots(n_vars,
                             n_vars,
                             figsize=(n_vars * 3.7, n_vars * 2.5))
    for i, j in np.ndindex(axes.shape):
        ax = axes[i, j]
        x_var, y_var = axis_vars[j], axis_vars[i]
        if i == j:
            sns.histplot(data=plot_df,
                         x=x_var,
                         ax=ax,
                         legend=False,
                         log_scale=(x_var in log_scale_vars),
                         color="#555555")
        else:
            sns.scatterplot(data=plot_df,
                            x=x_var,
                            y=y_var,
                            hue=hue_var,
                            style=style_var,
                            ax=ax,
                            legend=False,
                            alpha=0.8,
                            s=40,
                            palette=hue_map,
                            markers=marker_map)

        if x_var in axis_limits:
            xmin, xmax = axis_limits[x_var]
            ax.set_xlim(xmin, xmax)
            if i == n_vars - 1:
                ax.set_xticks([xmin, xmax])
        if y_var in axis_limits:
            ymin, ymax = axis_limits[y_var]
            if i != j:
                ax.set_ylim(ymin, ymax)
            if j == 0 and i != j:
                ax.set_yticks([ymin, ymax])

        if i != j and x_var in log_scale_vars:
            ax.set_xscale('log')
        if i != j and y_var in log_scale_vars:
            ax.set_yscale('log')
        if j == 0:
            ax.set_ylabel(y_var)
        else:
            ax.set_ylabel('')
        if i == n_vars - 1:
            ax.set_xlabel(x_var)
        else:
            ax.set_xlabel('')

    legend_elements, legend_title = [], ", ".join(ordered_cat_vars)
    if hue_var and style_var:
        for style_cat, marker in marker_map.items():
            for hue_cat, color in hue_map.items():
                label = f"{str(hue_cat)}, {str(style_cat)}"
                legend_elements.append(
                    Line2D([0], [0],
                           marker=marker,
                           color='w',
                           label=label,
                           markerfacecolor=color,
                           markersize=8))
    elif hue_var:
        for cat, color in hue_map.items():
            legend_elements.append(
                Line2D([0], [0],
                       marker='o',
                       color='w',
                       label=cat,
                       markerfacecolor=color,
                       markersize=8))

    fig.subplots_adjust(right=0.85)
    fig.legend(handles=legend_elements,
               loc='center left',
               bbox_to_anchor=(0.88, 0.5),
               title=legend_title)
    fig.suptitle('Pair Plot of Sampled Parameters', fontsize=16)

    base_folder = config['base_folder']
    base_filename = config['base_filename']
    os.makedirs(base_folder, exist_ok=True)
    base_path = Path(os.path.join(base_folder, base_filename))

    # base_path = Path(config['base_filename'])
    png_path = base_path.with_suffix('.png')
    pdf_path = base_path.with_suffix('.pdf')
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(
        f"Main pair plot successfully saved to '{png_path}' and '{pdf_path}'")
    plt.close()


def plot_category_counts(df: pd.DataFrame, category_vars: List, config: Dict):
    """Generates a separate frequency bar chart for each categorical variable."""
    if not category_vars:
        return
    print("\n--- Generating Category Frequency Plots ---")
    for cat_var in category_vars:
        plt.figure(figsize=(8, 5))
        ax = sns.countplot(data=df,
                           x=cat_var,
                           palette="viridis",
                           order=sorted(df[cat_var].unique()))
        ax.set_title(f"Frequency of Categories for '{cat_var}'")
        ax.set_ylabel("Count")
        ax.set_xlabel("Category")
        plt.xticks(rotation=15, ha='right')
        plt.tight_layout()
        safe_cat_var = "".join(c for c in cat_var
                               if c.isalnum() or c in (' ', '_')).rstrip()

        base_folder = config['base_folder']
        base_filename = f"counts_of_{safe_cat_var.replace(' ', '_')}"
        os.makedirs(base_folder, exist_ok=True)
        base_path = Path(os.path.join(base_folder, base_filename))

        # base_filename = f"counts_of_{safe_cat_var.replace(' ', '_')}"
        png_path = Path(f"{base_path}.png")
        pdf_path = Path(f"{base_path}.pdf")
        plt.savefig(png_path, dpi=150)
        plt.savefig(pdf_path)
        print(f"Category count plot saved to '{png_path}' and '{pdf_path}'")
        plt.close()


def generate_samples_for_sobol(parameter_config: List[Dict],
                               n_base_samples: int) -> Dict[str, Any]:
    """Generates a sample set using the Saltelli method, required for Sobol analysis.

    This function defines the 'problem' for SALib, generates the specific sample
    matrix using saltelli.sample, and then transforms the uniform [0,1] samples
    to their respective target distributions using the inverse transform method.

    Args:
        parameter_config (List[Dict]): The parameter space definition.
        n_base_samples (int): The base number of samples (N). The total number
            of simulations will be N * (2 * D + 2) for second-order indices,
            where D is the number of continuous parameters.

    Returns:
        Dict[str, Any]: A dictionary of the generated samples.
    """
    # 1. Construct the 'problem' dictionary for SALib's continuous variables
    continuous_params = [
        p for p in parameter_config if p['type'] == 'continuous'
    ]
    problem = {
        'num_vars': len(continuous_params),
        'names': [p['name'] for p in continuous_params],
        # Bounds are temporarily [0,1] as we will transform the samples later
        'bounds': [[0, 1]] * len(continuous_params)
    }

    # 2. Generate the Saltelli sample matrix in the unit hypercube [0,1]
    # Using calc_second_order=True is standard for a comprehensive analysis.
    # SALib will generate N * (2D + 2) samples.
    unit_samples = saltelli.sample(problem,
                                   n_base_samples,
                                   calc_second_order=True)

    total_samples_required = len(unit_samples)
    print(
        f"Saltelli sampling for {problem['num_vars']} variables with N={n_base_samples} requires a total of {total_samples_required} simulation runs."
    )

    samples = {}
    param_map = {p['name']: p for p in continuous_params}

    # 3. Transform each column of the unit_samples matrix to its target distribution
    for i, name in enumerate(problem['names']):
        param = param_map[name]
        distribution = param.get('distribution', 'linear')
        y = unit_samples[:, i]  # The uniform [0, 1] samples for this parameter

        if distribution == 'linear':
            min_val, max_val = param['range']
            samples[name] = min_val + y * (max_val - min_val)
        elif distribution == 'log':
            min_val, max_val = param['range']
            if min_val == max_val:
                samples[name] = np.full(total_samples_required, min_val)
            else:
                samples[name] = np.exp(
                    np.log(min_val) + y * (np.log(max_val) - np.log(min_val)))
        elif distribution == 'truncnorm':
            mean, std_dev, (
                min_val,
                max_val) = param['mean'], param['std_dev'], param['range']
            a, b = (min_val - mean) / std_dev, (max_val - mean) / std_dev
            dist = truncnorm(a, b, loc=mean, scale=std_dev)
            samples[name] = dist.ppf(y)
        elif distribution == 'powerlaw':
            exponent, (min_val, max_val) = param['exponent'], param['range']
            k = exponent + 1.0
            samples[name] = (y * (max_val**k - min_val**k) +
                             min_val**k)**(1.0 / k)
        elif distribution == 'trunclognorm':
            log_mean, log_std_dev, (
                min_val, max_val
            ) = param['log_mean'], param['log_std_dev'], param['range']
            log_min, log_max = np.log(min_val), np.log(max_val)
            a, b = (log_min - log_mean) / \
                log_std_dev, (log_max - log_mean) / log_std_dev
            dist_log = truncnorm(a, b, loc=log_mean, scale=log_std_dev)
            log_samples = dist_log.ppf(y)
            samples[name] = np.exp(log_samples)

    # 4. Handle discrete parameters: these are added on top of the continuous sample
    #    by simple random sampling, as they are not part of the Sobol analysis itself.
    discrete_params = [p for p in parameter_config if p['type'] == 'discrete']
    for param in discrete_params:
        values = param['values']
        weights = param.get('weights', None)

        # Create a probability distribution if weights are provided
        probabilities = np.array(weights) / sum(weights) if weights else None

        # Generate random choices for the entire sample size
        sampled_values = np.random.choice(len(values),
                                          size=total_samples_required,
                                          p=probabilities)

        if param.get('unpack', False):
            for i, col_name in enumerate(param['names']):
                samples[col_name] = [values[idx][i] for idx in sampled_values]
        else:
            samples[param['name']] = [values[idx] for idx in sampled_values]

    return samples


##########################################
# MAIN FUNCTION
##########################################


def main():
    """Main script to generate, save, and visualize a parameter set for simulations."""
    # N_SAMPLES is now the *base* number for the Saltelli sampler.
    # The total number of samples will be N * (2*D + 2), where D is the number of continuous params.
    # Powers of 2 are often recommended for N.

    try:
        with open(CONFIG_FILE, 'r') as f:
            parameter_config = json.load(f)
        print(f"Parameter configuration successfully loaded from '{CONFIG_FILE}'")
    except FileNotFoundError:
        print(f"ERROR: Configuration file '{CONFIG_FILE}' not found.")
        return  # Esce dallo script se il file non esiste
    except json.JSONDecodeError:
        print(f"ERROR: Could not decode JSON from '{CONFIG_FILE}'. Check for syntax errors.")
        return

    # Save the parameter configuration to a JSON file for the analysis script
    os.makedirs(CSV_FOLDER, exist_ok=True)

    config_path = Path(os.path.join(CSV_FOLDER, 'parameters.json'))
    with config_path.open('w') as f:
        json.dump(parameter_config, f, indent=4)
    print(f"Parameter configuration saved to '{config_path}'")

    # Use the new Saltelli-based sampling function
    samples = generate_samples_for_sobol(parameter_config, N_BASE_SAMPLES)
    df = pd.DataFrame(samples)
    df.insert(0, 'sample_id', df.index)

    output_file = os.path.join(CSV_FOLDER, OUTPUT_CSV_FILE)

    df.to_csv(output_file, index=False)
    print(
        f"\nSuccessfully generated '{OUTPUT_CSV_FILE}' with {len(df)} samples."
    )

    if GENERATE_PLOT:
        # Prepare data for plotting (this part remains the same)
        plot_df = df.copy()
        axis_vars, category_vars = [], []
        log_scale_vars = set()
        axis_limits = {}
        for param in parameter_config:
            if param.get('plot_as') == 'category':
                if param.get('unpack'):
                    source_name = f"({param['names'][0]},{param['names'][1]})"
                    plot_df[source_name] = plot_df.apply(
                        lambda r:
                        f"P({r[param['names'][0]]},{r[param['names'][1]]})",
                        axis=1)
                    category_vars.append(source_name)
                else:
                    category_vars.append(param['name'])
            else:
                if param.get('unpack'):
                    axis_vars.extend(param['names'])
                else:
                    axis_vars.append(param['name'])
                    if 'range' in param:
                        axis_limits[param['name']] = param['range']
            if param.get('distribution') in [
                    'log', 'powerlaw', 'trunclognorm'
            ]:
                log_scale_vars.add(param['name'])

        plot_config = {
            'base_filename': OUTPUT_PLOT_BASE_FILENAME,
            'base_folder': OUTPUT_PLOT_FOLDER
        }
        plot_pair_grid(plot_df, axis_vars, category_vars, log_scale_vars,
                       axis_limits, plot_config)
        plot_category_counts(plot_df, category_vars, plot_config)

    print("\nScript finished.")


if __name__ == '__main__':
    main()
