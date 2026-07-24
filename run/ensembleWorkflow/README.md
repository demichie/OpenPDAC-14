# Framework for Parametric Study Automation

## 1. Introduction

This folder contains a two-script Python framework designed to automate the setup and execution of parametric studies and simulation ensembles. The workflow is composed of two primary components: **(1) a sophisticated parameter sampling script** that generates a set of input parameters using the Latin Hypercube Sampling (LHS) method, and **(2) an ensemble generation script** that uses these parameters to create a complete directory structure for a batch of simulations. This framework supports a wide variety of statistical distributions for parameter definition, provides rich visual feedback, and includes robust safety features to prevent accidental data loss and aid in debugging.

## 2. Workflow Overview

The end-to-end process follows a clear and logical sequence:

1.  **Configuration**: The user defines the entire parameter space within the sampling script, specifying parameter names, types, and statistical distributions.
2.  **Sampling (Script 1)**: The sampling script is executed. It interprets the configuration and generates a `samples.csv` file, where each row represents a complete, unique set of parameters for one simulation. It also produces a series of plots to visualize the generated sample space.
3.  **Template Creation**: The user prepares a `templatedir` directory containing all necessary template files for a single simulation run. Placeholders are used where sampled values should be inserted.
4.  **Ensemble Generation (Script 2)**: The ensemble generation script is executed. It reads `samples.csv`, and for each row, it creates a new simulation directory (e.g., `ensemble.00000`), copies the contents of `templatedir`, and recursively replaces all placeholders with the corresponding parameter values.

This decouples the statistical sampling from the file system manipulation, making the framework modular and easy to manage.

## 3. Script 1: Parameter Sampling

This script is the heart of the framework, responsible for the statistical definition of the parameter space.

### 3.1. Purpose

The primary goal of this script is to generate a `samples.csv` file containing `N` sets of parameters, sampled efficiently using Latin Hypercube Sampling (LHS). LHS ensures that each parameter's range is evenly stratified, providing excellent coverage of the parameter space with a relatively small number of samples compared to pure Monte Carlo methods.

### 3.2. Configuration: The `parameter_config` List

The entire parameter space is defined in a single Python list named `parameter_config`. Each element of the list is a dictionary that defines one parameter (or a related group of parameters).

#### 3.2.1. Common Keys

-   `name` (str) or `names` (List[str]): The name of the parameter, which will become a column header in `samples.csv`. Use `names` for discrete parameters that unpack into multiple columns (e.g., coordinates).
-   `type` (str): The type of parameter. Must be either `'continuous'` or `'discrete'`.

#### 3.2.2. Continuous Distributions

For `'type': 'continuous'`, the following keys are used:

-   `distribution` (str): The statistical distribution to sample from.
-   `range` (List[float]): A list of two numbers, `[min_val, max_val]`, defining the inclusive bounds of the distribution.

The supported distributions are:

**a) Linear Uniform (`linear`)**
Samples are drawn with equal probability from anywhere within the specified range.

-   **Use Case**: For parameters where any value in the range is equally likely.
-   **Configuration**:
    ```python
    {'name': 'mu', 'type': 'continuous', 'distribution': 'linear', 'range': [0.1, 0.3]}
    ```

**b) Log Uniform (`log`)**
Samples are drawn uniformly from a logarithmic scale. This means that each order of magnitude within the range has an equal probability of being sampled.

-   **Use Case**: Excellent for parameters that span several orders of magnitude (e.g., 10 to 10,000), where the exact order of magnitude is the primary uncertainty.
-   **Configuration**:
    ```python
    {'name': 'Volume', 'type': 'continuous', 'distribution': 'log', 'range': [1000, 1e6]}
    ```

**c) Truncated Normal (`truncnorm`)**
Samples are drawn from a normal (Gaussian) distribution defined by a mean and standard deviation, but are strictly bounded by the `range`.

-   **Use Case**: For parameters that have a known or expected average value but are subject to physical or logical constraints.
-   **Required Keys**:
    -   `mean` (float): The mean (μ) of the underlying normal distribution.
    -   `std_dev` (float): The standard deviation (σ) of the underlying normal distribution.
    -   `range` (List[float]): The hard truncation limits `[min_val, max_val]`.
-   **Configuration**:
    ```python
    {'name': 'Pressure', 'type': 'continuous', 'distribution': 'truncnorm',
     'mean': 50.0, 'std_dev': 15.0, 'range': [30.0, 90.0]}
    ```

**d) Power Law (`powerlaw`)**
Samples are drawn from a distribution where the probability of a value `x` is proportional to `x^k`, where `k` is the specified exponent.

-   **Use Case**: For scale-free phenomena common in nature (e.g., earthquake magnitudes, particle sizes). A negative exponent means small values are vastly more probable than large ones.
-   **Required Keys**:
    -   `exponent` (float): The exponent `k` of the power law.
    -   `range` (List[float]): The inclusive bounds of the distribution. Must be positive.
-   **Configuration**:
    ```python
    {'name': 'GrainSize', 'type': 'continuous', 'distribution': 'powerlaw',
     'exponent': -1.5, 'range': [0.1, 100.0]}
    ```

**e) Truncated Log-Normal (`trunclognorm`)**
Samples are drawn from a log-normal distribution, but are strictly bounded by the `range`. A variable `X` is log-normal if `ln(X)` is normally distributed. This distribution is defined by the mean and standard deviation of the *logarithm* of the variable.

-   **Use Case**: For parameters that are known to follow a log-normal pattern (e.g., permeability, hydraulic conductivity, pollutant concentrations) but are constrained within a physically possible range.
-   **Required Keys**:
    -   `log_mean` (float): The mean (μ) of the *logarithm* of the variable.
    -   `log_std_dev` (float): The standard deviation (σ) of the *logarithm* of the variable.
    -   `range` (List[float]): The hard truncation limits `[min_val, max_val]` for the variable itself (not its logarithm). Must be positive.
-   **Configuration**:
    ```python
    {'name': 'Permeability', 'type': 'continuous', 'distribution': 'trunclognorm',
     'log_mean': -12.0, 'log_std_dev': 1.5, 'range': [1e-7, 1e-4]}
    ```

#### 3.2.3. Discrete Distributions

For `'type': 'discrete'`, the following keys are used:

-   `values` (List): A list of the possible discrete values the parameter can take. These can be numbers, strings, or even lists (for coordinates).
-   `weights` (List[float], optional): A list of relative weights corresponding to each item in `values`. If omitted, all values are assumed to have equal probability.

**a) Simple Discrete Variable**

-   **Configuration**:
    ```python
    {'name': 'SolverType', 'type': 'discrete', 'values': ['explicit', 'implicit'], 'weights':}
    ```
    *(In this example, 'explicit' is three times more likely to be chosen than 'implicit').*

**b) Paired Discrete Variables (e.g., Coordinates)**
This is a special case for handling related discrete values.

-   **Required Keys**:
    -   `names` (List[str]): The names of the output columns (e.g., `['x', 'y']`).
    -   `values` (List[List]): A list of pairs (or tuples).
    -   `unpack` (bool): Must be set to `True`.
-   **Configuration**:
    ```python
    {'names': ['x', 'y'], 'type': 'discrete', 'values': [,],
     'weights':, 'unpack': True, 'plot_as': 'category'}
    ```

### 3.3. Plotting Configuration

The script automatically generates plots to visualize the sample space. You can control how variables are represented using an optional key:

-   `plot_as` (str, optional): If set to `'category'`, the variable will be used for styling (color or marker shape) in the plots instead of being plotted on an axis.

The script intelligently assigns styles: the categorical variable with the most unique values is used for color (`hue`), and the one with the second-most is used for marker shape (`style`).

### 3.4. Output

1.  **`samples.csv`**: A CSV file where each row is a sample and each column is a parameter.
2.  **`samples_plot.png` / `.pdf`**: A pair plot matrix showing the distribution of each parameter (on the diagonal) and the relationships between parameter pairs (on the off-diagonals).
3.  **`counts_of_*.png` / `.pdf`**: For each variable marked as a category, a separate bar chart is generated showing the frequency of each discrete value.

## 4. Script 2: Ensemble Generation

This script uses the output of the sampling script to build the simulation directories.

### 4.1. Purpose

The script automates the creation of an ensemble of simulation directories. For each row in `samples.csv`, it creates a dedicated directory and populates it with the necessary input files, replacing placeholder strings with the specific parameter values for that run.

### 4.2. Prerequisites

1.  **`samples.csv`**: Must be present in the same directory.
2.  **`templatedir`**: A directory containing all the template files and subdirectories needed for a single simulation.

### 4.3. Placeholder Syntax

Within any text file in `templatedir`, you must use a specific syntax for placeholders that will be replaced. The syntax is:
`ENSEMBLE_ParameterName`

Where `ParameterName` must exactly match a column header in `samples.csv`.

**Example `parameters.txt` in `templatedir`**:
```
# Simulation Parameters
Volume = ENSEMBLE_Volume
Pressure = ENSEMBLE_Pressure
```

### 4.4. Operation

The script iterates through each row of `samples.csv`. For row `i`, it:
1.  Creates a directory named `ensemble.XXXXX` (e.g., `ensemble.00000` for row 0).
2.  Copies the entire contents of `templatedir` into this new directory.
3.  **Recursively** scans every file in the new directory and its subdirectories.
4.  In each file, it replaces every instance of a placeholder (e.g., `ENSEMBLE_Volume`) with the corresponding value from row `i`.

### 4.5. Safety and Debugging Features

-   **Data Loss Prevention**: The script will check if a target directory (e.g., `ensemble.00000`) already exists. If it does, the script will **stop with an error message**, preventing the accidental deletion of previous simulation results.
-   **Unused Parameter Warning**: After processing a directory, the script will check if any parameters from the CSV row were not used (i.e., their `ENSEMBLE_` placeholder was not found in any file). If so, it prints a warning to the console, helping to catch typos or forgotten placeholders in the template files.
-   **Detailed Logging**: The script prints a detailed log of every file it modifies and every specific replacement it makes, allowing for easy verification.

## 5. Quick Start Guide

1.  **Configure Script 1**: Open the sampling script (`sampling_script.py`) and carefully define all your parameters in the `parameter_config` list.
2.  **Create `templatedir`**: In the same location as the scripts, create a directory named `templatedir`. Populate it with all necessary input files, using the `ENSEMBLE_ParameterName` syntax for placeholders.
3.  **Make Scripts Executable**: Before the first run, make the utility scripts executable:
    ```bash
    chmod +x Allrun
    chmod +x Allclean
    ```
4.  **Execute the Workflow**: Run the main execution script. This will run both Python scripts in the correct sequence.
    ```bash
    ./Allrun
    ```
5.  **Review Output**: Check the generated `samples.csv` to ensure the values are reasonable. Examine the plot files (`.png`, `.pdf`) to visually confirm the distributions and correlations. Inspect a few `ensemble.XXXXX` directories to confirm that placeholders have been replaced.
6.  You are now ready to run your simulations.
7.  **Cleanup (Optional)**: When you are finished and have backed up your results, run the cleanup script to remove all generated files and directories.
    ```bash
    ./Allclean
    ```

## 6. Workflow Execution Scripts (`Allrun` and `Allclean`)

To streamline the execution of the entire workflow, two Bash utility scripts are provided: `Allrun` and `Allclean`. These scripts automate the process of running the Python scripts and cleaning up the generated files.

### 6.1. The `Allrun` Script

This script serves as the main entry point to execute the entire workflow. It runs the sampling and ensemble generation scripts in the correct sequence.

**Key Features:**
-   **Sequential Execution**: Ensures `sampling_script.py` is run before `ensemble_generator.py`.
-   **Pre-run Checks**: Verifies that both Python scripts and the `templatedir` exist before starting.
-   **Error Handling**: Uses `set -e` to immediately stop the entire process if any step fails, preventing partial or corrupt outputs.

**Usage:**
```bash
./Allrun
```

### 6.2. The `Allclean` Script

This script provides a safe and convenient way to reset your project directory by removing all files and folders generated by the `Allrun` script.

**Key Features:**
-   **Targeted Cleaning**: Specifically removes `ensemble.*` directories, `samples.csv`, and all generated plot files (`.png`, `.pdf`).
-   **User Confirmation**: **Crucially, it first lists all items that will be deleted and prompts the user for confirmation (`y/N`) before proceeding.** This is a critical safety feature to prevent accidental data loss.

**Usage:**
```bash
./Allclean
```
