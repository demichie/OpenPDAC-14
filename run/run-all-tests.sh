#!/bin/bash

# =============================================================================
# CI test runner for OpenPDAC test cases
# -----------------------------------------------------------------------------
# Purpose
#   This script is intended to be used in Continuous Integration (CI) to run
#   OpenPDAC test cases in a controlled and reproducible way.
#
# What it does
#   - Iterates over the selected test-case directories
#   - Skips test cases marked with a .noTest file
#   - Optionally shortens selected runs when executed in CI
#   - Runs the case workflow through Allrun
#   - Optionally validates results through Alltest
#   - Cleans the case with Allclean
#   - Restores any control dictionary modified for CI
#
# CI-specific behavior
#   If the environment variable CI=true and the marker file .shorten_in_ci is
#   present in a test-case directory, the script temporarily shortens the run by
#   modifying:
#     1. system/controlDict.run, if it exists
#     2. otherwise system/controlDict
#
# Notes
#   - The target dictionary is backed up before modification and restored at the
#     end of the test-case execution.
#   - Each test case is run in a subshell so that directory changes do not
#     affect the main script.
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# shorten_run
# -----------------------------------------------------------------------------
# Intelligently shortens a simulation for fast CI execution.
#
# Priority:
#   1. Modify 'system/controlDict.run' if it exists
#   2. Otherwise modify 'system/controlDict' if it exists
#
# Arguments:
#   $1 -> new endTime
#   $2 -> new writeInterval
# -----------------------------------------------------------------------------
shorten_run() {
  local new_end_time=$1
  local new_write_interval=$2
  local target_dict=""

  echo "CI environment detected. Shortening this test case."

  if [ -f "system/controlDict.run" ]; then
    target_dict="system/controlDict.run"
    echo "Found 'system/controlDict.run'. This file will be modified."
  elif [ -f "system/controlDict" ]; then
    target_dict="system/controlDict"
    echo "Found 'system/controlDict'. This file will be modified."
  else
    echo "ERROR: Could not find 'system/controlDict.run' or 'system/controlDict'."
    exit 1
  fi

  echo "--> Updating ${target_dict}: endTime=${new_end_time}, writeInterval=${new_write_interval}"

  # Create a backup before editing
  cp "$target_dict" "${target_dict}.bak"

  # Update endTime and writeInterval
  sed -i "s/^endTime .*/endTime         ${new_end_time};/" "$target_dict"
  sed -i "s/^writeInterval .*/writeInterval   ${new_write_interval};/" "$target_dict"
}

# -----------------------------------------------------------------------------
# Main loop over test-case directories
# -----------------------------------------------------------------------------
for test_case_dir in synthTopo2D/; do
  test_case_dir=${test_case_dir%/}

  if [ -f "${test_case_dir}/.noTest" ]; then
    echo "----------------------------------------------------"
    echo "Skipping test case: ${test_case_dir} (.noTest file found)"
    echo "----------------------------------------------------"
    continue
  fi

  echo "----------------------------------------------------"
  echo "Processing test case: ${test_case_dir}"
  echo "----------------------------------------------------"

  (
    cd "${test_case_dir}"

    # Shorten the run only in CI and only if explicitly requested
    if [ "$CI" = "true" ] && [ -f .shorten_in_ci ]; then
      # Use a short final time and a write interval that still guarantees output
      shorten_run "0.01" "2"
    fi

    # Run the simulation workflow
    if [ -f Allrun ]; then
      echo "Running case workflow (Allrun)..."
      chmod +x ./Allrun
      ./Allrun
    else
      echo "ERROR: Allrun script not found."
      exit 1
    fi

    # Run validation if available
    if [ -f Alltest ]; then
      echo "Validating results (Alltest)..."
      chmod +x ./Alltest
      ./Alltest

      if [ -f Allclean ]; then
        echo "Cleaning case (Allclean)..."
        chmod +x ./Allclean
        ./Allclean
      fi
    else
      echo "WARNING: Alltest script not found for ${test_case_dir}. Validation skipped."
    fi

    # Restore any modified dictionary from backup
    for bak_file in system/*.bak; do
      if [ -f "$bak_file" ]; then
        original_file="${bak_file%.bak}"
        mv "$bak_file" "$original_file"
      fi
    done
  )
done

echo "----------------------------------------------------"
echo "All test cases completed successfully."
echo "----------------------------------------------------"
