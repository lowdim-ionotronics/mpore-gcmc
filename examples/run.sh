#!/bin/sh
# Example voltage-sweep run using a JSON config (see run_config.json).
# Equivalent to passing every parameter as a --flag; see run_gcmc.py --help
# for the full option list. Output (.count/.restart/analysis) goes into
# ./output/ (created automatically), per run_config.json's "output_dir".
cd "$(dirname "$0")"
python3 ../run_gcmc.py -C run_config.json

# Analyze the run's output (charge, capacitance, density, energy vs.
# voltage); see analyze_gcmc.py --help for the full option list.
python3 ../analyze_gcmc.py -C run_config.json
