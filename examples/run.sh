#!/bin/sh
# Example single-voltage run using a JSON config (see run_config.json).
# Equivalent to passing every parameter as a --flag; see run_gcmc.py --help
# for the full option list.
cd "$(dirname "$0")"
python3 ../run_gcmc.py -C run_config.json
