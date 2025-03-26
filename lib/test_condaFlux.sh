#!/bin/bash

# Test if Flux is installed by trying to load it with Perl.
if perl -e "use Flux;" ; then
    echo "\nCondaFlux is installed and working!\n"
else
    echo "\nCondaFlux is NOT installed!\n"
fi

# Abort if the correct conda environment isn't active
if [[ "$CONDA_PREFIX" != "$FL_PREFIX" ]]; then
  echo "Error: This script must be run from within the 'fluxenv' conda environment." >&2
  echo "       Run 'conda activate fluxenv' and try again.\n" >&2
  exit 1
fi