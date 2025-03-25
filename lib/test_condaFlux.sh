#!/bin/bash

# Source conda's initialization script to enable conda commands.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the conda environment 'fluxenv'
conda activate fluxenv

# Test if Flux is installed by trying to load it with Perl.
if perl -e "use Flux;" ; then
    echo "\nFlux is installed and working!\n"
else
    echo "\nConda is NOT installed!\n"
fi