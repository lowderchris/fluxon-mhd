#!/bin/bash

echo "Uninstalling condaFLUX..."
# Source conda's initialization script to enable conda commands.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the conda environment 'fluxenv'
conda activate fluxenv

# Remove specific FLUX files and directories.
rm -rf "${FL_PREFIX}/lib/libflux.a"
rm -rf "${FL_PREFIX}/include/flux"

# Remove the first directory in PDLLIB if it contains "flux".
d_pdllib=$(echo "$PDLLIB" | awk -F: '{print $1}')
if echo "$d_pdllib" | grep -qi flux; then
    echo "Deleting $d_pdllib"
    rm -rf "$d_pdllib"
else
    echo "Skipping deletion of PDLLIB directory: $d_pdllib (does not contain 'flux')"
fi

# Remove the first directory in PERL5LIB if it contains "flux".
d_perl5lib=$(echo "$PERL5LIB" | awk -F: '{print $1}')
if echo "$d_perl5lib" | grep -qi flux; then
    rm -rf "$d_perl5lib"
else
    echo "Skipping deletion of PERL5LIB directory: $d_perl5lib (does not contain 'flux')\n"
fi

# Run test_conda.sh from the lib directory if it exists.
if [ -e "test_condaFlux.sh" ]; then
    (source test_condaFlux.sh)
fi

echo "\ncondaFlux uninstall complete. \n"