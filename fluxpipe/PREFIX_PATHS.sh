#!/bin/zsh

# Definitions: Adjust these values to your system configuration
PL_PATH="/Users/cgilbert/perl5/perlbrew/perls/perl-5.32.0"
FL_PATH="/Users/cgilbert/vscode/fluxons/fluxon-mhd"
SHELL_RC="~/.zshrc"


## Do not adjust this code --------------------------------
# Check if PL_PREFIX is not set
if [ -z "$PL_PREFIX" ]; then
    export PL_PREFIX=$PL_PATH
    echo "export PL_PREFIX=$PL_PREFIX" >> $SHELL_RC
    echo "PL_PREFIX set to $PL_PREFIX"
else
    echo "PL_PREFIX is already set to $PL_PREFIX"
fi

# Check if FL_PREFIX is not set
if [ -z "$FL_PREFIX" ]; then
    export FL_PREFIX=$FL_PATH
    echo "export FL_PREFIX=$FL_PREFIX" >> $SHELL_RC
    echo "FL_PREFIX set to $FL_PREFIX"
else
    echo "FL_PREFIX is already set to $FL_PREFIX"
fi

echo "\nFLUX Paths set!"

