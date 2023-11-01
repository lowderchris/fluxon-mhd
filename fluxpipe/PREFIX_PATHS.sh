#!/bin/zsh

# Definitions: Adjust these values to your system configuration
PL_PATH="/Users/cgilbert/perl5/perlbrew/perls/perl-5.32.0"
FL_PATH="/Users/cgilbert/vscode/fluxons/fluxon-local"
FL_MHDLIB="/Users/cgilbert/vscode/fluxons/fluxon-mhd"
SHELL_RC_CUSTOM="$HOME/.zshrc_custom"

# Function to update or set the environment variables
update_var() {
    var_name=$1
    new_path=$2
    current_path=$(eval echo \$$var_name)

    # Check and update the shell rc custom file
    if grep -q "^export $var_name=" $SHELL_RC_CUSTOM; then
        # If the line exists, replace it
        awk -v var="$var_name" -v path="$new_path" '
            $0 ~ "^export " var "=" { print "export " var "=" path; next }
            { print }
        ' $SHELL_RC_CUSTOM > temp_file && mv temp_file $SHELL_RC_CUSTOM
        echo "Updated $var_name to $new_path in $SHELL_RC_CUSTOM"
    else
        # If the line doesn't exist, append it
        echo "export $var_name=$new_path" >> $SHELL_RC_CUSTOM
        echo "Added $var_name to $new_path in $SHELL_RC_CUSTOM"
    fi
    # Export the new value
    export $var_name=$new_path
}

# Update or set PL_PREFIX
update_var PL_PREFIX $PL_PATH

# Update or set FL_PREFIX
update_var FL_PREFIX $FL_PATH

# Update or set FL_PREFIX
update_var FL_MHDLIB $FL_MHDLIB

echo "\nFLUX Paths set!"
