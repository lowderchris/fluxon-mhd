#!/bin/bash

# Source conda's initialization script to enable conda commands.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the conda environment 'fluxenv'
conda activate fluxenv

# Required content for the ~/.perldlrc file
required_content=$(cat <<'EOF'
require(q|PDL/default.perldlrc|);
use PDL::AutoLoader;
$PDL::AutoLoader::Rescan=1;
1;
EOF
)

# Path to the ~/.perldlrc file
perldlrc="$HOME/.perldlrc"

# Check if the required content is already in the ~/.perldlrc file
if [ -e "$perldlrc" ]; then
  file_content=$(cat "$perldlrc")
  if [[ "$file_content" == *"$required_content"* ]]; then
    echo "Content already exists in ~/.perldlrc"
  else
    echo "$required_content" >> "$perldlrc"
    echo "Appended required content to ~/.perldlrc"
  fi
else
  echo "$required_content" > "$perldlrc"
  echo "Created ~/.perldlrc and added required content"
fi

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: \$CONDA_PREFIX is not set" >&2
    return 1
fi

# Create activate.d directory if it doesn't exist
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"

# Create the environment variable setup script (activation script)
cat > "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh" <<EOL
#!/bin/bash

# Dynamically locate Flux.pm and set PERL5LIB
export PERL5LIB=\$(find "\$CONDA_PREFIX" -name Flux.pm -exec dirname {} \; | head -n 1)
# echo "Found PERL5LIB: \$PERL5LIB"

# Dynamically locate the Flux module directory and set PDLLIB
export PDLLIB=\$(find "\$CONDA_PREFIX" -type d -name "Flux" 2>/dev/null | grep "/share/dist/" | awk 'NR==1')
# echo "Found PDLLIB: \$PDLLIB"

# Set the FLUX prefixes for future recompiling
export PL_PREFIX="\$CONDA_PREFIX/lib/perl5/site_perl"
export FL_PREFIX="\$CONDA_PREFIX"

# Print the environment variables for verification
echo "FL_PREFIX is set to: \$FL_PREFIX"
echo "PL_PREFIX is set to: \$PL_PREFIX"
echo "PERL5LIB  is set to: \$PERL5LIB"
echo "PDLLIB    is set to: \$PDLLIB" \n
EOL

# Make the activation script executable
chmod +x "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh"

# Create the deactivation script directory
mkdir -p "$CONDA_PREFIX/etc/conda/deactivate.d"

# Create deactivation script for restoring original environment variables
cat > "$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh" <<'EOL'
#!/bin/bash

# Deactivation script for restoring original environment variables

  unset PERL5LIB
  unset PDLLIB
  unset PL_PREFIX
  unset FL_PREFIX

# Inform the user
echo "Deactivation: Environment variables have been unset. \n"
EOL

# Make the deactivation script executable
chmod +x "$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh"


# Confirmation message
echo "Activation and deactivation scripts created successfully. Environment variables will be managed upon conda activation/deactivation. \n"
