#!/bin/bash

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
  if grep -Fxq "$required_content" "$perldlrc"; then
    echo "Content already exists in ~/.perldlrc"
  else
    echo "$required_content" >> "$perldlrc"
    echo "Appended required content to ~/.perldlrc"
  fi
else
  echo "$required_content" > "$perldlrc"
  echo "Created ~/.perldlrc and added required content"
fi

# Create activate.d directory if it doesn't exist
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"

# Locate Flux.pm and set PERL5LIB
export PERL5LIB=$(find "$CONDA_PREFIX" -name Flux.pm -exec dirname {} \; | tr '\n' ':')$PERL5LIB

# Locate Flux module directory and set PDLLIB
export PDLLIB=$(find "$CONDA_PREFIX" -type d -name "Flux" 2>/dev/null | grep "/share/dist/" | awk 'NR==1'):$PDLLIB

# Create the environment variable setup script
cat > "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh" <<EOL
#!/bin/bash

# Set PERL5LIB to the directory containing Flux.pm
export PERL5LIB=$PERL5LIB

# Set PDLLIB to the Flux module directory
export PDLLIB=$PDLLIB

# Print the environment variables for verification
echo "PERL5LIB is set to: \$PERL5LIB"
echo "PDLLIB is set to: \$PDLLIB"
EOL

# Make the activation script executable
chmod +x "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh"

# Confirmation message
echo "Activation script created successfully. Environment variables will be set upon activating the conda environment."
conda activate fluxenv
