#!/bin/zsh
# ==============================================================================
# Fluxon-MHD Modeling Framework Installation Script for macOS using perlbrew
# ==============================================================================

export PL_PREFIX="$HOME/Library/perl"
export FL_PREFIX="$HOME/Library/flux"
export FL_MHDLIB="$HOME/flux"
export PROFILE_FILE="$HOME/.zshenv"  # Use .zshenv for all profile-related changes
export DESIRED_PERL_VERSION="perl-5.36.0"

# ==============================================================================
# Save environment variables to the profile file if they have changed
# ==============================================================================
ENV_VARS=(
    "export PL_PREFIX=\"$PL_PREFIX\""
    "export FL_PREFIX=\"$FL_PREFIX\""
    "export FL_MHDLIB=\"$FL_MHDLIB\""
    "export PROFILE_FILE=\"$PROFILE_FILE\""
    "export TARGET_LIB_DIR=\"$PL_PREFIX/lib/perl5\""
    "export DESIRED_PERL_VERSION=\"$DESIRED_PERL_VERSION\""
)

# Define color codes for the script
orange="\e[32m"
reset="\e[0m"

# Wrapper function for colored echo
colored_echo() {
    echo -e "${orange}$1${reset}"
}

# Function to update profile and source it right after
update_profile_file_and_source() {
    local var_to_check="$1"
    local file="$2"
    local var_name=$(echo "$var_to_check" | cut -d'=' -f1)

    # Append safely only if the variable doesn't exist, ensure proper formatting
    if ! grep -q "^$var_name=" "$file"; then
        echo "$var_to_check" >> "$file"
        colored_echo "Added $var_name to $file"
    else
        # Replace existing line safely
        sed -i '' "s|^.*$var_name=.*|$var_to_check|" "$file"
        colored_echo "Updated $var_name in $file"
    fi

    # Source profile after every update to apply the changes
    source "$file"
    colored_echo "Sourced $file to apply changes."
}

# Function to add a path to an environment variable if it doesn't already exist (zsh compatible)
add_to_env_var_if_not_exists() {
    local env_var_name="$1"
    local new_value="$2"
    local file="$3"

    local current_value=$(printenv "$env_var_name")

    if [[ ":$current_value:" != *":$new_value:"* ]]; then
        export "$env_var_name"="${current_value:+$current_value:}$new_value"
        update_profile_file_and_source "export $env_var_name=\"${current_value:+$current_value:}$new_value\"" "$file"
    fi
}

colored_echo "Sourcing profile to apply changes..."
source "$PROFILE_FILE"

# ==============================================================================
# Helper Functions
# ==============================================================================
function check_error {
    if [ $? -ne 0 ]; then
        colored_echo "Error encountered. Running diagnostics..."
        DIAGNOSTICS_DIR="$FL_PREFIX/fluxon-mhd/doc"
        if [ -f "$DIAGNOSTICS_DIR/flux_diagnostics.pl" ]; then
            perl "$DIAGNOSTICS_DIR/flux_diagnostics.pl"
        else
            colored_echo "Diagnostics script flux_diagnostics.pl not found in $DIAGNOSTICS_DIR"
        fi
        exit 1
    fi
}

# ==============================================================================
# Xcode Command Line Tools Installation
# ==============================================================================
if ! xcode-select --print-path &>/dev/null; then
    colored_echo "Xcode Command Line Tools not found. Installing them..."
    xcode-select --install
else
    colored_echo "Xcode Command Line Tools already installed."
fi

# ==============================================================================
# Homebrew Installation (for non-Perl dependencies)
# ==============================================================================
if ! command -v brew &> /dev/null; then
    colored_echo "Homebrew not found. Installing Homebrew..."
    BREW_OUTPUT=$(mktemp)
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" 2>&1 | tee "$BREW_OUTPUT"
    check_error

    colored_echo "Homebrew installed successfully."

    SHELL_SETUP_CMD=$(grep 'eval' "$BREW_OUTPUT" | head -n 1)

    if [ -n "$SHELL_SETUP_CMD" ]; then
        if ! grep -q "$SHELL_SETUP_CMD" "$PROFILE_FILE"; then
            colored_echo "Adding Homebrew shell setup to $PROFILE_FILE"
            echo "$SHELL_SETUP_CMD" >> "$PROFILE_FILE"
            source $PROFILE_FILE
            check_error
        fi
    else
        colored_echo "Error: Could not find the Homebrew shell setup command."
        exit 1
    fi
    rm "$BREW_OUTPUT"
else
    colored_echo "Homebrew is already installed."
fi

source "$PROFILE_FILE"

colored_echo "Installing gnuplot and fftw via Homebrew..."
brew install gnuplot fftw cpanm qt
check_error

# ==============================================================================
# Perlbrew Installation and Perl Version Setup
# ==============================================================================
if ! command -v perlbrew &> /dev/null; then
    colored_echo "Perlbrew not found. Installing perlbrew..."
    PERLBREW_OUTPUT=$(mktemp)
    curl -L https://install.perlbrew.pl | bash 2>&1 | tee "$PERLBREW_OUTPUT"
    check_error

    SHELL_SETUP_CMD=$(grep 'source ~' "$PERLBREW_OUTPUT" | head -n 1)

    if [ -n "$SHELL_SETUP_CMD" ]; then
        if ! grep -q "$SHELL_SETUP_CMD" "$PROFILE_FILE"; then
            colored_echo "Adding Perlbrew shell setup to $PROFILE_FILE"
            echo "$SHELL_SETUP_CMD" >> "$PROFILE_FILE"
            source "$PROFILE_FILE"
            check_error
        fi
    else
        colored_echo "Error: Could not find the Perlbrew shell setup command."
        exit 1
    fi
    rm "$PERLBREW_OUTPUT"
else
    colored_echo "Perlbrew is already installed."
    perlbrew switch "$DESIRED_PERL_VERSION"
    perlbrew --version
    check_error
fi

if ! perlbrew list | grep -q "$DESIRED_PERL_VERSION"; then
    colored_echo "$DESIRED_PERL_VERSION not found. Installing $DESIRED_PERL_VERSION..."
    perlbrew install --notest $DESIRED_PERL_VERSION &
    install_pid=$!
    sleep 2
    osascript <<EOF
tell application "Terminal"
    do script "sleep 2 && tail -f ~/perl5/perlbrew/build.$DESIRED_PERL_VERSION.log"
end tell
EOF
    wait $install_pid
    check_error
    perlbrew switch $DESIRED_PERL_VERSION
    check_error
else
    colored_echo "$DESIRED_PERL_VERSION is already installed."
    perlbrew switch $DESIRED_PERL_VERSION
    check_error
fi

colored_echo "Setting up local::lib with perlbrew..."
PERL_MM_OPT="INSTALL_BASE=$PL_PREFIX" cpanm --local-lib="$PL_PREFIX" local::lib
check_error

if ! grep -q "PERL5LIB=.*$TARGET_LIB_DIR" "$PROFILE_FILE"; then
    echo "export PERL5LIB=\"$TARGET_LIB_DIR:\$PERL5LIB\"" >> "$PROFILE_FILE"
    echo "eval \"\$(perl -I$TARGET_LIB_DIR -Mlocal::lib=$PL_PREFIX)\"" >> "$PROFILE_FILE"
    check_error
fi

colored_echo "Installing required Perl modules via cpanm..."
cpanm PDL
cpanm File::ShareDir
cpanm PDL::Graphics::Gnuplot
cpanm Math::RungeKutta
cpanm Term::ReadKey
check_error

PERL_VERSION=$(perl -e 'print $^V;')
colored_echo "Library instantiated at $TARGET_LIB_DIR, Perl version: $PERL_VERSION"

PERLDLRC_FILE="$HOME/.perldlrc"
if [ ! -f "$PERLDLRC_FILE" ]; then
    touch "$PERLDLRC_FILE"
fi

if ! grep -q "$TARGET_LIB_DIR" "$PERLDLRC_FILE"; then
    cat <<EOL >> "$PERLDLRC_FILE"
push(\@INC, "$TARGET_LIB_DIR");

require('PDL/default.perldlrc');
use PDL::AutoLoader;
\$PDL::AutoLoader::Rescan=1;
1;
EOL
    colored_echo "Configuration added to ~/.perldlrc"
fi

# ==============================================================================
# Fluxon-MHD Installation
# ==============================================================================
if [ ! -d "$FL_MHDLIB" ]; then
    colored_echo "Cloning the fluxon-mhd repository..."
    git clone https://github.com/lowderchris/fluxon-mhd.git "$FL_MHDLIB"
    check_error
fi

colored_echo "Navigating to the flux directory and compiling the project..."
cd "$FL_MHDLIB" || exit
PDLINSTALL_OUTPUT=$(mktemp)
make libbuild libinstall pdlbuild pdlinstall 2>&1 | tee "$PDLINSTALL_OUTPUT"
check_error

PDLLIB_DIR=$(grep -o '/Users[^ ]*auto/share/dist/Flux' "$PDLINSTALL_OUTPUT" | head -n 1)

if [ -n "$PDLLIB_DIR" ]; then
    colored_echo "Found Flux autoload directory: $PDLLIB_DIR "
    add_to_env_var_if_not_exists "PDLLIB" "+$PDLLIB_DIR" "$PROFILE_FILE"
else
    colored_echo "Error: Could not find the Flux autoload directory in the output."
    exit 1
fi

rm "$PDLINSTALL_OUTPUT"

colored_echo "Fluxon-MHD setup complete!"