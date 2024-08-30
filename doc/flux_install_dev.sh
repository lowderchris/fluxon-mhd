#!/bin/zsh

# ==============================================================================
# Fluxon-MHD Modeling Framework Installation Script for macOS
#
# This script automates the installation of the Fluxon-MHD modeling framework
# and its dependencies on macOS. The process involves three main steps:
#
# 1. Perl Installation: Installs Perl using Homebrew, sets up local Perl libraries
#    using local::lib, and configures cpanminus for easy module management.
#
# 2. Perl Data Language (PDL) Installation: Installs PDL and additional required
#    Perl modules using cpanminus.
#
# 3. Fluxon-MHD Installation: Clones the Fluxon-MHD repository, sets up necessary
#    environment variables, and compiles the software.
#
# Some steps require admin/root access, which will be indicated where necessary.
# ==============================================================================

# ==============================================================================
# Environment Setup
# ==============================================================================
# Set environment variables for installation directories
# PL_PREFIX: Directory for local Perl libraries
# FL_PREFIX: Directory where Fluxon-MHD will be installed
export PL_PREFIX="$HOME/Library/perl5"
export FL_PREFIX="$HOME/Library/flux"

# ==============================================================================
# Perl Installation
# ==============================================================================
# Homebrew is used to install Perl, as it provides an easy way to manage packages
# on macOS. Admin access is required to install Homebrew and Perl.
if ! command -v brew &> /dev/null; then
    echo "Homebrew not found. Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    echo "Homebrew installed. Please run the commands displayed by the installer to add Homebrew to your PATH."
    exit 1
fi

# Install Perl using Homebrew (Admin required)
echo "Installing Perl via Homebrew..."
brew install perl

# Optionally, pin Perl to prevent updates (Admin required)
# Pinning the Perl version ensures that it won't be automatically updated,
# which could introduce compatibility issues with your Perl environment.
echo "Pinning Perl to prevent automatic updates..."
brew pin perl

# Install system dependencies (Admin required)
# gnuplot: Required for generating plots
# fftw: Required for performing Fast Fourier Transforms
echo "Installing gnuplot and fftw via Homebrew..."
brew install gnuplot fftw

# Install cpanminus (Admin required)
# cpanminus is a lightweight Perl module installer that simplifies the process of
# managing Perl modules, making it easier to install required dependencies.
echo "Installing cpanminus..."
curl -L https://cpanmin.us | perl - --sudo App::cpanminus

# ==============================================================================
# Perl Local Library Setup
# ==============================================================================
# Set up local::lib to manage Perl modules locally within the user's home directory.
# This allows Perl modules to be installed without requiring admin access.
echo "Setting up local::lib..."
PERL_MM_OPT="INSTALL_BASE=$PL_PREFIX" cpanm --local-lib="$PL_PREFIX" local::lib

# Add local::lib setup to .zprofile for Zsh
# This ensures that the local::lib configuration is applied automatically whenever
# the user opens a new terminal session.
echo 'eval "$(perl -I$HOME/Library/perl5/lib/perl5 -Mlocal::lib=$HOME/Library/perl5)"' >> ~/.zprofile
source ~/.zprofile

# ==============================================================================
# Configure Perl Data Language (PDL)
# ==============================================================================
# Define the Perl version and local library path
# The Perl version is dynamically determined to ensure the correct path is used
# when configuring the PDL environment.
PERL_VERSION=$(perl -e 'print $^V;')
LOCAL_LIB_DIR="$PL_PREFIX/lib/perl5/$PERL_VERSION"

# Create ~/.perldlrc if it doesn't exist and configure it
# The .perldlrc file is used to configure the Perl Data Language (PDL) environment.
# Here, we ensure that PDL is correctly set up to use the local::lib directory.
echo "Configuring ~/.perldlrc..."
if [ ! -f "$HOME/.perldlrc" ]; then
    touch "$HOME/.perldlrc"
fi

cat <<EOL >> "$HOME/.perldlrc"
push(\@INC, "$LOCAL_LIB_DIR");

require('PDL/default.perldlrc');

use PDL::AutoLoader;
\$PDL::AutoLoader::Rescan=1;

1;
EOL

echo "Configuration added to ~/.perldlrc"

# ==============================================================================
# Install Perl Data Language (PDL) and Dependencies
# ==============================================================================
# Install required Perl modules using cpanminus
# PDL: The Perl Data Language, essential for numerical calculations in Perl
# File::ShareDir: Manages shared files in Perl distributions
# PDL::Graphics::Gnuplot: Interface between PDL and Gnuplot for plotting
# Math::RungeKutta: Provides Runge-Kutta methods for numerical integration
# Term::ReadKey: Allows reading of keystrokes in terminal applications
echo "Installing required Perl modules..."
cpanm PDL
cpanm File::ShareDir
cpanm PDL::Graphics::Gnuplot
cpanm Math::RungeKutta
cpanm Term::ReadKey

# ==============================================================================
# Fluxon-MHD Installation
# ==============================================================================
# Set environment variables for the flux installation
# These variables define the installation directories for Perl libraries and
# the Fluxon-MHD software.
echo "Setting up environment variables for FLUX installation..."
export PL_PREFIX="$HOME/Library/perl5"
export FL_PREFIX="$HOME/Library/flux"

# Clone the fluxon-mhd repository if not already cloned
# If the repository is not already present, clone it from GitHub into the
# specified directory. This ensures that the latest version of the software
# is available for installation.
if [ ! -d "$FL_PREFIX" ]; then
    echo "Cloning the fluxon-mhd repository..."
    git clone https://github.com/yourusername/fluxon-mhd.git "$FL_PREFIX"
fi

# Navigate to the flux directory and compile
# Run the make command to compile the Fluxon-MHD software. The 'make everything'
# command compiles the entire project, ensuring that all components are ready for use.
echo "Navigating to the flux directory and compiling the project..."
cd "$FL_PREFIX" || exit
make everything

# ==============================================================================
# Completion Message
# ==============================================================================
echo "Fluxon-mhd setup complete!"