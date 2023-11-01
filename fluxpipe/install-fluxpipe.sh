#!/bin/zsh

# Ensure the script stops on first error
set -e

[[ ! $(pwd) =~ fluxpipe ]] && cd fluxpipe

echo "\n\nInstalling FluxPipe..."

echo "Setting up PREFIX_PATHS.sh"
# Source the PREFIX_PATHS.sh script
chmod +x PREFIX_PATHS.sh
./PREFIX_PATHS.sh


# Check if Homebrew is installed
if ! command -v brew &>/dev/null; then
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    sudo brew install gnuplot
    sudo brew install fftw
    sudo brew install qt

else
    echo "\tHomebrew already installed!"
fi



# Check if perlbrew is installed
if ! command -v perlbrew &>/dev/null; then
    \curl -L https://install.perlbrew.pl | bash
    perlbrew install perl-5.36.1
    perlbrew switch perl-5.36.1
    curl -L https://cpanmin.us | perl - --sudo App::cpanminus || cpan App::cpanminus
else
    echo "\tperlbrew already installed!"
fi



# Check if conda command exists
if command -v conda &>/dev/null; then
    # Check if conda is initialized for zsh
    if [ ! -f ~/.zshrc ] || ! grep -q "conda.sh" ~/.zshrc; then
        conda init zsh
        conda init bash
        source ~/.zshrc
        # You can add conda-specific operations here if needed
    else
        echo "\tAnaconda already initialized in zshrc!"
    fi

    # Check if the fluxenv conda environment exists
    conda_envs=$(conda info --envs)
    if ! echo "$conda_envs" | grep -q "fluxenv"; then
        # Your code to create the fluxenv if it doesn't exist
        yes | conda create --name fluxenv --file requirements-conda.txt
        yes | conda create -n fluxenv
        conda activate fluxenv
    fi

else
    echo "\tConda is not installed. Falling back to virtualenv and pip."

    # Check if virtualenv is installed, and if not, install it
    if ! command -v virtualenv &>/dev/null; then
        pip install virtualenv
    fi

    # Create a new virtual environment
    virtualenv fluxenv_pip

    # Activate the virtual environment
    source fluxenv_pip/bin/activate

    # Install packages from the two requirements.txt files
    pip install -r requirements-pip.txt
    # Install packages from requirements-conda.txt and continue even if it fails
    pip install -r requirements-conda.txt || true
fi

echo "Python Environment Creation:"
# cd fluxpipe
(
conda activate fluxenv
yes | conda install pip
yes | pip install -e .
yes | pip install -r requirements-pip.txt
)
source ~/.zshrc
conda activate fluxenv

# make it so zshenv calls zshrc
grep -qF "source ~/.zshrc" ~/.zshenv || echo "source ~/.zshrc" >> ~/.zshenv

# echo the current directory


# Check and install PDL and other dependencies
cd $PL_PREFIX
if ! perl -MPDL -e 1 &>/dev/null; then
    cpanm PDL
fi
export PERL_MM_OPT="INSTALL_BASE=$PL_PREFIX/lib/perl5"
cpanm local::lib
echo 'eval "$(perl -I$PL_PREFIX/lib/perl5/lib/perl5 -Mlocal::lib=$PL_PREFIX/lib/perl5)"' >> ~/.zshrc
cpanm File::ShareDir
cpanm PDL::Graphics::Gnuplot
cpanm Math::RungeKutta
# ... Add other dependencies as needed ...

# Install FLUXpipe pdl and python
(
cd $FL_MHDLIB/fluxpipe
cpanm --installdeps .
cpanm --notest --local-lib=$FL_PREFIX .
if ! grep -q 'export PERL5LIB=$FL_PREFIX/lib/perl5:$PERL5LIB' ~/.zshrc_custom; then
    echo 'export PERL5LIB=$FL_PREFIX/lib/perl5:$PERL5LIB' >> ~/.zshrc_custom
fi
rm -rf blib\
)


echo "\n\nFluxpipe installation complete!\n\n"

if command -v conda &>/dev/null; then
    echo "Remember to run 'conda activate fluxenv' before running 'python fluxpipe/fluxpipe/runners/config_runner.py'\n\n"
else
    echo "Remember to run 'source fluxenv_pip/bin/activate' before running 'python fluxpipe/fluxpipe/runners/config_runner.py'\n\n"
fi

echo -n "Do you want to test the config_runner? ([yes]/no): "
read response
if [[ "$response" == "yes" || "$response" == "y" || "$response" == "" ]]; then
    echo "Running the file..."
    cd $FL_MHDLIB/fluxpipe/fluxpipe
    conda activate fluxenv
    python runners/config_runner.py
fi
