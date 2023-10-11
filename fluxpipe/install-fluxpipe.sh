#!/bin/zsh

# Ensure the script stops on first error
set -e

echo "\n\nInstalling FluxPipe..."

# Check if Homebrew is installed
if ! command -v brew &>/dev/null; then
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install gnuplot
    brew install fftw
    brew install qt

else
    echo "Homebrew already installed!"
fi



# Check if perlbrew is installed
if ! command -v perlbrew &>/dev/null; then
    \curl -L https://install.perlbrew.pl | bash
    perlbrew install perl-5.36.1
    perlbrew switch perl-5.36.1
    curl -L https://cpanmin.us | perl - --sudo App::cpanminus || cpan App::cpanminus
else
    echo "perlbrew already installed!"
fi

# Check if conda is initialized for zsh
(
if [ ! -f ~/.zshrc ] || ! grep -q "conda.sh" ~/.zshrc; then
    conda init zsh
    conda init bash
    source ~/.zshrc
else
    echo "Anaconda already initialized in zshrc!"
fi
)

# make it so zshenv calls zshrc
grep -qF "source ~/.zshrc" ~/.zshenv || echo "source ~/.zshrc" >> ~/.zshenv



# Check if the fluxenv conda environment exists

conda_envs=$(conda info --envs)
if ! echo "$conda_envs" | grep -q "fluxenv"; then
    # Your code to create the fluxenv if it doesn't exist
    yes | conda create --name fluxenv --file requirements-conda.txt
    yes | conda create -n fluxenv
    conda activate fluxenv

fi

# Source the PREFIX_PATHS.sh script
chmod +x PREFIX_PATHS.sh
./PREFIX_PATHS.sh

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
cd $FL_PREFIX/fluxpipe
cpanm --installdeps .
cpanm --notest --local-lib=local/ .
echo 'export PERL5LIB=$FL_PREFIX/fluxpipe/local/lib/perl5:$PERL5LIB' >> ~/.zshrc

(
    conda activate fluxenv
    yes | conda install pip
    yes | pip install -e .
    yes | pip install -r requirements-pip.txt
)
source ~/.zshrc
conda activate fluxenv

echo "Fluxpipe installation complete!\n\n"
echo "Remember to run 'conda activate fluxenv' before running 'python fluxpipe/fluxpipe/runners/config_runner.py'\n\n"

read "response?Do you want to test the config_runner? ([yes]/no): "
if [[ "$response" == "yes" || "$response" == "y" || "$response" == "" ]]; then
    echo "Running the file..."
    cd $FL_PREFIX/fluxpipe/fluxpipe
    python runners/config_runner.py
fi

