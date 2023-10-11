# Macintosh Installation

Here's the path to installing FLUX. Best to run everything as Admin if possible.
All of these commands should be typed into the terminal.
If you don't use zsh, you'll need to substitute the shell you use.

## Makefile Installation -----------------------------------

The whole installation process is being integrated into the makefile. See install_flux.sh.

#### Put the PL and FL_PREFIX paths into PREFIX_PATHS.sh

 Examples:
  ``PL_PATH="/Users/cgilbert/perl5/perlbrew/perls/perl-5.32.0"``
  ``FL_PATH="/Users/cgilbert/vscode/fluxons/fluxon-mhd"``

#### **If you don't have python installed, get anaconda:**

- Download the installer from [https://www.anaconda.com/products/individual](https://www.anaconda.com/products/individual)
- Run the installer

#### Build and Run the Makefile

``perl Makefile.PL``
``cd ..``
``make everything`` (The main fluxon-mhd one)

## Manual Installation --------------------------------------

This is already depricated because all of these steps have been moved into the makefile

#### Install Homebrew and Packages: ([https://brew.sh/](https://brew.sh/))) ---------------

- **To get Homebrew, use:**
  ``% /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"``
- **Install packages**
  ``% brew install gnuplot``
  ``% brew install fftw``

#### Install Perl -------------------------------------------

- **Get perlbrew (recommended, [https://perlbrew.pl/](https://perlbrew.pl/))**
  ``% \curl -L https://install.perlbrew.pl | bash``
- **Install the latest version of perl**
  ``% perlbrew install perl-5.36.1``
- **Set the system perl to this version**
  ``% perlbrew switch perl-5.36.1``
- **Install cpanminus, the perl package manager**
  ``% curl -L https://cpanmin.us | perl - --sudo App::cpanminus``
  **~ OR ~**
  ``% cpan App::cpanminus``

#### Install Python ---------------------------------------

- **If you don't have python installed, start by installing anaconda:**

  - Download the installer from [https://www.anaconda.com/products/individual](https://www.anaconda.com/products/individual)
  - Run the installer
- **Prepare the Environment:**

  - initialize conda (once ever)
    ``% conda init zsh``
  - restart your terminal (after init)
    ``% source ~/.zshrc``  ~ OR ~   %``zsh``
  - make a new environment using conda and activate it
    ``conda create -n fluxenv``
    ``conda activate fluxenv``
    - If you are using an Apple Silicon M1/M2 architecture, please see the footnotes.
  - navigate to the fluxpipe directory and
    install the package by running
    ``pip install -e .``

#### Set the environment variables ---------------------

* **Set the two prefixes**
  ``export FL_PREFIX='_absolute/_path/_to/fluxon-mhd' ``
  ``export PL_PREFIX=_absolute/_path/_to/perl ``

  * example: ``export FL_PREFIX=$HOME/vscode/fluxons/fluxon-mhd ``
    example: ``export PL_PREFIX=/$HOME/perl5/perlbrew/perls/perl-5.32.0 ``
* **Add them to your environment permanently by running the following commands:**
  ``echo 'export PL_PREFIX='$PL_PREFIX >> ~/.zshrc ``
  ``echo 'export FL_PREFIX='$FL_PREFIX >> ~/.zshrc ``
* **make it so zshenv calls zshrc**
  ``grep -qF "source ~/.zshrc" ~/.zshenv || echo "source ~/.zshrc" >> ~/.zshenv``

#### Install PDL and other dependencies --------------

- **Install PDL and its dependencies**
  - Go to Perl directory
    ``% cd $PL_PREFIX ``
  - Install PDL
    ``% cpanm PDL ``
  - Configure 'local::lib' to set the local build path (leave this as-is).
    ``% PERL_MM_OPT="INSTALL_BASE=$PL_PREFIX/lib/perl5"``
    ``% cpanm local::lib ``
    ``% echo 'eval "$(perl -I$PL_PREFIX/lib/perl5/lib/perl5 -Mlocal::lib=$PL_PREFIX/lib/perl5)"' >> ~/.zshrc ``
- **Install other dependencies of the project**
  ``% cpanm File::ShareDir ``
  ``% cpanm PDL::Graphics::Gnuplot ``
  ``% cpanm Math::RungeKutta ``
  and others...

#### Install FLUX -----------------------------------------

- **Get Flux**
  ``% cd $FL_PREFIX/.. ``
  ``% git clone https://github.com/lowderchris/fluxon-mhd.git ``
- **Compile Flux**
  ``% cd $FL_PREFIX ``
  ``% sudo make everything ``

### Install FLUXpipe -------------------------------

- Install the fluxpipe pdl package locally
  ``% cd $FL_PREFIX/fluxpipe ``
  ``% cpanm --installdeps . ``
  ``% cpanm --notest --local-lib=local/ .``
- Set the perl path to the local build and reload the terminal
  ``echo 'export PERL5LIB=$FL_PREFIX/fluxpipe/local/lib/perl5:$PERL5LIB' >> ~/.zshrc ``
  ``source ~/.zshrc``
- - Unfortunately this build must be recompiled after every change. To do this, run
    ``cd $FL_PREFIX/fluxpipe ``
    ``cpanm --notest --local-lib=local/ .``
  - This can be configured to happen automatically in VSCode. You'll need the File Watcher extension.
    If you place the following lines in your settings.json file in vscode, it will keep your pdl build up to date whenever you save a file.

    ``"filewatcher.commands": [ { "match": ".*\\.(pdl|pm|json)$", "isAsync": true, "event": "onFileChange", "cmd": "cd $FL_PREFIX/fluxpipe && cpanm --notest --local-lib=local/ . && echo $(date) " } ],``

### Footnotes --------------------------------------

[1] Functions to Create Conda Envs for Specific Architectures

**Create x86 conda environment**
``echo 'conda_create_x86 () { CONDA_SUBDIR=osx-64 conda create -n $@; conda activate $1; }' >> ~/.zshrc``

**Create ARM conda environment**
``echo 'conda_create_ARM () { CONDA_SUBDIR=osx-arm64 conda create -n $@; conda activate $1; }' >> ~/.zshrc ``

**Run any command with a specified architecture**
``arch -arm64 <command> <args>``

examples:
    ``conda_create_ARM myenv_ARM python=3.9``
    ``conda_create_x86 myenv_x86 python=3.9``
    ``arch -arm64 brew install qemu``
