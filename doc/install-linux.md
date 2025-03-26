# Linux Installation Instructions

These instructions were created by installing FLUX on two different linux servers running AlmaLinux. Please let us know if anything breaks! Extra troubleshooting steps are listed at the bottom.

## 1. Get Conda Environment Manager

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompts:
- Hold down the arrow key and then accept terms.
- Accept Conda activation at startup.
- Exit and reconnect to the SSH session/ shell session.

## 2. Create Conda Environment and Get Dependencies

```sh
conda create -n fluxenv
conda activate fluxenv
conda install -c conda-forge gcc gxx binutils c-compiler gsl perl perl-app-cpanminus perl-extutils-makemaker make cmake automake autoconf libtool m4 patch libxcrypt gnuplot cairo pango qt pfsspy sqlite
```

Confirm installation by typing `y` when prompted.

## 3. Install the Perl Data Language (PDL) and Dependencies

```sh
cpanm Alien::Build::Plugin::Gather::Dino Capture::Tiny Chart::Gnuplot Config::IniFiles Devel::CheckLib File::HomeDir File::Map File::ShareDir File::ShareDir::Install File::Which Inline Inline::C Inline::Python List::MoreUtils Math::GSL Math::GSL::Alien Math::Interpolate Math::Interpolator Math::RungeKutta Moo::Role Net::SSLeay PDL PDL::GSL PDL::GSL::INTEG PDL::Graphics::Gnuplot PDL::Graphics::Simple Parallel::ForkManager Term::ReadKey Test::Builder Text::CSV local::lib
```

If the system struggles with `File::Map`, try:

```sh
export LDFLAGS=""
export CFLAGS="-O2"
export PERL_LDFLAGS=""
export PERL_CFLAGS="-O2"
export LD=x86_64-conda-linux-gnu-gcc
export CC=x86_64-conda-linux-gnu-gcc
```

If `Inline::C` isn't found, run:

```sh
cpanm Inline::C
```

## 4. Get FLUX

Clone and install Flux:

```sh
git clone https://github.com/lowderchris/fluxon-mhd.git
cd fluxon-mhd
make condaflux
```
## 5. Get Fluxpype Wrapper (To run FLUX with)

```sh
cd ..
git clone https://github.com/GillySpace27/fluxpype.git
cd fluxpype
pip install -e .
```

If the GCC version is too old, fix it with:

```sh
export CC=$(which gcc)
export CXX=$(which g++)
```

## 6. Configure the `config.ini` File

```sh
nano fluxpype/config.ini
```

Set the code to do your desired run using the documentation (TBD)

Run the configuration run script:

```sh
python fluxpype/config_runner.py
```
- For the paths to work out, the config runner must be called from one directory above the file.


## FAQs & Troubleshooting
- If gnuplot is not detected correctly in PDL, try the following, and if it works, put this line into either the conda activation script or the bashrc file:
    - ```export GNUPLOT_BINARY=$(realpath $(which gnuplot))```

- If GSL is not found, use
    - ```conda install gsl```
    - ```cpanm --notest PDL::GSL```

- If Flux.pm isn't found at the end of the flux install, try:

```sh
export PERL5LIB=$(dirname $(find $CONDA_PREFIX -name Flux.pm)):$PERL5LIB
make everything
```

- To manually set the prefixes for a conda installation before using the makefile:

```sh
export PL_PREFIX="$CONDA_PREFIX/lib/perl5/site_perl"
export FL_PREFIX="$CONDA_PREFIX"
```
- To confirm installation success:

```sh
perl -e "use Flux;"
```