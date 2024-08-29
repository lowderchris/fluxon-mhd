# macOS Installation

To install the fluxon-mhd modeling framework, there are three main requirements - a perl installation with libraries, the perl data language toolset, and finally compiling the fluxon-mhd software itself. Some of the background perl installation requires root / admin access, but beyond that packages can be installed locally for ease of editing and recompilation.

## Perl installation

First, install the homebrew package, which will require sudo access.
```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Next, use homebrew as an admin user to install perl.
```shell
brew install perl
```

If you plan on using homebrew for other packages, and would like to keep perl from updating alongside other packages, the installation can be frozen.
```shell
brew pin perl
```

Two more homebrew packages will be needed later, which can also be installed via an admin user.
```shell
brew install gnuplot
brew install fftw
```

The final suggested step as an admin user is to setup the cpanminus tool for installing perl packages.
```shell
curl -L https://cpanmin.us | perl - --sudo App::cpanminus
```

Set perl to install packages to a local library directory. Here we've chosen to use a perl5 directory inside the ~/Library folder, but this could be anywhere your user has permissions. This is then stored to automatically update in a .zprofile startup script for Zsh - modify if using a different shell.
```shell
PERL_MM_OPT="INSTALL_BASE=~/Library/perl5" cpanm --local-lib=~/Library/perl5 local::lib
echo 'eval `perl -I ~/Library/perl5/lib/perl5 -Mlocal::lib=~/Library/perl5`' >> ~/.zprofile
```

Write the following to ~/.perldlrc, making sure to reference the appropriate location local perl libraries will be installed to.
```shell
push(@INC, "/Users/username/Library/perl5/lib/perl5/5.38.2");

require('PDL/default.perldlrc');

use PDL::AutoLoader;
$PDL::AutoLoader::Rescan=1;

1;

```

## Perl Data Language (PDL) installation

Using cpanm, install PDL and required depedencies.
```shell
cpanm PDL
```

A few other packages of note will be needed for the later fluxon-mhd install.
```shell
cpanm File::ShareDir
cpanm PDL::Graphics::Gnuplot
cpanm Math::RungeKutta
cpanm Term::ReadKey
```

## FLUX installation

Finally, for the first install, and for future recompilations, set the perl and flux prefix paths.
```shell
export PL_PREFIX='/Users/username/Library/perl5'
export FL_PREFIX='/Users/username/Library/flux'
```

If you haven't already, clone the fluxon-mhd repository and then make everything from the base directory.
```shell
make everything
```
