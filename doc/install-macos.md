# macOS Installation

To install the fluxon-mhd modeling framework, there are three main requirements - a perl installation with libraries, the perl data language toolset, and finally compiling the fluxon-mhd software itself.

Homebrew and the base perl installation will be installed as a user that has full system access, with further packages and flux itself installed locally. This provides both a buffer if the system perl is updated, and also allows regular access for users who do not have regular administrative privileges. Note that if you have an admin user account, you can also install everything within the homebrew perl library directory, updating paths below accordingly.

## Perl installation

First, install the homebrew package, which will require being logged in as an admin user. Follow through the installation process.
```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Make sure at the end of this to follow the instructions for adding homebrew to your path
```shell
(echo; echo 'eval "$(/opt/homebrew/bin/brew shellenv)"') >> /Users/<USERNAME>/.zprofile
eval "$(/opt/homebrew)
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
brew install cpanminus
```

Set perl to install packages to a local library directory. Here we've chosen to use a perl5 directory inside the ~/Library folder, but this could be anywhere your user has permissions. This is then stored to automatically update in a .zprofile startup script for Zsh - modify if using a different shell.
```shell
PERL_MM_OPT="INSTALL_BASE=~/Library/perl5" cpanm --local-lib=~/Library/perl5 local::lib
echo 'eval `perl -I ~/Library/perl5/lib/perl5 -Mlocal::lib=~/Library/perl5`' >> ~/.zprofile
source ~/.zprofile
```

Write the following to ~/.perldlrc, making sure to reference the appropriate location local perl libraries will be installed to.
```shell
push(@INC, "/Users/<USERNAME>/Library/perl5/lib/perl5/5.38.2");

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
export PL_PREFIX='/Users/<USERNAME>/Library/perl5'
export FL_PREFIX='/Users/<USERNAME>/Library/flux'
```

If you haven't already, clone the fluxon-mhd repository and then make everything from the base directory.
```shell
git clone https://github.com/lowderchris/fluxon-mhd.git
cd fluxon-mhd

make everything
```
