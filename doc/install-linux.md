# Linux Installation

To install the fluxon-mhd modeling framework, there are three main requirements - a perl installation with libraries, the perl data language toolset, and finally compiling the fluxon-mhd software itself.

Here we'll outline steps for installing on a Linux system - notably tested on an Ubuntu system. The broad steps can be applied for other flavors of Linux.

## Perl installation

First, use apt-get to install cpanminus, an easier tool to use for installing perl packages.
```shell
sudo apt-get install cpanminus
```

Install PDL using cpanminus, forcing to bypass any dependency warnings.
```shell
sudo cpanm --force PDL
```

If problems are encountered with that PDL installation route, one can install an older version of PDL using apt-get. This may run into code compilation issues down the road.
```shell
sudo apt-get install pdl
```

Write the following to ~/.perldlrc, making sure to reference the appropriate location local perl libraries will be installed to.
```shell
push(@INC, "/usr/local/lib/aarch64-linux-gnu/perl/5.38.2");

require('PDL/default.perldlrc');

use PDL::AutoLoader;
$PDL::AutoLoader::Rescan=1;

1;

```

A few other packages of note will be needed for the later fluxon-mhd install.
```shell
sudo cpanm File::ShareDir
sudo cpanm Math::RungeKutta
sudo cpanm PDL::GSL::INTERP
sudo cpanm Term::ReadKey
```

Optionally, if you need graphics, install gnuplot
```shell
sudo apt-get install gnuplot
sudo cpanm PDL::Graphics::Gnuplot
```

## FLUX installation

To build some of the C code, make sure you have a compiler installed.
```shell
sudo apt-get install gcc
```

If you haven't already, clone the fluxon-mhd repository and then make everything from the base directory.
```shell
git clone https://github.com/lowderchris/fluxon-mhd.git
cd fluxon-mhd

make everything
```
