The Field Line Universal Relaxer (FLUX) is a prototype MHD solver that is currently able to find nonlinear force-free fields with an exactly specified topology in a quasi-static, infinitely conductive, low-beta medium. FLUX accepts as input a collection of field lines that define the topology that you want to study. The FLUX engine relaxes the position of those field lines until they approximate the physical force-free state with the same topology that you specified.

FLUX consists of two parts: a relaxation engine, written in C, and control interfaces in Python and Perl. Most user interfaces are provided in Python.

For more details see [Lowder, Gilly, and DeForest 2024](https://ui.adsabs.harvard.edu/abs/2024ApJ...965...13L/abstract) or the FLUX [wiki](https://github.com/lowderchris/fluxon-mhd/wiki).


## Installation

Instructions for installing FLUX on Linux or macOS machines can be found in the [documentation](https://github.com/lowderchris/fluxon-mhd/tree/main/doc) folder.


## Code structure

- doc -- install guide and some documentation (see the [wiki](https://github.com/lowderchris/fluxon-mhd/wiki) for more up-to-date information.)

- include -- 'c' header files (these get installed into $FL_PREFIX/include/flux/)

- lib -- the central libflux.a core simulation library (this gets installed into $FL_PREFIX/lib/)

- pdl -- PDL glue code and perl modules

  - pdl/PDL -- PDL autoload subroutines

- py -- Python visualization and analysis tools

- t -- Test suite


## Acknowledgements

Authors:
- Craig DeForest (SwRI)
- Charles Kankelborg (MSU)
- Derek Lamb (SwRI)
- Chris Lowder (SwRI)
- Chris Gilly (then at SwRI)
- Laurel Rachmeler (then at CU Boulder & SwRI)
- Alisdair Davey (then at SwRI)

Contributors (gratefully acknowledged):
- Nathan Schwadron
- Teddy Walls
- Christopher Plumberg

The authors gratefully acknowledge support from NASA's Living With A Star and Solar/Heliospheric Physics programs, the AFOSR, and from the Southwest Research Institute.
