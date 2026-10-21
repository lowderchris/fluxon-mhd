# FLUX Python Tools

Here find some assembled Python tools for interfacing with and interpreting results from the FLUX core PDL code.

## Plotting tools
- plt-world.py - Plots a single snapshot of a fluxon world
- plt-relax.py - Plots a sequence of relaxation steps for a fluxon world


## Scripts
- rdworld.py - Read a .flux world output file into a Python class object. This allows for easy manipulation of fluxon data, and for visualization using...
- vizworld.py - Render an interactive 3D visualization of the fluxon world. For now, this renders a pseudo-photospheric magnetic field, and interpolates the solar wind speed distribution at the defined outer boundary.

## Usage

For a more friendly working environment, iPython and the pylab modules can be loaded,

    ipython --pylab

A saved fluxon world can be restored to life using,

    from rdworld import rdworld
    world = rdworld('filename-goes-here.flux')

This resulting world class contains some (for now) of the information contained within the original PDL .flux world file.

Flux concentrations:

    world.fc.x(,y,z)    - Cartesian position
    world.fc.fl         - Magnetic flux
    world.fc.id         - Flux concentration ID

Fluxons:

    world.fx.x(,y,z)    - Cartesian positions of vertex points
    world.fx.fc0(,fc1)  - Associated flux concentration IDs
    world.fx.fl         - Magnetic flux
    world.fx.id         - Fluxon ID

More details (connectivity, etc) can be siphoned from the original .flux world file as required in future.

## Package dependencies

Python 3, NumPy, Scipy are the standard set that these routines rely upon. iPython is optional, but provides a friendly environment in which to work.

For the visualization routines, the package MayaVi is required, which provides a more user-friendly interface for VTK. The pip3 version of this package weaves a tangled web of dependencies, but at the time of writing this, the bleeding-edge Github version of [MayaVi](https://github.com/enthought/mayavi) version is working.

See installation notes for a full description of getting an environment set up.
