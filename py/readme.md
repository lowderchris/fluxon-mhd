# FLUXON Python Tools

## Notes
Just a few quick notes about the Python tools being developed here for FLUXON. This is a fairly rough migration of odds-and-ends scripts to a more centralized location, and things will be tidied up during development.

## Scripts
- rdworld.py - Read a .flux world output file into a Python class object. This allows for easy manipulation of fluxon data, and for visualization using...

- vizworld.py - Render an interactive 3D visualization of the fluxon world. For now, this renders a pseudo-photospheric magnetic field, and interpolates the solar wind speed distribution at the defined outer boundary.

## Setup
A few variables will need to be defined within vizworld.py, depending on installation route and desired parameters.

## Usage... for now

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

For the visualization routines, the package MayaVi is required, which provides a more user-friendly interface for VTK. Installation has in the past proven a tangled web of dependencies and abandoned packages, but as of recent versions installation should be possible via:

    brew install vtk --with-python --with-qt --without-python@2
    pip3 install (--user) PyQt5
    pip3 install (--user) mayavi

Note that as of late, the installation / compilation of new(er) versions of VTK appears to be problematic, with a fix requiring installing from the most updated source. This causes problems down the road with installing mayavi, as it is unable to see a specific release version number of VTK.

Good luck.
