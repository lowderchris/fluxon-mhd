# Field Line Universal relaXer

The Field Line Universal relaXer (**FLUX**) is a fluxon-based model that follows magnetic field evolution in a conductive atmosphere, with exactly prescribed field topology. It has been used to solve plasma problems that are quasi-static and involve low &beta;. As of the 1-May-2009 release, FLUX can also handle the &beta;&approx;1 case and plasma dynamics, though these features are currently considered experimental. FLUX can be used to study time-dependent 3-D field evolution with a desktop workstation or laptop computer, and is intended for studying CME onset and related eruptive events in the solar corona.

FLUX is [free software](http://en.wikipedia.org/wiki/Free_software), available under the [GPLv2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt) license. It depends only on other free software, notably [GCC](http://gcc.gnu.org), [Perl](http://www.perl.org) and [PDL](http://pdl.perl.org). Rendering uses [Gnuplot](http://www.gnuplot.info) and [PDL::Graphics::Gnuplot](http://search.cpan.org/dist/PDL-Graphics-Gnuplot/).

What is FLUX?
=============

FLUX is a prototype physics simulation code intended for studying MHD evolution with slow evolution interrupted by periodic disruptions. It has been used extensively as a time-dependent nonlinear force free field solver with fixed topology and limited reconnection loci. Recent additions to the code make it useful for solving cases where the plasma pressure and dynamics are nontrivial. The v2.3 release contains plasma pressure and dynamic forces, allowing study of events such as CME onset that require inertia. The code is not (yet) suitable for studying wave phenomena.

FLUX is first-stage parallelized: it can make use of multiple CPUs that share a single workstation. This allows timely simulation of systems with up to a few times 10<sup>3</sup> fluxons (a few times 10<sup>4</sup> vertices). FLUX may also be used on a laptop or desktop workstation. Typical idealized plasma systems require only a few MB to represent, and small-but-nontrivial simulations can be constructed that fit entirely with the L2 cache of a modern 64-bit workstation.

FLUX uses the "fluxon" approach to MHD, modeling the field as a collection of curvilinear 1-D manifolds that approximate field lines. Fluxons differ from field lines primarily in that each fluxon carries a finite amount of magnetic flux, while a field line carries only an infinitesimal amount of flux. FLUX stores each fluxon as an ordered list of vertices, each of which contains spatial coordinates and ancillary information of use to the simulation engine.

Unlike more conventional codes, FLUX preserves magnetic topology exactly: unless specifically triggered by the modeling code, no reconnection at all occurs in the modeled magnetic flux system. FLUX is a coordinate-free code: each vertex is free to wander throughout all space, unless you define specific boundary conditions that limit the size of the simulation.

How it works
------------

FLUX is a force-balance relaxation engine, with additional pseudoforces -- it relaxes to solve the equation "F - ma = 0". At each relaxation step, the code calculates several vector forces (and, if desired, the inertial pseudoforce) acting at each vertex, and then moves the vertex in the direction of the forces. After a sufficient number of small steps, the system reaches equilibrium (if one exists). Typically, three forces are used: the magnetic tension ("curvature") and magnetic pressure ("field gradient") components of the Lorenz force, and a pseudoforce that keeps the vertices placed optimally along the length of each fluxon. Additional pseudoforce terms can be added to track inertial effects; these forces include a damping term to prevent interaction of the grid speed and wave speed.

The magnetic pressure and curvature forces are calculated using an analytic geometry technique called "Voronoi analysis": the code derives the locus that is closer to a particular fluxon segment ("flux element" or "fluxel") than to any other, and calculates the magnetic field by considering that all the flux associated with the fluxel is compressed into that locus. The cross-sectional area of the locus gives an estimate of the local magnetic field strength, and the asymmetry of the locus gives an estimate of the magnetic field gradient. The fluxons, in a sense, discretize the topology of the field: topological changes require reconnection of a pair of fluxons, either by exchange or by passing through one another. For time-dependent modeling, each time step consists of an update to the boundary condition, followed by a relaxation of the simulation field to match the new boundaries. Thus we distinguish between "relaxation steps" which take no physical time, and "time steps". A typical time step might require a thousand relaxation steps.

How it performs
---------------

Because FLUX grid size need only be large enough to resolve the structures being studied, grids are much smaller than in conventional 3-D magnetic field codes. Typical ideal ("toy") problems use only a few thousand fluxels in less than 100 fluxons. Problems of this size require minutes to relax on large multicore workstations. Larger problems like semi-empirical studies of solar active regions might require up to 10<sup>4</sup> fluxons, and require hours to relax. Most solar magnetic field configurations are scale invariant: the total amount of curvature of each closed field line is roughly constant against total length of the field line. This yields a major scaling advantage for FLUX, which seems to scale more like a conventional 2-D code than 3-D code: the number of vertices in the simulation appears to be proportional to the complexity of the boundary for a large range of systems.

In a direct comparison with a more conventional MHD model (ARMS) in early 2008, FLUX obtained similar answers in its appropriate computational domain, but used under 1% of the computing power as ARMS.

System requirements
-------------------

### Software environment

FLUX needs:

-   a POSIX system under which to run. It has been tested under Fedora Linux, Ubuntu Linux, and MacOS X. It is thought to work under Cygwin (but why bother?)

-   Perl 5.8.x or 5.10.x (standard on most POSIX systems these days)

-   PDL 2.4.x (standard on most Linuxes and available through Fink; source is available from [<http://pdl.perl.org>](http://pdl.perl.org)

-   gcc & make (standard essentially everywhere). (Other C compilers may work; nobody has yet tried.)

-   NetPBM (standard on most Linuxes and Macs; required for rendering; used by the PDL::IO::Pic module).

-   OpenGL (required for 3-D rendering; used by the PDL::Graphics::TriD module)

### Architecture

FLUX is thought to work on all major architectures. It has been compiled on and for 64-bit Mac ARM chips, 32- and 64-bit Intel chips, 64-bit AMD Athlon architecture, and G4/G5 architectures. FLUX is first-stage parallelized using direct POSIX system calls (fork() and pipe()).

Release Notes for the current release
-------------------------------------

TODO - update / organize / link this section

v2.2 (released 22-Nov-2008) includes many new features, both bug fix and performance boost:

-   **Mass** is actively tracked in the data structures
-   **Binary save files** are supported
-   **First-stage parallelization** of the relaxation process seems functional. ( the code can spawn on the local machine, making use of multiple processor cores that share a memory bus ).
-   **Interpolation** of all scientific values onto arbitrary locations is now possible; this should enable data export to other codes.
-   Several new reconnection criteria make reconnection "in the wrong direction" much less common
-   Multiple bug fixes in the perl interface code
-   Data structure labels are now editable via the perl hash interface
-   Better support for motion and interaction of flux concentrations

v2.1 includes several minor bug fixes, improvements to the rendering code, and Perl convenience functions for the geometry library.

v2.0 was the first revision for which all the basic field manipulation and Perl interface infrastructure seems to be 100% present. Notable improvements from v1.2 are:

-   The "under-the-hood" Perl interface is much more robust and no longer leaks memory or dumps core on Perl quit.
-   All tied-hash fields are implemented (so, e.g., you can explore the FLUX tree structures from within your Perl script).
-   Magnetic reconnection is fully tested and robust.
-   Grid de-refinement is fully tested and robust.
-   Open boundaries are fully supported.
-   U-loop and omega-loop traversals of the open boundary are supported.
-   Plasmoids are fully supported, including reconnection that involves plasmoids.
-   Flux emergence and submergence are fully supported.
-   New force laws provide better convergence and higher accuracy
-   Some improvements to the rendering code
-   Addition of footpoint-motion tools to the main tree (the "duction" infrastructure)


Theory
======

The analogy between [field lines](http://en.wikipedia.org/wiki/Field_line) and [vector field values](http://en.wikipedia.org/wiki/Vector_field) is pretty well known. It's generally used for visualization of the magnetic field. One develops a vector field, either analytically or numerically, and shoots a collection of field lines through the calculated field. Since it's hard to render a 3-D vector field directly, one then renders the field lines instead.

Fluxon codes use the field line analogy in reverse, treating the field as a collection of discretized field lines. They are discretized in two senses: (A) each line (a ["fluxon"](http://en.wikipedia.org/wiki/Fluxon)) carries a finite, rather than infinitesimal, amount of [magnetic flux](http://en.wikipedia.org/wiki/Magnetic_flux); and (B) the field lines are piecewise linear, composed of vertices and line segments, rather than smooth curves. Approximations to the magnetic pressure and curvature forces can be calculated from the local geometry around each vertex and/or line segment.

FLUX uses 2-D analytic geometry to construct the [Voronoi cell](http://en.wikipedia.org/wiki/Voronoi_cell) around each line segment - the locus (in the plane perpendicular to the fluxon line segment) that is closer to the current line segment than to any other. The Voronoi cell is calculated by a simple [hull algorithm](User_Manual:_hull_algorithm) acting in the plane perpendicular to each fluxel in turn. The local field strength (hence magnetic pressure) and gradient may be inferred directly from the geometry of this Voronoi cell. The direction is determined from the direction of the line segment itself.

FLUX is a coordinate-free code: while the internal representation of space is Cartesian, the physical simulation interacts with that representation only through vector operations such as inner and cross product. Hence there are no inherent boundaries to the location, and no preferred directions: the code works equally well in spherical geometries as planar ones. In that respect, fluxon codes in general are like smooth-particle hydrodynamic codes: they represent a continuous variable (the vector magnetic field for FLUX, or the fluid density field for a smooth-particle code) as a collection of discrete objects (fluxons vs. particles) in an unbounded space.

For more detail, read [DeForest & Kankelborg 2007](http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2007JATP...69..116D&link_type=PREPRINT&db_key=PHY).

Why bother?
-----------

Fluxons have the advantage that they discretize the topology of the magnetic field. Because the differential topology of a true field line map is represented using the discrete topology of individual field lines, it is not possible for fluxons to reconnect "accidentally", and numerical reconnection is avoided. Other advantages include possible scaling and stability benefits compared to conventional codes. FLUX is believed to be the only code capable of running useful 3-D (quasi-)time-dependent MHD simulations on a desktop workstation.

FLUX is currently the only MHD simulation code that can maintain and evolve small open-field regions embedded in the closed "quiet" solar corona. These regions are the leading candidate source region for the slow solar wind that permeates the Solar System near the plane of the ecliptic. It is also the only MHD code that can evolve and maintain indefinitely the large metastable systems of magnetic flux that are thought to cause coronal mass ejections.

Publications
============

Papers
------

-   Rachmeler, DeForest, Pariat, Antiochos, & Kankelborg 2009: demonstration that reconnection \*is\* required for certain symmetric systems to erupt.

<!-- -->

-   [Rachmeler, DeForest & Kankelborg, APJ 2009](http://data.boulder.swri.edu/derek/FLUX_files/media/Rachmeler2009.pdf): use of the fluxon code to demonstrate an ideal MHD instability that can in principle cause CME eruptions even without reconnection.

<!-- -->

-   [DeForest & Kankelborg, JASTP 2007](http://data.boulder.swri.edu/derek/FLUX_files/media/Flux.pdf): Description paper of the Voronoi method for estimating magnetic field forces.

Presentations
-------------

There have been many. Mediawiki doesn't like to deal with huge files, so it doesn't make good sense to upload whole presentations here into the wiki -- but I'd like to put them on the web somewhere. This is a partial catalog starting in late 2006:

### 2009

-   Rachmeler et al. 2009: dissertation talk at AAS/SPD (Boulder): description of recent transition to global solar modeling

### 2008

-   [[Rachmeler, DeForest, Pariat, & Antiochos|SPD_08_poster]], poster at AAS/SPD (Ft. Lauderdale): overview of progress, demo of the recent comparison between FLUX and ARMS. Winner: best student poster.

<!-- -->

-   DeForest et al. 2008: poster at the Evershed centenary meeting (Bangalore, India)

### 2007

-   DeForest & Rachmeler, talk at AAS/SPD (Honolulu): overview of progress, and demo of reconnection
-   Rachmeler & DeForest: poster at AAS/SPD (Honolulu): herniation result. Runner-up for best student poster.
-   Bruder: Final report, 2007 MSU REU program (Bozeman): [archived at MSU](http://solar.physics.montana.edu/reu/2007/dbruer/presentation/hello.html). Overview of debugging/validation efforts.
-   DeForest & Rachmeler: two (!) invited talks at SHINE (Whistler, BC): overview of fluxon modeling and results; presentation of the herniation result.

### 2006

-   Rachmeler & DeForest: poster at AGU (San Francisco), describing the herniation result. Won an award for best student poster.

Useful movies & plots
---------------------

TODO - update / organize / link this section

Current projects
================

TODO - update / organize / link this section

User's Manual
=============

TODO - update / organize / link this section

It is largely a programmer's reference, containing descriptions of the routines and calling conventions in the libraries and code, and some of the design ideas that went into them. It is being expanded to include more examples on how to use the code.

Beefs/bugs
==========

Please submit any issues to the GitHub Issues tracker at the top of this page.

Author(s)
=========

Fluxons are a joint invention of [Craig DeForest](http://www.boulder.swri.edu/~deforest) and [Charles Kankelborg](http://solar.physics.montana.edu/kankel). The code contains contributions from many individuals, including:

-   Craig DeForest
-   Charles Kankelborg
-   Laurel Rachmeler
-   Alisdair Davey
-   Daniel Bruder
-   Teddy Walls
-   Derek Lamb
-   Chris Lowder

FLUX is currently (2017) maintained by Chris Lowder, Derek Lamb, and Craig DeForest.

Acknowledgements
================

The initial development and a 2015--2017 revival were funded internally by the [Southwest Research Institute](http://www.swri.org).  Significant additional funding has also come from NASA's Science Mission Directorate, Heliophysics Division, through the [Living With a Star](http://lws.gsfc.nasa.gov) and Supporting Research programs.