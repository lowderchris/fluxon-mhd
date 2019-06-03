# A quick script to generate some fluxon relaxation animations in Mayavi

# Define any plotting options
keynote = 0                 # 0 for LaTeX, 1 for Keynote
wsize = [1024, 1024]        # window size
if keynote:
    bgcol = (0.0, 0.0, 0.0)     # background color rgb tuple
    fgcol = (1.0, 1.0, 1.0)     # foreground color rgb tuple
else:
    bgcol = (1.0, 1.0, 1.0)     # background color rgb tuple
    fgcol = (0.0, 0.0, 0.0)     # foreground color rgb tuple

# Define input file information
preffile = 'dat/frtest-rlx'    # flux world file [.flux]
prewfile = 'dat/frtest-rlx'     # fluxon wind file [.dat]

# Define the photospheric and outer boundary radii
r0 = 1.0
r1 = 20.0

# Define any wind parameters for plotting
wv0 = 0.0
wv1 = 1000.0

# Import any libraries
import sys
#sys.path.insert(0, mayavilib)
from mayavi import mlab
import numpy as np
import scipy.interpolate

# Define a function to infer an idealized magnetogram from fluxon magnetic sources
def fluxon_magnetogram(oph, oth, ovals, nph, nth):
    fcwph = 2 * np.pi / 360
    fcwth = 2 * np.pi / 360
    #padth = 10
    #radth = padth * (np.pi / nth)
    padph = 10
    radph = padph * (2*np.pi / nph)
    padbr = np.zeros([nth, nph + 2*padph])
    gph = np.linspace(0-radph, 2*np.pi+radph, num=nph+2*padph)
    gth = np.linspace(0, np.pi, num=nth)
    mgph, mgth = meshgrid(gph, gth)
    #fcwph = (sin(mgth+np.pi)+1)*0.05+0.05
    #fcwth = sin(mgth)*0.05+0.01
    for i in np.arange(len(oph)):
        padbr += ovals[i] * exp(-1*( ((mgph-oph[i])**2)/(2*fcwph**2) + ((mgth-oth[i])**2)/(2*fcwth**2) ))
    br = copy(padbr[:, padph:nph+padph])
    br[:,0:padph] += padbr[:,-padph:]
    br[:,-padph:] += padbr[:,0:padph]

    br = np.roll(br, int(nph/2), axis=1)

    return br

#Quick counter
pnum = 0

# Loop through the files
for rlx in np.append(np.arange(0,100), np.arange(1,31)*100):

    ffile = preffile + str(rlx) + '.flux'
    wfile = prewfile + str(rlx) + '.dat'

    # Read the world
    from rdworld import rdworld
    world = rdworld(ffile)

    # Generate a representative photospheric Br map from fluxon magnetic concentrations
    nph = 360
    nth = 180
    fcx = np.array(world.fc.x)
    fcy = np.array(world.fc.y)
    fcz = np.array(world.fc.z)
    fcth = np.arccos(fcz / np.sqrt(fcx**2 + fcy**2 + fcz**2))
    fcph = np.arctan2(fcy, fcx) + np.pi
    del fcx, fcy, fcz
    fcvals = np.array(world.fc.fl)
    br = fluxon_magnetogram(fcph, fcth, fcvals, nph, nth)

    # Plotting coordinates
    sph = np.linspace(0, 2*np.pi, num=nph)
    sth = np.linspace(0, np.pi, num=nth)

    # Clear up some space
    if pnum != 0:
            fig.scene.disable_render = True
            mlab.clf()
            mlab.close(all=True)
    
    # Create a blank canvas
    fig = mlab.figure(1, bgcolor = bgcol, fgcolor = fgcol, size=wsize)
    fig.scene.disable_render = False
    mlab.view(0, 70, 22, (0,0,0))
    #fig.scene.light_manager.light_mode = 'vtk'

    # Create a sun
    r = r0
    phi, theta = np.meshgrid(sph, sth)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    #sol = mlab.mesh(x, y, z, scalars=np.zeros((180,360)), colormap='Greys')
    sol = mlab.mesh(x, y, z, scalars=np.clip(br,-1,1), colormap='Greys')
    sol.module_manager.scalar_lut_manager.reverse_lut = True

    # Plot the fluxon fieldlines
    #for i in np.unique(world.fx.id)[0::5]:
    #wdat = np.loadtxt(wfile)
    fig.scene.disable_render = True
    for i in np.unique(world.fx.id):
        flw = world.fx.id.index(i)
        # Alternatively, plot each fieldline with the corresponding vrend value
        #svel = wdat[where(wdat[:,7] == i)[0],0]
        #if (len(svel) == 0): svel = 0
        slen = np.sqrt((world.fx.x[flw])**2 + (world.fx.y[flw])**2 + (world.fx.z[flw])**2)
        if (slen[0]<1.05) & (slen[-1]<1.05):
            cmap = 'Blues'
            rmax = r1/4.
            rcmap = 1 if keynote else 1
        else:
            cmap = 'magma'
            rmax = r1
            rcmap = 1 if keynote else 0
        tfl = mlab.plot3d(world.fx.x[flw], world.fx.y[flw], world.fx.z[flw], slen, tube_radius=None, colormap=cmap, vmin=r0, vmax=rmax, reset_zoom=False)
        #tfl = mlab.plot3d(world.fx.x[flw], world.fx.y[flw], world.fx.z[flw], np.repeat(svel, len(world.fx.x[flw])), tube_radius=None, colormap='inferno', vmin=0.0, vmax=140, reset_zoom=False)
        if rcmap : tfl.module_manager.scalar_lut_manager.reverse_lut = True

    # Resume normal rendering
    fig.scene.disable_render = False

    # Orient the world
    mlab.view(0, 70, 22, (0,0,0))
    # mlab.view(-35,90,5, (0,0,0))

    # Stop, annotate
    mlab.text(0.01,0.01,'%05.f'%rlx,opacity=0.5,width=0.07)

    # Save the figure
    mlab.savefig('./frm/frm'+'%05.f'%pnum+'.png')

    # Advance the counter
    pnum = pnum + 1

# After all that finishes up...
# ffmpeg -pattern_type glob -i 'frm/frm*.png' -r 10 -vcodec libx264 -pix_fmt yuv420p -q 0 plt-relax.mp4

# To vingnette output images...
# convert fr-snapshot.png -background white -vignette 0x20+10+10 fr-snapshot2.png
