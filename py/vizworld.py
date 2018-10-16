# A quick script to generate some fluxon plots in Mayavi

# Define any plotting options
wsize = [1024, 1024]        # window size
bgcol = (1.0, 1.0, 1.0)     # background color rgb tuple
fgcol = (0.0, 0.0, 0.0)     # foreground color rgb tuple

# Define input file information
ffile = 'chtest-rlx3000.flux'    # flux world file
wfile = 'chtest-rlx3000.dat'     # fluxon wind file

# Define the photospheric and outer boundary radii
r0 = 1.0
r1 = 20.0

# Define any wind parameters for plotting
wv0 = 0.0
wv1 = 1000.0

# Import any libraries
from mayavi import mlab
import numpy as np
import scipy.interpolate

# Define a function to interpolate from a series of points onto a spherical map
def sphint(ocrds, ovals, iph, ith):
    # Boundary at 0 longitude
    mcrds = np.copy(ocrds)
    mcrds[:,0] = mcrds[:,0] - 2*np.pi
    icrds = np.append(ocrds, mcrds, 0)
    ivals = np.append(ovals, ovals)
    # Boundary at 360 longitude
    mcrds = np.copy(ocrds)
    mcrds[:,0] = mcrds[:,0] + 2*np.pi
    icrds = np.append(icrds, mcrds, 0)
    ivals = np.append(ivals, ovals)
    # Boundary at 0 latitude
    mcrds = np.copy(ocrds)
    mcrds[:,0] = mcrds[:,0] + np.pi
    mcrds[where(mcrds[:,0] > 2*np.pi),0] = mcrds[where(mcrds[:,0] > 2*np.pi),0] - 2*np.pi
    mcrds[:,1] = ((-1) * mcrds[:,1])
    icrds = np.append(icrds, mcrds, 0)
    ivals = np.append(ivals, ovals)
    # Boundary at 180 latitude
    mcrds = np.copy(ocrds)
    mcrds[:,0] = mcrds[:,0] + np.pi
    mcrds[where(mcrds[:,0] > 2*np.pi),0] = mcrds[where(mcrds[:,0] > 2*np.pi),0] - 2*np.pi
    mcrds[:,1] = ((-1) * mcrds[:,1]) + 2*np.pi
    icrds = np.append(icrds, mcrds, 0)
    ivals = np.append(ivals, ovals)
    igph, igth = meshgrid(iph, ith)
    iws = scipy.interpolate.griddata(icrds, ivals, (igph, igth), method='cubic')

    return iws

# Read the world
from rdworld import rdworld
world = rdworld(ffile)

# Create a fake photosphere for testing
#br = np.zeros((180,360),dtype=np.double)
#lat = np.linspace(0,np.pi,180)
#lon = np.linspace(0,2*np.pi,360)
#llon, llat = np.meshgrid(lon, lat)
#br += (-1) * np.cos(llat)**7

# Generate a representative photospheric Br map by interpolating flux concentration points and a scattering of zero points
# There's probably a better way to do this, but this is just for visualization for now
fcx = np.array(world.fc.x)
fcy = np.array(world.fc.y)
fcz = np.array(world.fc.z)
fcth = np.arccos(fcz / np.sqrt(fcx**2 + fcy**2 + fcz**2))
fcph = np.arctan2(fcy, fcx)
del fcx, fcy, fcz
fcvals = np.array(world.fc.fl)
iph = np.linspace(0, 2*np.pi, num=360)
ith = np.linspace(0, np.pi, num=180)
fcph = np.append(fcph, np.random.rand(1000) * 2 * np.pi - np.pi)
fcth = np.append(fcth, np.random.rand(1000) * (110*np.pi/180) + (35*np.pi/180))
fcvals = np.append(fcvals, np.zeros(1000))
br = sphint(np.stack((fcph, fcth), axis=1), fcvals, iph, ith)

# Alternatively, generate a representative photospheric Br map from flux concentrations and synthetic distributions of flux
#fcx = np.array(world.fc.x)
#fcy = np.array(world.fc.y)
#fcz = np.array(world.fc.z)
#fcth = np.arccos(fcz / np.sqrt(fcx**2 + fcy**2 + fcz**2))
#fcph = np.arctan2(fcy, fcx)
#del fcx, fcy, fcz
#fcvals = np.array(world.fc.fl)
#iph = np.linspace(0, 2*np.pi, num=360)
#ith = np.linspace(0, np.pi, num=180)
#br = numpy.zeros((180,360))
# ...

# Create a blank canvas
fig = mlab.figure(1, bgcolor = bgcol, fgcolor = fgcol, size=wsize)

# Create a sun
r = r0
phi, theta = np.meshgrid(iph, ith)
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
z = r * np.cos(theta)

sol = mlab.mesh(x, y, z, scalars=np.clip(br,-1,1), colormap='Greys')
sol.module_manager.scalar_lut_manager.reverse_lut = True

# Plot the fluxon fieldlines
#for i in np.unique(world.fx.id)[0::5]:
#wdat = np.loadtxt(wfile)
for i in np.unique(world.fx.id):
    flw = world.fx.id.index(i)
    # Alternatively, plot each fieldline with the corresponding vrend value
    #svel = wdat[where(wdat[:,7] == i)[0],0]
    #if (len(svel) == 0): svel = 0
    slen = np.sqrt((world.fx.x[flw])**2 + (world.fx.y[flw])**2 + (world.fx.z[flw])**2)
    tfl = mlab.plot3d(world.fx.x[flw], world.fx.y[flw], world.fx.z[flw], slen, tube_radius=None, colormap='YlGnBu', vmin=r0, vmax=r1, reset_zoom=False)
    #tfl = mlab.plot3d(world.fx.x[flw], world.fx.y[flw], world.fx.z[flw], np.repeat(svel, len(world.fx.x[flw])), tube_radius=None, colormap='inferno', vmin=0.0, vmax=140, reset_zoom=False)
    tfl.module_manager.scalar_lut_manager.reverse_lut = True

# Map interpolated solar wind values at an outer boundary
wdat = np.loadtxt(wfile)
gdr = where((wdat[:,6]>2.5) & (wdat[:,0]<5000))[0]
ocrds = np.rollaxis(np.array([wdat[gdr,5], wdat[gdr,3]]), 1)
ovals = wdat[gdr,0]

iph = np.linspace(0, 2*np.pi, num=360)
ith = np.linspace(0, np.pi, num=180)
iws = sphint(ocrds, ovals, iph, ith)

r = r1
phi, theta = np.meshgrid(iph, ith)
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
z = r * np.cos(theta)

ssol = mlab.mesh(x, y, z, scalars=np.clip(iws,wv0,wv1), colormap='inferno', opacity=0.75)
ssol.enable_contours = True
ssol.actor.property.lighting = False
ssol.actor.property.line_width = 2.0

cb = mlab.colorbar(object=ssol, title='Vr', orientation='horizontal')
cb.data_range = array([wv0,wv1])
cb.scalar_bar_representation.position = array([0.1,0.01])
cb.scalar_bar_representation.position2 = array([0.8, 0.1])
#cb.scalar_bar.bar_ratio = 0.2
cb.scalar_bar.title_ratio = 0.6
