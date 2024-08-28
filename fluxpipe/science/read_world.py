# A quick function to read output FLUX world data

# Import libraries
import numpy

def rdworld(filename):

    # Some quick testing on storing this data in a class object
    class world(object):

        def __init__(self):
            self.fc = self.fc()
            self.fx = self.fx()

        # Flux concentrations
        class fc:
            def __init__(self):
                self.id = []
                self.x = []
                self.y = []
                self.z = []
                self.fl = []

        def add_fc(self, id, x, y, z, fl):
            self.fc.id.append(id)
            self.fc.x.append(x)
            self.fc.y.append(y)
            self.fc.z.append(z)
            self.fc.fl.append(fl)

        # Fluxons
        class fx:
            def __init__(self):
                self.id = []
                self.fc0 = []
                self.fc1 = []
                self.fl = []
                self.x = []
                self.y = []
                self.z = []

        def add_fx(self, id, fc0, fc1, fl, r0, r1, x, y, z):
            self.fx.id.append(id)
            self.fx.fc0.append(fc0)
            self.fx.fc1.append(fc1)
            self.fx.fl.append(fl)
            self.fx.x.append(numpy.array([r0[0]] + x + [r1[0]]))
            self.fx.y.append(numpy.array([r0[1]] + y + [r1[1]]))
            self.fx.z.append(numpy.array([r0[2]] + z + [r1[2]]))

    # Let there be light
    w = world()

    # Now actually go about reading the data
    df = open(filename, 'r')

    # Initialize data storage
    ln_id = []
    ln_fc0 = []
    ln_fc1 = []
    ln_fl = []
    ln_r0 = []
    ln_r1 = []

    vx_lid = []
    vx_id = []
    vx_pos = []
    vx_x = []
    vx_y = []
    vx_z = []

    for rl in df:
        sl = rl.strip().split()

        # Skip ahead for blank lines
        if len(sl) == 0:
            continue

        # Read out flux concentrations
        if sl[0] == 'NEW':
            w.add_fc(int(sl[1]), numpy.double(sl[2]), numpy.double(sl[3]), numpy.double(sl[4]), numpy.double(sl[5]))

        # Read out fluxon data
        if sl[0] == 'LINE':
            if int(sl[1])>0:
                ln_id.append(int(sl[1]))
                ln_fc0.append(int(sl[4]))
                ln_fc1.append(int(sl[5]))
                ln_fl.append(numpy.double(sl[6]))
                ln_r0.append([numpy.double(sl[7]), numpy.double(sl[8]), numpy.double(sl[9])])
                ln_r1.append([numpy.double(sl[10]), numpy.double(sl[11]), numpy.double(sl[12])])

        # Read out vertex points
        if sl[0] == 'VERTEX':
            vx_lid.append(int(sl[1]))
            vx_id.append(int(sl[2]))
            vx_pos.append(int(sl[3]))
            vx_x.append(numpy.double(sl[4]))
            vx_y.append(numpy.double(sl[5]))
            vx_z.append(numpy.double(sl[6]))

        # Exit on neighbor information
        if 'VNEIGHBOR' in sl[0]:
            break

    df.close()

    # Convert lists into numpy arrays
    ln_id  = numpy.array(ln_id)
    ln_fc0 = numpy.array(ln_fc0)
    ln_fc1 = numpy.array(ln_fc1)
    ln_fl = numpy.array(ln_fl)
    ln_r0 = numpy.array(ln_r0)
    ln_r1 = numpy.array(ln_r1)

    vx_lid = numpy.array(vx_lid)
    vx_id = numpy.array(vx_id)
    vx_pos = numpy.array(vx_pos)
    vx_x = numpy.array(vx_x)
    vx_y = numpy.array(vx_y)
    vx_z = numpy.array(vx_z)

    # Parse the line and vertex lists and create fluxon objects
    for lid in ln_id:
        wl = numpy.where(ln_id == lid)[0]
        wv = numpy.where(vx_lid == lid)[0]

        w.add_fx(lid, ln_fc0[wl][0], ln_fc1[wl][0], ln_fl[wl][0], ln_r0[wl,:][0].tolist(), ln_r1[wl,:][0].tolist(), vx_x[wv].tolist(), vx_y[wv].tolist(), vx_z[wv].tolist())

    return w

# A quick comment on reading FLUX output world files:

# Flux concentrations
# NEW	100	0.314000	-0.059000	-0.948000	-1	(null)	0
# FC    ID      X               Y               Z               Flux    BND PT  Rad

# Line
# LINE	101	-302	-301	-1	100	1            7.61951 -2.51477 -7.52456   0.314 -0.059 -0.948
# Line  ID      Start   End     FC0     FC1     Flux        StX     StY     StZ     EndX    EndY    EndZ

# Vertex
# VERTEX	101	10650	1	7.203456	-2.377456	-7.113698
# Vertex        LineID  VtxID   Pos     X               Y               Z
