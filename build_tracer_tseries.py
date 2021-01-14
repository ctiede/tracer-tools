import sys
import h5py 
import numpy as np 
from tracers import TracerData_t

class TracerTimeseries:
    def __init__(self, ID):
        self.id = ID
        self.x  = []
        self.y  = []
        self.vx = []
        self.vy = []
        self.density  = []
        self.pressure = []

    def append(self, new):
        self.x.append(       new[1])
        self.y.append(       new[2])
        self.vx.append(      new[3])
        self.vy.append(      new[4])
        self.density.append( new[5])
        self.pressure.append(new[6])



def unpack_tracer_data(tracer):
    data = 0
    velo = 1
    return (tracer[data][0], # id       (0)
            tracer[data][1], # x        (1)
            tracer[data][2], # y        (2)
            tracer[velo][0], # vx       (3)
            tracer[velo][1], # vy       (4)
            tracer[2],       # density  (5)
            tracer[3])       # pressure (6)

def append_timeseries(fname, tseries):
    h5f     = h5py.File(fname, 'r')
    unpack  = np.vectorize(lambda t: unpack_tracer_data(t))
    tracer_data = np.column_stack(unpack(h5f['tracers'][...]))
    for t in tracer_data:
        ID = int(t[0])
        tseries[ID].append(t)

def write_tracer_timeseries(group, tseries):
    group.create_dataset('x', data=tseries.x)
    group.create_dataset('y', data=tseries.y)
    group.create_dataset('vx', data=tseries.vx)
    group.create_dataset('vy', data=tseries.vy)
    group.create_dataset('density', data=tseries.density)
    group.create_dataset('pressure', data=tseries.pressure)

if __name__ == '__main__':

    files    = sys.argv[1:]
    ntracers = np.max(TracerData_t(files[0]).ids()) + 1
    tseries  = [TracerTimeseries(i) for i in range(int(ntracers))]
    time     = []

    for f in files:
        print(f)
        time.append(h5py.File(f, 'r')['time'][...])
        append_timeseries(f, tseries)

    fname = 'tracer_tseries.h5'
    h5f = h5py.File(fname, 'w')
    h5f.create_dataset('time', data=time)
    print("   Saving {}...".format(fname))

    for t in tseries:
        g = h5f.create_group("{}".format(t.id))
        write_tracer_timeseries(g, t)

    h5f.close()
