import sys
import h5py 
import time
import numpy as np 
from argparse import ArgumentParser
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

def select_ids(fname, id_file=None):
    if id_file is None:
        return TracerData_t(fname).ids()
    else:
        return np.loadtxt(id_file)

def append_timeseries(fname, ids, tseries):
    tracer_data = TracerData_t(fname).stack()
    for t in tracer_data:
        ID = int(t[0])
        if ID not in ids:
            continue
        ix, = np.where(ids == ID)[0]
        tseries[ix].append(t)

def write_tracer_timeseries(group, tseries):
    group.create_dataset('x', data=tseries.x)
    group.create_dataset('y', data=tseries.y)
    group.create_dataset('vx', data=tseries.vx)
    group.create_dataset('vy', data=tseries.vy)
    group.create_dataset('density', data=tseries.density)
    group.create_dataset('pressure', data=tseries.pressure)

if __name__ == '__main__':

    # files   = sys.argv[1:]
    parser = ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('--id-file', default=None)
    args = parser.parse_args()

    # ===============================================================
    ids     = select_ids(args.files[0], id_file=args.id_file)
    tseries = np.array([TracerTimeseries(i) for i in ids])
    ts      = []

    # ===============================================================
    for f in args.files:
        start = time.time()
        ts.append(h5py.File(f, 'r')['time'][...])
        append_timeseries(f, ids, tseries)
        end = time.time()
        print(f, "[{:.2f} s]".format(end - start))

    # ===============================================================
    fname = 'tracer_tseries.h5'
    h5f = h5py.File(fname, 'w')
    h5f.create_dataset('time', data=ts)
    print("   Saving {}...".format(fname))

    # ===============================================================
    for t in tseries:
        g = h5f.create_group("{:d}".format(int(t.id)))
        write_tracer_timeseries(g, t)

    # ===============================================================
    h5f.close()
