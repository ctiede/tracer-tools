import sys
import h5py 
import time
import numpy as np 
from argparse import ArgumentParser
from tracers import TracerData_t


def get_tracer_data(fname, n):
    print(fname)
    tracer_heap  = TracerData_t(fname).stack()
    tracer_array = np.zeros((n, tracer_heap.shape[1]))
    tracer_ids   = tracer_heap[:,0].astype(int)
    tracer_array[:,0] = np.arange(n)
    tracer_array[:,5] = -1 #set density of lost tracers to -1
    tracer_array[tracer_ids] = tracer_heap
    return tracer_array


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('--out-file', default='tracer_tseries.h5')
    parser.add_argument('--tof', default=None, type=int)
    args = parser.parse_args()  

    fname = args.out_file
    h5f   = h5py.File(fname, 'w')
    n_tracers_max = int(np.max(TracerData_t(args.files[0]).ids()) + 1)
    if args.tof is not None:
        times = []
        files = [args.file[i:i+args.tof] for i in range(0, len(args.files), args.tof)]
        for i, fs in enumerate(files):
            tseries_i = np.dstack([get_tracer_data(f, n_tracers_max) for f in fs])
            ts = np.array([h5py.File(f, 'r')['time'][()] for f in fs]) 
            times.append(ts)
            g = h5f.create_group('{:}'.format(i))
            g.create_dataset('time', data=ts)
            g.create_dataset('tseries', data=tseries_i)
            g.close()
        h5f.create_dataset(np.concatenate(times))
        h5f.close()
    else:
        tracer_series = np.dstack([get_tracer_data(f, n_tracers_max) for f in args.files])
        times = np.array([h5py.File(f, 'r')['time'][()] for f in args.files])
        print("   Saving {}...".format(fname))
        h5f.create_dataset('time', data=times)
        h5f.create_dataset('tracer_tseries', data=tracer_series)
        h5f.close()


    