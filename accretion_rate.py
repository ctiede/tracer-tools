import sys
import h5py 
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from argparse import ArgumentParser
from tracers import TracerData_t


def split_tracers(file, rcut=1.0):
    data  = TracerData_t(file)
    radii = data.radii()
    return (data.time(), len(radii[radii < rcut]), len(radii[radii >= rcut]))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    if len(args.files) < 2:
        print("Require more files to calcualte tracer accretion")

    time   = []
    n_accr = []
    n_disk = []
    for f in args.files:
        (t, na, nd) = split_tracers(f)
        n_accr.append(na)
        n_disk.append(nd)
        time.append(t)

    np.save('tracer_counts.npy', np.column_stack([time, n_accr, n_disk]))
    # n_accr = np.array(n_accr)
    # n_disk = np.array(n_disk)
    # time   = np.array(time)

    # dt    = np.diff(time)
    # dorbs = dt / (2 * np.pi)
    # dn_accrete = np.diff(n_accr)
    # dn_percent = dn_accr / n_disk[:-1]