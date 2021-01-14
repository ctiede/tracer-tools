from argparse import ArgumentParser
import sys
import h5py 
import numpy as np
from tracers import TracerData_t


def tracers_in_outer_disk(fname):
    data = TracerData_t(fname)
    ids  = data.ids()
    r    = data.radii()
    return ids[r > 2.5]

def find_accreted_tracers(f_start, f_finish):
    id_list  = tracers_in_outer_disk(f_start)
    data     = TracerData_t(f_finish)
    ids      = data.ids()
    r        = data.radii()
    accreted = ids[r < 0.75]
    return np.intersect1d(id_list, accreted)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('start' , help='file to start search for accreted tracers')
    parser.add_argument('finish', help='file to end search for accreted tracers')
    args = parser.parse_args()

    ids = find_accreted_tracers(args.start, args.finish)

    print("Found {} accreted tracers".format(len(ids)))
    np.savetxt('accreted_traceres.txt', ids.astype(int))
