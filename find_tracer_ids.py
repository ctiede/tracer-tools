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

def find_tracers_with_spam(file, range):
    data = TracerData_t(file)
    ids  = data.ids()
    j    = data.specific_angular_momentum()
    return ids[(range[0] < j) & (j < range[1])]

def find_tracers_with_radius(file, range):
    data = TracerData_t(file)
    ids  = data.ids()
    r    = data.radii()
    return ids[(range[0] < r) & (r < range[1])]

def find_switched_minidisk(f_start, f_finish):
    rh     = 0.4
    start  = TracerData_t(f_start)
    finish = TracerData_t(f_finish)
    start_on_bh1  = start.ids()[start.distance_to_component_1() < rh]
    start_on_bh2  = start.ids()[start.distance_to_component_2() < rh]
    finish_on_bh1 = finish.ids()[finish.distance_to_component_1() < rh]
    finish_on_bh2 = finish.ids()[finish.distance_to_component_2() < rh]
    switch1 = np.intersect1d(start_on_bh1, finish_on_bh2)
    switch2 = np.intersect1d(start_on_bh2, finish_on_bh1)
    return np.concatenate([switch1, switch2])


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('file', 
        help='File of tracers')
    parser.add_argument('--accreted-by', '-ab', default=None, 
        help='File to find tracers that accreted since `init-file`')
    parser.add_argument('--specific-ang-mom-in', nargs=2, type=float, default=None, 
        help='Find tracers with specific ang mom between [j0, j1]')
    parser.add_argument('--radii-in', nargs=2, type=float, default=None,
        help='Find tracers with radii between [r0, r1]')
    parser.add_argument('--switched-mini-by', default=None,
        help='File to find tracers that switched minidisks since `init-file`')
    args = parser.parse_args()
    print(args)

    if args.accreted_by is not None:
        ids = find_accreted_tracers(args.file, args.accreted_by)
        print("Found {} accreted tracers".format(len(ids)))
        np.savetxt('accreted_tracers.txt', ids.astype(int), fmt='%i')

    if args.specific_ang_mom_in is not None:
        ids = find_tracers_with_spam(args.file, args.specific_ang_mom_in)
        print("Found {} tracers with specific ang mom between {}".format(len(ids), args.specific_ang_mom_in))
        np.savetxt('tracers_specific_ang_mom.txt', ids.astype(int), fmt='%i')

    if args.radii_in is not None:
        ids = find_tracers_with_radius(args.file, args.radii_in)
        print("Found {} tracers with radius between {}".format(len(ids), args.radii_in))
        np.savetxt('tracers_radius.txt', ids.astype(int), fmt='%i')

    if args.switched_mini_by is not None:
        ids = find_switched_minidisk(args.file, args.switched_mini_by)
        print("Found {} tracers that swithced minidisks".format(len(ids)))
        np.savetxt('tracers_switched.txt', ids.astype(int), fmt='%i')