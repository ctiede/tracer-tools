import sys
import h5py
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
from tracers import TracerData_t
from figures import get_tracer_tseries
import figures as fs


def inv_gauss(x, mu=0.0, sigma=1.0, weight=1.0):
    return 1 - weight * np.exp(-(x - mu)**2 / sigma**2)


def gausspace(x0, xf, n, mu=0.0, sigma=1.0, w=1.0):
    xs = np.linspace(x0, xf, n)
    return xs * inv_gauss(xs, mu=0.0, sigma=5.0, weight=0.7)


def make_figure_djdt_tracers(fname, dt=0.1, nbins=25, toi=0.1):
    fig, ax = plt.subplots(1, figsize=[fs.column_width, fs.column_width * fs.gold_ratio])
    plt.subplots_adjust(top=0.99, bottom=0.18, left=0.15, right=0.99)
    bins = gausspace(-7, 6, nbins, sigma=2, mu=0.0, w=0.7)
    # tof = int(1 / h5py.File(fname, 'r')['modle']['toi'][()])
    tof = int(1 / toi)
    rs = np.array([1.0, 2.5, 4.5, 10.0])
    cs = mpl.cm.magma(np.linspace(0, 0.9, len(rs) - 1))
    r  = get_tracer_tseries(fname, 'r')
    j  = get_tracer_tseries(fname, 'j')
    r0 = r[:,:-int(dt*tof)].flatten()
    jdot = (j[:,int(dt*tof):] - j[:,:-int(dt*tof)]).flatten() / dt
    for (a, b, c) in zip(rs[:-1], rs[1:], cs):
        ann = np.where((a < r0) & (r0 < b))
        sze = len(jdot[ann])
        ws  = np.ones(sze) / sze
        lab = r'{:.1f} $<$ r $<$ {:.1f}'.format(a, b)
        ax.hist(jdot[ann], weights=ws, bins=bins, histtype='step', lw=1.5, color=c, label=lab)
    ax.set_yscale('log')
    ax.set_xlabel(r'Change in specific angular momentum $\partial j / \partial t$')
    ax.set_ylabel('Percent of total')
    ax.legend(loc='upper left')
    return fig


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('file', help='a tracer_tseries.h5 files')
    parser.add_argument('--hardcopy', action='store_true')
    args = parser.parse_args()
    
    file = args.file
    fs.configure_matplotlib()
    fig = make_figure_djdt_tracers(file)
    plt.show()


