import sys
import h5py 
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from argparse import ArgumentParser
from tracers import TracerData_t
from figures import get_tracer_tseries, get_dataset
import figures as fs



def vr_dispersion_fluid(fname, ax):
    r  = get_dataset(fname, 'r')
    vr = get_dataset(fname, 'vr')

    rs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    cs = mpl.cm.gist_earth(np.linspace(0, 1, len(rs) + 1))
    for i, (a, b) in enumerate(zip(rs[:-1], rs[1:])):
        ann = (a < r) & (r < b)
        lab = r'{} < r < {}'.format(a, b)
        ax.hist(vr[ann], bins=20, density=True, histtype='step', lw=2.0, color=cs[i], alpha=0.8, label=lab)
    ax.set_ylim([0.0, 12.5])
    ax.set_xlabel('Radial velocity', fontsize=11)
    ax.set_ylabel('Probability density', fontsize=11)
    ax.axvline(0.0, color='grey', alpha=0.7, lw=0.75)
    ax.legend(loc='best', frameon=False, fontsize=8)

    dsp = []
    rc  = [np.mean([b, b + 0.1]) for b in np.arange(1.5, 8.0, 0.1)]
    for x in np.arange(1.5, 8, 0.1):
        ann = np.where((b < r) & (r < b + 0.1))
        dsp.append(np.std(vr[ann].flat))
    return dsp


def make_figure_vr_dispersion_tracers(fname, savefile='vr_dispersion.h5'):
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=[fs.text_width, fs.column_width])
    h5f  = h5py.File(savefile, 'w')
    prof = h5f.create_group('profiles')
    disp = h5f.create_group('dispersion')

    print('Loading data...')
    r  = get_tracer_tseries(fname, 'r')
    vr = get_tracer_tseries(fname, 'vr')
    rs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    prof.create_dataset('annuli', data=rs)

    cs = mpl.cm.cubehelix(np.linspace(0, 1, len(rs) + 1))
    print('Building histograms...')
    for i, (a, b) in enumerate(zip(rs[:-1], rs[1:])):
        ann = np.where((a < r) & (r < b))
        lab = r'{} $<$ r $<$ {}'.format(a, b)
        ws  = np.ones(len(vr[ann].flat)) / len(vr[ann].flat)
        bins = np.linspace(-0.75, 0.75, 100)
        n, e, p = ax0.hist(vr[ann].flat, bins=bins, weights=ws, histtype='step', lw=1.5, color=cs[i], alpha=0.8, label=lab)
        g = prof.create_group('{:.1f}:{:.1f}'.format(a, b))
        g.create_dataset('values', data=n)
        g.create_dataset('edges' , data=e)
    
    ax0.set_ylim([0, 0.45])
    ax0.set_xlim([-0.5, 0.5])
    ax0.set_xlabel('Radial velocity')
    ax0.set_ylabel('Probability density')
    ax0.axvline(0.0, color='grey', alpha=0.7, ls='--', lw=0.75)
    ax0.legend(loc='upper left', frameon=False)
    
    print('Calculating dispersion curve...')
    dsp  = []
    rdsp = np.arange(1.5, 8.0, 0.1)
    rc  = [np.mean([b, b + 0.1]) for b in rdsp]
    disp.create_dataset('radii', data=rc)
    for b in rdsp:
        ann = np.where((b < r) & (r < b + 0.1))
        dsp.append(np.std(vr[ann].flat))
    disp.create_dataset('dispersion', data=dsp)
    ax1.plot(rc, dsp, lw=2.0, color='salmon', label=r'$STD[v_r]$')
    ax1.set_xlabel('Radius')
    ax1.legend(loc='upper right', frameon=False)
    h5f.close()
    return fig


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('file', help='a tracer_tseries.h5 files')
    parser.add_argument('--hardcopy', action='store_true')
    args = parser.parse_args()

    file = args.file
    fs.configure_matplotlib()
    fig = make_figure_vr_dispersion_tracers(file)
    # vr_dispersion_fluid(chkpt, range=[1,8], dr=0.5)
    
    if args.hardcopy is True:
        print('   Saving...')
        plt.savefig('tracer_vr_dispersion.pdf', dpi=400, bbox_inches='tight', pad_inches=0.05)
    else:
        plt.show()
