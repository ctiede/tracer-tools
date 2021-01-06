import sys
import h5py 
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 


red     = [237/255, 102/255,  93/255]
blue    = [114/255, 158/255, 206/255]
purp    = [123/255, 102/255, 210/255]
green   = [105/255, 183/255, 100/255]
orange  = [255/255, 187/255, 120/255]


# 
# Returns tuple of data for a given tracer:
#   
#   (id, x, y, vx, vy, density, pressure)
# 
def unpack_tracer_data(tracer):
    data = 0
    velo = 1
    return (tracer[data][0], 
            tracer[data][1], 
            tracer[data][2], 
            tracer[velo][0], 
            tracer[velo][1], 
            tracer[2], 
            tracer[3])


def plot_block_tracers(ax, tracers, **kwargs):
    f = np.vectorize(lambda t: unpack_tracer_data(t))
    (iD, x, y, vx, vy, rho, p) = f(tracers)

    cmap = mpl.cm.magma
    ax.scatter(x, y, s=5.0, c=rho, cmap=cmap)
    # ax.quiver(x, y, vx, vy, scale=85.0)


def plot_tracers(fname):
    h5f = h5py.File(fname, 'r')
    
    fig, ax = plt.subplots(1, figsize=[8,8])
    plot_block_tracers(ax, h5f['tracers'][...])

    # TODO: output config options in tracer.h5 files? Or just load a checkpoint from same run?
    domain_radius = 6.0
    ax.set_xlim([-domain_radius, domain_radius])
    ax.set_ylim([-domain_radius, domain_radius])

    file_num = fname.split('.')[-2]
    filename = 'tracers_' + file_num + '.png'
    
    print('Saving ' + filename)
    plt.savefig(filename)

    # plt.show()
    plt.close()


if __name__ == '__main__':
    files = sys.argv[1:]
    for f in files:
        plot_tracers(f)
