import sys
import h5py 
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from argparse import ArgumentParser
from tracers import TracerData_t



red     = [237/255, 102/255,  93/255]
blue    = [114/255, 158/255, 206/255]
purp    = [123/255, 102/255, 210/255]
green   = [105/255, 183/255, 100/255]
orange  = [255/255, 187/255, 120/255]


Omega = 1.0

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


def jacobi_constant(x, y, vx, vy, t, corot=False):
    # Add flag to calculate binary location or use corot frame
    # 
    # Cj_crit = 3.6409381618949714 ~ 3.64
    # 
    n   = len(x)
    r2  = x * x + y * y
    v2  = vx * vx + vy * vy
    r_0 = np.sqrt((x - 0.5)**2 + y * y) # corot frame
    r_1 = np.sqrt((x + 0.5)**2 + y * y)
    U   = 0.5 / r_0 + 0.5 / r_1 + 0.5 * Omega * Omega * r2
    return 2.0 * U - v2


def specific_angular_momentum(x, y, vx, vy):
    r = np.column_stack((x , y ))
    v = np.column_stack((vx, vy))
    return np.cross(r, v)


def plot_block_tracers(ax, file, q, corot=False, **kwargs):
    (iD, x, y, vx, vy, rho, p) = TracerData_t(file).unpack()
    
    time = h5py.File(file, 'r')['time'][()]
    if corot:
        r     = np.sqrt(x * x + y * y)
        phi   = np.arctan2(y, x)
        theta = Omega * time
        x  = r * np.cos(phi - theta)
        y  = r * np.sin(phi - theta)

    if q is 'density':
        cmap = mpl.cm.magma
        ax.scatter(x, y, s=5.0, c=np.log10(rho), cmap=cmap)

    if q is 'pressure':
        cmap = mpl.gist_heat
        ax.scatter(x, y, s=5, c=np.log10(p), cmap=cmap)

    if q is 'veclocity':
        cmap = mpl.cm.viridis
        ax.scatter(x, y, s=5.0, c=np.sqrt(vx**2 + vy**2))
        ax.quiver(x, y, vx, vy, scale=85.0)
        
    elif q is 'jacobi_constant':
        cmap = mpl.cm.coolwarm_r
        cj = jacobi_constant(x, y, vx, vy, time, corot=corot)
        ax.scatter(x, y, s=5.0, c=cj, cmap=cmap, vmin=3.58, vmax=3.7)
    
    elif q is 'specific_angular_momentum':
        cmap = mpl.cm.twilight
        spam = specific_angular_momentum(x, y, vx, vy)
        ax.scatter(x, y, s=5.0, c=spam, cmap=cmap, vmin=-2.5, vmax=2.5)

    elif q is 'streams':
        r   = np.sqrt( x *  x +  y *  y)
        c   = jacobi_constant(x, y, vx, vy, time, corot=corot)
        d1  = TracerData_t(file).distance_to_component_1()
        d2  = TracerData_t(file).distance_to_component_2()
        rh  = 0.25        
        ck  = 3.64
        eps = 0.1
        # stream = np.where((ck - c > eps * ck)&(d2 > 2 * rh)&(d1 < 6 * rh))
        stream = np.where((ck - c > eps * ck)&(x > -0.5)&(y < 1.07)&(d2 > rh)&(d1 > rh))
        # stream = np.where((ck - c > eps * ck))
        ax.scatter(x[stream], y[stream], s=5.0, c='r')
        ax.scatter(x, y, s=3.0, c='grey', alpha=0.1)

        phi = np.arctan2(y, x)

        print(len(x[stream]), np.sum(j[stream]) / len(x[stream]))
        trunc = plt.Circle((0.0, 0.0), 1.07, fill=False, color='lightblue', lw=0.75, label=r'$r_{\rm trunc}$')
        ax.add_artist(trunc)
        # ax.legend(loc='best')

    else:
        print("Enter valid quantity for plotting")


def plot_tracers(fname, q='density', domain_radius=10.0, corot=False, binary=False):
    fig, ax = plt.subplots(1, figsize=[8,8])
    plot_block_tracers(ax, fname, q, corot=corot)

    time = h5py.File(fname, 'r')['time'][...]
    ax.set_xlim([-domain_radius, domain_radius])
    ax.set_ylim([-domain_radius, domain_radius])
    ax.text(0.83, 1.01, "{:.2f} orbs".format(time / (2 * np.pi)), transform=ax.transAxes)

    if binary is True:
        theta = Omega * time
        x = 0.5 * np.cos(theta)
        y = 0.5 * np.sin(theta)
        ax.plot( x,  y, marker='+', ms=6, mew=1.5, color='k')
        ax.plot(-x, -y, marker='+', ms=6, mew=1.5, color='k')

    file_num = fname.split('.')[-2]
    filename = 'tracers_' + file_num + '.png'
    
    print('Saving ' + filename)
    plt.savefig(filename)
    # plt.show()
    plt.close()


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('--all-tracers', '-a', action='store_true')
    parser.add_argument('--domain-radius', '-r', type=float, default=None)
    parser.add_argument('--show-binary', '-b', action='store_true')
    parser.add_argument('--corotate', '-c', action='store_true')
    parser.add_argument('--other', '-o', action='store_true') # change for later
    
    args = parser.parse_args()
    radius = args.domain_radius
    corot  = args.corotate
    binary = args.show_binary
    if args.all_tracers is True:
        for f in args.files:
            plot_tracers(f, q='streams', domain_radius=radius, corot=corot, binary=binary)

    if args.other is True:
        #todo
        quit()


