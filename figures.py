import h5py 
import numpy as np 
import matplotlib.pyplot as plt 

text_width   = 7.1
column_width = 3.35225
gold_ratio   = 0.618

Omega = 1.0

def configure_matplotlib():
    plt.rc('xtick' , labelsize=8)
    plt.rc('ytick' , labelsize=8)
    plt.rc('axes'  , labelsize=8)
    plt.rc('legend', fontsize=8)
    plt.rc('font', family='DejaVu Sans', size=8)
    plt.rc('text', usetex=True)


def get_tracer_tseries(fname, key):
    h5f = h5py.File(fname, 'r')
    ts  = h5f['tracer_tseries'][()]

    if key == 't':
        return h5f['time'][()]

    if key == 'x':
        return ts[:,1,:]

    if key == 'y':
        return ts[:,2,:]

    if key == 'r':
        return np.sqrt(ts[:,1,:]**2 + ts[:,2,:]**2)

    if key == 'vx':
        return ts[:,3,:]

    if key == 'vy':
        return ts[:,4,:]

    if key == 'v':
        return np.sqrt(ts[:,3,:]**2 + ts[:,4,:]**2)

    if key == 'vr':
        x  = ts[:,1,:]
        y  = ts[:,2,:]
        vx = ts[:,3,:]
        vy = ts[:,4,:]
        r  = np.sqrt(x * x + y * y)
        return (x * vx + y * vy) / r

    if key == 'density':
        return ts[:,5,:]

    if key == 'pressure':
        return ts[:,6,:]

    if key == 'j':
        r = np.dstack((get_tracer_tseries(fname, 'x') , get_tracer_tseries(fname, 'y')))
        v = np.dstack((get_tracer_tseries(fname, 'vx'), get_tracer_tseries(fname, 'vy')))
        return np.cross(r, v)

    if key == 'J':
        try:
            return get_tracer_tseries(fname, 'j') * ts[:, 7, :]
        except:
            print("Timeseries doesn't have weights")


def get_tracer_tslice(fname, t):
    h5f  = h5py.File(fname, 'r')
    ts   = h5f['tracer_tseries'][()]
    time = h5f['time'][()]
    return ts[:,:,np.where(np.abs(time - t) is np.min(np.abs(time - t)))]
    # return ts[:, :, t]


def get_dataset(chkpt, key):
    f5 = h5.File(chkpt, 'r')
    nb = f5['model']['num_blocks'][()]
    bs = f5['model']['block_size'][()]
    dl = f5['model']['domain_radius'][()]
    xc = midpoint(np.linspace(-dl, dl, bs * nb + 1))

    if key is 'x':
        x, _ = [a.T for a in np.meshgrid(xc, xc)]
        return x

    if key is 'y':
        _, y = [a.T for a in np.meshgrid(xc, xc)]
        return y

    if key is 'r':
        x, y = [a.T for a in np.meshgrid(xc, xc)]
        return np.sqrt(x * x + y * y)

    if key in ['vx', 'vy', 'rho']:
        q = np.zeros([bs * nb, bs * nb])
        for block in f5['state']['solution']:
            level, index = block.split(":")
            (i, j) = [int(i) for i in index.split('-')]
            (m, n) = (i * bs, j * bs)  
            u = f5['state']['solution'][block]['conserved']
            d = u[()][str(0)]
            if key is 'rho':
                q[m:m+bs, n:n+bs] = d
            elif key is 'vx':
                q[m:m+bs, n:n+bs] = u[()][str(1)] / d
            elif key is 'vy':
                q[m:m+bs, n:n+bs] = u[()][str(2)] / d
        return q

    if key is 'vr': 
        x  = get_dataset(fname, 'x')
        y  = get_dataset(fname, 'y')
        vx = get_dataset(fname, 'vx')
        vy = get_dataset(fname, 'vy')
        return (x * vx + y * vy) / get_dataset(fname, 'r')


def midpoint(a):
    return (a[:-1] + a[1:]) / 2.