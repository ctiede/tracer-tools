import h5py
import numpy as np

def unpack_tracer_data(tracer):
    data = 0
    velo = 1
    return (tracer[data][0], # id       (0)
            tracer[data][1], # x        (1)
            tracer[data][2], # y        (2)
            tracer[velo][0], # vx       (3)
            tracer[velo][1], # vy       (4)
            tracer[2],       # density  (5)
            tracer[3])       # pressure (6)


class TracerData_t:
    def __init__(self, fname):
        self.file = fname

    def time(self):
        return h5py.File(self.file, 'r')['time'][...]

    def unpack(self):
        h5f    = h5py.File(self.file, 'r')
        unpack = np.vectorize(lambda t: unpack_tracer_data(t))
        return unpack(h5f['tracers'][...])

    def stack(self):
        return np.column_stack(self.unpack())

    def time(self):
        h5f = h5py.File(self.file, 'r')
        return h5f['time'][()]

    def ids(self):
        return self.unpack()[0]

    def radii(self):
        (ID, x, y, vx, vy, rho, p) = self.unpack()
        return np.sqrt(x * x + y * y)

    def velocities(self):
        (ID, x, y, vx, vy, rho, p) = self.unpack()
        return np.sqrt(vx * vx + vy * vy)

    def densities(self):
        return self.unpack()[5]

    def angular_momentum_density(self):
        (ID, x, y, vx, vy, rho, p) = self.unpack()
        r = np.column_stack((x , y))
        v = np.column_stack((vx, vy))
        return np.cross(r, v) * rho

    def specific_angular_momentum(self):
        (ID, x, y, vx, vy, rho, p) = self.unpack()
        r = np.column_stack((x , y))
        v = np.column_stack((vx, vy))
        return np.cross(r, v)

    def distance_component_1(self):
        (ID, x, y, vx, vy, rho, p) = self.unpack()
        Omega = 1.0
        time  = self.time()
        xh = 0.5 * np.cos(Omega * time)
        yh = 0.5 * np.sin(Omega * time)
        return np.sqrt((x - xh)**2 + (y - yh)**2)

    def distance_component_2(self):
        (ID, x, y, vx, vy, rho, p) = self.unpack()
        Omega = 1.0
        time  = self.time()
        xh = - 0.5 * np.cos(Omega * time)
        yh = - 0.5 * np.sin(Omega * time)
        return np.sqrt((x - xh)**2 + (y - yh)**2)


