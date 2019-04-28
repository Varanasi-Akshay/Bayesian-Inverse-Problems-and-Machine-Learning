from pyStructModel import calstress
import numpy as np
from torch import from_numpy
from scipy.linalg import sqrtm
import IPython
from math import pi

# default structural model setting
parameters = np.array([100, 0.0, 1.0 , 1.0, 100000, 0.0, 25./180*pi, 1.155, 0.03, 1.0, 1.24, 1.0, 1000, 1.0, 0.0, 0.0, 1.0, 0.0001, 1000, 1.0, 0.0, 0.0, 1.0])


def sampleGenerator(Erange = (0.115, 0.31, 0.0), n = 21, crosslinked=True, withoutmatrix=True):
# default sample range:
#     0 <= E_11 <= 0.15
#     0 <= E_22 <= 0.25
#     - 0.15 <= E_12 <= 0.15
# return [[E_11, E_22, E12], ...], [[S_11, S_22, S_12], ...]

    E_11_range = Erange[0]
    E_22_range = Erange[1]
    E_12_range = Erange[2]

    if isinstance(n, int):
        n = (n, n, 1)
    
    num_of_samples = n[0]*n[1]*n[2]

    # # log scaled
    # E11 = E_11_range - E_11_range/np.exp(2)*(
    #     np.exp(np.linspace(0, 2, n[0])) - 1
    # )
    # E11 = np.flip(E11)
    # E22 = E_22_range - E_22_range/np.exp(2)*(
    #     np.exp(np.linspace(0, 2, n[0])) - 1
    # )
    # E22 = np.flip(E22)
    # E12 = np.linspace(-E_12_range, E_12_range, n[2])

    # linear scaled
    E11 = np.linspace(0, E_11_range, n[0])
    E22 = np.linspace(0, E_22_range, n[0])
    E12 = np.linspace(-E_12_range, E_12_range, n[2])

    xv, yv, zv = np.meshgrid(E11, E22, E12)
    xv = xv.flatten()
    yv = yv.flatten()
    zv = zv.flatten()

    E_samples = np.array([xv, zv, zv, yv]).T

    S_samples = np.zeros_like(E_samples)
    F_samples = np.zeros_like(E_samples)

    for i in range(num_of_samples):
        e = np.copy(E_samples[i])
        c = np.zeros_like(e)
        s = np.zeros_like(e)

        c = 2*e.reshape((2,2))+np.identity(2)

        f = sqrtm(c)
        f = f.flatten()

        calstress(parameters, f, s, crosslinked, withoutmatrix)

        S_samples[i] = s
        F_samples[i] = f


    # IPython.embed()
    return E_samples[:, [0,3,1]], S_samples[:, [0,3,1]]
