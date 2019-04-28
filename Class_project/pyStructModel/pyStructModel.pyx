# distutils: language = c++

cimport numpy as np
import numpy as np
from libcpp cimport bool

from cython.view cimport array as cvarray

cdef extern from "BetaDistribution.hpp":
    pass

cdef extern from "ModelCollagenTissueExp.hpp":
    pass

cdef extern from "ModelPSstruc.hpp":
    pass


cdef extern from "stress.hpp":
    void stress(double parameters[23], double strain[4], double res[4], bool crosslinked, bool withoutmatrix)

def calstress( np.ndarray[np.float_t, ndim=1, mode="c"] parameters not None,\
               np.ndarray[np.float_t, ndim=1, mode="c"] strain not None, \
               np.ndarray[np.float_t, ndim=1, mode="c"] res not None, \
               bool crosslinked, \
               bool withoutmatrix):
    """ wrap  """

    cdef double[:] viewparameters = parameters
    cdef double[:] viewstrain = strain
    cdef double[:] viewres = res

    stress(&viewparameters[0], &viewstrain[0], &viewres[0], crosslinked, withoutmatrix)
    return None
