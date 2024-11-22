import cython
cimport numpy as cnp
cnp.import_array()

def gauss_seidel_iteration(cnp.ndarray[cnp.double_t, ndim=2] p,
                           cnp.ndarray[cnp.double_t, ndim=2] b,
                           double pcoef,
                           double dy,
                           double dx):
    cdef int i,j

    for i in range(1,p.shape[0]-1):
        for j in range(1,p.shape[1]-1):
            p[i,j] = pcoef * ((p[i,j+1] + p[i,j-1]) * dy**2
                                + (p[i+1,j] + p[i-1,j]) * dx**2) - b[i-1,j-1]

    return p