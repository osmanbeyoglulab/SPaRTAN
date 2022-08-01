import cython
import numpy as np
cimport scipy.linalg.cython_blas as blas


def kron(double[:, ::1] a, double[:, ::1] b):
    cdef int i = int(a.shape[0])
    cdef int j = int(a.shape[1])
    cdef int k = int(b.shape[0])
    cdef int l = int(b.shape[1])
    cdef int onei = 1
    cdef double oned = 1
    cdef int m, n
    result = np.zeros((i*k, j*l), float)
    cdef double[:, ::1] result_v = result
    for n in range(i):
        for m in range(k):
            blas.dger(&l, &j, &oned, &b[m, 0], &onei, &a[n, 0], &onei, &result_v[m+k*n, 0], &l)
    return result


def removeDiagC(double[:, ::1] L, double[::1] Yv, int[::1] diag):
    cdef:
        int Yrows = int(L.shape[0])
        int Ycols = int(L.shape[1])
        int Yaxis = int(diag.shape[0])
        int i, j, moveStart, moveEnd

    for i in range(1,Yaxis):
        moveStart = diag[i-1] + 1
        moveEnd = diag[i] - 1
        for j in range(moveStart,moveEnd+1):
            Yv[j-i] = Yv[j]
            L[j-i,] = L[j]
    return( np.asarray(L[:-Yaxis,]), np.asarray(Yv[:-Yaxis]))
