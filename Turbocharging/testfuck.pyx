import cython
import sys
sys.setrecursionlimit(5000)


def powerOfCompressor(T0, pik, k=1.4, R=287.15, m_a=1, etak=1):
    return T0


def p():
    print("hello")

cpdef int f(int n):
    cdef int a=1
    cdef int b=1
    cdef int c
    for i in range(n-1):
        c=a
        a=b
        b=c+b
    return b
