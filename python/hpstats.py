#!/usr/bin/python3
import os
import cffi

from cffi import FFI
ffibuilder = FFI()

with open('hpstats.hpp') as h:
    header = h.read()
    ffibuilder.cdef(header)

    
ffibuilder.set_source("PascalX_hpstats",header,sources=["hpstats.cpp"], include_dirs=["../build/include/"],libraries=["quadmath"],library_dirs=["../build/lib/"],extra_compile_args=["-O3"])

if __name__ == "__main__":
    ffibuilder.compile(verbose=False)
