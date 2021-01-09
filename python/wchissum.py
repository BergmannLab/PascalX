#!/usr/bin/python3
import os
import cffi

from cffi import FFI
ffibuilder = FFI()

with open('wchissum.hpp') as h:
    header = h.read()
    ffibuilder.cdef(header)

    
ffibuilder.set_source("PascalX_core",header,sources=["wchissum.cpp"], include_dirs=["../build/include/"],libraries=["ruben","davies","quadmath"],library_dirs=["../build/lib/"],extra_compile_args=["-O2"])

if __name__ == "__main__":
    ffibuilder.compile(verbose=False)
