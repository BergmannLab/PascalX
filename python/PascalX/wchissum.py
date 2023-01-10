#    PascalX - A python3 library for high precision gene and pathway scoring for 
#              GWAS summary statistics with C++ backend.
#              https://github.com/BergmannLab/PascalX
#
#    Copyright (C) 2021 Bergmann lab and contributors
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from PascalX_core import lib,ffi
from PascalX import hpstats

import time

def onemin_cdf_davies(X, lb, lim=100000, acc=1e-16, mode=''):
    """
    Calculates tail probability for linear combination of chi2 distributed random variables (1-cdf(X))
    via Davies algorithm 
    
    X: Point to evaluate
    lb: weights 
    lim: Max # integration terms
    acc: Requested accuracy
    mode: '','128b','100d' the internal precision to use
    """
    _L = np.ascontiguousarray(lb, dtype='float64')
    _trace = np.zeros(7,dtype='float64')
    _ifault = np.zeros(1,dtype='int32')

    tic = time.time()

    if mode == '128b':
        res = lib.oneminwchissum_m1nc0_davies_128b(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-32),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
    elif mode == '100d':
    
        res = lib.oneminwchissum_m1nc0_davies_100d(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
    elif mode == 'auto':
        res = lib.oneminwchissum_m1nc0_davies_auto(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
        
    else:
        res = lib.oneminwchissum_m1nc0_davies(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-16),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )

    toc = time.time()
    
    return [res,_ifault[0],_trace, round(toc-tic,5)]

def fconstmin_cdf_davies(F,c, X, lb, lim=100000, acc=1e-16, mode=''):
    """
    Calculates tail probability for linear combination of chi2 distributed random variables (c-cdf(X))
    via Davies algorithm 
    
    c: Constant
    X: Point to evaluate
    lb: weights 
    lim: Max # integration terms
    acc: Requested accuracy
    mode: '','128b','100d' the internal precision to use
    """
    _L = np.ascontiguousarray(lb, dtype='float64')
    _trace = np.zeros(7,dtype='float64')
    _ifault = np.zeros(1,dtype='int32')

    tic = time.time()

    if mode == '128b':
        res = F*lib.constminwchissum_m1nc0_davies_128b(
            c,
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-32),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
    elif mode == '100d':
    
        res = F*lib.constminwchissum_m1nc0_davies_100d(
            c,
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
    elif mode == 'auto':
        res = lib.fconstminwchissum_m1nc0_davies_auto(
            F,
            c,
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            min(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
        
    else:
        res = F*lib.constminwchissum_m1nc0_davies(
            c,
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-16),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )

    toc = time.time()
    
    return [res,_ifault[0],_trace, round(toc-tic,5)]
    

    
def onemin_cdf_davies_nc(X, lb, nc, lim=100000, acc=1e-16, mode=''):
    """
    Calculates tail probability for linear combination of non-central chi2 distributed random variables (1-cdf(X))
    via Davies algorithm 
    
    X: Point to evaluate
    lb: weights 
    nc: non-centrality parameters
    lim: Max # integration terms
    acc: Requested accuracy
    mode: '','128b','100d' the internal precision to use
    """
    _L = np.ascontiguousarray(lb, dtype='float64')
    _nc = np.ascontiguousarray(nc, dtype='float64')
    _trace = np.zeros(7,dtype='float64')
    _ifault = np.zeros(1,dtype='int32')

    tic = time.time()

   
    if mode == '128b':
        res = lib.oneminwchissum_m1_davies_128b(
            ffi.cast("double *",_L.ctypes.data),
            ffi.cast("double *",_nc.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-32),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
    elif mode == '100d':
    
        res = lib.oneminwchissum_m1_davies_100d(
            ffi.cast("double *",_L.ctypes.data),
            ffi.cast("double *",_nc.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
    elif mode == 'auto':
        res = lib.oneminwchissum_m1_davies_auto(
            ffi.cast("double *",_L.ctypes.data),
            ffi.cast("double *",_nc.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )
        
    else:
        res = lib.oneminwchissum_m1_davies(
            ffi.cast("double *",_L.ctypes.data),
            ffi.cast("double *",_nc.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-16),
            ffi.cast("int *",_ifault.ctypes.data),
            ffi.cast("double *",_trace.ctypes.data)
        )

    toc = time.time()
    
    return [res,_ifault[0],_trace, round(toc-tic,5)]    
    
    
    
def onemin_cdf_ruben(X, lb, lim=100000, acc=1e-16, mode=''):
    """
    Calculates tail probability for linear combination of chi2 distributed random variables (1-cdf(X))
    via Ruben's algorithm 
    
    X: Point to evaluate
    lb: weights 
    lim: Max # integration terms
    acc: Requested accuracy
    mode: '','128b','100d' the internal precision to use
    """
    _L = np.ascontiguousarray(lb, dtype='float64')
    _ifault = np.zeros(1,dtype='int32')

    tic = time.time()

    if mode == '128b':
        res = lib.oneminwchissum_m1nc0_ruben_128b(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-32),
            ffi.cast("int *",_ifault.ctypes.data)
        )
    elif mode == '100d':
    
        res = lib.oneminwchissum_m1nc0_ruben_100d(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-100),
            ffi.cast("int *",_ifault.ctypes.data)
        )
    elif mode == '200d':
    
        res = lib.oneminwchissum_m1nc0_ruben_200d(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-200),
            ffi.cast("int *",_ifault.ctypes.data)
        )
        
    else:
        res = lib.oneminwchissum_m1nc0_ruben(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            max(acc,1e-16),
            ffi.cast("int *",_ifault.ctypes.data)
        )

    toc = time.time()
    
    return [res,_ifault[0], round(toc-tic,5)]
   

def onemin_cdf_auto(X, lb, lim=1000000, acc=1e-100, mode=''):
    
    _L = np.ascontiguousarray(lb, dtype='float64')
    _ifault = np.zeros(1,dtype='int32')

    tic = time.time()
    res = lib.oneminwchissum_m1nc0_auto(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X,
            int(lim),
            acc,
            ffi.cast("int *",_ifault.ctypes.data)
        )
    toc = time.time()
    
    return [res, _ifault[0], round(toc-tic,5)]


def onemin_cdf_satterthwaite(X, lb, mode='auto'):
    tic = time.time()
    
    _L = np.ascontiguousarray(lb, dtype='float64')
    
    if mode == 'auto':
        res = lib.oneminwchissum_m1nc0_satterthwaite_auto(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-300
    elif mode == '200d':
        res = lib.oneminwchissum_m1nc0_satterthwaite_200d(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-200
    elif mode == '100d':
        res = lib.oneminwchissum_m1nc0_satterthwaite_100d(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-100
    elif mode == '128b':
        res = lib.oneminwchissum_m1nc0_satterthwaite_float128(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-32    
    else:
        res = lib.oneminwchissum_m1nc0_satterthwaite(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-15
            
    toc = time.time()
    
    return [max(rf,res), 0, round(toc-tic,5)]
    
    
    
def onemin_cdf_pearson(X, lb, mode='auto'):
    tic = time.time()

    _L = np.ascontiguousarray(lb, dtype='float64')

    if mode == 'auto':
        res = lib.oneminwchissum_m1nc0_pearson_auto(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-300
    elif mode == '200d':
        res = lib.oneminwchissum_m1nc0_pearson_200d(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-200
    elif mode == '100d':
        res = lib.oneminwchissum_m1nc0_pearson_100d(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-100
    elif mode == '128b':
        res = lib.oneminwchissum_m1nc0_pearson_float128(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-32    
    else:
        res = lib.oneminwchissum_m1nc0_pearson(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-15
    
    toc = time.time()
     
    return [max(rf,res), 0, round(toc-tic,5)]
    
def onemin_cdf_saddle(X, lb, mode='auto'):
    tic = time.time()

    _L = np.ascontiguousarray(lb, dtype='float64')

    if mode == 'auto':
        res = lib.oneminwchissum_m1nc0_saddle_auto(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-300
    elif mode == '200d':
        res = lib.oneminwchissum_m1nc0_saddle_200d(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-200
    elif mode == '100d':
        res = lib.oneminwchissum_m1nc0_saddle_100d(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-100
    elif mode == '128b':
        res = lib.oneminwchissum_m1nc0_saddle_float128(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-32    
    else:
        res = lib.oneminwchissum_m1nc0_saddle(ffi.cast("double *",_L.ctypes.data),len(_L),X)
        rf = 1e-15
    
    toc = time.time()
    
    if res > 0:
        return [max(rf,res), 0, round(toc-tic,5)]
    else:
        return [-1,1,round(toc-tic,5)]