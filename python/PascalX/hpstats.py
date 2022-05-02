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

from PascalX_hpstats import lib,ffi
import numpy as np

def chi2_invcdf_1mx(x, dof=1):
    """
    Inverse chi2 cdf of 1-x
    
    Args:
        x(float) : Point to evaluate
        dof(int) : Degrees of freedom
    """
    
    assert dof > 0,"Degrees of freedom need to be > 0"
    
    if x > 1e-15:
        return lib.invchi2cdf_1mx(x,dof)
    elif x > 1e-32:
        return lib.invchi2cdf_1mx_128b(x,dof)
    elif x > 1e-100:
        return lib.invchi2cdf_1mx_100d(x,dof)
    
    raise Exception("x=",str(x),"has to be > 1e-100")
    
def onemin_chi2_cdf(x, dof=1):
    """
    1-chi2 cdf
    
    x  : Point to evaluate
    dof: Degrees of freedom
    """
    
    assert dof > 0,"Degrees of freedom need to be > 0"
    
    return lib.onemin_chi2cdf_100d(x,dof)
    
def onemin_cdf_satterthwaite_100d(X, lb):
    
    _L = np.ascontiguousarray(lb, dtype='float64')

    res = lib.oneminwchissum_m1nc0_satterthwaite_100d(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X)
  
    return res

def onemin_cdf_satterthwaite_200d(X, lb):
    
    _L = np.ascontiguousarray(lb, dtype='float64')

    res = lib.oneminwchissum_m1nc0_satterthwaite_200d(
            ffi.cast("double *",_L.ctypes.data),
            len(_L),
            X)
  
    return res

def norm_cdf_100d(x, m, s):
    """
    normal cdf
    
    x : Point to evaluate
    m : Mean
    s : Std
    """
    
    return lib.normcdf_100d(x,m,s)

def onemin_norm_cdf_100d(x, m, s):
    """
    normal cdf
    
    x : Point to evaluate
    m : Mean
    s : Std
    """
    
    return lib.onemin_normcdf_100d(x,m,s)

    