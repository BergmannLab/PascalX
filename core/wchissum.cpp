/*
    PascalX - A python3 library for high precision gene and pathway scoring for 
              GWAS summary statistics with C++ backend.
              https://github.com/BergmannLab/PascalX

    Copyright (C) 2021 Bergmann lab and contributors

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <ruben.hpp>
#include <davies.hpp>

using namespace boost::multiprecision;
using namespace boost::math;

double onemin_smart_algo(double* lambda, int n, double X, int* ifault) {
    
    std::cout << std::scientific;
   
    // Dummy centrality parameters and multiplicities
    double* nc = new double[n];
    int* ml = new int[n];
    
    // find min and max lambda
    double lmax = lambda[0];
    double lmin = lambda[0];
    double mean = 0;
    
    for(int i = 0; i < n; i++) {
        nc[i] = 0.0;
        ml[i] = 1;
        mean += lambda[i];
        if(lambda[i] < lmin) lmin = lambda[i];
        if(lambda[i] > lmax) lmax = lambda[i];
    }
    
    mean /= n;
    
    double R = lmax/lmin;
    double RESULT;
    std::cout<<"\n"<<std::endl;
    
    if(R < 100 || (n < 20 && R < 1000)) {
        
        Ruben:
        
        // Use Ruben
        std::cout << "Using Ruben (N="<< n <<  " ,M=" << mean << " ,R=" << R << ", X=" << X << ")" <<std::endl;
  
        double dnsty;
       
        // Small precision run
        RESULT = 1 - ruben<double>(lambda, ml, nc, n, X, 1, 100000, 1e-16, &dnsty, ifault);
        
        if(RESULT < 1e-14 || ifault[0] != 0) {
            float128 dn;
            float128 tmp;
            
            tmp = 1 - ruben<float128>(lambda, ml, nc, n, X, 1, 1000000, 1e-32, &dn, ifault);
                
            if(tmp < 1e-30 || ifault[0] != 0) {
                
                cpp_bin_float_100 dn;
                cpp_bin_float_100 tmp;
                
                tmp = 1 - ruben<cpp_bin_float_100>(lambda, ml, nc, n, X, 1, 10000000, 1e-100, &dn, ifault);
                RESULT = tmp.convert_to<double>();
                
             } else {
                RESULT = tmp.convert_to<double>();
            }
            
        } 
        
    } else {
        // Use Davies
        std::cout << "Using Davies (N="<< n <<  " ,M=" << mean << " ,R=" << R << ", X=" << X << ")" <<std::endl;
  
        double* trace = new double[7];
        
        RESULT = 1 - davies<double>(lambda, ml, nc, n, X, 20000, 1e-12, ifault, trace);
        
        std::cout << "RES1: " << RESULT << " [" << ifault[0] << "]"<< std::endl;
        
        if(RESULT < 1e-14 || ifault[0] != 0) {
            float128 tmp = 1 - davies<float128>(lambda, ml, nc, n, X, 100000, 1e-32, ifault, trace);
            
            if(tmp < 1e-30 || ifault[0] != 0) {
                   
                cpp_bin_float_100 tmp = 1 - davies<cpp_bin_float_100>(lambda, ml, nc, n, X, 1000000, 1e-100, ifault, trace);
                
                 RESULT = tmp.convert_to<double>();
            } else {
                 RESULT = tmp.convert_to<double>();
            }
                
        }
        
        
        if(ifault[0]==1) goto Ruben;
        
        delete(trace);
    }
    
    
    // Cleanup
    delete(nc);
    delete(ml);
 
    return RESULT;
}


extern "C"
double onemin_wchissum_m1nc0(double* lambda, int n, double X, int* ifault) {
 
    return onemin_smart_algo(lambda, n, X, ifault);
}