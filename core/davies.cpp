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

#include "davies.hpp"

double onemin_davies(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace) {
   
    return  1 - davies<double>(lb, n, nc, r, c, lim, acc, ifault, trace);
}

double onemin_davies_128b(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace) {
    float128 tmp =  1 - davies<float128>(lb, n, nc, r, c, lim, acc, ifault, trace);
    
    return tmp.convert_to<double>();
}

double onemin_davies_100d(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace) {
    cpp_bin_float_100 tmp =  1 - davies<cpp_bin_float_100>(lb, n, nc, r, c, lim, acc, ifault, trace);
    
    return tmp.convert_to<double>();
}


extern "C"
void onemin_davies_str(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace, char* result) {
   
    double tmp = 1 - davies<double>(lb, n, nc, r, c, lim, acc, ifault, trace);

    std::ostringstream out;
    out.precision(16);
    out << std::fixed << tmp;
    
    strcpy(result, out.str().c_str());      
}

extern "C"
void onemin_davies_128b_str(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace, char* result) {
   
    float128 tmp = 1 - davies<float128>(lb, n, nc, r, c, lim, acc, ifault, trace);

    std::string s = tmp.convert_to<std::string>(); //d::to_string(tmp);
    
    strcpy(result, s.c_str());      
}

extern "C"
void onemin_davies_100d_str(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace, char* result) {
   
    cpp_bin_float_100 tmp = 1 - davies<cpp_bin_float_100>(lb, n, nc, r, c, lim, acc, ifault, trace);

    std::string s = tmp.convert_to<std::string>(); //d::to_string(tmp);
    
    strcpy(result, s.c_str());      
}
