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
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>

using namespace boost::multiprecision;
using namespace boost::math;

extern "C"
double invchi2cdf_1mx(double x, int dof) {
    chi_squared chi2 = chi_squared(dof);
    
    return quantile(chi2,1.-x);
}

extern "C"
double invchi2cdf_1mx_128b(double x, int dof) {
    float128 tmp = float128(1.) - float128(x);
    
    chi_squared_distribution<float128> chi2(dof);
    
    return quantile(chi2,tmp).convert_to<double>();
}

extern "C"
double invchi2cdf_1mx_100d(double x, int dof) {
    cpp_bin_float_100 tmp = cpp_bin_float_100(1.) - cpp_bin_float_100(x);
    
    chi_squared_distribution<cpp_bin_float_100> chi2(dof);
    
    return quantile(chi2,tmp).convert_to<double>();
}

extern "C"
double onemin_chi2cdf(double x, int dof) {
    chi_squared chi2 = chi_squared(dof);
    
    return 1. - cdf(chi2,x);
}

extern "C"
double onemin_chi2cdf_128b(double x, int dof) {
    float128 tmp =  float128(x);
    
    chi_squared_distribution<float128> chi2(dof);
    
    tmp = float128(1.) - cdf(chi2,tmp);
    
    return tmp.convert_to<double>();
}

extern "C"
double onemin_chi2cdf_100d(double x, int dof) {
    cpp_bin_float_100 tmp =  cpp_bin_float_100(x);
    
    chi_squared_distribution<cpp_bin_float_100> chi2(dof);
    
    tmp = cpp_bin_float_100(1.) - cdf(chi2,tmp);
    
    return tmp.convert_to<double>();
}


/*
    Satterthwaite-Welch approximation
    
    following Box 1956: 
    "Some Theorems on Quadratic Forms Applied in the Study of Analysis of Variance Problems, I. Effect of Inequality of Variance in the One-Way Classification "
    
    Ann. Math. Statist.
    Volume 25, Number 2 (1954), 290-302.
    
     ~ g*chi(2) ~ Gamma(k=h/2,theta=2g)
*/
extern "C"
double oneminwchissum_m1nc0_satterthwaite_100d(double* lambda, int N, double X) {
    
    // Calc g,h
    double sum1 = 0;
    double sum2 = 0;
    for(int i =0; i < N; i++) {
        sum1 += lambda[i];
        sum2 += lambda[i]*lambda[i];
    }
    
    double g = sum2/sum1;
    double h = sum1*sum1/sum2;
   
    // Calc cdf
    cpp_bin_float_100 x = cpp_bin_float_100(X);
    
    gamma_distribution<cpp_bin_float_100> gamma(h/2,2*g);
   
    x = cpp_bin_float_100(1.) - cdf(gamma,x);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_satterthwaite_200d(double* lambda, int N, double X) {
    
    // Calc g,h
    double sum1 = 0;
    double sum2 = 0;
    for(int i =0; i < N; i++) {
        sum1 += lambda[i];
        sum2 += lambda[i]*lambda[i];
    }
    
    double g = sum2/sum1;
    double h = sum1*sum1/sum2;
   
    // Calc cdf
    number<cpp_bin_float<200>> x = number<cpp_bin_float<200>>(X);
    
    gamma_distribution<number<cpp_bin_float<200>>> gamma(h/2,2*g);
   
    x = number<cpp_bin_float<200>>(1.) - cdf(gamma,x);
    
    return x.convert_to<double>(); 
}