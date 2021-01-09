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

#include <iostream>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <ruben.hpp>
#include <cmath>
#include <chrono> 

using namespace boost::multiprecision;
using namespace std::chrono; 

int main()
{
    /*
    std::cout << std::setprecision(std::numeric_limits<cpp_bin_float_100>::max_digits10);
    cpp_bin_float_100 pi = boost::math::constants::pi<cpp_bin_float_100>();
    std::cout << pi << std::endl;
    
    cpp_bin_float_100 L = log(cpp_bin_float_100(2.0))/8.0;
    std::cout << L << std::endl;
    */
   
/*
    std::cout << std::setprecision(std::numeric_limits<cpp_bin_float_100>::max_digits10);
    cpp_bin_float_100 L = 5e-35Q;
    std::cout << atan(L) << std::endl;
    L = 5e-50Q;
    std::cout << atan(L) << std::endl;
*/
    
   
    
    /*
    float128 T = 1e-1000Q;
    std::cout << -log(T) << std::endl;
    
     */
   
    int ifault[1];
    float128 dnsty[1];
    
    
   /*
    double lambda[2] = {0.5,0.25};
    int mult[2] = {1,1};
    double delta[2] = {0,0};
    
    // Generate Random values 
    auto f = []() -> double { return 1*(double)rand() / RAND_MAX; }; 
    
    std::vector<double> values(10); 
    std::generate(values.begin(), values.end(), f); 
    
    for(int i = 0; i < 10; i++) {
        auto start = high_resolution_clock::now(); 

        double res = ruben(lambda, mult, delta, 2, values[i], 1, 10000, 1e-10, dnsty, ifault);

        auto stop = high_resolution_clock::now(); 
        auto duration = duration_cast<milliseconds>(stop - start); 

        std::cout << "in: " << values[i] << " prob: " << res << " ifault: " << ifault[0] << " ["<< duration.count() << "ms]" << std::endl;
    }
    */
    
    double lambda[3][3] = {{6.,3.,1.},{6.,3.,1.},{6.,3.,1.}};
    int  mult[3][3] = {{1,1,1},{1,1,1},{1,1,1}};
    double delta[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double c[3] = {1.,7., 20.};
    
    double T[3] = {0.0542,0.4936,0.8760};
   
    std::cout << "Starting RUBEN algorithm tests:" << std::endl;
    //std::cout << std::setprecision(std::numeric_limits<float128>::max_digits10 );
   // std::cout << std::scientific;
    
    for(int i = 0; i < 3; i++) {
        auto start = high_resolution_clock::now(); 

        float128 res = ruben<float128>(lambda[i], mult[i], delta[i], 3, c[i], -1, 10000, 1e-50, dnsty, ifault);

        auto stop = high_resolution_clock::now(); 
        auto duration = duration_cast<milliseconds>(stop - start); 

        double diff = round( res.convert_to<double>() * 10000.0 ) / 10000.0 - T[i];
        
        if(diff == 0 && ifault[0] == 0) {
            std::cout << i <<": "<< res << " | "<< T[i] << " | ifault: " << ifault[0] << " ["<< duration.count() << "ms] OK!" << std::endl;
        } else {
            std::cout << i <<": "<< res << " | "<< T[i] << " | ifault: " << ifault[0] << " ["<< duration.count() << "ms] FAIL!" << std::endl;
        }
    }
}
