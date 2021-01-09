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
/*
 * CDF and PDF of positive linear combination of chi2 random variables.
 *
 * Implementation based on:
 * 
 * Probability content of regions under spherical normal distributions IV
 * Ann. Math. Statist. (1962), Vol. 33, No.2 
 * H. Ruben
 *
 * Algorithm AS 106, Appl. Statist. (1977) Vol. 26, No.1
 * J. Sheil and I. O'Muircheartaigh
 *
 * Algorithm AS 204, Appl. Statist. (1984) Vol. 33, No.3 
 * R.W. Farebrother
 * 
 * Changes made to AS 204:
 *
 *   - Multi-precision
 *  
 * ifaults:
 *  -i: Constraints not satisfied
 *   1: non-fatal underflow of a0
 *   2: constraints not satisfied
 *   3: probability < -1
 *   4: required accuracy not obtained
 *   5: not in range [0,1]
 *   6: dnsty is negative
 *   9: 4 + 5
 *  10: 4 + 6 
 *   0: ok
 * 
 */

#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/distributions/normal.hpp>
#include <string>

using namespace boost::multiprecision;
using namespace boost::math;

/*
   Ruben
   
   in:
   
   lambda: Weights of linear combination
   mult  : Multiplicities
   delta : Non-centrality parameters
   n     : # terms
   c     : Evaluation point
   mode  : Expansion parameter. Mode >0: mode*EVmin, else 2/(1/EVmin+1/Evmax) (cf. Ruben 1962)
   maxit : Maximum number of terms to sum
   eps   : Requested accuracy
   
   
   out:
   
   ruben : CDF(c)
   dnsty : PDF(c)
   ifault: Error code
*/
template <typename REAL> REAL ruben(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, REAL* dnsty, int* ifault) 
{
    REAL LN_SQRT_PIH = log(sqrt(constants::pi<REAL>() / REAL(2.0)));
    
    REAL TOL;
    if (std::is_same<REAL, double>::value) { 
        TOL = -200;
    }  else {
        TOL = -1000;
    }
    
	if(n < 1 || c <= 0.0 || maxit < 1 || eps <= 0.0) {
		// Consistency checks failed
		ifault[0] = 2;
        
		return -2.0;	
        
	} else {
		// Init
		int i;
        
		// More consistency checks
        // Find Lmin and Lmax
		REAL beta = lambda[0];
		REAL sum = lambda[0];
        REAL hold,hold2;
        REAL prbty;
        
        for(i = 0; i < n; i++) {
			hold = lambda[i];
			if(hold <=0 || mult[i] < 1 || delta[i] < 0.0) {
				ifault[0] = -i;
                
				return -7.0;
			}
            
			if(beta > hold) beta = hold;
			if(sum < hold) sum = hold;
		}

		// Set beta
        if(mode > 0) {
            beta *= mode;
        } else {
            beta = 2.0/(1.0/beta+1.0/sum);
        }
        
        int k = 0;
        sum = 1.0;
        REAL sum1 = 0.0;
        
        REAL* theta = new REAL[n];
        REAL* gamma = new REAL[n];
        
        for(i=0; i < n; i++) {
            hold = beta/lambda[i];
            gamma[i] = 1.0 - hold;
            sum *= pow(hold,mult[i]);
            sum1 += delta[i];
            k += mult[i];
            theta[i] = 1.0;
        }
	
        REAL ao = exp(0.5*(log(sum)-sum1));
        
        if(ao <= 0.0) {
            ifault[0] = 1;
            dnsty[0] = 0.0;
            delete(theta);
            delete(gamma);
            
            return 0.0;
            
        } else {
            
            // Main part of algo
         
            REAL z = c / beta;
            REAL lans,dans,pans;
        
            if(k % 2 == 0) {
                
                i = 2;
                lans = -0.5*z;
                dans = exp(lans);
                pans = 1.0 - dans;
                
            } else {
                i = 1;
                lans = -0.5 *(z+log(z)) - LN_SQRT_PIH;
                dans = exp(lans);
                
                normal_distribution<REAL> NDIST(0,1);
                
                pans = cdf(NDIST,sqrt(z)) - cdf(NDIST,-sqrt(z));
            }
            
            k -= 2;
            
            for(i; i <=k; i+=2) {
                
                if(lans < TOL) {
                    lans += log(z/i);
                    dans = exp(lans);
                } else {
                    dans *= z/i;
                }
                
                pans -= dans;
            }
            
                      
            // Calc remaining terms
            REAL eps2 = eps/ao;
            REAL aoinv = 1.0/ao;
            sum = aoinv - 1.0;
            dnsty[0] = dans;
            prbty = pans;
            
            REAL* a = new REAL[maxit];
            REAL* b = new REAL[maxit];
            
            // Main loop (iterations)
            int m;
            for(m = 1; m <= maxit; m++) {
                sum1 = 0.0;
                  
                for(i = 0; i < n; i++) {
                   
                    hold2 = theta[i]*gamma[i];
                    sum1 += hold2*mult[i] + m*delta[i]*(theta[i]-hold2);
                
                    theta[i] = hold2;
                }
                
                sum1 *= 0.5;
                b[m-1] = sum1;
                
                for(i = m-1; i>=1; i--) {
                   sum1 += b[i-1]*a[m-i-1]; 
                }
                
                sum1 /= m;
                a[m-1] = sum1;
                
                k+=2;
                
                if(lans < TOL) {
                    lans += log(z/k);
                    dans = exp(lans);
                } else {
                    dans *= z/k;
                }
                
                pans -= dans;
                sum -= sum1;
                dnsty[0] += dans*sum1;
                sum1 *= pans;
                prbty += sum1;
               
                if(prbty < -aoinv) {
                    ifault[0] = 3;
                    
                    delete(a);
                    delete(b);
                    delete(theta);
                    delete(gamma);
                    
                    return -3.0;
                }
                
                
                if(abs(pans*sum) < eps2 && abs(sum1) < eps2) {
                    // Done
                    ifault[0] = 0;
                   
                    // Cleanup
                    delete(a);
                    delete(b);
                
                    goto OKEXIT;
                }
                
                
            }
            
 
            // Cleanup
            delete(a);
            delete(b);
        }
        
        ifault[0] = 4;
        
        OKEXIT:
        
        // Cleanup
        delete(theta);
        delete(gamma);
      
        dnsty[0] *= ao/(2*beta);
        prbty *= ao;
        
        if(prbty < 0.0 || prbty > 1.0 ) {
            ifault[0] += 5;
        } else {
            
            if(dnsty[0] < 0.0) {
                ifault[0] += 6;
            }
        }
        
        
        return prbty;
 	}
}


/*
    Ruben, double precision
*/
double onemin_ruben(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault) {
    double density;
    
    return 1 - ruben<double>(lb, n, nc, r, c, 1, lim, acc, &density, ifault);
}

/*
    Ruben, quad precision
*/
double onemin_ruben_128b(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault) {
    float128 density;
    float128 tmp =  1 - ruben<float128>(lb, n, nc, r, c, 1, lim, acc, &density, ifault);
    
    return tmp.convert_to<double>();
}

/*
    Ruben, 100d precision
*/
double onemin_ruben_100d(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault) {
    cpp_bin_float_100 density;
    cpp_bin_float_100  tmp =  1 - ruben<cpp_bin_float_100>(lb, n, nc, r, c, 1, lim, acc, &density, ifault);
    
    return tmp.convert_to<double>();
}

/*
    Ruben, 100d precision
*/
double onemin_ruben_200d(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault) {
    number<cpp_bin_float<200>> density;
    number<cpp_bin_float<200>> tmp =  1 - ruben<number<cpp_bin_float<200>>>(lb, n, nc, r, c, 1, lim, acc, &density, ifault);
    
    return tmp.convert_to<double>();
}



/*
   Ruben, double precision, char* return
*/
extern "C"      
void ruben_128b_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault, char* result) {
    
    float128* density = new float128[1];
    float128 tmp = ruben<float128>(lambda,mult,delta,n, c, mode, maxit, eps, density, ifault);
    
    std::string s = density[0].convert_to<std::string>();
    
    delete(density);
    strcpy(dnsty, s.c_str()); 
    
    s = tmp.convert_to<std::string>();
    strcpy(result, s.c_str()); 
}


/*
    1 - Ruben, double precision, char* return
*/
extern "C"      
void onemin_ruben_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault, char* result) {
    
    double* density = new double[1];
    double tmp = 1 - ruben<double>(lambda,mult,delta,n, c, mode, maxit, eps, density, ifault);

    std::string s = std::to_string(density[0]);
    
    delete(density);
    
    std::ostringstream out;
    out.precision(16);
    out << std::fixed << tmp;
    
    strcpy(result, out.str().c_str());       
}


/*
    1 - Ruben, 100d precision, char* return
*/
extern "C"      
void onemin_ruben_100d_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault, char* result) {
    
    cpp_bin_float_100* density = new cpp_bin_float_100[1];
    cpp_bin_float_100 tmp = 1 - ruben<cpp_bin_float_100>(lambda,mult,delta,n, c, mode, maxit, eps, density, ifault);

    std::string s = density[0].convert_to<std::string>();
    
    delete(density);
    strcpy(dnsty, s.c_str()); 
    
    s = tmp.convert_to<std::string>();
    strcpy(result, s.c_str());     
}

/*
    1 - Ruben, quad precision, char* return
*/
extern "C"      
void onemin_ruben_128b_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault, char* result) {
    
    float128* density = new float128[1];
    float128 tmp = 1 - ruben<float128>(lambda,mult,delta,n, c, mode, maxit, eps, density, ifault);
    
    std::string s = density[0].convert_to<std::string>();
    
    delete(density);
    strcpy(dnsty, s.c_str()); 
    
    s = tmp.convert_to<std::string>();
    strcpy(result, s.c_str());     
}
