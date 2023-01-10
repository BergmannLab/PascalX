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
#define _USE_MATH_DEFINES

#include "ruben.hpp"
#include "davies.hpp"
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>

extern "C"
double oneminwchissum_m1_davies(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_davies(lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    
    return ret;
}


extern "C"
double oneminwchissum_m1nc0_davies(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    double ret = oneminwchissum_m1_davies(lambda, nc, N, X, lim, acc, ifault, trace); 

    free(nc); 
    
    return ret;
}




extern "C"
double oneminwchissum_m1_davies_128b(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_davies_128b(lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    
    return ret;
}

extern "C"
double oneminwchissum_m1nc0_davies_128b(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
  
    double ret = oneminwchissum_m1_davies_128b(lambda,nc,N,X,lim,acc,ifault,trace);

    free(nc); 
    
    return ret;
}

extern "C"
double oneminwchissum_m1_davies_100d(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_davies_100d(lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    
    return ret;
}


extern "C"
double oneminwchissum_m1nc0_davies_100d(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
    
    double ret = oneminwchissum_m1_davies_100d(lambda,nc,N,X,lim,acc,ifault,trace);

    free(nc);
    
    return ret;
}

extern "C"
double oneminwchissum_m1_davies_auto(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace) {
    
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double prec = 1e-6;
    int iterms = (N < 10) ? 100000000 : 1000000;
    
    int counter = 0;
    
    
    START:
     
    // Try first davies @ low precision
    double ret;
    if(prec > 1e-16) {
        ret = onemin_davies(lambda,mu,nc,N,X,iterms,prec,ifault,trace);
    } else {
        if(prec > 1e-32) {
            ret = onemin_davies_128b(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
        } else {
            ret = onemin_davies_100d(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
        }
    }
   
    counter++;
    
    if(counter > 8) return ret;
    
    switch(ifault[0]) {
        case 0: {
           
            if(ret > 0 && ret > prec*1e3) {
                
                // All ok
                return ret;
            
            } else {
                // Not enough precision
                prec *= 1e-6;
                
                goto START;
            }
        }
            
        case 1: {
            // Not enough int terms
           
            iterms += 1000000;
            
            goto START;
        }
            
        case 2: {
            // Not accurate, increase interal precision
            
           if(prec > 1e-32) {
               ret = onemin_davies_128b(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
            } else {
               ret = onemin_davies_100d(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
            }
            
            if(ifault[0]==0 && ret > 0 && ret > prec*1e3) {
                // All ok
                return ret;
            } 
           
            prec *= 1e-6;
            
            goto START;
        }  
            
        case 4: {
            // Increase # integration terms
            iterms *= 2;
            
            goto START;
        }      
            
    }
    
    
    
    free(mu); 

    return ret;
    
}

extern "C"
double oneminwchissum_m1nc0_davies_auto(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
    
    double ret = oneminwchissum_m1_davies_auto(lambda, nc, N, X, lim, acc, ifault, trace); 
    
    free(nc);
  
    return ret;
    
}




extern "C"
double constminwchissum_m1_davies(double x,double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = constmin_davies(x,lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    
    return ret;
}


extern "C"
double constminwchissum_m1nc0_davies(double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    double ret = constminwchissum_m1_davies(x, lambda, nc, N, X, lim, acc, ifault, trace); 

    free(nc); 
    
    return ret;
}


extern "C"
double constminwchissum_m1nc0_davies_128b(double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = constmin_davies_128b(x,lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double constminwchissum_m1nc0_davies_100d(double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = constmin_davies_100d(x,lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double fconstminwchissum_m1nc0_davies_auto(double F,double x,double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
    
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double prec = 1e-6;
    int iterms = (N < 10) ? 100000000 : 1000000;
    
    int counter = 0;
    
    
    START:
     
    // Try first davies @ low precision
    double ret;
    if(prec > 1e-16) {
        ret = F*constmin_davies(x,lambda,mu,nc,N,X,iterms,prec,ifault,trace);
    } else {
        if(prec > 1e-32) {
            ret = F*constmin_davies_128b(x,lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
        } else {
            ret = F*constmin_davies_100d(x,lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
        }
    }
   
    counter++;
    
    if(counter > 8) return ret;
    
    switch(ifault[0]) {
        case 0: {
           
            if(ret > 0 && ret > prec*1e3) {
                
                // All ok
                return ret;
            
            } else {
                // Not enough precision
                prec *= 1e-6;
                
                goto START;
            }
        }
            
        case 1: {
            // Not enough int terms
           
            iterms += 1000000;
            
            goto START;
        }
            
        case 2: {
            // Not accurate, increase interal precision
            
           if(prec > 1e-32) {
               ret = F*constmin_davies_128b(x,lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
            } else {
               ret = F*constmin_davies_100d(x,lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
            }
            
            if(ifault[0]==0 && ret > 0 && ret > prec*1e3) {
                // All ok
                return ret;
            } 
           
            prec *= 1e-6;
            
            goto START;
        }  
            
        case 4: {
            // Increase # integration terms
            iterms *= 2;
            
            goto START;
        }      
            
    }
    
    
    
    free(mu); 
    free(nc);
    free(trace);
    
    return ret;
    
}




extern "C"
double oneminwchissum_m1nc0_ruben(double* lambda, int N, double X, int lim, double acc, int* ifault) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_ruben(lambda,mu,nc,N,X,lim,acc,ifault);

    free(mu); 
    free(nc);
    
    return ret;
}


extern "C"
double oneminwchissum_m1nc0_ruben_128b(double* lambda, int N, double X, int lim, double acc, int* ifault) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_ruben_128b(lambda,mu,nc,N,X,lim,acc,ifault);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double oneminwchissum_m1nc0_ruben_100d(double* lambda, int N, double X, int lim, double acc, int* ifault) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_ruben_100d(lambda,mu,nc,N,X,lim,acc,ifault);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double oneminwchissum_m1nc0_ruben_200d(double* lambda, int N, double X, int lim, double acc, int* ifault) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_ruben_200d(lambda,mu,nc,N,X,lim,acc,ifault);

    free(mu); 
    free(nc);
    
    return ret;
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
double oneminwchissum_m1nc0_satterthwaite(double* lambda, int N, double X) {
    
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
    
    gamma_distribution<double> gamma(h/2,2*g);
   
    return 1. - cdf(gamma,X);
    
}

extern "C"
double oneminwchissum_m1nc0_satterthwaite_float128(double* lambda, int N, double X) {
    
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
    float128 x = float128(X);
    
    gamma_distribution<float128> gamma(h/2,2*g);
   
    x = 1. - cdf(gamma,x);
    
    return x.convert_to<double>(); 
}

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

extern "C"
double oneminwchissum_m1nc0_satterthwaite_300d(double* lambda, int N, double X) {
    
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
    number<cpp_bin_float<300>> x = number<cpp_bin_float<300>>(X);
    
    gamma_distribution<number<cpp_bin_float<300>>> gamma(h/2,2*g);
   
    x = number<cpp_bin_float<300>>(1.) - cdf(gamma,x);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_satterthwaite_500d(double* lambda, int N, double X) {
    
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
    number<cpp_bin_float<500>> x = number<cpp_bin_float<500>>(X);
    
    gamma_distribution<number<cpp_bin_float<500>>> gamma(h/2,2*g);
   
    x = number<cpp_bin_float<500>>(1.) - cdf(gamma,x);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_satterthwaite_1000d(double* lambda, int N, double X) {
    
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
    number<cpp_bin_float<1000>> x = number<cpp_bin_float<1000>>(X);
    
    gamma_distribution<number<cpp_bin_float<1000>>> gamma(h/2,2*g);
   
    x = number<cpp_bin_float<1000>>(1.) - cdf(gamma,x);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_satterthwaite_auto(double* lambda, int N, double X) {
    double res = oneminwchissum_m1nc0_satterthwaite(lambda,N,X);
    
    if (res<1e-15) {
        res = oneminwchissum_m1nc0_satterthwaite_float128(lambda,N,X);
        
        if(res < 1e-32) {
            res = oneminwchissum_m1nc0_satterthwaite_100d(lambda,N,X);
        
            if(res < 1e-98) {
                res = oneminwchissum_m1nc0_satterthwaite_200d(lambda,N,X);
                
                if(res < 1e-195) {
                    res = oneminwchissum_m1nc0_satterthwaite_300d(lambda,N,X);
                }
            }
        }
    }
    
    return res;
}

/*
    Pearson's approximation
    (central case)
    
    following J. P. Imhof:
    Biometrika , Dec., 1961, Vol. 48, No. 3/4 (Dec., 1961), pp. 419-426
*/
    
extern "C"
double oneminwchissum_m1nc0_pearson(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    
    // Calc cdf
    
    chi_squared chisq(h);
   
    double x = 1. - cdf(chisq,y);
    
    return x; 
}

extern "C"
double oneminwchissum_m1nc0_pearson_float128(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    
    float128 Y(y);
    
    // Calc cdf
    
    chi_squared_distribution<float128> chisq(h);
   
    float128 x = 1. - cdf(chisq,Y);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_pearson_100d(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    number<cpp_bin_float<100>> Y(y);
    
    // Calc cdf  
    chi_squared_distribution<number<cpp_bin_float<100>>> chisq(h);
   
    number<cpp_bin_float<100>> x = 1. - cdf(chisq,Y);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_pearson_200d(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    number<cpp_bin_float<200>> Y(y);
    
    // Calc cdf  
    chi_squared_distribution<number<cpp_bin_float<200>>> chisq(h);
   
    number<cpp_bin_float<200>> x = 1. - cdf(chisq,Y);
    
    return x.convert_to<double>(); 
}


extern "C"
double oneminwchissum_m1nc0_pearson_300d(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    number<cpp_bin_float<300>> Y(y);
    
    // Calc cdf  
    chi_squared_distribution<number<cpp_bin_float<300>>> chisq(h);
   
    number<cpp_bin_float<300>> x = 1. - cdf(chisq,Y);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_pearson_500d(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    number<cpp_bin_float<500>> Y(y);
    
    // Calc cdf  
    chi_squared_distribution<number<cpp_bin_float<500>>> chisq(h);
   
    number<cpp_bin_float<500>> x = 1. - cdf(chisq,Y);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_pearson_1000d(double* lambda, int N, double X) {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    
    for(int i = 0; i < N; i++) {
        c1 += lambda[i];
        c2 += lambda[i]*lambda[i];
        c3 += lambda[i]*lambda[i]*lambda[i];
    }
    
    double h = c2*c2*c2/(c3*c3);
    double y = std::max(0.,(X-c1)*sqrt(h/c2)+h);
    number<cpp_bin_float<1000>> Y(y);
    
    // Calc cdf  
    chi_squared_distribution<number<cpp_bin_float<1000>>> chisq(h);
   
    number<cpp_bin_float<1000>> x = 1. - cdf(chisq,Y);
    
    return x.convert_to<double>(); 
}

extern "C"
double oneminwchissum_m1nc0_pearson_auto(double* lambda, int N, double X) {
    double res = oneminwchissum_m1nc0_pearson(lambda,N,X);
    
    if (res < 1e-15) {
        res = oneminwchissum_m1nc0_pearson_float128(lambda,N,X);
        
        if(res < 1e-32) {
            res = oneminwchissum_m1nc0_pearson_100d(lambda,N,X);
        
            if(res < 1e-98) {
                res = oneminwchissum_m1nc0_pearson_200d(lambda,N,X);
                
                if(res < 1e-195) {
                    res = oneminwchissum_m1nc0_pearson_300d(lambda,N,X);
                }
            }
        }
    }
    
    return res;
}

/*
    Saddle-point approximation
    (central case)
    
    following Kuonen:
    Biometrika (1999) 86, 4, pp. 929-935
    
    returns -1 if error occurs
*/
    
double K(double* lambda, int N, double zeta) {
    double r = 0;
    for(int i =0; i < N; i++) {
        r += log(1-2*zeta*lambda[i]);
    }
    
    return -0.5*r;
}

double K1(double* lambda, int N, double zeta) {
    double r = 0;
    for(int i =0; i < N; i++) {
        r += lambda[i]/(1-2*zeta*lambda[i]);
    }
    
    return r;
}

double K2(double* lambda, int N, double zeta) {
    double r = 0;
    for(int i =0; i < N; i++) {
        r += lambda[i]*lambda[i]/pow(1-2*zeta*lambda[i],2);
    }
    
    return 2*r;
}

extern "C"
double oneminwchissum_m1nc0_saddle(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        double Z = (w + log(v/w)/w)/sqrt(2.);
        
        // Calc and return solution
        return (0.5 - 0.5*erf( Z ));
            
    } catch(...) {
        return -1;
    }   
}

extern "C"
double oneminwchissum_m1nc0_saddle_float128(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        
        double z = (w + log(v/w)/w)/sqrt(2.);
        
        float128 Z(z);
        
        // Calc and return solution
        return (0.5 - 0.5*erf( Z )).convert_to<double>();
        
    } catch(...) {
        return -1;
    }
}


extern "C"
double oneminwchissum_m1nc0_saddle_100d(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        double z = (w + log(v/w)/w)/sqrt(2.);
        
        number<cpp_bin_float<100>> Z(z);
        
        // Calc and return solution
        return (0.5 - 0.5*erf( Z )).convert_to<double>();
        
    } catch(...) {
        return -1;
    }   
}



extern "C"
double oneminwchissum_m1nc0_saddle_200d(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        double z = (w + log(v/w)/w)/sqrt(2.);
        
        number<cpp_bin_float<200>> Z(z);
        
        // Calc and return solution
        return (0.5 - 0.5*erf( Z )).convert_to<double>();
        
    } catch(...) {
        return -1;
    }   
}

extern "C"
double oneminwchissum_m1nc0_saddle_300d(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        double z = (w + log(v/w)/w)/sqrt(2.);
        
        number<cpp_bin_float<300>> Z(z);
        
        // Calc and return solution
        return (0.5 - 0.5*erf( Z )).convert_to<double>();
        
    } catch(...) {
        return -1;
    }   
}

extern "C"
double oneminwchissum_m1nc0_saddle_500d(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        double z = (w + log(v/w)/w)/sqrt(2.);
        
        number<cpp_bin_float<500>> Z(z);
        
        // Calc and return solution
        return (0.5 - 0.5*erf( Z )).convert_to<double>();
        
    } catch(...) {
        return -1;
    }   
}

/*
extern "C"
double oneminwchissum_m1nc0_saddle_auto(double* lambda, int N, double X) {
    double res = oneminwchissum_m1nc0_saddle(lambda,N,X);
    
    if ((res < 1e-15) && (res >= 0)) {
        res = oneminwchissum_m1nc0_saddle_float128(lambda,N,X);
        
        if ((res < 1e-32)  && (res >= 0)) {
            res = oneminwchissum_m1nc0_saddle_100d(lambda,N,X);
        
            if ((res < 1e-98)  && (res >= 0)) {
                res = oneminwchissum_m1nc0_saddle_200d(lambda,N,X);
                
                if ((res < 1e-195)  && (res >= 0)) {
                    res = oneminwchissum_m1nc0_saddle_300d(lambda,N,X);
                }
            }
        }
    }
    
    return res;
}
*/


extern "C"
double oneminwchissum_m1nc0_saddle_auto(double* lambda, int N, double X) {
    double sum = lambda[0];
        
    // find maxb
    double ma = 1./lambda[0];
    for(int i=1;i<N;i++) {
        double tmp = 1./lambda[i];
        sum += lambda[i];
        
        if(tmp < ma) {
            ma = tmp;
        }
    }
    ma *= 0.5;
    
    // Do not use in unstable regime
    if (abs(sum-X)/X < 1e-5) {
        return -1;
    }
    
    const int digits = std::numeric_limits<double>::digits; 
    int get_digits = static_cast<int>(digits * 0.6);
    
    // Solve for zeta 
    const boost::uintmax_t maxit = 10000;
    boost::uintmax_t it = maxit;
    try {
        double zeta = boost::math::tools::newton_raphson_iterate(
        [lambda,N,X](double x) {
            return std::make_pair(K1(lambda,N,x)-X,K2(lambda,N,x));
        },
        -0.5*N/X,-0.5*N/X-1,ma,
        get_digits, it
        );
        
        // Calc parameters    
        double v = zeta*sqrt(K2(lambda,N,zeta));
        double w = sign(zeta)*sqrt(2*(zeta*X - K(lambda,N,zeta)));
        double z = (w + log(v/w)/w)/sqrt(2.);
        
        // Calc cdf using error function
        double Z = z;
        
        double res = ( 0.5*(1 - erf( Z )) );
        
        if (res < 1e-15) {
            float128 Z(z);
            res = ( 0.5*(1 - erf( Z )) ).convert_to<double>();
        
            if (res < 1e-32) {
                number<cpp_bin_float<100>> Z(z);
                res = ( 0.5*(1 - erf( Z )) ).convert_to<double>();

                if (res < 1e-98) {
                    number<cpp_bin_float<200>> Z(z);
                    res = ( 0.5*(1 - erf( Z )) ).convert_to<double>();

                    if (res < 1e-195) {
                        number<cpp_bin_float<300>> Z(z);
                        res = ( 0.5*(1 - erf( Z )) ).convert_to<double>();
                    }
                }
            }
        }
        
        // Calc and return solution
        return res;
        
    } catch(...) {
        return -1;
    }   
}



extern "C"
double oneminwchissum_m1nc0_auto(double* lambda, int N, double X, int lim, double acc, int* ifault) {
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    

    // Init trace for davies
    double* trace = (double*) calloc(7,sizeof(double));
    
    // Calc required terms for Davies (w. prec)
    double min = lambda[0];
    double max = lambda[0];
    double reqt = lambda[0];
    for(int i = 1; i < N; i++) {
        if (lambda[i] < min) min = lambda[i];
        if (lambda[i] < max) max = lambda[i];
        reqt += lambda[i];
    }
    reqt = pow(20/M_PI,2./N)*abs(X - reqt)/min/4./M_PI;
    
    // Estimate precision needed
    double tmp = oneminwchissum_m1nc0_pearson_auto(lambda, N, X);
    
    /* Seems to slow like this ...
    double tmp = oneminwchissum_m1nc0_saddle_auto(lambda, N, X);
    
    if( tmp < 0) {
        tmp = oneminwchissum_m1nc0_pearson_auto(lambda, N, X);
    }
    */
    
    double prec = tmp*1e-3;
    double iprec = prec;
    
    int iterms = 20000;
    if (prec < 1e-16) {
        iterms = 100000;
    } 
    
    if(prec < 1e-32) {
        iterms = 1000000;
    }
    
    if(prec < 1e-95) {
        iterms = 10000000;
    }
    
    int counter = 0;
    double ret;
    double et;
    
    START:
    
    prec = std::max(prec,1e-100);
    
    et = reqt*pow(prec,-2./N)+1000;
    if (max/min > 2 && N > 5 && et < iterms) {
        // Davies
        if (iprec > 1e-16 && prec > 1e-16) {
            ret = onemin_davies(lambda,mu,nc,N,X,et,prec,ifault,trace);
        } else {
            if (iprec > 1e-32 && prec > 1e-32) {
                ret = onemin_davies_128b(lambda,mu,nc,N,X,et,prec,ifault,trace);
            } else {
                ret = onemin_davies_100d(lambda,mu,nc,N,X,et,prec,ifault,trace);
            }
        }
    } else {
        // Ruben
        if (iprec > 1e-16 && prec > 1e-16) {
            ret = onemin_ruben(lambda,mu,nc,N,X,iterms,prec,ifault);
        } else { 
            if (iprec > 1e-32 && prec > 1e-32) {
                ret = onemin_ruben_128b(lambda,mu,nc,N,X,iterms,prec,ifault);
            } else {
                ret = onemin_ruben_100d(lambda,mu,nc,N,X,iterms,prec,ifault);
            }
        }
    }
    
    switch(ifault[0]) {
        case 0: {
           
            if(ret > 0 && ret < 1 && (ret > prec*1e3 || prec < 1e-99)) {
                // All ok
                goto END;
            
            } else {
                // Not enough precision
                prec *= 1e-6;
                
                break;
            }
        }
            
        case 1: {
            // Not enough int terms
            iterms += 20000;
            break;
        }
            
        case 2: {
            // Not accurate, increase interal precision
            if (iprec > 1e-16) {
                iprec = 1e-31;
            } else {
                iprec = 1e-99;
            }
            
            prec *= 1e-6;
            break;
        }  
            
        case 7: {
            // Not accurate, increase interal precision
            if (iprec > 1e-16) {
                iprec = 1e-31;
            } else {
                iprec = 1e-99;
            }
            
            prec *= 1e-6;
            break;
        }
            
        case 4: {
            // Increase # integration terms
            iterms *= 2;
            break;
        } 
            
        case 5: {
            // Out of bounds
            
            // Max prec reached
            if(iprec <= 1e-99) {
                ret = abs(ret);
                goto END;
            }
            
            // Increase precision
            if (iprec > 1e-16) {
                iprec = 1e-31;
            } else {
                iprec = 1e-99;
            }
            
            prec *= 1e-6;
            break;
        }
            
    }
    
    counter += 1;
    if (counter > 5 ) { 
        goto END;
    }
    
    goto START;
    
            
    //HP_RUBEN:
    //ret = onemin_ruben_100d(lambda,mu,nc,N,X,10000000,1e-100,ifault); 
    //ret = oneminwchissum_m1nc0_satterthwaite_100d(lambda,N,X);
    
    END: 
    free(mu); 
    free(nc);
    free(trace);
    
    return ret;    
}



/*

extern "C"
double oneminwchissum_m1nc0_auto(double* lambda, int N, double X, int lim, double acc, int* ifault) {
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    // Init trace for davies
    double* trace = (double*) calloc(7,sizeof(double));
    
    // Calc required terms for Davies (w. prec)
    double min = lambda[0];
    double max = lambda[0];
    double reqt = lambda[0];
    for(int i = 1; i < N; i++) {
        if (lambda[i] < min) min = lambda[i];
        if (lambda[i] < max) max = lambda[i];
        reqt += lambda[i];
    }
    reqt = pow(20/M_PI,2./N)*abs(X - reqt)/min/4./M_PI;
    
    double prec = 1e-6;
    double iprec = 1e-15;
    int iterms = 20000;
    
    int counter = 0;
    double ret;
    double et;
    
    START:
    counter += 1;
    if (counter > 20) { 
        goto HP_RUBEN;
    }
    
    
    et = reqt*pow(prec,-2./N)+1000;
    if (max/min > 2 && N > 5 && et < iterms) {
        // Davies
        if (iprec > 1e-16 && prec > 1e-16) {
            ret = onemin_davies(lambda,mu,nc,N,X,et,prec,ifault,trace);
        } else {
            if (iprec > 1e-32 && prec > 1e-32) {
                ret = onemin_davies_128b(lambda,mu,nc,N,X,et,prec,ifault,trace);
            } else {
                ret = onemin_davies_100d(lambda,mu,nc,N,X,et,prec,ifault,trace);
            }
        }
    } else {
        // Ruben
        if (iprec > 1e-16 && prec > 1e-16) {
            ret = onemin_ruben(lambda,mu,nc,N,X,iterms,prec,ifault);
        } else { 
            if (iprec > 1e-32 && prec > 1e-32) {
                ret = onemin_ruben_128b(lambda,mu,nc,N,X,iterms,prec,ifault);
            } else {
                ret = onemin_ruben_100d(lambda,mu,nc,N,X,iterms,prec,ifault);
            }
        }
    }
    
    switch(ifault[0]) {
        case 0: {
           
            if(ret > 0 && ret > prec*1e3) {
                // All ok
                goto END;
            
            } else {
                // Not enough precision
                prec *= 1e-6;
                
                goto START;
            }
        }
            
        case 1: {
            // Not enough int terms
            iterms += 20000;
            goto START;
        }
            
        case 2: {
            // Not accurate, increase interal precision
            if (iprec > 1e-16) {
                iprec = 1e-31;
            } else {
                iprec = 1e-99;
            }
            
            prec *= 1e-6;
            
            goto START;
        }  
            
        case 4: {
            // Increase # integration terms
            iterms *= 2;
            
            goto START;
        } 
            
    }
    
    goto START;
    
            
    HP_RUBEN:
    //ret = onemin_ruben_100d(lambda,mu,nc,N,X,10000000,1e-100,ifault); 
    //ret = oneminwchissum_m1nc0_satterthwaite_100d(lambda,N,X);
    
    END: 
    free(mu); 
    free(nc);
    free(trace);
    
    return ret;
    
    
}

*/

/*
extern "C"
double oneminwchissum_m1nc0_auto(double* lambda, int N, double X, int lim, double acc, int* ifault) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    // Init trace for davies
    double* trace = (double*) calloc(7,sizeof(double));
    
    double prec = 1e-6;
    int iterms = 10000;
    
    int counter = 0;
    double ret;
    
    START:
    counter++;
    
    if(counter > 5) goto HP_RUBEN; 
    
    // Try first davies @ low precision
   
    if(prec > 1e-16) {
        ret = onemin_davies(lambda,mu,nc,N,X,iterms,prec,ifault,trace);
    } else {
        if(prec > 1e-32) {
            ret = onemin_davies_128b(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
        } else {
            ret = onemin_davies_100d(lambda,mu,nc,N,X,iterms,prec,ifault,trace);            
        }
    }
   
    switch(ifault[0]) {
        case 0: {
           
            if(ret > 0 && ret > prec*1e3) {
                // All ok
                goto END;
            
            } else {
                // Not enough precision
                prec *= 1e-6;
                
                goto START;
            }
        }
            
        case 1: {
            // Not enough int terms
            // Try Ruben and increase iterms
            
            iterms += 20000;
             
            if(prec < 1e-14) {
                ret = onemin_ruben(lambda,mu,nc,N,X,iterms,prec,ifault);
                if(ifault[0]==7) {
                    ret = onemin_ruben_128b(lambda,mu,nc,N,X,iterms,prec,ifault);
                }
            } else { 
                if(prec < 1e-32) {
                    ret = onemin_ruben_128b(lambda,mu,nc,N,X,iterms,prec,ifault);
                } else {
                    ret = onemin_ruben_100d(lambda,mu,nc,N,X,iterms,prec,ifault);
                }
            }
            
            if(ifault[0]==0 && ret > 0 && ret > prec*1e3) {
                // All ok
                goto END;
            }
            
            iterms += 70000;
            
            goto START;
        }
            
        case 2: {
            // Not accurate, increase interal precision
            
           if(prec > 1e-32) {
               ret = onemin_davies_128b(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
            } else {
               ret = onemin_davies_100d(lambda,mu,nc,N,X,iterms,prec,ifault,trace); 
            }
            
            if(ifault[0]==0 && ret > 0 && ret > prec*1e3) {
                // All ok
                goto END;
            } 
           
            prec *= 1e-6;
            
            goto START;
        }  
            
        case 4: {
            // Increase # integration terms
            iterms *= 2;
            
            goto START;
        }        
    }
    
    goto START;
    
    HP_RUBEN:
    //ret = onemin_davies_100d(lambda,mu,nc,N,X,10000000,1e-100,ifault,trace); 
    
    END:
    
    free(mu); 
    free(nc);
    free(trace);
    
    return ret;
}
*/