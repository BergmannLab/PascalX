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

#include "ruben.hpp"
#include "davies.hpp"


extern "C"
double oneminwchissum_m1nc0_davies(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_davies(lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double oneminwchissum_m1nc0_davies_128b(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_davies_128b(lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double oneminwchissum_m1nc0_davies_100d(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
	
    // Init non-centralities to 0
    double* nc = (double*) calloc(N, sizeof(double));
   
    // Init multiplicities to 1
    int* mu = (int*) malloc(N*sizeof(int));
    
    for(int i = 0; i < N; i++) {
        mu[i] = 1;
    }
    
    double ret = onemin_davies_100d(lambda,mu,nc,N,X,lim,acc,ifault,trace);

    free(mu); 
    free(nc);
    
    return ret;
}

extern "C"
double oneminwchissum_m1nc0_davies_auto(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace) {
    
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
    
    if(counter > 5) return ret;
    
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
            // Try Ruben and increase iterms
            
            iterms += 20000;
                
            ret = onemin_ruben(lambda,mu,nc,N,X,iterms,prec,ifault);
            
            if(ifault[0]==0 && ret > 0 && ret > prec*1e3) {
                // All ok
                return ret;
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