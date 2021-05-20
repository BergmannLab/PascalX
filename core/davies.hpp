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
    The code in davies.hpp is based on Davies original C implementation (qfc.c, MIT licence)  
    of Algorithm AS 155: The Distribution of a Linear Combination of chi^2 Random Variables 
    
    The original C code from Robert Davies is available at http://robertnz.net/
    
    Changes made:
        - Multi-precision
        - Class structure 
        - Other minor mods
*/


#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <setjmp.h>

#define TRUE  1
#define FALSE 0
typedef int BOOL;

using namespace boost::multiprecision;
using namespace boost::math;


template <class REAL>
class DaviesAlgo {
    const REAL pi = constants::pi<REAL>();
    const REAL log28 = log(REAL(2.))/8.;
    
    REAL divis[4] = {2.0, 1.4, 1.2, 1.1};
    int rats[4]   = {1, 2, 4, 8};
 
    static REAL square(REAL x)  { return  pow(x,2); }
    static REAL cube(REAL x)  { return  pow(x,3); }
    jmp_buf env;

    REAL sigsq, lmax, lmin, mean, c;
    REAL intl, ersm;
    int count, r, lim; 
    BOOL ndtsrt, fail;
    int *n,*th; 
    
    REAL *lb,*nc;
    
    REAL CUTOFF;
    
    void counter(void) /*  count number of calls to errbd, truncation, cfe */
    {
        count = count + 1;
        if ( count > lim ) longjmp(env, 1);
    }

    REAL exp1(REAL x) { 
		return x < CUTOFF ? 0.0 : exp(x);
 	}
    
    static REAL log1(REAL x, BOOL first) /* if (first) log(1 + x) ; else  log(1 + x) - x */
    { 
        if (fabs(x) > 0.1) {
        
            return (first ? log(1.0 + x) : (log(1.0 + x) - x));
        
        } else {
            REAL s, s1, term, y, k;
            y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
            s = (first ? 2.0 : - x) * y;
            y = square(y);
            
            for (s1 = s + term / k; s1 != s; s1 = s + term / k)
            { 
                k = k + 2.0; term = term * y; s = s1; 
            }
            
            return s;
        }
        
    }
    
    void order() /* find order of absolute values of lb */
    {
        int j, k; 
        REAL lj;

        for ( j = 0; j < r; j++ )
        {
            lj = fabs(lb[j]);
            for (k = j - 1; k >= 0; k--)
            {
                if ( lj > fabs(lb[th[k]]) ) th[k + 1] = th[k];
                else goto l1;
            }
            
            k = -1;
            
            l1:
            
            th[k + 1] = j;
        }
        
        ndtsrt = FALSE;
    }
    
    REAL errbd(REAL u, REAL* cx) /*  find bound on tail probability using mgf, cutoff point returned to *cx */
    {
        REAL sum1, lj, ncj, x, y, xconst; 
        int j, nj;
        
        counter();
    
        xconst = u * sigsq;  
        sum1 = u * xconst;  
        u = 2.0 * u;
        
        for (j = r - 1; j >= 0; j--)
        {
             nj = n[j]; lj = lb[j]; ncj = nc[j];
             x = u * lj; y = 1.0 - x;
             xconst = xconst + lj * (ncj / y + nj) / y;
             sum1 = sum1 + ncj * square(x / y) + nj * (square(x) / y + log1(-x, FALSE )); 
        }
        
        *cx = xconst; 
        
        return exp1(-0.5 * sum1); 
    }
    
    
    REAL ctff(REAL accx, REAL* upn) /*  find ctff so that p(qf > ctff) < accx  if (upn > 0, p(qf < ctff) < accx otherwise */
    {
        REAL u1, u2, u, rb, xconst, c1, c2;
    
        u2 = *upn;   u1 = 0.0;  c1 = mean;
        rb = 2.0 * ((u2 > 0.0) ? lmax : lmin);
        
        for (u = u2 / (1.0 + u2 * rb); errbd(u, &c2) > accx; u = u2 / (1.0 + u2 * rb))
        {
            u1 = u2;  c1 = c2;  u2 = 2.0 * u2;
        }
        
        for (u = (c1 - mean) / (c2 - mean); u < 0.9; u = (c1 - mean) / (c2 - mean))
        {
            u = (u1 + u2) / 2.0;
            if (errbd(u / (1.0 + u * rb), &xconst) > accx)
            {  
                u1 = u; c1 = xconst;  
            } else
            {  
                u2 = u;  c2 = xconst; 
            }
        }
        
        *upn = u2; 
        
        return c2;
    }   
    
    REAL truncation(REAL u, REAL tausq) /* bound integration error due to truncation at u */
    {
        REAL sum1, sum2, prod1, prod2, prod3, lj, ncj, x, y, err1, err2;
        int j, nj, s;

        counter();

        sum1  = 0.0; prod2 = 0.0;  prod3 = 0.0;  s = 0;
        sum2 = (sigsq + tausq) * square(u); prod1 = 2.0 * sum2;
        
        u = 2.0 * u;
        
        for (j = 0; j < r; j++ )
        {
            lj = lb[j];  ncj = nc[j]; nj = n[j];
            x = square(u * lj);
            sum1 = sum1 + ncj * x / (1.0 + x);
            
            if (x > 1.0)
            {
                prod2 = prod2 + nj * log(x);
                prod3 = prod3 + nj * log1(x, TRUE );
                s = s + nj;

            } else {
                prod1 = prod1 + nj * log1(x, TRUE );
            }  
        }
        
        sum1 = 0.5 * sum1;
        prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
        
        x = exp1(-sum1 - 0.25 * prod2) / pi; 
        y = exp1(-sum1 - 0.25 * prod3) / pi; 
        
        err1 =  ( s  ==  0 )  ? 1.0 : x * 2.0 / s;
        err2 =  ( prod3 > 1.0 )  ? 2.5 * y : 1.0;
        
        if (err2 < err1) err1 = err2;
        
        x = 0.5 * sum2;
        err2 =  ( x  <=  y )  ? 1.0  : y / x;
        
        return  ( err1 < err2 )  ? err1  :  err2;
    }
    
    void findu(REAL* utx, REAL accx) /*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
    {
        REAL u, ut; int i;
        
        ut = *utx; u = ut / 4.0;
        if ( truncation(u, 0.0) > accx )
        {
            for ( u = ut; truncation(u, 0.0) > accx; u = ut) ut = ut * 4.0;
        } else {
        
            ut = u;
            for ( u = u / 4.0; truncation(u, 0.0) <=  accx; u = u / 4.0 )
            ut = u;
        }
        
        for ( i = 0; i < 4; i++)
        { 
            u = ut / divis[i]; 
            if ( truncation(u, 0.0)  <=  accx )  ut = u; 
        }
        
        *utx = ut;
    }

    REAL cfe(REAL x) /*  coef of tausq in error when convergence factor of exp1(-0.5 * tausq * u ^ 2) is used when df is evaluated at x */
    {
        REAL axl, axl1, axl2, sxl, sum1, lj; int j, k, t;

        counter();
        
        if (ndtsrt) order();
        
        axl = fabs(x);  sxl = (x > 0.0) ? 1.0 : -1.0;  sum1 = 0.0;
        
        for ( j = r - 1; j >= 0; j-- ) { 
            t = th[j];
            if ( lb[t] * sxl > 0.0 )
            {
                lj = fabs(lb[t]);
                axl1 = axl - lj * (n[t] + nc[t]);  axl2 = lj / log28;
                if ( axl1 > axl2 )  axl = axl1; 
                else
                {
                    if ( axl > axl2 )  axl = axl2;
                        sum1 = (axl - axl1) / lj;
                        for ( k = j - 1; k >= 0; k--)
                            sum1 = sum1 + (n[th[k]] + nc[th[k]]);
                        goto  l;
                    }
            }
        }
        
        l:
        
        if (sum1 > 100.0)
        { 
            fail = TRUE; 
            return 1.0; 
        } 
        
        return pow(REAL(2.0), (sum1 / 4.0)) / (pi * square(axl));
    }
    
    void integrate(int nterm, REAL interv, REAL tausq, BOOL mainx) 
    /*  carry out integration with nterm terms, at stepsize
          interv.  if (! mainx) multiply integrand by
             1.0 - exp(-0.5 * tausq * u ^ 2) 
    */
    {
        REAL inpi, u, sum1, sum2, sum3, x, y, z;
        int k, j, nj;

        inpi = interv / pi;
        //std::cout << "** DEBUG: " << interv << " | " << inpi << std::endl;
        
        for ( k = nterm; k >= 0; k--)
        {
            u = (k + 0.5) * interv;
            sum1 = - 2.0 * u * c;  
            sum2 = fabs(sum1);
            sum3 = - 0.5 * sigsq * square(u);

            for ( j = r - 1; j >= 0; j--)
            {
                nj = n[j];  
                x = 2.0 * lb[j] * u;  
                y = square(x);
                sum3 = sum3 - 0.25 * nj * log1(y, TRUE );
                y = nc[j] * x / (1.0 + y);
                z = nj * atan(x) + y;
                sum1 = sum1 + z;   
                sum2 = sum2 + fabs(z);
                sum3 = sum3 - 0.5 * x * y;
            }

            x = inpi * exp1(sum3) / u; 

            if ( !mainx ) {
                x = x * (1.0 - exp1(-0.5 * tausq * square(u))); 
            }

            sum1 = sin(0.5 * sum1) * x;  
            sum2 = 0.5 * sum2 * x;
            
            intl = intl + sum1; 
            
            //std::cout << "DEBUG: " << (intl) << " | " << ersm << " | " << x << " | " << z << std::endl;
            
            ersm = ersm + sum2;
        }
    }
    
public:
    DaviesAlgo () {
    
        if (std::is_same<REAL, double>::value) {
            CUTOFF = -50;
        } else if (std::is_same<REAL, float128>::value)  {
            CUTOFF = -100;
        } else if (std::is_same<REAL, cpp_bin_float_100>::value ) {
            CUTOFF = -300;
        }
        
    }
    
    REAL cdf(double* lb1, double* nc1, int* n1, int r1, double sigma, double c1, int lim1, double acc, double* trace, int* ifault)

    /*  distribution function of a linear combination of non-central
       chi-squared random variables :
    input:
       lb[j]            coefficient of j-th chi-squared variable
       nc[j]            non-centrality parameter
       n[j]             degrees of freedom
       j = 0, 2 ... r-1 # coefficients
       sigma            coefficient of standard normal variable
       c                point at which df is to be evaluated
       lim              maximum number of terms in integration
       acc              maximum error
    output:
       ifault = 1       required accuracy NOT achieved
                2       round-off error possibly significant
                3       invalid parameters
                4       unable to locate integration parameters
                5       out of memory
       trace[0]         absolute sum
       trace[1]         total number of integration terms
       trace[2]         number of integrations
       trace[3]         integration interval in final integration
       trace[4]         truncation point in initial integration
       trace[5]         s.d. of initial convergence factor
       trace[6]         cycles to locate integration parameters     */

    {
        int j, nj, nt, ntm;  
        REAL acc1, almx, xlim, xnt, xntm;
        REAL utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj, ncj;

        REAL qfval = -1.0;

        if (setjmp(env) != 0) { 
          *ifault = 4; goto endofproc; 
        }

        r = r1; lim = lim1; c = REAL(c1);
        n = n1; 
        
        lb = new REAL[r];
        nc = new REAL[r];
        
        for( j = 0; j < r; j++) {
            lb[j] = lb1[j];
            nc[j] = nc1[j];
        }
        
        for ( j = 0; j < 7; j++ ) {
            trace[j] = 0.0;
        }

        *ifault = 0; count = 0;
        intl = 0.0; ersm = 0.0;
        qfval = -1.0; 
        acc1 = acc; 
        ndtsrt = TRUE;  
        fail = FALSE;
        xlim = REAL(lim);

        th = (int*)malloc(r*(sizeof(int)));

        if (!th) { 
            *ifault = 5;  
            goto  endofproc; 
        } 

        /* find mean, sd, max and min of lb, check that parameter values are valid */
        sigsq = square(sigma); 
        sd = sigsq;
       
        lmax = 0.0; 
        lmin = 0.0; 
        mean = 0.0;

        for (j = 0; j < r; j++ )
        {
            nj = n[j];  lj = lb[j];  ncj = nc[j];
            if ( nj < 0  ||  ncj < 0.0 ) { 
                *ifault = 3;  
                goto  endofproc;  
            }
            sd  = sd  + square(lj) * (REAL(2 * nj) + 4.0 * ncj);
            mean = mean + lj * (nj + ncj);

            if (lmax < lj) {
                lmax = lj;
            }
            else if (lmin > lj) {
                lmin = lj;
            }
        }

        //std::cout << "DEBUG: " << lmin << " | " << lmax << " | " <<  mean << " | " << sd << std::endl;

        if ( sd == 0.0  ) {  
            qfval = (c > 0.0) ? 1.0 : 0.0; 
            goto  endofproc;  
        }

        if ( lmin == 0.0 && lmax == 0.0 && sigma == 0.0 ){
            *ifault = 3;  
            goto  endofproc;  
        }

        sd = sqrt(sd);
        almx = (lmax < - lmin) ? - lmin : lmax;

        /* starting values for findu, ctff */
        utx = 16.0 / sd;  up = 4.5 / sd;  un = - up;

        /* truncation point with no convergence factor */
        findu(&utx, .5 * acc1);

        /* does convergence factor help */
        if (c != 0.0  && (almx > 0.07 * sd))
        {
            tausq = .25 * acc1 / cfe(c);
            if (fail) fail = FALSE ;
            else if (truncation(utx, tausq) < .2 * acc1)
            {
                sigsq = sigsq + tausq;
                findu(&utx, .25 * acc1);
                trace[5] = (double)sqrt(tausq);
            }
        }
        
        trace[4] = (double)utx;  
        acc1 = 0.5 * acc1;
        
        //std::cout << "[DAVIES DEBUG](utx): " << utx << std::endl;

        /* find RANGE of distribution, quit if outside this */

        l1:

        d1 = ctff(acc1, &up) - REAL(c);
        if (d1 < 0.0) { 
            qfval = 1.0; 
            goto endofproc; 
        }

        d2 = REAL(c) - ctff(acc1, &un);
        if (d2 < 0.0) { 
            qfval = 0.0; 
            goto endofproc; 
        }

        /* find integration interval */
        intv = 2.0 * pi / ((d1 > d2) ? d1 : d2);
        
        //std::cout << "[DEBUG]: " << intv << " | " << utx << " | "<< d1 << " | " << d2 << " | " << std::endl;
        
        /* calculate number of terms required for main and auxillary integrations */
        xnt = utx / intv;  xntm = 3.0 / sqrt(acc1);

        //std::cout << "[DAVIES DEBUG](xnt): " << xnt << std::endl;

        if (xnt > xntm * 1.5)
        {
            //std::cout << "[DEBUG]: Entering aux int" << std::endl;
            
            /* parameters for auxillary integration */
            if (xntm > xlim) { 
                *ifault = 1; 
                goto endofproc; 
            }

            ntm = (int)floor(xntm + 0.5);
            intv1 = utx / ntm;  
            x = 2.0 * pi / intv1;

            if (x <= fabs(c)) goto l2;

            /* calculate convergence factor */
            tausq = .33 * acc1 / (1.1 * (cfe(c - x) + cfe(c + x)));

            if (fail) goto l2;

            acc1 = .67 * acc1;

            /* auxillary integration */
            integrate(ntm, intv1, tausq, FALSE );
            xlim = xlim - xntm;  sigsq = sigsq + tausq;
            trace[2] = trace[2] + 1; trace[1] = trace[1] + ntm + 1;

            /* find truncation point with new convergence factor */
            findu(&utx, .25 * acc1);  acc1 = 0.75 * acc1;

            goto l1;
        }

        /* main integration */
        l2:

        trace[3] = (double)intv;

        if (xnt > xlim) { 
            *ifault = 1; 
            goto endofproc; 
        }

        nt = (int)floor(xnt + 0.5);

        integrate(nt, intv, 0.0, TRUE );

        trace[2] = trace[2] + 1; 
        trace[1] = trace[1] + nt + 1;
        
        qfval = 0.5 - intl;
        
        trace[0] = (double)ersm;

        /* test whether round-off error could be significant
         allow for radix 8 or 16 machines */
        up = ersm; 
        x = up + REAL(acc) / 10.0;

        for (j = 0; j < 4; j++) { 
            if (rats[j] * x == rats[j] * up) *ifault = 2; 
        }

        endofproc:

        delete(lb);
        delete(nc);
        
        free(th); 
        trace[6] = (double)count;

        return qfval;
    }
};


template <typename REAL> REAL davies(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace) {
    
    DaviesAlgo<REAL> D;
    
    return D.cdf(lb, nc, n, r, 0.0, c, lim, acc, trace, ifault);
}

double onemin_davies(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace);
double onemin_davies_128b(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace);
double onemin_davies_100d(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace);

double constmin_davies(double x,double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace);
double constmin_davies_128b(double x,double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace);
double constmin_davies_100d(double x,double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault, double* trace);
