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

extern double oneminwchissum_m1nc0_davies(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double oneminwchissum_m1nc0_davies_128b(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double oneminwchissum_m1nc0_davies_100d(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double oneminwchissum_m1nc0_davies_auto(double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);

extern double constminwchissum_m1nc0_davies(double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double constminwchissum_m1nc0_davies_128b(double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double constminwchissum_m1nc0_davies_100d(double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double fconstminwchissum_m1nc0_davies_auto(double F, double x, double* lambda, int N, double X, int lim, double acc, int* ifault, double* trace);

extern double oneminwchissum_m1nc0_ruben(double* lambda, int N, double X, int lim, double acc, int* ifault);
extern double oneminwchissum_m1nc0_ruben_128b(double* lambda, int N, double X, int lim, double acc, int* ifault);
extern double oneminwchissum_m1nc0_ruben_100d(double* lambda, int N, double X, int lim, double acc, int* ifault);
extern double oneminwchissum_m1nc0_ruben_200d(double* lambda, int N, double X, int lim, double acc, int* ifault);

extern double oneminwchissum_m1nc0_auto(double* lambda, int N, double X, int lim, double acc, int* ifault);

extern double oneminwchissum_m1_davies(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double oneminwchissum_m1_davies_128b(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double oneminwchissum_m1_davies_100d(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace);
extern double oneminwchissum_m1_davies_auto(double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace);


extern double constminwchissum_m1_davies(double x,double* lambda, double* nc, int N, double X, int lim, double acc, int* ifault, double* trace);
