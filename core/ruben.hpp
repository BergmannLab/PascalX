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

template <typename REAL> REAL ruben(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, REAL* dnsty, int* ifault);

double onemin_ruben(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault);
double onemin_ruben_128b(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault);
double onemin_ruben_100d(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault);
double onemin_ruben_200d(double* lb, int* n, double* nc, int r, double c, int lim, double acc, int* ifault);

extern void ruben_128b_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault);
extern void onemin_ruben_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault);
extern void onemin_ruben_128b_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault);
extern void onemin_ruben_100d_str(double* lambda, int* mult, double* delta, int n, double c, double mode, int maxit, double eps, char* dnsty, int* ifault);

