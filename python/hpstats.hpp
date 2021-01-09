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

extern double invchi2cdf_1mx(double x, int dof);
extern double invchi2cdf_1mx_128b(double x, int dof);
extern double invchi2cdf_1mx_100d(double x, int dof);

extern double onemin_chi2cdf(double x, int dof);
extern double onemin_chi2cdf_128b(double x, int dof);
extern double onemin_chi2cdf_100d(double x, int dof);

extern double oneminwchissum_m1nc0_satterthwaite_100d(double* lambda, int N, double X);
extern double oneminwchissum_m1nc0_satterthwaite_200d(double* lambda, int N, double X);
