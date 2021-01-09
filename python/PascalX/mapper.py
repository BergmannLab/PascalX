#    PascalX - A python3 library for high precision gene and pathway scoring for 
#              GWAS summary statistics with C++ backend.
#              https://github.com/BergmannLab/PascalX
#
#    Copyright (C) 2021 Bergmann lab and contributors
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import gzip

class mapper:
    
    def __init__(self):
        self._GENEIDtoSNP = {}
        
    
    def load_mapping(self,file,gcol=0,rcol=1,ocol=3,splitchr="\t",header=False):
        """
        Loads a SNP to gene mapping
        
        Args:
            
            file(string): File to load
            gcol(int): Column with gene id
            rcol(int): Column with SNP id
            ocol(int): Column with weight
            splitchr(string): Character used to separate columns
            header(bool): Header present
            
        """
        self._GENEIDtoSNP = {}
        with gzip.open(file,'rt') as f:
            if header:
                f.readline()
            
            c = 0
            for line in f:
                line = line.rstrip().split(splitchr)
                
                gid = line[gcol]
                
                if gid not in self._GENEIDtoSNP:
                    self._GENEIDtoSNP[gid] = [[],[]]
                    
                    c = c + 1
                    
                self._GENEIDtoSNP[gid][0].append(line[rcol])
                self._GENEIDtoSNP[gid][1].append(int(line[ocol]))
                
            print(c,"gene to SNP mappings loaded")
        