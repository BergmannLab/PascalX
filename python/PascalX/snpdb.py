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

import sys
import zlib
import gzip
import os.path
import pickle
from sortedcontainers import SortedList

class db:
    """
    Class for handling storage of the raw genotype data. The data is indexed and stored for each chromosome individually as zlib compressed pickle. The indexing allows fast random access via SNP ids or positions.
    
    """
    
    def __init__(self):
        self._modified = False;
        self._idx = [{},{}]
    
        pass
        
    def open(self,filename):
        """
        Opens storage file. A new file is created if not exists.
        
        Args:
            
            filename(string): Name to use for the storage file
        """
        # Load index
        self._filename = filename
        
        if os.path.isfile(filename+".idx.gz"):
            fp = gzip.open(filename+'.idx.gz','rb')
            self._idx = pickle.load(fp)
            fp.close()
        else:
            self._idx = [{},{}]
        
        # open file
        self._datafile = open(filename+".db","a+b")
   
    
    def insert(self,data):
        """
        Stores set of rows into the file storage
        
        Args:
            
            data(dict): First storage key is the key in the dictionary. Second storage key is the first element of the inner list.
            
        Warning:
            If all insert calls are done, the close function has to be called once to make the index persistent.
        """
        self._modified = True;
        
        self._datafile.seek(0,2)
        
        for D in data:
            I = [0,0]
            I[0] = self._datafile.tell()
            self._datafile.write(zlib.compress(pickle.dumps(data[D],protocol=pickle.HIGHEST_PROTOCOL)))
            I[1] = self._datafile.tell()
            
            # Pos based index
            if D not in self._idx[0]:
                self._idx[0][D] = [ I ]
            else:
                self._idx[0][D].append( I )
            
            # Snp based index
            rid = data[D][0].split(";")
            for r in rid:
                if r not in self._idx[1]:
                    self._idx[1][r] = [ I ]
                else:
                    self._idx[1][r].append( I )
            
            
    def get(self,pos):
        """
        Returns all stored data for a set of SNPs indexed via positions
        
        Args:
            
            pos(list): Positions of SNPs to retrieve 
            
        """
        E = []
        for R in pos:
            if R in self._idx[0]:
                for p in self._idx[0][R]:
                    self._datafile.seek(p[0])
                    data = self._datafile.read(p[1]-p[0])

                    E.append( pickle.loads(zlib.decompress(data) ) )
            else:
                #print("Error:",R,"not found in index")
                E.append(None)
        
        return E
    
    def getSNPatPos(self,pos):
        """
        Returns SNP id at position
        
        Args:
            
            pos(list): Positions of SNPs to retrieve 
            
        """
        E = []
        for R in pos:
            if R in self._idx[0]:
                for p in self._idx[0][R]:
                    self._datafile.seek(p[0])
                    data = self._datafile.read(p[1]-p[0])

                    E.append( pickle.loads(zlib.decompress(data) )[0] )
            else:
                #print("Error:",R,"not found in index")
                E.append(None)
        
        return E
    
    def getSNPs(self,snps):
        """
        Returns all stored data for a set of SNPs indexed via SNP ids
        
        Args:
            
            snp(list): ids of SNPs to retrieve 
            
        """
        E = []
        for R in snps:
            if R in self._idx[1]:
                for p in self._idx[1][R]:
                    self._datafile.seek(p[0])
                    data = self._datafile.read(p[1]-p[0])

                    E.append( pickle.loads(zlib.decompress(data) ) )
            else:
                #print("Error:",R,"not found in index")
                E.append(None)
        
        return E
    
    
    def getSNPKeys(self):
        """
        Returns the SNP ids in storage 
        """
        return self._idx[1].keys()
   
    def getKeys(self):
        """
        Returns SNP positions in storage
        """
        return self._idx[0].keys()
   
    def getSortedKeys(self):
        """
        Returns a sorted list of SNP positions in storage
        """
        return SortedList(self._idx[0].keys())
    
    def getSNPpos(self,snpid):
        """
        Returns the position corresponding to a snpid
        
        Warning: 
            The current implementation is inefficient.
        """
        if snpid in self._idx[1]:
            fseek = self._idx[1][snpid][0]
        
            for pos, fs in self._idx[0].items(): 
                for i in range(0,len(fs[0])):
                    if fseek[0] == fs[0][i]:
                        return pos

        return None
    
    def close(self):
        """
        Closes open storage file. 
        
        Warning:
            After all inserts are done this function has to be called once to re-generate the index and close the storage file.
        """
        if self._modified:
            fp = gzip.open(self._filename+'.idx.gz','wb')
            pickle.dump(self._idx, fp, protocol=pickle.HIGHEST_PROTOCOL)
            fp.close()
        
        self._idx = [{},{}]
        
        # close 
        self._datafile.close()
        