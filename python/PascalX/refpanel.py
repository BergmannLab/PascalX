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

from PascalX import snpdb

from tqdm import tqdm

import multiprocessing as mp

import os.path
import sys

import gzip
import numpy as np

class refpanel:
    
    def __init__(self):
        pass
    
    def load_pos_reference(self,cr):
        """
        Returns a snpdb object for a chromosome and a sorted list of SNP positions on the chromosome
        
        Args:
        
            cr(int): Chromosome number
        
        """
        db = snpdb.db()
        db.open(self._refData+'.chr'+str(cr))
        
        return [db,db.getSortedKeys()]
    
    def load_snp_reference(self,cr):
        """
        Returns a snpdb object for a chromosome
        
        Args:
            
            cr(int): Chromosome number
        
        """
        db = snpdb.db()
        db.open(self._refData+'.chr'+str(cr))
        
        return db
    
    def set_refpanel(self,filename, parallel=1):
        """
        Sets the reference panel to use
        
        Args:
            
            filename(string): /path/filename (without .chr#.db ending)
            parallel(int): Number of cores to use for parallel import of reference panel
            
        Note:
            One file per chromosome with ending .chr#.db required (#: 1-22). If imported reference panel is not present, he will automatically try to import from .chr#.tped.gz files.        
        """
        self._refData = filename
        
        NF = []
        for i in range(1,23):
            if not os.path.isfile(filename+".chr"+str(i)+".idx.gz") or not os.path.isfile(filename+".chr"+str(i)+".db"):
                NF.append(i)
        
        # Import if missing
        if len(NF) > 0:
            print("Reference panel data not imported. Trying to import...")
            self._import_reference(chrs=NF,parallel=parallel)
            
    def _import_reference_thread(self,i):
        
        # Load
        with gzip.open(self._refData+'.chr'+str(i)+'.tped.gz','rt') as f:
            
            db = snpdb.db()
            db.open(self._refData+'.chr'+str(i))
            
            for line in f:
                L = line.split() # [chr,rid,irrelevant,pos,genotype]

                # Get genotype and dephase
                if L[1][0:2] == 'rs':
                    dephased = (np.array(L[4:],dtype='int32') - 1)
                    genotype = np.sum(dephased.reshape((int(len(dephased)/2),2)),axis=1)
                    m = np.mean(genotype)
                    s = np.std(genotype)

                    if s!=0:
                        # Compute MAF
                        MAF = m/2.
                        if (MAF > 0.5):
                            MAF = 1.0 - MAF;

                        T = [L[1],MAF,(genotype-m)/s]

                        # Store                             
                        db.insert({int(L[3]):T})
                                  
            db.close()
            
        return True
        
    def _import_reference(self,chrs=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],parallel=1):
        """
        Imports reference data from .tped.gz files.
        (Has only to be run once. The imported data is stored on disk for later usage.)
        
        chrs    : List of chromosomes to import
        parallel: # of cores to use (WARNING: Take care that you have sufficient memory!)
        
        """
        
        # Check if ref exist
        for i in range(1,23):
            if not os.path.isfile(self._refData+".chr"+str(i)+".tped.gz"):
                print("ERROR: ", self._refData+".chr"+str(i)+".tped.gz", "not found")   
                return

            
        # Start import    
        pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
           
        with tqdm(total=len(chrs), desc="Importing reference panel", bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]",file=sys.stdout) as pbar:
        
            def update(*a):
                pbar.update(1)

            res = []
            for i in chrs:
                res.append(pool.apply_async(self._import_reference_thread, args=(i,), callback=update))
                
            # Wait to finish
            for r in res:
                r.get()

        pool.close()
        
        
    def getSNPtoChrMap(self):
        """
        Returns a dictionary mapping SNP id to corresponding chromosome number
        """
        MAP = {}
        for cr in range(1,23):
            db = self.load_snp_reference(cr)
            snps = db.getSNPKeys()
            
            for S in snps:
                MAP[S] = cr
    
        return MAP
    
    def getChrSNPs(self,cr):
        """
        Returns SNP ids for a chromosome
        
        Args:
            
            cr(int): Chromosome number
            
        """
        db = self.load_snp_reference(cr)
       
        return db.getSNPKeys()
       