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

from fastnumbers import int

import re

class refpanel:
    
    def __init__(self):
        pass
    
    def load_pos_reference(self,cr,keep_idx = None):
        """
        Returns a snpdb object for a chromosome and a sorted list of SNP positions on the chromosome
        
        Args:
        
            cr(int): Chromosome number
        
        """
        if keep_idx is None:
            db = snpdb.db()
        else:
            db = snpdb.db_subset(keep_idx)
            
        db.open(self._refData+'.chr'+str(cr))
        
        return [db,db.getSortedKeys()]
    
    def load_snp_reference(self,cr,keep_idx = None):
        """
        Returns a snpdb object for a chromosome
        
        Args:
            
            cr(int): Chromosome number
        
        """
        if keep_idx is None:
            db = snpdb.db()
        else:
            db = snpdb.db_subset(keep_idx)
            
        db.open(self._refData+'.chr'+str(cr))
        
        return db
    
    def set_refpanel(self,filename, parallel=1, keepfile=None, qualityT=100, SNPonly=False, chrlist=None, sourcefilename=None,regEx=None,nobar=True):
        """
        Sets the reference panel to use
        
        Args:
            
            filename(string): /path/filename (without .chr#.db ending)
            parallel(int): Number of cores to use for parallel import of reference panel
            
            keepfile(string): [only for .vcf] File with sample ids (one per line) to keep.  None to keep all.
            qualityT(int): [only for .vcf] Quality threshold for variant to keep (None to ignore)
            SNPonly(bool): [only for .vcf] Load only SNPs 
            chrlist(list): List of chromosomes to import. (None to import 1-22)
            sourcefilename(string): /path/filename (without .chr#. ending) of .tped | .vcf files. None to use same as filename
            regEx(string): Regular expression to filter sample ids. First capture group is kept. [only for .vcf]
            nobar(bool): Show progress bar (updates only if a chromosome finished)
            
        Note:
        
            One file per chromosome with ending .chr#.db required (#: 1-22). If imported reference panel is not present, PascalX will automatically try to import from .chr#.tped.gz or .chr#.vcf.gz files.
            
        Note:
        
            Alleles (under .vcf import) are stored internally in the order [ALT,REF].
            
        """
        self._refData = filename
        self._srcData = sourcefilename
        
        if chrlist is None:
            NF = []
            for i in range(1,23):
                if not os.path.isfile(filename+".chr"+str(i)+".idx.gz") or not os.path.isfile(filename+".chr"+str(i)+".db"):
                    NF.append(i)
        else: 
            NF = []
            for i in chrlist:
                if not os.path.isfile(filename+".chr"+str(i)+".idx.gz") or not os.path.isfile(filename+".chr"+str(i)+".db"):
                    NF.append(i)
            
        # Import if missing
        if len(NF) > 0:
            print("Reference panel data not imported. Trying to import...")
            self._import_reference(chrs=NF,parallel=parallel,keepfile=keepfile,qualityT=qualityT,SNPonly=SNPonly,regEx=regEx,nobar=nobar)
            
    def _import_reference_thread_tped(self,i):
        
        # Load
        if self._srcData is None:
            fn = self._refData
        else:
            fn = self._srcData
            
        with gzip.open(fn+'.chr'+str(i)+'.tped.gz','rt') as f:
            
            db = snpdb.db()
            db.open(self._refData+'.chr'+str(i))
            
            for line in f:
               
                L = line.split() # [chr,rid,irrelevant,pos,genotype]
            
                # Get genotype and dephase
                if L[1][0:2] == 'rs':
                    
                    dephased = (np.array(L[4:],dtype='b') - 1)
                    genotype = np.sum(dephased.reshape((int(len(dephased)/2),2)),axis=1)
                    
                    # PLINK uses 1 for minor allele -> Convert to minor allele count
                    genotype[genotype==2] = -1
                    genotype[genotype==0] = 2
                    genotype[genotype==-1] = 0
                    
                    m = np.mean(genotype)
                    s = np.std(genotype)
                    
                    if s!=0:
                        # Compute MAF
                        MAF = m/2.
                        if (MAF > 0.5):
                            MAF = 1.0 - MAF;

                        T = [L[1],round(MAF,3),genotype]
                    
                        # Store                             
                        db.insert({int(L[3]):T})
                                  
            db.close()
            
        return True
        
    def _import_reference_thread_vcf(self,i,keepfile,qualityT,SNPonly,regEx=None):
        # Pre-compile genotype splitter
        Gsplitter = re.compile('/|\|')
      
        # Load filter info
        keep = set([])
        if keepfile is not None:
            f = open(keepfile,'r')
            for line in f:
                S = line.split("\t")[0]
                keep.add(S)

            f.close()
        
        sampleMap = {}
        
        # Load
        if self._srcData is None:
            fn = self._refData
        else:
            fn = self._srcData
            
        with gzip.open(fn+'.chr'+str(i)+'.vcf.gz','rt') as f:
            
            # Find header
            for line in f:
                # Detect infos and headers
                if line[:2] == "##":
                    continue

                # Detect sample names
                if line[:2] == "#C":
                    data = line.split("\t")
                    tmp = data[9:]
                    
                    # RegEx processing of sample names
                    if regEx is not None:
                        for j in range(0,len(tmp)):
                            m = re.search(regEx,tmp[j])
                            try:
                                tmp[j] = m.group(1)
                            except:
                                continue
                            
                    for j in range(0,len(tmp)):
                        if (keepfile is None) or (tmp[j] in keep):
                            sampleMap[j] = tmp[j]

                    break

            sampleKeys = list(sampleMap.keys())
            
            # Store sample keys
            with open(self._refData+'.sampleIds.txt','wt') as g:
                s = ""
                for j in range(0,len(sampleKeys)-1):
                    s += str(sampleMap[sampleKeys[j]])+"\t"
                
                s += str(sampleMap[sampleKeys[len(sampleKeys)-1]])
                
                g.write(s+'\n')
                
                
            db = snpdb.db()
            db.open(self._refData+'.chr'+str(i))
            
            # Main data import loop
            for line in f:
                try:
                    # Data line
                    data = line.split("\t")
    
                    # Get GT pos
                    tmp = data[8].split(":")
                    GT = -1
                    for j in range(0,len(tmp)):
                        if tmp[j] == 'GT':
                            GT = j
                            break
    
                    # Checks
                    if (GT == -1) or (data[2][:2] != 'rs') or (data[6] != 'PASS' and qualityT is not None and (int(data[5]) < qualityT)):
                        continue
    
                    # Infer alternate alleles (pos 0: ref allele)
                    alleles = [data[3]]
                    alleles.extend(data[4].split(","))
                    
                    if SNPonly and (len(data[3]) > 1):
                        continue
                                   
                    counter = np.zeros(len(alleles),dtype='int') 
                    
                    # Read genotype
                    genotypes = data[9:]

                    # Multi-allelic flag
                    mallelic = False
                    
                    # Only read samples in sampleMap
                    genomap = {}
                    for j in range(0,len(sampleKeys)):
                        
                        geno = Gsplitter.split( genotypes[sampleKeys[j]].split(":")[GT] )

                        if len(geno) > 2:
                            mallelic = True
                            break
                            
                        # Ignore half-calls
                        if geno[0] != "." and geno[1] != ".":
                            counter[int(geno[0])] += 1
                            counter[int(geno[1])] += 1
                        
                        genomap[sampleKeys[j]] = geno

                    if not mallelic: # Skip if more than two alleles 
                        # Reference allele
                        refp = 0
                        
                        SC = np.argsort(counter) # Sort alleles count
                        for p in SC:
                            
                            if p != refp:
                                if SNPonly and len(alleles[p]) > 1:
                                    continue
        
                                minp = str(p)
        
                                gd = np.zeros(len(sampleKeys),dtype='B')
        
                                for j in range(0,len(sampleKeys)):
                                    geno = genomap[sampleKeys[j]]
        
                                    # Ignore half-calls
                                    if geno[0] != '.' and geno[1] != '.':
        
                                        if geno[0] == minp:
                                            gd[j] += 1
        
                                        if geno[1] == minp:
                                            gd[j] += 1
        
                                # Compute MAF
                                MAF = np.mean(gd)/2.
                                if (MAF > 0.5):
                                    MAF = 1.0 - MAF;
                                
                                T = [data[2],MAF,gd,alleles[p],alleles[refp]] # Stores alt and ref allele
        
                                db.insert({int(data[1]):T})
                
                except Exception as e:
                    print("[WARNING] Data line broken:",line)   
                    
            f.close()
            db.close()
            
        return True
        

    def _import_reference(self,chrs=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],parallel=1,keepfile=None,qualityT=100,SNPonly=False,regEx=None,nobar=True):
        """
        Imports reference data from .tped.gz or .vcf.gz files.
        (Has only to be run once. The imported data is stored on disk for later usage.)
        
        chrs    : List of chromosomes to import
        parallel: # of cores to use (WARNING: Take care that you have sufficient memory!)
        keepfile: File with sample ids (one per line) to keep (only for .vcf) 
        qualityT: Quality threshold for variant to keep (only for .vcf) (None to ignore)
        
        Warning: 
            Direct .vcf import is currently only experimental !
        
        """
        
        # Check if ref exist
        if self._srcData is None:
            fn = self._refData
        else:
            fn = self._srcData
            
        for i in chrs:
            if not os.path.isfile(fn+".chr"+str(i)+".tped.gz") and not os.path.isfile(fn+".chr"+str(i)+".vcf.gz"):
                print("ERROR: ", fn+".chr"+str(i)+".(tped|vcf).gz", "not found")   
                return

            if os.path.isfile(fn+".chr"+str(i)+".tped.gz"):
                cmd = 'tped'
            else:
                cmd = 'vcf'
            
        # Start import    
        pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
           
        with tqdm(total=len(chrs), desc="Importing reference panel", bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]",file=sys.stdout,disable=nobar) as pbar:
        
            def update(*a):
                pbar.update(1)

            res = []
            for i in chrs:
                if cmd == 'tped':
                    res.append(pool.apply_async(self._import_reference_thread_tped, args=(i,), callback=update))
                elif cmd == 'vcf':
                    res.append(pool.apply_async(self._import_reference_thread_vcf, args=(i,keepfile,qualityT,SNPonly,regEx), callback=update))
            
                
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
       