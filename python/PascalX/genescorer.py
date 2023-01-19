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

from PascalX import wchissum,snpdb,tools,refpanel,genome,hpstats
from PascalX.mapper import mapper

import gzip
import pickle

import numpy as np

#from tqdm import tqdm
from tqdm.auto import tqdm

import multiprocessing as mp

import matplotlib.pyplot as plt

import seaborn as sns

import sys
import os.path

from scipy.stats import norm

from sortedcontainers import SortedSet

from abc import ABC

import time

try:
    import cupy as cp
    pool = cp.cuda.MemoryPool(cp.cuda.malloc_managed)
    cp.cuda.set_allocator(pool.malloc)
except ModuleNotFoundError:
    cp = None

    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class genescorer(ABC):
    """
    
    Genescorer base class
    
    """
    
    # Mapper vars are static
    _MAP = None
    _iMAP = None
   
    
    def __init__(self):
        pass
    
    
    def load_refpanel(self, filename, parallel=1,keepfile=None,qualityT=100,SNPonly=False,chrlist=None):
        """
        Sets the reference panel to use
        
        Args:
        
            filename(string): /path/filename (without .chr#.db ending)
            parallel(int): Number of cores to use for parallel import of reference panel
            
            keepfile: File with sample ids (one per line) to keep (only for .vcf) 
            qualityT: Quality threshold for variant to keep (only for .vcf)
            SNPonly : Import only SNPs (only for .vcf)
            chrlist(list): List of chromosomes to import. (None to import 1-22)
            
        Note:
        
            One file per chromosome with ending .chr#.db required (#: 1-22). If imported reference panel is not present, PascalX will automatically try to import from .chr#.tped.gz or .chr#.vcf.gz files.
            
        """
        self._ref = refpanel.refpanel()
        self._ref.set_refpanel(filename=filename,parallel=parallel,keepfile=keepfile,qualityT=qualityT,SNPonly=SNPonly,chrlist=chrlist)

    
    def load_genome(self,file,ccol=1,cid=0,csymb=5,cstx=2,cetx=3,cs=4,cb=None,chrStart=0,splitchr='\t',NAgeneid='n/a',useNAgenes=False,header=False):
        """
        Imports gene annotation from text file

        Args:
            file(text): File to load
            ccol(int): Column containing chromosome number
            cid(int): Column containing gene id
            csymb(int): Column containing gene symbol
            cstx(int): Column containing transcription start
            cetx(int): Column containing transcription end
            cs(int): Column containing strand
            cb(int): Column containing band (None if not supplied)
            chrStart(int): Number of leading characters to skip in ccol
            splitchr(string): Character used to separate columns in text file
            NAgeneid(string): Identifier for not available gene id
            useNAgenes(bool): Import genes without gene id
            header(bool): First line is header


        **Internal:**
            * ``_GENEID`` (dict) - Contains the raw data with gene id (cid) as key
            * ``_GENESYMB`` (dict) - Mapping from gene symbols (csymb) to gene ids (cid)
            * ``_GENEIDtoSYMB`` (dict) - Mapping from gene ids (cid) to gene symbols (csymb)
            * ``_CHR`` (dict) - Mapping from chromosomes to list of gene symbols
            * ``_BAND`` (dict) - Mapping from band to list of gene symbols
            * ``_SKIPPED`` (dict) - Genes (cid) which could not be imported

        Note:
           An unique gene id is automatically generated for n/a gene ids if ``useNAgenes=true``.

        Note:
           Identical gene ids occuring in more than one row are merged to a single gene, if on same chromosome and positional gap is < 1Mb. The strand is taken from the longer segment.    

        """
        GEN = genome.genome()
        GEN.load_genome(file,ccol,cid,csymb,cstx,cetx,cs,cb,chrStart,splitchr,NAgeneid,useNAgenes,header)
        
        self._GENOME = GEN
        self._GENEID = GEN._GENEID
        self._GENESYMB = GEN._GENESYMB
        self._GENEIDtoSYMB = GEN._GENEIDtoSYMB
        self._CHR = GEN._CHR
        self._BAND = GEN._BAND
        
        self._SKIPPED = GEN._SKIPPED

    
    def load_mapping(self,file,gcol=0,rcol=1,wcol=None,delimiter="\t",a1col=None,a2col=None,bcol=None,pcol=None,pfilter=1,header=False,joint=True,symbol=False):
        """
        Loads a SNP to gene mapping
        
        Args:
            
            file(string): File to load
            gcol(int): Column with gene id
            rcol(int): Column with SNP id
            wcol(int): Column with weight
            a1col(int): Column of alternate allele (None for ignoring alleles)
            a2col(int): Column of reference allele (None for ignoring alleles)
            bcol(int): Column with additional weight
            delimiter(string): Character used to separate columns
            pfilter(float): Only include rows with wcol < pfilter
            header(bool): Header present
            joint(bool): Use mapping SNPs and gene window based SNPs
            symbol(bool): True: Gene id is a gene symbol; False: Gene id is an ensembl gene id
            
        Note:
        
            * A loaded mapping takes precedence over a loaded positional gene annotation
            * The mapping data is stored statically (same for all class initializations)
            
        """
        if symbol and self._GENOME is None:
            print("For symbol==True a genome has to be loaded first (use .load_genome)")
        else:
            M = mapper(self._GENOME)
            M.load_mapping(file,gcol,rcol,wcol,a1col,a2col,bcol,pcol,delimiter,pfilter,header,symbol)
            self._MAP = M._GENEIDtoSNP
            self._iMAP = M._SNPtoGENEID
            self._joint = joint

    def load_GWAS(self,file,rscol=0,pcol=1,bcol=None,a1col=None,a2col=None,delimiter=None,header=False,NAid='NA',log10p=False,cutoff=1e-300):
        """
        Load GWAS summary statistics p-values
        
        Args:
        
            file(string): File containing the GWAS summary statistics data. Either as textfile or gzip compressed with ending .gz
            rscol(int): Column of SNP ids
            pcol(int) : Column of p-values
            bcol(int) : Column of betas (optional)
            a1col(int): Column of alternate allele (None for ignoring) 
            a2col(int): Column of reference allele (None for ignoring)
            delimiter(String): Split character 
            header(bool): Header present
            NAid(String): Code for not available (rows are ignored)
            log10p(bool): p-values are given -log10 transformed 
            cutoff(float): Cutoff value for p-values (None for no cutoff)
        """
        self._GWAS = {}
        self._GWAS_beta = {}
        self._GWAS_alleles = {}
        
        if file[-3:] == '.gz':
            f = gzip.open(file,'rt')
        else:
            f = open(file,'rt')
            
        if header:
            f.readline()
        
        c = 0
        for line in f:
            if delimiter is None:
                L = line.split()
            else:
                L = line.split(delimiter)

            if L[pcol] != NAid and L[rscol] != NAid:
                p = float(L[pcol])
                
                # Threshold very small SNPs (1e-300 recommended for numerical stability)
                if cutoff is not None and p < cutoff:
                    p = cutoff
                    c +=1
                    
                if log10p:
                    p = 10**(-p)
                    
                if p > 0 and p < 1 and L[rscol][:2]=='rs':
                    self._GWAS[L[rscol]] = p
                else:
                    continue
            else:
                continue
                    
            if bcol is not None:
                b = float(L[bcol])
                self._GWAS_beta[L[rscol]] = b
              
            if a1col is not None and a2col is not None:
                self._GWAS_alleles[L[rscol]] = [L[a1col].upper(),L[a2col].upper()]
            
        f.close()

        if c > 0:
            print(c,"SNPs cutoff to",cutoff)
                
        print(len(self._GWAS),"SNPs loaded")
        
    def matchAlleles(self, SNPonly=False):
        """
        Matches alleles between loaded GWAS and reference panel 
        (SNPs with non matching alleles are removed)
        
        Args:
        
            SNPonly(bool) : Keep only SNPs
            
        """
        db = {}
        for i in range(1,23):
            db[i] = self._ref.load_snp_reference(i)
    
        todel = []
        Nr = 0
        Ns = 0
        N = 0
        for x in self._GWAS_alleles:
            N += 1
            found = False
            
            A = self._GWAS_alleles[x]
            if SNPonly == False or (len(A[0])==1 and len(A[1])==1):  
            
                for c in db:
                    snp = db[c].getSNPs([x])
                    if len(snp) > 0:
                        for s in snp:
                            if s is not None and [s[3],s[4]] == A: 
                                found = True
                                break

                    if found:
                        break
            else:
                Ns += 1
                
            if not found:
                # Remove if alleles do not match to ref panel
                Nr += 1
                todel.append(x)
            
        
        del db
    
        for x in todel:
                
            del self._GWAS_alleles[x]
            
            if x in self._GWAS:
                del self._GWAS[x]
                
            if x in self._GWAS_beta:
                del self._GWAS_beta[x]

            
        print(N,"GWAS SNPs")
        
        if SNPonly:
            print(Ns,"Non-SNPs removed")
            
        if (N-Ns) != 0:
            print(round(Nr/(N-Ns)*100,2),"% non-matching with ref panel -> ",Nr,"SNPs removed")    
       
    def rank(self):
        """
        QQ normalizes the p-values of loaded GWAS
        
        """
        SNPs = list(self._GWAS.keys())
     
        pA = np.ones(len(SNPs))
        
        for i in range(0,len(SNPs)):
            pA[i] = self._GWAS[SNPs[i]]
         
        self._GWAS = {}
        
        # Rank
        p = np.argsort(pA)
        wr = np.zeros(len(p))

        for i in range(0,len(p)):
            wr[p[i]] = (i+1.)/(len(p)+1.) 

        for i in range(0,len(SNPs)):
            self._GWAS[SNPs[i]] = wr[i]

        print(len(SNPs),"SNPs ( min p:", f'{1./(len(p)+1):.2e}',")")
    
    
    def rank_mapper(self):
        """
        QQ normalizes the p-values of loaded Mapper
        
        """
        SNPs = np.array(list(self._iMAP.keys()))
        
        # Build up data from mapper
        pA = []
        for i in range(0,len(SNPs)):
            G = self._iMAP[SNPs[i]]
            for j in range(0,len(G)):
                R = self._MAP[G[j]]    
                if SNPs[i] in R.keys():
                    pA.append(R[SNPs[i]][0])
                    
        # Rank
        rA = np.argsort(pA)
        map_A = {}
        for i in range(0,len(rA)):
            map_A[rA[i]] = (i+1.)/(len(rA)+1.) 
        
        # Set data in mapper
        # Generate new dataset based on hashmaps
        c = 0
        for i in range(0,len(SNPs)):
            G = self._iMAP[SNPs[i]]
            for j in range(0,len(G)):
                R = self._MAP[G[j]]
                if SNPs[i] in R.keys():     
                    R[SNPs[i]][0] = map_A[c]
                    c += 1
                      
        print(len(SNPs),"shared SNPs ( min p:", f'{1./(len(rA)+1):.2e}',")")
       
    
    
    
    def save_GWAS(self,file):
        """
        Save GWAS p-values
        
        Args:
        
            file(string): Filename to store data
        """
        f = open(file,"w")
        
        for K in self._GWAS.keys():
            f.write(str(K)+"\t"+str(self._GWAS[K])+"\n")
            
        f.close()
    
    def save_scores(self,file):
        """
        Save computed gene scores 
        
        Args:
            
            file(string): Filename to store data
        
        """
        f = open(file,"w")
        
        for K in self._SCORES.keys():
            f.write(str(K)+"\t"+str(self._SCORES[K])+"\n")
            
        f.close()
    
    def load_scores(self,file,gcol=0,pcol=1,header=False):
        """
        Load computed gene scores
        
        Args:
            
            file(string): Filename of data to load
            gcol(int): Column with gene symbol
            pcol(int): Column with p-value
            header(bool): File contains a header (True|False)
        """
        self._SCORES = {}
        
        f = open(file,"r")
        if header:
            f.readline()
            
        for line in f:
            L = line.rstrip('\n').split("\t")
            self._SCORES[L[gcol]] = float(L[pcol])
            
        f.close()
        
        print(len(self._SCORES),"scores loaded")
    
    
    
    def get_topscores(self,N=10):
        """
        Prints and returns the top gene scores
        
        Args:
        
            N(int): # to show
            
        Returns:
        
            list: Ordered list of top scores
            
        """
        K = []
        V = []
        
        for key, value in self._SCORES.items():
            K.append(key)
            V.append(value)
    
        I = np.argsort(V)[:N]
        
        R = []
        
        for i in I:
            print(K[i]," ",V[i])
            R.append([K[i],V[i]])
        
        return R
    
    def get_geneinfo(self,gene):
        """
        Shows details of the loaded annotation for a gene
        
        Args:
        
            gene(String): Gene symbol
            
        """
        if gene in self._GENESYMB:
            G = self._GENEID[self._GENESYMB[gene]]
            return G
        else:
            print(gene,"not in annotation")
            return None
        
    
    def plot_Manhattan(self,ScoringResult=None,region=None,sigLine=0,logsigThreshold=0,labelSig=True,labelList=[],style='colorful'):
        """
        Produces a Manhattan plot
        
        Args:
        
            ScoringResult(list): List of gene,p-value pairs [['Gene',p-value],...] (if None, generate from internal _SCORES)
            region(str): Band region to plot
            sigLine(float): Draws a horizontal line at p-value (if > 0)
            logsigThreshold(float): Significance threshold above which to label genes (log(p) values)
            labelSig(bool): Label genes above significance threshold
            labelList(list): List of gene names to label
            style(string): Design of the plot ('classic' | 'colorful')
            
        """
        
        if ScoringResult is None:
            ScoringResult = []

            for K in self._SCORES:
                ScoringResult.append([K,self._SCORES[K]])
        
        
        A = np.array(ScoringResult)
        
        G = []
        for g in range(0,len(A[:,0])):
            if A[g,0] in self._GENESYMB:
                
                if region is None:
                    G.append(self._GENESYMB[A[g,0]])
                else:
                    if A[g,0] in self._BAND[region]:
                        G.append(self._GENESYMB[A[g,0]])
                    else:
                        G.append(None)
            else:
                print("[WARNING]: "+A[g,0]+" not in annotation -> ignoring")
                G.append(None)
                
        p = -np.log10(A[:,1].astype('float64'))
        
        # Calc positions
        Lc1 = []
        Lp1 = []
        
        Lc2 = []
        Lp2 = []
        
        Lc3 = []
        Lp3 = []
        Ln3 = []
        chrs = []
       
        for i in range(0,len(G)):
            if G[i] is not None:
                cr = self._GENEID[G[i]][0]
                
                pos = self._CHR[cr][3] - self._CHR[cr][1] + 0.5*(self._GENEID[G[i]][1] + self._GENEID[G[i]][2]) # <- ok   
                
                #print("[DEBUG]:",G[i],self._GENEID[G[i]],pos, self._CHR[cr][3],self._CHR[cr][1], self._CHR[cr][3] - self._CHR[cr][1], (self._GENEID[G[i]][1] + self._GENEID[G[i]][2])   )
                
                if self._GENEID[G[i]][0] not in chrs:
                    chrs.append(self._GENEID[G[i]][0])
                
                symb = self._GENEIDtoSYMB[G[i]]
                
                if( (logsigThreshold > 0 and p[i] > logsigThreshold) or (symb in labelList) ):
                        Lc3.append(pos)
                        Lp3.append(p[i])
                        Ln3.append(symb)
                else:
                    if self._GENEID[G[i]][0] != 'X' and self._GENEID[G[i]][0] != 'Y' and self._GENEID[G[i]][0] != 'M':
                        if int(self._GENEID[G[i]][0]) % 2 == 0:
                            Lc1.append(pos)
                            Lp1.append(p[i])
                        else:
                            Lc2.append(pos)
                            Lp2.append(p[i])
                    else:
                        if self._GENEID[G[i]][0] != 'M' or self._GENEID[G[i]][0] != 'Y':
                            Lc1.append(pos)
                            Lp1.append(p[i])    
                        else:
                            Lc2.append(pos)
                            Lp2.append(p[i])

        # Calc xticks
        chrs = np.sort(chrs) 
        chrp = np.zeros(len(chrs))
        for i in range(0,len(chrs)):
            chrp[i] = self._CHR[chrs[i]][3] - self._CHR[chrs[i]][1]  + 0.5*(self._CHR[chrs[i]][1] + self._CHR[chrs[i]][2]) # Tick seems ok: chr start pos + (last tx - first tx)/2
            
       
        # Plot
        plt.ylabel("$-log_{10}(p)$")
        
        if not region:
            plt.xlabel("chr")
        else:
            plt.xlabel("chr"+region)
            
        plt.xticks(chrp, chrs) 
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        if style=='classic':
            plt.scatter(Lc1,Lp1,color='silver')
            plt.scatter(Lc2,Lp2,color='gray')
            plt.scatter(Lc3,Lp3,color='red', marker='v', s=30)
        else:
            Lx = []
            Ly = []
            
            if(len(Lc1) > 0):
                Lx.extend(Lc1)
                Ly.extend(Lp1)
            
            if(len(Lc2) > 0):
                Lx.extend(Lc2)
                Ly.extend(Lp2)
            
            if(len(Lc3) > 0):
                Lx.extend(Lc3)
                Ly.extend(Lp3)
        
            plt.scatter(Lx,Ly,c=Ly, cmap="coolwarm")
            #plt.colorbar()
            
        if labelSig:
            ax = plt.gca()  
       
            for i in range(0,len(Lc3)):
                ax.annotate(Ln3[i], (Lc3[i]+0.1, Lp3[i]+0.1))

    
        if(sigLine > 0):
            plt.axhline(y=-np.log10(sigLine), color='r', linestyle='dotted')
    

    
    def clean(self):
        """
        Removes scores obtained from previous runs
        
        """
        self._SCORES = {}
        self._SKIPPED = {}
        
               
    def score_chr(self,chrs,unloadRef=False,method='saddle',mode='auto',reqacc=1e-100,intlimit=100000,parallel=1,nobar=False,autorescore=False,keep_idx=None):
        """
        Perform gene scoring for full chromosomes
        
        Args:
        
            chrs(list): List of chromosomes to score.
            unloadRef(bool): Keep only reference data for one chromosome in memory (True, False) per core
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            parallel(int) : # of cores to use
            nobar(bool): Do not show progress bar
            autorescore(bool): Automatically try to re-score failed genes via Pearson's algorithm
        
        """
        tic = time.time()
        
        S = np.array(chrs)
        
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        lock = mp.Manager().Lock()
        
        # Build list of genes for chromosomes
        G = []
        for c in S:
            G.extend(self._CHR[str(c)][0])
        
        res = self.score(G,parallel,unloadRef,method,mode,reqacc,intlimit,nobar,autorescore,keep_idx)
        
        toc = time.time()
        
        if res is not None:
            print("[time]:",str(round(toc-tic,1))+"s;",round(len(G)/(toc-tic),2),"genes/s")
        
        return res
        

    def score_all(self,parallel=1,method='saddle',mode='auto',reqacc=1e-100,intlimit=100000,nobar=False,autorescore=False,keep_idx=None):
        """
        Perform full gene scoring
        
        Args:
        
            parallel(int) : # of cores to use
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            nobar(bool): Do not show progress bar
            autorescore(bool): Automatically try to re-score failed genes via Pearson's algorithm
        
        """
        
        self._SCORES = {}
        
        return self.score_chr([i for i in range(1,23)],True,method,mode,reqacc,intlimit,parallel,nobar,autorescore,keep_idx)
        
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
class chi2sum(genescorer):
    """
    
    Implementation of chi2 sum based genescorer
    
    """
    
    def __init__(self,window=50000,varcutoff=0.99,MAF=0.05,genome=None,gpu=False):
        """
        Gene scoring via sum of chi2
        
        Args:
        
            window(int): Window size around gene tx start and end
            varcutoff(float): Variance to keep
            MAF(double): MAF cutoff 
            genome(Genome): Set gene annotation
            gpu(bool): Use GPU for linear algebra operations (requires cupy library)

        """
        
        self._window = window
        self._varcutoff = varcutoff
        self._MAF = MAF

        self._GWAS = {}
        self._GWAS_beta = {}
        self._GWAS_alleles = {}
        
        self._GENES = {}
        self._CHR = {}
        self._SKIPPED = {}
        self._SCORES = {}
        
        self._joint = False
        self._WEIGHT = {}
                        
        # Set annotation
        if genome is not None:
            self._GENEID = genome._GENEID
            self._GENESYMB = genome._GENESYMB
            self._GENEIDtoSYMB = genome._GENEIDtoSYMB
            self._CHR = genome._CHR

            self._SKIPPED = genome._SKIPPED

        # Set GPU
        if gpu and cp is not None:
            self._useGPU = True
        else:
            self._useGPU = False
            if gpu and cp is None:
                print("Error: Cupy library not detected => Using CPUs")

       
    def _calcGeneSNPcorr(self,cr,gene,REF,useAll=False):
        
        if self._joint and self._MAP is not None:
            G = self._GENEID[gene]
            P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

            if gene in self._MAP:
                P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                #P = list(set(P))
     
        elif self._MAP is None:
            G = self._GENEID[gene]

            P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
        else:
            if gene in self._MAP:
                P = set(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys())))
            else:
                P = []
        
        DATA = REF[str(cr)][0].get(list(P))
      
        filtered = {}
        
        #use = []
        #RID = []
        
        # Sort out
        for D in DATA:
            # Select
            if (D[0] in self._GWAS or useAll) and (D[1] > self._MAF) and (D[0] not in filtered or D[1] < filtered[D[0]][0]):
                filtered[D[0]] = [D[1],D[2]]
                #use.append(D[2])
                #RID.append(s)
                

        # Calc corr
        RID = list(filtered.keys())
        use = []
        for i in range(0,len(RID)):
            use.append(filtered[RID[i]][1])
            
        use = np.array(use)
        
        if len(use) > 1:
            if self._useGPU:
                C = cp.asnumpy(cp.corrcoef(cp.asarray(use)))
            else:
                C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        
        return C,np.array(RID)

    
    def _calcGeneSNPcorr_wAlleles(self,cr,gene,REF,useAll=False):
        
        if self._joint and self._MAP is not None:
            G = self._GENEID[gene]
            P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

            if gene in self._MAP:
                P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                #P = list(set(P))
          
        elif self._MAP is None:
            G = self._GENEID[gene]

            P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
        else:
            if gene in self._MAP:
                P = set(REF[str(cr)][0].getSNPsPos( list(self._MAP[gene].keys()) ))
            else:
                P = []
        
        DATA = REF[str(cr)][0].get(list(P))
            
        filtered = {}
        
        #use = []
        #RID = []
        
        # Sort out
        for D in DATA:
            # Select
            if (D[0] in self._GWAS or useAll) and D[1] > self._MAF and (D[0] not in filtered or D[1] < filtered[D[0]][0]) and self._GWAS_alleles[D[0]][0] == D[3] and self._GWAS_alleles[D[0]][1] == D[4]:
                          
                filtered[D[0]] = [D[1],D[2]]

                #use.append(D[2])
                #RID.append(s)

        # Calc corr
        RID = list(filtered.keys())
        use = []
        for i in range(0,len(RID)):
            use.append(filtered[RID[i]][1])
            
        use = np.array(use)
        
        if len(use) > 1:
            if self._useGPU:
                C = cp.asnumpy(cp.corrcoef(cp.asarray(use)))
            else:
                C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        
        return C,np.array(RID)

    
    def _getChi2Sum_mapper(self,RIDs,gene):
        #print([GWAS[x] for x in RIDs])
        ps = np.zeros(len(RIDs))
        for i in range(0,len(ps)):
            #ps = chi2.ppf(1- np.array([GWAS[x] for x in RIDs]),1)
            if RIDs[i] in self._MAP[gene] and self._MAP[gene][RIDs[i]][4] is not None:
                ps[i] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._MAP[gene][RIDs[i]][4])
            else:
                ps[i] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._GWAS[RIDs[i]])
                
        return np.sum(ps)
        
    def _getChi2Sum(self,RIDs):
        #print([GWAS[x] for x in RIDs])
        ps = np.zeros(len(RIDs))
        for i in range(0,len(ps)):
            #ps = chi2.ppf(1- np.array([GWAS[x] for x in RIDs]),1)
            ps[i] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._GWAS[RIDs[i]])
        return np.sum(ps)
    
    def _calcAndFilterEV(self,C):
        try:
            if self._useGPU:
                L = cp.asnumpy(cp.linalg.eigvalsh(cp.asarray(C)))
            else:
                L = np.linalg.eigvalsh(C)
        except: 
            return None
        
        F = L > 0
        L = L[F][::-1]
        
        if len(L) > 0:
            N_L = []

            # Leading EV
            c = L[0]
            N_L.append(L[0])
            
            T = np.sum(L)
            
            # Cutoff variance for remaining EVs
            for i in range(1,len(L)):
                c = c + L[i]
                if c < self._varcutoff*T:
                    N_L.append(L[i])

            return N_L
        
        else:
            return None
    
    def _scoreThread(self,N_L,S,g,method,mode,reqacc,intlimit):
        
        if N_L is not None:
            if method=='davies':
                RESULT = [g,wchissum.onemin_cdf_davies(S,N_L,acc=reqacc,mode=mode,lim=intlimit)]
            elif method=='ruben':
                RESULT = [g,wchissum.onemin_cdf_ruben(S,N_L,acc=reqacc,mode=mode,lim=intlimit)]
            elif method=='satterthwaite':
                RESULT = [g,wchissum.onemin_cdf_satterthwaite(S,N_L,mode=mode)]
            elif method=='pearson':
                RESULT = [g,wchissum.onemin_cdf_pearson(S,N_L,mode=mode)]
            elif method=='saddle':
                RESULT = [g,wchissum.onemin_cdf_saddle(S,N_L,mode=mode)]        
            else:
                RESULT = [g,wchissum.onemin_cdf_auto(S,N_L,acc=reqacc,mode=mode,lim=intlimit)]

            return RESULT
        else:
            return None
        
    def _scoremain(self,gene,unloadRef=False,method='saddle',mode='auto',reqacc=1e-100,intlimit=100000,label='',baroffset=0,nobar=False,lock=None,keep_idx=None):
        
        G = np.array(gene)
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        #cores = max(1,min(parallel,mp.cpu_count()))
        #pool = mp.Pool(cores)
        REF = {}
        
        #print("# cores:",max(1,min(parallel,mp.cpu_count())))
        #with tqdm(total=len(G), bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]", file=sys.stdout, position=baroffset, leave=True,disable=nobar) as pbar:
       
        if not nobar:
            print(' ', end='', flush=True) # Hack to work with jupyter notebook 
       
        with lock:
            pbar = tqdm(total=len(G), bar_format="{l_bar}{bar} [ estimated time left: {remaining} ] {postfix}", position=baroffset, leave=True,disable=nobar)
            #pbar.set_description(label.rjust(15,"_"))#+"("+str(self._GENEID[G[i]][4]).ljust(15)+")")

        for i in range(pbar.total):
            #print(i)
            if G[i] in self._GENEID:
                
                with lock:
                    pbar.set_postfix_str(str(self._GENEID[G[i]][4]).ljust(15))
   
                    
                cr = self._GENEID[G[i]][0]

                if not cr in REF:
                    #with lock:
                        #pbar.set_description(label+"(loading       )")

                    if unloadRef:
                        REF = {}

                    REF[cr] = self._ref.load_pos_reference(cr,keep_idx)

                    
                if len(self._GWAS_alleles)==0:
                    C,R = self._calcGeneSNPcorr(cr,G[i],REF)
                else:
                    C,R = self._calcGeneSNPcorr_wAlleles(cr,G[i],REF)

                if len(R) > 1:
                    # Score
                    if self._MAP is not None and self._joint == False:
                        S = self._getChi2Sum_mapper(R,G[i]) 
                    else:
                        S = self._getChi2Sum(R)

                    RES = self._scoreThread(self._calcAndFilterEV(C),S,G[i],method,mode,reqacc,intlimit)

                    if RES is not None and (RES[1][1]==0 or RES[1][1]==5) and RES[1][0] > 0 and RES[1][0] <= 1 and (RES[1][0] > reqacc*1e3 or ( (method=='auto' or method=='satterthwaite' or method=='pearson' or method=='saddle')  )):
                        RESULT.append( [self._GENEIDtoSYMB[RES[0]],float(RES[1][0]),len(R)])
                    elif RES is not None:
                        FAIL.append([self._GENEIDtoSYMB[RES[0]],len(R),RES[1]])
                    else:
                        TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"Singular covariance matrix"])
                        
                        
                    with lock:
                        pbar.update(1)

                else:

                    if len(R) == 1:
                        if self._MAP is not None and self._joint == False:
                            if R[0] in self._MAP[G[i]] and self._MAP[G[i]][R[0]][0] is not None:
                                RESULT.append( [self._GENEIDtoSYMB[G[i]],float(self._MAP[G[i]][R[0]][0]),1] )
                            else:
                                RESULT.append( [self._GENEIDtoSYMB[G[i]],float(self._GWAS[R[0]]),1] )
                                
                        else:
                            RESULT.append( [self._GENEIDtoSYMB[G[i]],float(self._GWAS[R[0]]),1] )
                    else:
                        TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"No SNPs"])

                    with lock:
                        pbar.update(1)
            else:
                TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"Not in annotation"])

                with lock:
                    pbar.update(1)

                    
        with lock:  
            pbar.set_postfix_str("done".ljust(15))
            pbar.close()
            
        return RESULT,FAIL,TOTALFAIL
    
    def score(self,gene,parallel=1,unloadRef=False,method='saddle',mode='auto',reqacc=1e-100,intlimit=1000000,nobar=False,autorescore=False,keep_idx=None):
        """
        Performs gene scoring for a given list of gene symbols
        
        Args:
        
            gene(list): gene symbols to score.
            parallel(int) : # of cores to use
            unloadRef(bool): Keep only reference data for one chromosome in memory (True, False) per core
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            nobar(bool): Do not show progress bar
            autorescore(bool): Automatically try to re-score failed genes via Pearson's algorithm
        
        """
        
        # Check method
        methods = ['auto','saddle','pearson','satterthwaite','ruben','davies']
        
        if method not in methods:
            print("No valid scoring method set. Available methods:",methods)
            return None
        else:
            print("Scoring with method",method)
        
        if type(autorescore) is str and autorescore not in methods:
            print("No valid scoring method set for rescoring. Available methods:",methods)
            return None
            
                  
        G = []
        
        for i in range(0,len(gene)):
            if gene[i] in self._GENESYMB:
                G.append(self._GENESYMB[gene[i]])
            elif gene[i] in self._GENEID:
                G.append(gene[i])        
            else:
                print("[WARNING]: "+gene[i]+" not in annotation -> ignoring")
        
        lock = mp.Manager().Lock()
        
        if parallel <= 1:
            R = self._scoremain(G,unloadRef,method,mode,reqacc,intlimit,'',0,nobar,lock,keep_idx)
        else:
            R = [[],[],[]]
            S = np.array_split(G,parallel)
            
            result_objs = []
            pool = mp.Pool(max(1,min(len(S),mp.cpu_count())))
                
            for i in range(0,len(S)): 
                
                if len(S[i]) > 0:
                    result = pool.apply_async(self._scoremain, (S[i],True,method,mode,reqacc,intlimit,'',i,nobar,lock,keep_idx))
                    result_objs.append(result)

            results = [result.get() for result in result_objs]    

            for r in results:
                R[0].extend(r[0])
                R[1].extend(r[1])
                R[2].extend(r[2])
            
            pool.close()
     
    
        print(len(R[0]),"genes scored")
        if len(R[1])>0:
            print(len(R[1]),"genes failed (try to .rescore with other settings)")
        if len(R[2])>0:
            print(len(R[2]),"genes can not be scored (check annotation)")
        
        # Store in _SCORES:
        for X in R[0]:
            self._SCORES[X[0]] = float(X[1])
               
        if (len(R[1]) > 0 and 
            (
                (type(autorescore)==bool and autorescore)
                or autorescore in methods
            )
           ):
            
            if type(autorescore)==bool:
                method = 'pearson'
            else:
                method = autorescore
            
            print("Rescoring failed genes with method",method)
                
            R = self.rescore(R,method=method,mode='auto',reqacc=1e-100,intlimit=10000000,parallel=parallel,nobar=nobar,keep_idx=keep_idx)
            if len(R[1])>0:
                print(len(R[1]),"genes failed to be scored")
                
        return R
     
    def activateFails(self,RESULT):
        """
        Helper method to force activate failed genes to success genes 
        
        Args:
            
            Result(list): Return of scoring function
        
        Warning:
            Only use if you know what you are doing !
        """
        for F in RESULT[1]:
            RESULT[0].append([F[0],F[2][0],F[1]])
            self._SCORES[F[0]] = F[2][0]
            
        RESULT[1].clear()
        
    def rescore(self,RESULT,method='pearson',mode='auto',reqacc=1e-100,intlimit=100000,parallel=1,nobar=False,keep_idx=None):
        """
        Function to re-score only the failed gene scorings of a previous scoring run with different scorer settings. 
       
       Args:
       
            RESULT(list): Return of one of the gene scorring methods
            parallel(int) : # of cores to use
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            nobar(bool): Do not show progress bar
        
        Warning:
        
            Does NOT deep copy the input RESULT
        
        """
        GENES = np.array(RESULT[1],dtype=object)[:,0]
        G = []
        for i in range(0,len(GENES)):
            if GENES[i] in self._GENESYMB:
                G.append(self._GENESYMB[GENES[i]])
            else:
                print("[WARNING]: "+GENES[i]+" not in annotation -> ignoring")
        
        lock = mp.Manager().Lock()
        
        if parallel <= 1:
            RES = self._scoremain(G,True,method,mode,reqacc,intlimit,'',i,nobar,lock,keep_idx)
        else:
            RES = [[],[],[]]
            S = np.array_split(G,parallel)
            
            pool = mp.Pool(max(1,min(len(S),mp.cpu_count())))
            
            result_objs = []
                
            for i in range(0,len(S)): 

                if len(S[i]) > 0:
                    result = pool.apply_async(self._scoremain, (S[i],True,method,mode,reqacc,intlimit,'',i,nobar,lock,keep_idx))
                    result_objs.append(result)

            results = [result.get() for result in result_objs]    

            for r in results:
                RES[0].extend(r[0])
                RES[1].extend(r[1])
                RES[2].extend(r[2])
            
            pool.close()
            
        RESULT[0].extend(RES[0])
       
        
        # Store in _SCORES:
        for X in RES[0]:
            self._SCORES[X[0]] = float(X[1])
    
        print(len(RES[0]),"genes successfully rescored.")
        if len(RES[1]) > 0:
            print(len(RES[1]),"genes still failed.")
    
        RESULT[2].extend(RES[2])
        RESULT[1].clear()
        RESULT[1].extend(RES[1])
        
        return RESULT
    
    
    
    def score_gene_bulk_chr(self,chrs,gene,data,unloadRef=False,method='saddle',mode='auto',reqacc=1e-100,intlimit=100000,autorescore=False):
        """
        Perform scoring in bulk for supplied set of SNPs
        
        Args:
            chrs(int) : Chromosome number the supplied SNPs are located on
            gene(string): Gene symbol for the SNPs
            data(list} : List of SNP data in format [ [rsid1,rsid2,...], [GWASid1, GWASid2,...], M ] with M a pvalue matrix (rows: GWAS, cols: rsid)
            unloadRef(bool): Remove loaded reference data from memory
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            autorescore(bool): Automatically try to re-score failed genes via Pearson's algorithm
        
        """
        
        # Prepare ref panel
        if not hasattr(self, '_REF') or not chrs in self._REF:
                
            if not hasattr(self, '_REF') or unloadRef:
                self._REF = {}

            self._REF[chrs] = self._ref.load_pos_reference(chrs)

        # Init GWAS dummy data
        self._GWAS = {}
        for rsid in data[0]:
            self._GWAS[rsid]= None
        
        
        # Calc SNP-SNP correlation
        C,R = self._calcGeneSNPcorr(chrs,self._GENESYMB[gene],self._REF)
        
        # Calc and filter EV
        EVL = self._calcAndFilterEV(C)
        
        # Loop over GWAS
        RESULT = {}
        for i in range(0,len(data[1])):
            GID = data[1][i]
            RESULT[GID] = [[],[],[]]
            
            if gene not in self._GENESYMB:
                RESULT[GID][2].append([self._GENEIDtoSYMB[G[i]],"Not in annotation"])
                continue
                
            # Set SNP data
            for j in range(0,len(data[0])):
                self._GWAS[data[0][j]] = data[2][i,j]
            
            if len(R) > 1:
                S = self._getChi2Sum(R)

                RES = self._scoreThread(EVL,S,self._GENESYMB[gene],method,mode,reqacc,intlimit)

                if (RES[1][1]==0 or RES[1][1]==5) and RES[1][0] > 0 and RES[1][0] <= 1 and (RES[1][0] > reqacc*1e3 or ( (method=='auto' or method=='satterthwaite' or method=='pearson' or method=='saddle')  )):
                    RESULT[GID][0].append( [self._GENEIDtoSYMB[RES[0]],float(RES[1][0]),len(R)])
                else:
                    RESULT[GID][1].append([self._GENEIDtoSYMB[RES[0]],len(R),RES[1]])

            
            else:

                if len(R) == 1:
                    RESULT[GID][0].append( [self._GENEIDtoSYMB[G[i]],float(self._GWAS[R[0]]),1] )
                else:
                    RESULT[GID][2].append([self._GENEIDtoSYMB[G[i]],"No SNPs"])

        
                
        return RESULT
        
    
    
    
    
    def _calcMultiGeneSNPcorr(self,cr,genes,REF,wAlleles=True):
        
        filtered = {}
        
        use = []
        RID = []
        pos = []
        
        for gene in genes:
            DATA = []
            if self._joint and self._MAP is not None:
                G = self._GENEID[gene]
                P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

                if gene in self._MAP:
                    P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                    #P = list(set(P))
                
            elif self._MAP is None:
                G = self._GENEID[gene]

                P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
            else:
                if gene in self._MAP:
                    P = set(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys())))
                else:
                    P = []
                    
            DATA = REF[str(cr)][0].get(list(P))
    
            # Sort out
            for D in DATA:
                # Select
                if D[0] in self._GWAS and D[1] > self._MAF and (D[0] not in filtered or filtered[D[0]][0] < D[1]) and (not wAlleles or (self._GWAS_alleles[D[0]][0] == D[3] and self._GWAS_alleles[D[0]][1] == D[4])):
                    filtered[D[0]] = [D[1],D[2]]
                    #use.append(D[2])
                    #RID.append(s)

            pos.append(len(filtered))
            
        # Calc corr
        RID = list(filtered.keys())
        use = []
        for i in range(0,len(RID)):
            use.append(filtered[RID[i]][1])
            
        use = np.array(use)
        
        if len(use) > 1:
            if self._useGPU:
                C = cp.asnumpy(cp.corrcoef(cp.asarray(use)))
            else:
                C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        
        return C,np.array(RID),pos
    
                
    def plot_genesnps(self,G,show_correlation=False,mark_window=False,tickspacing=10,color='limegreen',corrcmap=None):
        """
        Plots the SNP p-values for a list of genes and the genotypic SNP-SNP correlation matrix
        
        Args:
        
            G(list): List of gene symbols
            show_correlation(bool): Plot the corresponding SNP-SNP correlation matrix 
            tickspacing(int): Spacing of ticks for correlation plot
            color(color): Color for SNP associations
            corrcmap(cmap): Colormap to use for correlation plot (None for default)
            
        """
        print("Window size:",self._window)

        if isinstance(G,list):
            if len(G) < 2:
                G = G[0]
                
        if not isinstance(G, list):
            if G not in self._GENESYMB:
                print(G,"not in annotation!")
                return
            
            D = self._GENEID[self._GENESYMB[G]]
            print("Chr:",D[0])

            REF = {}
            REF[D[0]] = self._ref.load_pos_reference(D[0])
            DB  = self._ref.load_snp_reference(D[0]) 
            
            if len(self._GWAS_alleles) == 0:
                corr,RID,pos = self._calcMultiGeneSNPcorr(D[0],[self._GENESYMB[G]],REF,False)
            else:
                corr,RID,pos = self._calcMultiGeneSNPcorr(D[0],[self._GENESYMB[G]],REF,True)
              
            if mark_window:
                P = list(DB.getSortedKeys().irange(self._GENEID[self._GENESYMB[G]][1],self._GENEID[self._GENESYMB[G]][2]))

                #print(P[0],"-",P[-1])
                DATA = np.array(DB.getSNPatPos(P))
                #print("[DEBUG]:",len(RID),len(DATA),self._GENEID[self._GENESYMB[G]][2]-self._GENEID[self._GENESYMB[G]][1],P[0],P[-1])
                # Find start
                sid = 0
                stop = False
                for i in range(0,len(DATA)):
                    if not stop:
                        for k in range(0,len(RID)):
                            if RID[k] == DATA[i]:
                                sid = k
                                stop = True
                                #print("First tx SNP:",DATA[i],"@",DB.getSNPpos(DATA[i]))
                                break
                    else:
                        break

                # Find end            
                eid = len(RID)-1
                stop = False
                for i in range(len(DATA)-1,0,-1):
                    if not stop:
                        for k in range(len(RID)-1,0,-1):
                            if RID[k] == DATA[i]:
                                eid = k
                                stop = True
                                #print("Last  tx SNP:",DATA[i],"@",DB.getSNPpos(DATA[i]))
                                break
                    else:
                        break

        else:
            
            ID = []
            D = []
             
            for g in G:
                if g not in self._GENESYMB:
                    print(g,"not in annotation!")
                    
                    return
                ID.append(self._GENESYMB[g])
                D.append(self._GENEID[self._GENESYMB[g]])
            
            cr = D[0][0]
            for i in range(1,len(D)):
                if D[i][0] != cr:
                    print("All genes need to be on the same chromosome!")
                    return
            
            print("Chr:",cr)

            REF = {}
            REF[cr] = self._ref.load_pos_reference(cr)

            corr,RID,pos = self._calcMultiGeneSNPcorr(cr,ID,REF)

            
        print("# SNP:",len(RID))
        
        
        
        DICT = {}
        x =  []
        h1 = []
        h2 = []
        c1 = []
        c2 = []
        for i in range(0,len(RID)):
            
            pA = self._GWAS[RID[i]]
            x.append(i)
            h1.append(-np.log10(pA))
            
            DICT[i] = [RID[i],pA]
                
        if show_correlation:
            plt.subplot(1, 2, 1)

            plt.bar(x,h1,color=color)    
            plt.ylabel("$-\log_{10}(p)$")
           
                
            for i in range(0,len(pos)-1):
                plt.axvline(x=pos[i], color='black', ls=':')
        
            if mark_window and not isinstance(G, list):
                plt.axvline(sid,color='black',linestyle='dotted')
                plt.axvline(eid,color='black',linestyle='dotted')
        
            plt.subplot(1, 2, 2)

            # Generate a custom diverging colormap
            if corrcmap is None:
                cmap = sns.diverging_palette(230, 20, as_cmap=True)
            else:
                cmap = corrcmap
            sns.heatmap(corr,cmap=cmap,square=True,vmin=-1,vmax=+1,xticklabels=tickspacing,yticklabels=tickspacing)
        
            for i in range(0,len(pos)-1):
                plt.axvline(x=pos[i], color='black', ls=':')
                plt.axhline(y=pos[i], color='black', ls=':')
        
            if mark_window and not isinstance(G, list):
                plt.axvline(x=sid, color='black', ls=':')
                plt.axhline(y=sid, color='black', ls=':')
                plt.axvline(x=eid, color='black', ls=':')
                plt.axhline(y=eid, color='black', ls=':')
        
        else:
            plt.bar(x,h1,color=color)    
            plt.ylabel("$-\log_{10}(p)$")
        
            if mark_window and not isinstance(G, list):
                plt.axvline(sid,color='black',linestyle='dotted')
                plt.axvline(eid,color='black',linestyle='dotted')
        
        return DICT,corr
        

    def test_gene_assocdir(self,gene,epsilon=1e-8):
        """
        Tests for directional association of the gene
        (Requires betas of GWAS)
        
        Args:
        
            gene(str): Gene symbol to test
            epsilon(float): Regularization parameter 
        """
        if gene not in self._GENESYMB:
            print(G,"not in annotation!")
            return
        
        # Get background
        D = self._GENEID[self._GENESYMB[gene]]
       
        REF = {}
        REF[D[0]] = self._ref.load_pos_reference(D[0])

        if len(self._GWAS_alleles)==0:
            C,RID = self._calcGeneSNPcorr(D[0],self._GENESYMB[gene],REF)
        else:
            C,RID = self._calcGeneSNPcorr_wAlleles(D[0],self._GENESYMB[gene],REF)
            
        # Regularize C
        A = (C + C.T)/2
        evalu, evec = np.linalg.eigh(A)
        evalu[evalu < epsilon] = epsilon

        # Decompose
        C = evec.dot(np.diag(evalu)).dot(evec.T)
        
        L = np.linalg.cholesky(C)
        
        F = np.linalg.norm(L)**2
        
        # Get and calc stats
        z = np.array([np.sign( self._GWAS_beta[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail( self._GWAS[x]  )) for x in RID])                   
        #print("Right tail:",hpstats.onemin_norm_cdf_100d(np.sum(z),0,F))
        #print("Left tail :",hpstats.norm_cdf_100d(np.sum(z),0,F))
        
        return hpstats.onemin_norm_cdf_100d(np.sum(z),0,F),hpstats.norm_cdf_100d(np.sum(z),0,F)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class wchi2sum(chi2sum):
    """
    
    Implementation of weighted chi2 sum based genescorer
    
    Note:
        
        SNP weights have to be supplied via Mapper
    
    """
    def _getChi2Sum_mapper(self,RIDs,gene):
        #print([GWAS[x] for x in RIDs])
        ps = np.zeros(len(RIDs))
        for i in range(0,len(ps)):
            
            if RIDs[i] in self._MAP[gene] and self._MAP[gene][RIDs[i]][0] is not None:
                w = self._MAP[gene][RIDs[i]][0] 
            else:
                w = 1
                
            if RIDs[i] in self._MAP[gene] and self._MAP[gene][RIDs[i]][4] is not None:
                ps[i] = w*tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._MAP[gene][RIDs[i]][4])
            else:
                ps[i] = w*tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._GWAS[RIDs[i]])
                
        return np.sum(ps)
    
    
    def _calcGeneSNPcorr(self,cr,gene,REF,useAll=False):
        
        if self._joint and self._MAP is not None:
            G = self._GENEID[gene]
            P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

            if gene in self._MAP:
                P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                #P = list(set(P))
     
        elif self._MAP is None:
            G = self._GENEID[gene]

            P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
        else:
            if gene in self._MAP:
                P = set(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys())))
            else:
                P = []
        
        DATA = REF[str(cr)][0].get(list(P))
      
        filtered = {}
        
        #use = []
        #RID = []
        
        # Sort out
        for D in DATA:
            # Select
            if (D[0] in self._GWAS or useAll) and (D[1] > self._MAF) and (D[0] not in filtered or D[1] < filtered[D[0]][0]):
                filtered[D[0]] = [D[1],D[2]]
                #use.append(D[2])
                #RID.append(s)
                

        # Calc corr
        RID = list(filtered.keys())
        use = []
        for i in range(0,len(RID)):
            use.append(filtered[RID[i]][1])
            
        use = np.array(use)
        
        # Get weights
        w = np.ones(len(RID))
        for i in range(0,len(RID)):
            if RID[i] in self._MAP[gene] and self._MAP[gene][RID[i]][0] is not None:
                w[i] = self._MAP[gene][RID[i]][0] 
        
        Wh = np.sqrt(np.diag(w))
        
        if len(use) > 1:
            if self._useGPU:
                C = cp.corrcoef(cp.asarray(use))
                Whc = cp.asarray(Wh)
                C = cp.asnumpy(Whc.dot(C.dot(Whc)))
            else:
                C = np.corrcoef(use)
                C = Wh.dot(C.dot(Wh))   
        else:
            C = np.ones((1,1))*Wh
        
        return C,np.array(RID)

    
    def _calcGeneSNPcorr_wAlleles(self,cr,gene,REF,useAll=False):
        
        if self._joint and self._MAP is not None:
            G = self._GENEID[gene]
            P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

            if gene in self._MAP:
                P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                #P = list(set(P))
          
        elif self._MAP is None:
            G = self._GENEID[gene]

            P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
        else:
            if gene in self._MAP:
                P = set(REF[str(cr)][0].getSNPsPos( list(self._MAP[gene].keys()) ))
            else:
                P = []
        
        DATA = REF[str(cr)][0].get(list(P))
            
        filtered = {}
        
        #use = []
        #RID = []
        
        # Sort out
        for D in DATA:
            # Select
            if (D[0] in self._GWAS or useAll) and D[1] > self._MAF and (D[0] not in filtered or D[1] < filtered[D[0]][0]) and self._GWAS_alleles[D[0]][0] == D[3] and self._GWAS_alleles[D[0]][1] == D[4]:
                          
                filtered[D[0]] = [D[1],D[2]]

                #use.append(D[2])
                #RID.append(s)

        # Calc corr
        RID = list(filtered.keys())
        use = []
        for i in range(0,len(RID)):
            use.append(filtered[RID[i]][1])
            
        use = np.array(use)
        
        # Calc corr
        RID = list(filtered.keys())
        use = []
        for i in range(0,len(RID)):
            use.append(filtered[RID[i]][1])
            
        use = np.array(use)
        
        # Get weights
        w = np.ones(len(RID))
        for i in range(0,len(RID)):
            if RID[i] in self._MAP[gene] and self._MAP[gene][RID[i]][0] is not None:
                w[i] = self._MAP[gene][RID[i]][0] 
        
        Wh = np.sqrt(np.diag(w))
        
        if len(use) > 1:
            if self._useGPU:
                C = cp.corrcoef(cp.asarray(use))
                Whc = cp.asarray(Wh)
                C = cp.asnumpy(Whc.dot(C.dot(Whc)))
            else:
                C = np.corrcoef(use)
                C = Wh.dot(C.dot(Wh))   
        else:
            C = np.ones((1,1))*Wh
        
        return C,np.array(RID)

    
    
    def score(self,gene,parallel=1,unloadRef=False,method='saddle',mode='auto',reqacc=1e-100,intlimit=1000000,nobar=False,autorescore=False,keep_idx=None):
        """
        Performs gene scoring for a given list of gene symbols
        
        Args:
        
            gene(list): gene symbols to score.
            parallel(int) : # of cores to use
            unloadRef(bool): Keep only reference data for one chromosome in memory (True, False) per core
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            nobar(bool): Do not show progress bar
            autorescore(bool): Automatically try to re-score failed genes via Pearson's algorithm
        
        """     
        if self._MAP is None:
            print("Weighted chi2sum scorer works only with gene-SNP associations set via Mapper")

            return None
                  
        else:
            return super().score(gene,parallel,unloadRef,method,mode,reqacc,intlimit,nobar,autorescore,keep_idx)
        
        
        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
class cauchysum(genescorer):
    """
    
    Implementation of cauchy distribution based genescorer
    
    following math of Yaowu Liua and Jun Xie,
    doi.org/10.1080/01621459.2018.1554485
    
    """
    
    def __init__(self,window=50000,varcutoff=0.99,MAF=0.05,genome=None,gpu=False):
        """
        Gene scoring via cauchy distribution
        
        Args:
        
            window(int): Window size around gene tx start and end
            varcutoff(float): Variance to keep
            MAF(double): MAF cutoff 
            genome(Genome): Set gene annotation

        """
        
        self._window = window
        self._varcutoff = varcutoff
        self._MAF = MAF

        self._GWAS = {}
        self._GWAS_beta = {}
        self._GWAS_alleles = {}
        
        self._GENES = {}
        self._CHR = {}
        self._SKIPPED = {}
        self._SCORES = {}
        
        self._joint = False
        self._WEIGHT = {}
                        
        # Set annotation
        if genome is not None:
            self._GENEID = genome._GENEID
            self._GENESYMB = genome._GENESYMB
            self._GENEIDtoSYMB = genome._GENEIDtoSYMB
            self._CHR = genome._CHR
            self._SKIPPED = genome._SKIPPED

            
            
    def _getSNPs(self,cr,gene,REF,useAll=False):
        
        if self._joint and self._MAP is not None:
            G = self._GENEID[gene]
            P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

            if gene in self._MAP:
                P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                #P = list(set(P))
     
        elif self._MAP is None:
            G = self._GENEID[gene]

            P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
        else:
            if gene in self._MAP:
                P = set(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys())))
            else:
                P = []
        
        DATA = REF[str(cr)][0].get(list(P))
           
        filtered = {}
         
        # Sort out
        for D in DATA:
            # Select
            if (D[0] in self._GWAS or useAll) and (D[1] > self._MAF) and (D[0] not in filtered or D[1] < filtered[D[0]][0]):
                filtered[D[0]] = [D[1],D[2]]
                #use.append(D[2])
                #RID.append(s)
                

        # Calc corr
        RID = list(filtered.keys())
        
        return np.array(RID)

    
    def _getSNPs_wAlleles(self,cr,gene,REF,useAll=False):
        
        if self._joint and self._MAP is not None:
            G = self._GENEID[gene]
            P = SortedSet(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

            if gene in self._MAP:
                P.update(list(REF[str(cr)][0].getSNPsPos(list(self._MAP[gene].keys()))))
                #P = list(set(P))
          
        elif self._MAP is None:
            G = self._GENEID[gene]

            P = REF[str(cr)][1].irange(G[1] - self._window, G[2] + self._window)
        else:
            if gene in self._MAP:
                P = set(REF[str(cr)][0].getSNPsPos( list(self._MAP[gene].keys()) ))
            else:
                P = []
        
        DATA = REF[str(cr)][0].get(list(P))
            
        filtered = {}
          
        # Sort out
        for D in DATA:
            # Select
            if (D[0] in self._GWAS or useAll) and D[1] > self._MAF and (D[0] not in filtered or D[1] < filtered[D[0]][0]) and self._GWAS_alleles[D[0]][0] == D[3] and self._GWAS_alleles[D[0]][1] == D[4]:
                          
                filtered[D[0]] = [D[1],D[2]]

                #use.append(D[2])
                #RID.append(s)

        # Calc corr
        RID = list(filtered.keys())
        
        return np.array(RID)
        
            
            
    def _getChi2Sum_mapper(self,RIDs,gene):
        #print([GWAS[x] for x in RIDs])
        ps = np.zeros(len(RIDs))
        for i in range(0,len(ps)):
            #ps = chi2.ppf(1- np.array([GWAS[x] for x in RIDs]),1)
            if RIDs[i] in self._MAP[gene] and self._MAP[gene][RIDs[i]][4] is not None:
                ps[i] = self._MAP[gene][RIDs[i]][4]
            else:
                ps[i] = self._GWAS[RIDs[i]]
                
        return ps
        
    def _getChi2Sum(self,RIDs):
        #print([GWAS[x] for x in RIDs])
        ps = np.zeros(len(RIDs))
        for i in range(0,len(ps)):
            #ps = chi2.ppf(1- np.array([GWAS[x] for x in RIDs]),1)
            ps[i] = self._GWAS[RIDs[i]]
        return ps      
        
    def _scoremain(self,gene,unloadRef,label='',baroffset=0,nobar=False,lock=None):
        
        G = np.array(gene)
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        REF = {}
        
        if not nobar:
            print(' ', end='', flush=True) # Hack to work with jupyter notebook 
       
        with lock:
            pbar = tqdm(total=len(G), bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]", position=baroffset, leave=True,disable=nobar)
            pbar.set_description(label)#+"("+str(self._GENEID[G[i]][4]).ljust(15)+")")

        for i in range(pbar.total):
            #print(i)
            if G[i] in self._GENEID:
                
                with lock:
                    pbar.set_postfix_str(str(self._GENEID[G[i]][4]).ljust(15))
   
                cr = self._GENEID[G[i]][0]
    
                # Load Reference panel
                if not cr in REF:
                   
                    if unloadRef:
                        REF = {}

                    REF[cr] = self._ref.load_pos_reference(cr)

                
                # Load SNPs
                if len(self._GWAS_alleles)==0:
                    R = self._getSNPs(cr,G[i],REF)
                else:
                    R = self._getSNPs_wAlleles(cr,G[i],REF)

                if len(R) > 0:
                                    
                    # Score
                    if self._MAP is not None and self._joint == False:
                        
                        S = self._getChi2Sum_mapper(R,G[i]) 
                    else:
                        S = self._getChi2Sum(R)

                    RES = [G[i],[hpstats.cauchytest_300d(S)]]

                    RESULT.append( [self._GENEIDtoSYMB[RES[0]],float(RES[1][0]),len(R)])
                    
                    with lock:
                        pbar.update(1)

                else:

                    TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"No SNPs"])

                    with lock:
                        pbar.update(1)
            else:
                TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"Not in annotation"])

                with lock:
                    pbar.update(1)

                    
        with lock:
            pbar.set_postfix_str("done".ljust(15))
            pbar.close()
        
        
        return RESULT,FAIL,TOTALFAIL
    
    
    
    
    def score(self,gene,parallel=1,unloadRef=False,method='saddle',mode='auto',reqacc=1e-100,intlimit=1000000,nobar=False,autorescore=False,keep_idx=None):
        """
        Performs gene scoring for a given list of gene symbols
        
        Args:
        
            gene(list): gene symbols to score.
            parallel(int) : # of cores to use
            unloadRef(bool): Keep only reference data for one chromosome in memory (True, False) per core
            nobar(bool): Do not show progress bar
           
        """
        
        G = []
        
        for i in range(0,len(gene)):
            if gene[i] in self._GENESYMB:
                G.append(self._GENESYMB[gene[i]])
            elif gene[i] in self._GENEID:
                G.append(gene[i])        
            else:
                print("[WARNING]: "+gene[i]+" not in annotation -> ignoring")
        
        lock = mp.Manager().Lock()
        
        if parallel <= 1:
            R = self._scoremain(G,unloadRef,'',0,nobar,lock)
        else:
            R = [[],[],[]]
            S = np.array_split(G,parallel)
            
            result_objs = []
            pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
                
            for i in range(0,len(S)): 

                result = pool.apply_async(self._scoremain, (S[i],True,'',i,nobar,lock))
                result_objs.append(result)

            results = [result.get() for result in result_objs]    

            for r in results:
                R[0].extend(r[0])
                R[1].extend(r[1])
                R[2].extend(r[2])
            
            pool.close()
     
    
        print(len(R[0]),"genes scored")
        if len(R[1])>0:
            print(len(R[1]),"genes failed (try to .rescore with other settings)")
        if len(R[2])>0:
            print(len(R[2]),"genes can not be scored (check annotation)")
        
        # Store in _SCORES:
        for X in R[0]:
            self._SCORES[X[0]] = float(X[1])
               
      
        return R