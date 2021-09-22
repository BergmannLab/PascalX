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

from PascalX import  wchissum,tools,refpanel,hpstats,genome

import numpy as np
import gzip

import multiprocessing as mp

from scipy.special import iti0k0
from scipy.stats import chi2

from tqdm import tqdm

import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#from numba import njit

import seaborn as sns
from abc import ABC

import time



class crosscorer(ABC):
    _ENTITIES_p = {}
    _ENTITIES_b = {}
    
    def __init__(self):
        pass
    
    def load_refpanel(self,filename, parallel=mp.cpu_count()):
        """
        Sets the reference panel to use
        
        Args:
        
            filename(string): /path/filename (without .chr#.db ending)
            parallel(int): Number of cores to use for parallel import of reference panel
            
        Note:
            One file per chromosome with ending .chr#.db required (#: 1-22). If imported reference panel is not present, he will automatically try to import from .chr#.tped.gz files.        
        """
        self._ref = refpanel.refpanel()
        self._ref.set_refpanel(filename, parallel)

        
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
       
        self._GENEID = GEN._GENEID
        self._GENESYMB = GEN._GENESYMB
        self._GENEIDtoSYMB = GEN._GENEIDtoSYMB
        self._CHR = GEN._CHR
        self._BAND = GEN._BAND
        
        self._SKIPPED = GEN._SKIPPED
        
        
    def load_GWAS(self,file,rscol=0,pcol=1,bcol=2,idcol=None,name='GWAS',delimiter=None,NAid='n/a',header=False,threshold=1,mincutoff=1e-1000,rank=False):
        """
        Load GWAS summary statistics p-values and betas

        Args:

            file(string): File containing the GWAS summary statistics data. Either as textfile or gzip compressed with ending .gz
            rscol(int): Column of SNP ids
            pcol(int) : Column of p-values
            bcol(int) : Column of betas
            idcol : Column of identifiers, if several different GWAS in one file
            name : Identifier code for GWAS (needs to be unique)
            delimiter(String): Split character 
            header(bool): Header present
            NAid(String): Code for not available (rows are ignored)
        """

        minp = 1

        if file[-3:] == '.gz':
            f = gzip.open(file,'rt')
        else:
            f = open(file,'rt')

        if header:
            f.readline()

        for line in f:
            if delimiter is None:
                L = line.split()
            else:
                L = line.split(delimiter)

            if L[bcol] != NAid and L[pcol] != NAid:
                b = float(L[bcol])
                p = float(L[pcol])

                if p > 0 and p < 1 and p < threshold:

                    if p < minp:
                        minp = p

                    if not idcol is None:
                        nid = L[idcol]
                    else:
                        nid = name

                    if not nid in crosscorer._ENTITIES_p:
                        crosscorer._ENTITIES_p[nid] = {}
                        crosscorer._ENTITIES_b[nid] = {}

                    #if L[rscol] in self._ENTITIES_p[nid]:
                    #    c_test += 1

                    crosscorer._ENTITIES_p[nid][L[rscol]] = max(p,mincutoff)
                    crosscorer._ENTITIES_b[nid][L[rscol]] = b


        f.close()

        # Rank
        if rank:
            rid = [x for x in crosscorer._ENTITIES_p[nid]]
            wd = [crosscorer._ENTITIES_p[nid][x] for x in crosscorer._ENTITIES_p[nid]]
            p = np.argsort(wd)
            wr = np.zeros(len(p))

            for i in range(0,len(p)):
                wr[p[i]] = (i+1.)/(len(p)+1.) 

            for i in range(0,len(rid)):
                crosscorer._ENTITIES_p[nid][rid[i]] = wr[i]


        print(name,len(crosscorer._ENTITIES_p[nid]),"SNPs loaded ( min p:",minp,")")



    def unload_entity(self, nid):
        del crosscorer._ENTITIES_p[nid]
        del crosscorer._ENTITIES_b[nid]
        self._SCORES={}
        self._last_EA = None
        self._last_EB = None
  

    def _calcSNPcorr(self,SNPs,D,MAF=0.05):
        DATA = D.getSNPs(SNPs)
        
        use = []
        RID = []
       
        # Sort out
        for D in DATA:
            # Select
            if D[1] > MAF and D[1] < 1-MAF:
                use.append(D[2])
                RID.append(D[0])

        # Calc corr
        use = np.array(use)
        
        if len(use) > 1:
            C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        return C,np.array(RID)

    def jointlyRank(self,E_A,E_B):
        """
        Jointly QQ normalizes the p-values of two GWAS
        
        Args:
            E_A(str) : Identifier of first GWAS
            E_B(str) : Identifier of second GWAS
        """
        SNPs = np.array(list(crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()))
     
        pA = np.ones(len(SNPs))
        pB = np.ones(len(SNPs))
        
        for i in range(0,len(SNPs)):
            pA[i] = crosscorer._ENTITIES_p[E_A][SNPs[i]]
            pB[i] = crosscorer._ENTITIES_p[E_B][SNPs[i]]
         
        crosscorer._ENTITIES_p[E_A] = {}
        crosscorer._ENTITIES_p[E_B] = {}
        
        # Rank
        p = np.argsort(pA)
        wr = np.zeros(len(p))

        for i in range(0,len(p)):
            wr[p[i]] = (i+1.)/(len(p)+1.) 

        for i in range(0,len(SNPs)):
            crosscorer._ENTITIES_p[E_A][SNPs[i]] = wr[i]

        p = np.argsort(pB)
        wr = np.zeros(len(p))

        for i in range(0,len(p)):
            wr[p[i]] = (i+1.)/(len(p)+1.) 

        for i in range(0,len(SNPs)):
            crosscorer._ENTITIES_p[E_B][SNPs[i]] = wr[i]
        
        print(len(SNPs),"shared SNPs ( min p:",1./(len(p)+1),")")
    
    def get_topscores(self,N=10):
        """
        Prints the top gene scores
        
        Args:
        
            N(int): # to show
        """
        K = []
        V = []
        
        for key, value in self._SCORES.items():
            K.append(key)
            V.append(value)
    
        I = np.argsort(V)[:N]
        
        for i in I:
            print(K[i]," ",V[i])
            
            
            
    def plot_genesnps(self,G,E_A,E_B,rank=False,zscore=False,show_correlation=False,mark_window=False):
        """
        Plots the SNP p-values for a list of genes and the genotypic SNP-SNP correlation matrix
        
        Args:
        
            G(list): List of gene symbols
            show_correlation(bool): Plot the corresponding SNP-SNP correlation matrix 
            mark_window(bool): Mark the gene transcription start and end positions
        """
        if G not in self._GENESYMB:
            print(G,"not in annotation!")
            return
        
        print("Window size:",self._window)
        
        D = self._GENEID[self._GENESYMB[G]]
        print("Chr:",D[0])
        
        DB  = self._ref.load_snp_reference(D[0])    
        REF = DB.getSortedKeys()
        
        # Load SNPs from ref panel
        P = list(REF.irange(self._GENEID[self._GENESYMB[G]][1]-self._window,self._GENEID[self._GENESYMB[G]][2]+self._window))
        #print(P[0],"-",P[-1])
            
        DATA = np.array(DB.getSNPatPos(P))
        #print(DATA)
        SNPs = frozenset(list(crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()))
        
        #UNIOND = np.intersect1d(SNPs,DATA)
    
        UNIOND = [x for x in DATA if x in SNPs]
        
        corr,RID = self._calcSNPcorr(UNIOND,DB,self._MAF)
        
        print("# SNP:",len(RID))
        #print("First SNP:",RID[0],"@",DB.getSNPpos(RID[0]))
        #print("Last  SNP:",RID[-1],"@",DB.getSNPpos(RID[-1]))

        if mark_window:
            P = list(REF.irange(self._GENEID[self._GENESYMB[G]][1],self._GENEID[self._GENESYMB[G]][2]))
            
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

        
        @ticker.FuncFormatter
        def major_formatter(x, pos):
            label = str(-x) if x < 0 else str(x)
            return label
        
        DICT = {}
        x =  []
        h1 = []
        h2 = []
        c1 = []
        c2 = []
        
        wd = [crosscorer._ENTITIES_p[E_A][x] for x in RID]
        zd = [crosscorer._ENTITIES_p[E_B][x] for x in RID]
        
        if rank:
            p = np.argsort(wd)
            wd = np.zeros(len(p))

            for i in range(0,len(p)):
                wd[p[i]] = (i+1.)/(len(p)+1.) 

            p = np.argsort(zd)
            zd = np.zeros(len(p))

            for i in range(0,len(p)):
                zd[p[i]] = (i+1.)/(len(p)+1.) 

        if zscore:
            ws = np.array([np.sign(crosscorer._ENTITIES_p[E_A][x]) for x in RID])
            zs = np.array([np.sign(crosscorer._ENTITIES_p[E_B][x] ) for x in RID])
            w = ws*np.sqrt([tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(x) for x in wd])
            z = zs*np.sqrt([tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(x) for x in zd])
            
            wt = ( (w-np.mean(w))/np.std(w) )**2
            zt = ( (z-np.mean(z))/np.std(z) )**2
            
            wd = [ 1-chi2.cdf(x,1) for x in wt]
            zd = [ 1-chi2.cdf(x,1) for x in zt]
            
        
        for i in range(0,len(RID)):

            #pA = self._ENTITIES_p[E_A][RID[i]]
            #pB = self._ENTITIES_p[E_B][RID[i]]
            pA = wd[i]
            pB = zd[i]

            x.append(i)
            h1.append(-np.log10(pA))
            c1.append( np.sign(crosscorer._ENTITIES_b[E_A][RID[i]]) + 1 )

            h2.append(np.log10(pB))
            c2.append( np.sign(crosscorer._ENTITIES_b[E_B][RID[i]]) + 1)

            DICT[i] = [RID[i],pA,pB,np.sign(crosscorer._ENTITIES_b[E_A][RID[i]]),np.sign(crosscorer._ENTITIES_b[E_B][RID[i]])]
                
        
        
        if show_correlation:
            plt.subplot(1, 2, 1)

            plt.bar(x,h1,color=plt.get_cmap("tab20")(c1))    
            plt.bar(x,h2,color=plt.get_cmap("tab20")(c2))    
            plt.ylabel("$-\log_{10}(p)$")
            plt.axhline(0,color='black',linestyle='dashed')
            ax = plt.gca()
            ax.yaxis.set_major_formatter(major_formatter)
            
            if mark_window:
                plt.axvline(sid,color='black',linestyle='dotted')
                plt.axvline(eid,color='black',linestyle='dotted')
            
            plt.subplot(1, 2, 2)

            # Generate a custom diverging colormap

            cmap = sns.diverging_palette(230, 20, as_cmap=True)

            sns.heatmap(corr,cmap=cmap,square=True)
            
            if mark_window:
                plt.axvline(x=sid, color='black', ls=':')
                plt.axhline(y=sid, color='black', ls=':')
                plt.axvline(x=eid, color='black', ls=':')
                plt.axhline(y=eid, color='black', ls=':')
        else:
            plt.bar(x,h1,color=plt.get_cmap("tab20")(c1))    
            plt.bar(x,h2,color=plt.get_cmap("tab20")(c2))    
            plt.ylabel("$-\log_{10}(p)$")
            plt.axhline(0,color='black',linestyle='dashed')
            ax = plt.gca()
            ax.yaxis.set_major_formatter(major_formatter)
             
            if mark_window:
                plt.axvline(sid,color='black',ls=':')
                plt.axvline(eid,color='black',ls=':')
           
        
        return DICT,corr
         

    
    def score(self,gene,E_A=None,E_B=None,threshold=1,parallel=1,method=None,mode=None,nobar=None,reqacc=None,autorescore=False):
        """
        Performs cross scoring for a given list of gene symbols
        
        Args:
        
            gene(list): gene symbols to score.
            E_A(str): First GWAS to use
            E_B(str): Second GWAS to use
            parallel(int) : # of cores to use
            unloadRef(bool): Keep only reference data for one chromosome in memory (True, False) per core
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite')
            mode(string): Precision mode to use ('','128b','100d')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            threshold(bool): Threshold p-value to reqacc
            nobar(bool): Show progress bar
            autorescore(bool): Automatically try to re-score failed genes
        
        """
        if E_A is None:
            if self._last_EA is not None:
                E_A = self._last_EA
            else:
                print("Entity A needs to be specified (parameter E_A='name')")
                return
        else:
            self._last_EA = E_A
            
        if E_B is None:
            if self._last_EA is not None:
                E_B = self._last_EB
            else:
                print("Entity B needs to be specified (parameter E_B='name')")
                return
        else:
            self._last_EB = E_B    
        
        
        #SNPs = np.array(list(crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()))
        
        
        # Check if in annotation and sort after chromosome
        G = []
        cr = []
        for i in range(0,len(gene)):
            if gene[i] in self._GENESYMB:
                G.append(self._GENESYMB[gene[i]])
                cr.append(self._GENEID[self._GENESYMB[gene[i]]][0])
            else:
                print("[WARNING]: "+gene[i]+" not in annotation -> ignoring")
        
        I = np.argsort(cr)
        G = np.array(G)[I]
                
        # Score gene-wise
        if parallel==1:
            R = self._score_gene_thread(G,E_A,E_B,0,nobar)
        else:
            R = [[],[],[]]
            S = np.array_split(G,parallel)
            
            pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
            
            result_objs = []

            for i in range(0,len(S)): 

                result = pool.apply_async(self._score_gene_thread, (S[i],E_A,E_B,i,nobar))
                result_objs.append(result)

            results = [r.get() for r in result_objs]    

            for r in results:
                R[0].extend(r[0])
                R[1].extend(r[1])
                R[2].extend(r[2])

            pool.close()
                
        # Store in _SCORES:
        for X in R[0]:
            self._SCORES[X[0]] = float(X[1])
                      
        print(len(R[0]),"genes scored")
        if len(R[1])>0:
            print(len(R[1]),"genes failed (try to .rescore with other settings)")
        if len(R[2])>0:
            print(len(R[2]),"genes can not be scored (check annotation)")
      
                
        return R
        
    
    def score_all(self,E_A=None,E_B=None,threshold=1,parallel=1,method=None,mode=None,nobar=None,reqacc=None):
        """
        Performs cross scoring for all gene symbols
        
        Args:
        
            E_A(str): First GWAS to use
            E_B(str): Second GWAS to use
            parallel(int) : # of cores to use
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite')
            mode(string): Precision mode to use ('','128b','100d')
            reqacc(float): requested accuracy 
            nobar(bool): Show progress bar
        
        """
        self._SCORES = {}
        
        if E_A is None:
            if self._last_EA is not None:
                E_A = self._last_EA
            else:
                print("Entity A needs to be specified (parameter E_A='name')")
                return
        else:
            self._last_EA = E_A
            
        if E_B is None:
            if self._last_EA is not None:
                E_B = self._last_EB
            else:
                print("Entity B needs to be specified (parameter E_B='name')")
                return
        else:
            self._last_EB = E_B    
        
        
        Clist = self.score_chr(E_A, E_B, parallel=parallel,threshold=threshold)
           
        return Clist

    def score_chr(self,E_A,E_B,chrs=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],threshold=1,parallel=1):
        """
        Performs cross scoring for gene symbols on given chromosomes
        
        Args:
        
            chr(list): Chromosomes to score.
            E_A(str): First GWAS to use
            E_B(str): Second GWAS to use
            parallel(int) : # of cores to use
          
        """
        if not E_A in crosscorer._ENTITIES_p:
            print("[ERROR]:",E_A," not loaded")
        
        if not E_B in crosscorer._ENTITIES_p:
            print("[ERROR]:",E_B," not loaded")
            
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        if parallel==1:
            
            with tqdm(total=len(chrs), desc="cross scoring ["+str(E_A).ljust(20)+" -> "+str(E_B).ljust(20)+"]", bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]",file=sys.stdout) as pbar:
                
                for i in chrs:
                    
                    C = self._score_chr_thread(i,E_A,E_B,threshold)
                    
                    RESULT.extend(C[0])
                    FAIL.extend(C[1])
                    TOTALFAIL.extend(C[2])
                    
                    pbar.update(1)
        else:
            # NOTE: PARALLEL SEEMS BROKEN ! -> Pool does not initialize !?!
            Clist = [None]*len(chrs)
       
            pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
         
            #with tqdm(total=len(chrs), desc="cross scoring ["+str(E_A).ljust(20)+" -> "+str(E_B).ljust(20)+"]", bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]",file=sys.stdout) as pbar:
                
                #def update(*a):
                #    pbar.update(1)
                
            res = []
           
            for i in chrs:
               
                #res.append(pool.apply_async(self._score_chr_thread, args=(i,E_A,E_B,threshold), callback=update))
                res.append(pool.apply_async(self._score_chr_thread, args=(i,E_A,E_B,threshold)))

            # Wait to finish
            for i in range(0,len(res)):
                C = res[i].get()

                RESULT.extend(C[0])
                FAIL.extend(C[1])
                TOTALFAIL.extend(C[2])
                    
            pool.close()
        
        
        # Store in _SCORES:
        for X in RESULT:
            self._SCORES[X[0]] = float(X[1])
        
        return RESULT,FAIL,TOTALFAIL
    
#@njit
def _zsum_EV_cutoff(S,VT,L):
    N_L = []

    # Leading EV
    c = L[0]
    N_L.append(0.5*L[0])
    N_L.append(-0.5*L[0])

    # Cutoff variance for remaining EVs
    for i in range(1,len(L)):
        c = c + L[i]
        N_L.append(0.5*L[i])
        N_L.append(-0.5*L[i])

        if c >= VT:
            break

    return N_L

class zsum(crosscorer):
    """This class implements the cross scorer based on SNP coherence over gene windows.
    """
   
    def __init__(self, window=50000, varcutoff=0.99, MAF=0.05, leftTail=False):
        """
        Initialization:
        
        Args:
        
            window(int): Window size around gene tx start and end
            varcutoff(float): Variance to keep
            MAF(double): MAF cutoff 
            leftTail(bool): Perform a test on the left or right tail (True|False)
        """ 
        self._CHR = {}
        
        self._window = window
        self._varcutoff = varcutoff
        self._MAF = MAF
        
        self._SKIPPED = {}
        self._SCORES = {} 
        
        self._last_EA = None
        self._last_EB = None
        
        self._leftTail = leftTail
        
    
    def _scoreThread(self,gi,C,S,g,mode,reqacc,intlimit,varcutoff,window,MAF):
      
        if len(C) > 1:
            L = np.linalg.eigvalsh(C)
            L = L[L>0][::-1]
            VT = varcutoff*np.sum(L)
            
            N_L = _zsum_EV_cutoff(S,VT,L)
          
            if not self._leftTail:
                return [g,wchissum.onemin_cdf_davies(S,N_L,acc=reqacc,mode=mode,lim=intlimit)] 
            else:
                return [g,wchissum.fconstmin_cdf_davies(-1,0,S,N_L,acc=reqacc,mode=mode,lim=intlimit)] 
            
        else:
           
            if not self._leftTail:
                if S > 0:
                    r = (0.5*np.pi - iti0k0(abs(S))[1])/np.pi # abs due to symmetry
                else:
                    r = (0.5*np.pi + iti0k0(abs(S))[1])/np.pi # abs due to symmetry
            else:
                if S > 0:
                    r = (0.5*np.pi + iti0k0(abs(S))[1])/np.pi # abs due to symmetry
                else:
                    r = (0.5*np.pi - iti0k0(abs(S))[1])/np.pi # abs due to symmetry
                
            # Precision cutoff
            if r == 0:
                r = 1e-16 # OK ?
           
            return [g,[r,0]]
   
    def _score_chr_thread(self, C, E_A, E_B, threshold):
        
       
        SNPs = crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()

        S = self._ref.getChrSNPs(C)
 
        # Intersect with available chromosome SNPs
        D = np.array(list(SNPs & S))
        
        del SNPs
        del S
        
        if threshold < 1:
            # Threshold SNPs according to E_A
            I = [False]*len(D)

            for i in range(len(D)):
                #if E_p_A[D[i]] < threshold:
                if crosscorer._ENTITIES_p[E_A][D[i]] < threshold:
                
                    I[i] = True

            D = D[I]

            
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        if len(D) > 0:
            
            db  = self._ref.load_snp_reference(C)
            REF = db.getSortedKeys()
            
            if len(self._CHR[C]) > 0:
                # Score on genes
              
                # ToDo: Move into own function
            
                # Loop over genes
                for gene in self._CHR[C][0]:
                  
                    G = self._GENEID[gene]
                    
                    P = list(REF.irange(G[1]-self._window,G[2]+self._window))
            
                    DATA = np.array(db.getSNPatPos(P))
                    
                    UNIOND = np.intersect1d(D,DATA)
                    
                    if len(UNIOND) > 0:
                        corr,RID = self._calcSNPcorr(UNIOND,db,self._MAF)
                       
                       
                            
                        w = np.array([np.sign( crosscorer._ENTITIES_b[E_A][x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail( crosscorer._ENTITIES_p[E_A][x]  )) for x in RID])
                        z = np.array([np.sign( crosscorer._ENTITIES_b[E_B][x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail( crosscorer._ENTITIES_p[E_B][x]  )) for x in RID])                   
                        
                        sumr = z.dot(w)

                        score = self._scoreThread(0,corr,sumr,E_B,'auto',1e-16,1000000,self._varcutoff,self._window,self._MAF)

                        if score[1][1] !=0 or score[1][0] <= 0.0:
                            #print("[WARNING]( chr",C,"):",score)
                            FAIL.append([G[4],score,len(w),sumr])
                            
                        else:
                            RESULT.append([G[4],score[1][0],len(w),np.sign(sumr)])
                            #RESULT.append([G[4],score[1][0],len(w),np.sign(sumr),sumr,score[1]])
                           
                    else:
                        TOTALFAIL.append([G[4],"No SNPs"])
                        
                return RESULT,FAIL,TOTALFAIL
            
            else:
                # Score on full chromosome
        
                corr,RID = self._calcSNPcorr(D,db)


                # Get w, z   
                w = np.array([np.sign(E_b_A[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_A[x])) for x in RID])
                z = np.array([np.sign(E_b_B[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_B[x])) for x in RID])
      
                
                sumr = z.dot(w)

                score = self._scoreThread(0,corr,sumr,E_B,'128b',1e-16,1000000,0.99,50000,0.05)

                if score[1][1] !=0 or score[1][0] <= 0.0:
                    return [],[[C,score]],[]
                else:
                    return [[C,score[1][0],len(w)]],[],[]
        
        else:
            return [],[],[C,"No SNPs"]
    

    def _score_gene_thread(self,G,E_A,E_B,baroffset=0,nobar=False):
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        REF = {}
        DB  = {}
        
        SNPs = np.array(list(crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()))
        
        E_p_A = crosscorer._ENTITIES_p[E_A]
        E_b_A = crosscorer._ENTITIES_b[E_A]
        E_p_B = crosscorer._ENTITIES_p[E_B]
        E_b_B = crosscorer._ENTITIES_b[E_B]
        
        with tqdm(total=len(G), bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]", file=sys.stdout, position=baroffset, leave=True,disable=nobar) as pbar:
            for g in G:

                cr = self._GENEID[g][0]

                if cr == 'X' or cr == 'Y': # Skip these
                    continue

                if not cr in REF:
                    DB[cr]  = self._ref.load_snp_reference(cr)    
                    REF[cr] = DB[cr].getSortedKeys()

                #print(g,cr,self._GENEID[g])

                P = list(REF[cr].irange(self._GENEID[g][1]-self._window,self._GENEID[g][2]+self._window))

                DATA = np.array(DB[cr].getSNPatPos(P))

                UNIOND = np.intersect1d(SNPs,DATA)

                if len(UNIOND) > 0:
                    corr,RID = self._calcSNPcorr(UNIOND,DB[cr])

                   
                    # Get w, z   
                    w = np.array([np.sign(E_b_A[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_A[x])) for x in RID])

                    z = np.array([np.sign(E_b_B[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_B[x])) for x in RID])

                  
                      
                    sumr = z.dot(w)

                    score = self._scoreThread(0,corr,sumr,E_B,'auto',1e-16,1000000,0.99,50000,0.05)

                    if score[1][1] !=0 or score[1][0] <= 0.0:
                        #print("[WARNING]( chr",C,"):",score)
                        FAIL.append([self._GENEID[g][4],score,len(w),sumr])

                    else:
                        RESULT.append([self._GENEID[g][4],score[1][0],len(w),np.sign(sumr)])
                    
                else:
                    TOTALFAIL.append([self._GENEID[g][4],"No SNPs"])


                pbar.update(1)

        return RESULT,FAIL,TOTALFAIL

    
    
            
#@njit        
def _rsum_EV_cutoff(S,VT,L):
    N_L = []

    # Leading EV
    c = L[0]
    N_L.append( L[0]*(1-S/np.sqrt((1+S**2)))/2*np.sqrt(1+S**2) )
    N_L.append(-L[0]*(1+S/np.sqrt((1+S**2)))/2*np.sqrt(1+S**2) )

    # Cutoff variance for remaining EVs
    for i in range(1,len(L)):
        c = c + L[i]
        N_L.append( L[i]*(1-S/np.sqrt((1+S**2)))/2*np.sqrt(1+S**2)  )
        N_L.append( -L[i]*(1+S/np.sqrt((1+S**2)))/2*np.sqrt(1+S**2) )

        if c >= VT:
            break

    return N_L

class rsum(crosscorer):
    """This class implements the ratio cross scorer based on SNP coherence/variance over gene windows.
    """

    def __init__(self, window=50000, varcutoff=0.99, MAF=0.05, leftTail=False):
        """
        Args:
        
            window(int): Window size around gene tx start and end
            varcutoff(float): Variance to keep
            MAF(double): MAF cutoff 
            leftTail(bool): Perform a test on the left or right tail (True|False)
        """
        self._CHR = {}
        
        self._window = window
        self._varcutoff = varcutoff
        self._MAF = MAF
        
        self._SKIPPED = {}
        self._SCORES = {} 
        
        self._last_EA = None
        self._last_EB = None
        
        self._leftTail = leftTail
        

    def _scoreThread(self,gi,C,S,g,mode,reqacc,intlimit,varcutoff,window,MAF):
      
        if len(C) > 1:
            L = np.linalg.eigvalsh(C)
            L = L[L>0][::-1]
            VT = varcutoff*np.sum(L)
            
            N_L = _rsum_EV_cutoff(S,VT,L)
                    
            if not self._leftTail:
                return [g,wchissum.onemin_cdf_davies(0,N_L,acc=reqacc,mode=mode,lim=intlimit)] 
            else:
                return [g,wchissum.fconstmin_cdf_davies(-1,0,0,N_L,acc=reqacc,mode=mode,lim=intlimit)] 
                
        else:
            # Use 1d case (Cauchy)
            
            if not self._leftTail:
                ret = 1-(0.5+(np.arctan( S )/np.pi))
            else:
                ret = (0.5+(np.arctan( S )/np.pi))
            
            # Precision cutoff
            if ret == 0:
                ret = 1e-16 # OK ?
           
            return [g,[ret,0]]

       
     
  
    def _score_chr_thread(self, C, E_A, E_B, threshold):
        
        SNPs = crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()

        S = self._ref.getChrSNPs(C)
 
        # Intersect with available chromosome SNPs
        D = np.array(list(SNPs & S))
         
        del SNPs
        del S
        
        if threshold < 1:
            # Threshold SNPs according to E_A
            I = [False]*len(D)

            for i in range(len(D)):
                #if E_p_A[D[i]] < threshold:
                if crosscorer._ENTITIES_p[E_A][D[i]] < threshold:
                
                    I[i] = True

            D = D[I]

            
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        if len(D) > 0:
            
            db  = self._ref.load_snp_reference(C)
            REF = db.getSortedKeys()
            
            if len(self._CHR[C]) > 0:
                # Score on genes
              
                # ToDo: Move into own function
            
                # Loop over genes
                for gene in self._CHR[C][0]:
                  
                    G = self._GENEID[gene]
                    
                    P = list(REF.irange(G[1]-self._window,G[2]+self._window))
            
                    DATA = np.array(db.getSNPatPos(P))
                    
                    UNIOND = np.intersect1d(D,DATA)
                    
                    if len(UNIOND) > 0:
                        corr,RID = self._calcSNPcorr(UNIOND,db,self._MAF)
                       
                        w = np.array([np.sign( crosscorer._ENTITIES_b[E_A][x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail( crosscorer._ENTITIES_p[E_A][x]  )) for x in RID])
                        z = np.array([np.sign( crosscorer._ENTITIES_b[E_B][x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail( crosscorer._ENTITIES_p[E_B][x]  )) for x in RID])

                        
                        norm = z.dot(z)
                        
                        if norm != 0:
                            sumr = z.dot(w)/norm

                            score = self._scoreThread(0,corr,sumr,E_B,'auto',1e-16,1000000,self._varcutoff,self._window,self._MAF)

                            if score[1][1] !=0 or score[1][0] <= 0.0:
                                #print("[WARNING]( chr",C,"):",score)
                                FAIL.append([G[4],score,len(w),sumr])

                            else:
                                RESULT.append([G[4],score[1][0],len(w),np.sign(sumr)])
                                #RESULT.append([G[4],score[1][0],len(w),np.sign(sumr),sumr,score[1]])
                        else:
                            TOTALFAIL.append([G[4],"NaN for denom"])
                    else:
                        TOTALFAIL.append([G[4],"No SNPs"])
                        
                return RESULT,FAIL,TOTALFAIL
            
            else:
                # Score on full chromosome
        
                corr,RID = self._calcSNPcorr(D,db)

              
                # Get w, z   
                w = np.array([np.sign(E_b_A[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_A[x])) for x in RID])
                z = np.array([np.sign(E_b_B[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_B[x])) for x in RID])

                # Nominator
                sumr = z.dot(w)
                
                # Denominator
                #sumn = z.dot(z)
                
                
                score = self._scoreThread(0,corr,sumr,E_B,'128b',1e-16,1000000,0.99,50000,0.05)

                if score[1][1] !=0 or score[1][0] <= 0.0:
                    return [],[[C,score]],[]
                else:
                    return [[C,score[1][0],len(w)]],[],[]
        
        else:
            return [],[],[C,"No SNPs"]
    

    def _score_gene_thread(self,G,E_A,E_B,baroffset=0,nobar=False):
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        REF = {}
        DB  = {}
        
        SNPs = np.array(list(crosscorer._ENTITIES_p[E_A].keys() & crosscorer._ENTITIES_p[E_B].keys()))
        
        E_p_A = crosscorer._ENTITIES_p[E_A]
        E_b_A = crosscorer._ENTITIES_b[E_A]
        E_p_B = crosscorer._ENTITIES_p[E_B]
        E_b_B = crosscorer._ENTITIES_b[E_B]
        
        with tqdm(total=len(G), bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]", file=sys.stdout, position=baroffset, leave=True,disable=nobar) as pbar:
            for g in G:

                cr = self._GENEID[g][0]

                if cr == 'X' or cr == 'Y': # Skip these
                    continue

                if not cr in REF:
                    DB[cr]  = self._ref.load_snp_reference(cr)    
                    REF[cr] = DB[cr].getSortedKeys()

                #print(g,cr,self._GENEID[g])

                P = list(REF[cr].irange(self._GENEID[g][1]-self._window,self._GENEID[g][2]+self._window))

                DATA = np.array(DB[cr].getSNPatPos(P))

                UNIOND = np.intersect1d(SNPs,DATA)

                if len(UNIOND) > 0:
                    corr,RID = self._calcSNPcorr(UNIOND,DB[cr])

                    # Get w, z   
                    w = np.array([np.sign(E_b_A[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_A[x])) for x in RID])

                    z = np.array([np.sign(E_b_B[x])*np.sqrt(tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(E_p_B[x])) for x in RID])
                      
                    norm = z.dot(z)
                    if norm != 0:
                        sumr = z.dot(w)/norm

                        score = self._scoreThread(0,corr,sumr,E_B,'auto',1e-16,1000000,0.99,50000,0.05)

                        if score[1][1] !=0 or score[1][0] <= 0.0:
                            #print("[WARNING]( chr",C,"):",score)
                            FAIL.append([self._GENEID[g][4],score,len(w),sumr])

                        else:
                            RESULT.append([self._GENEID[g][4],score[1][0],len(w),np.sign(sumr)])
                    else:
                        TOTALFAIL.append([self._GENEID[g][4],"NaN for denom"])
                else:
                    TOTALFAIL.append([self._GENEID[g][4],"No SNPs"])


                pbar.update(1)

        return RESULT,FAIL,TOTALFAIL

    
