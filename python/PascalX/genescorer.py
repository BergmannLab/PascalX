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

from tqdm import tqdm

import multiprocessing as mp

import matplotlib.pyplot as plt

import seaborn as sns

import sys
import os.path

from scipy.stats import norm

class chi2sum:
    
    def __init__(self,window=50000,varcutoff=0.99,MAF=0.05,genome=None):
        """
        Gene scoring via sum of chi2
        
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
        
        self._MAP = None
        
        # Set annotation
        if genome is not None:
            self._GENEID = genome._GENEID
            self._GENESYMB = genome._GENESYMB
            self._GENEIDtoSYMB = genome._GENEIDtoSYMB
            self._CHR = genome._CHR

            self._SKIPPED = genome._SKIPPED

            
    def load_refpanel(self, filename, parallel=1,keepfile=None,qualityT=100,SNPonly=False):
        """
        Sets the reference panel to use
        
        Args:
        
            filename(string): /path/filename (without .chr#.db ending)
            parallel(int): Number of cores to use for parallel import of reference panel
            
            keepfile: File with sample ids (one per line) to keep (only for .vcf) 
            qualityT: Quality threshold for variant to keep (only for .vcf)
            SNPonly : Import only SNPs (only for .vcf)
            
        Note:
        
            One file per chromosome with ending .chr#.db required (#: 1-22). If imported reference panel is not present, PascalX will automatically try to import from .chr#.tped.gz or .chr#.vcf.gz files.
            
        """
        self._ref = refpanel.refpanel()
        self._ref.set_refpanel(filename=filename,parallel=parallel,keepfile=keepfile,qualityT=qualityT,SNPonly=SNPonly)

        
  
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
        
        self._GENEID = GEN._GENEID
        self._GENESYMB = GEN._GENESYMB
        self._GENEIDtoSYMB = GEN._GENEIDtoSYMB
        self._CHR = GEN._CHR
        self._BAND = GEN._BAND
        
        self._SKIPPED = GEN._SKIPPED

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
            
        Note:
            A loaded mapping takes precedence over a loaded positional gene annotation
            
        """
        M = mapper()
        M.load_mapping(file,gcol,rcol,ocol,splitchr,header)
        self._MAP = M._GENEIDtoSNP
        
     
    def load_GWAS(self,file,rscol=0,pcol=1,bcol=None,a1col=None,a2col=None,delimiter=None,header=False,NAid='NA'):
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

        for line in f:
            if delimiter is None:
                L = line.split()
            else:
                L = line.split(delimiter)

            if L[pcol] != NAid and L[rscol] != NAid:
                p = float(L[pcol])
                if p > 0 and p < 1:
                    self._GWAS[L[rscol]] = p

            if bcol is not None:
                b = float(L[bcol])
                self._GWAS_beta[L[rscol]] = b
              
            if a1col is not None and a2col is not None:
                self._GWAS_alleles[L[rscol]] = [L[a1col],L[a2col]]
            
        f.close()
                    
        print(len(self._GWAS),"SNPs loaded")

    
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
        
    def _calcGeneSNPcorr(self,cr,gene,REF,useAll=False):
        
        if self._MAP is None:
            G = self._GENEID[gene]

            P = list(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

        else:
            if gene in self._MAP:
                P = self._MAP[gene][1]
            else:
                P = []
                
        DATA = REF[str(cr)][0].get(P)
           
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
            C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        
        return C,np.array(RID)

    
    def _calcGeneSNPcorr_wAlleles(self,cr,gene,REF,useAll=False):
        
        if self._MAP is None:
            G = self._GENEID[gene]

            P = list(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))

        else:
            if gene in self._MAP:
                P = self._MAP[gene][1]
            else:
                P = []
                
        DATA = REF[str(cr)][0].get(P)
            
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
            C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        
        return C,np.array(RID)

    
        
    def _getChi2Sum(self,RIDs):
        #print([GWAS[x] for x in RIDs])
        ps = np.zeros(len(RIDs))
        for i in range(0,len(ps)):
            #ps = chi2.ppf(1- np.array([GWAS[x] for x in RIDs]),1)
            ps[i] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._GWAS[RIDs[i]])
        return np.sum(ps)
    
    def _scoreThread(self,gi,C,S,g,method,mode,reqacc,intlimit):
        L = np.linalg.eigvalsh(C)
        L = L[L>0][::-1]
        N_L = []

        # Leading EV
        c = L[0]
        N_L.append(L[0])
        
        # Cutoff variance for remaining EVs
        for i in range(1,len(L)):
            c = c + L[i]
            if c < self._varcutoff*np.sum(L):
                N_L.append(L[i])

        if method=='davies':
            RESULT = [g,wchissum.onemin_cdf_davies(S,N_L,acc=reqacc,mode=mode,lim=intlimit)]
        elif method=='ruben':
            RESULT = [g,wchissum.onemin_cdf_ruben(S,N_L,acc=reqacc,mode=mode,lim=intlimit)]
        elif method=='satterthwaite':
            RESULT = [g,wchissum.onemin_cdf_satterthwaite(S,N_L,mode=mode)]
        else:
            RESULT = [g,wchissum.onemin_cdf_auto(S,N_L,acc=reqacc,mode=mode,lim=intlimit)]

        return RESULT
        
    def _scoremain(self,gene,unloadRef=False,method='auto',mode='',reqacc=1e-100,intlimit=100000,label='',baroffset=0,nobar=False,threshold=True):
        
        G = np.array(gene)
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        #cores = max(1,min(parallel,mp.cpu_count()))
        #pool = mp.Pool(cores)
        REF = {}
        
        #print("# cores:",max(1,min(parallel,mp.cpu_count())))
        with tqdm(total=len(G), bar_format="{l_bar}{bar} [ estimated time left: {remaining} ]", file=sys.stdout, position=baroffset, leave=True,disable=nobar) as pbar:
            for i in range(pbar.total):
                #print(i)
                if G[i] in self._GENEID:
                    cr = self._GENEID[G[i]][0]

                    if not cr in REF:
                        pbar.set_description(label+"(loading       )")
                        
                        if unloadRef:
                            REF = {}
                            
                        REF[cr] = self._ref.load_pos_reference(cr)
                    
                    pbar.set_description(label+"("+str(self._GENEID[G[i]][4]).ljust(15)+")")
                    
                    if len(self._GWAS_alleles)==0:
                        C,R = self._calcGeneSNPcorr(cr,G[i],REF)
                    else:
                        C,R = self._calcGeneSNPcorr_wAlleles(cr,G[i],REF)
                        
                    if len(R) > 1:
                        # Score
                        S = self._getChi2Sum(R)
                        
                        RES = self._scoreThread(i,C,S,G[i],method,mode,reqacc,intlimit)
                        
                        if RES[1][1]==0 and (RES[1][0] > reqacc*1e3 or ( (method=='auto' or method=='satterthwaite') and RES[1][0] > 0 )):
                            RESULT.append( [self._GENEIDtoSYMB[RES[0]],float(RES[1][0]),len(R)] )
                                
                        else:
                            if threshold and (RES[1][1] == 0 or RES[1][1] == 2 or RES[1][1] == 4 or RES[1][1] == 5) and abs(RES[1][0]) <= reqacc*1e3:
                                # If threshold, set to reqacc
                                RESULT.append( [self._GENEIDtoSYMB[RES[0]],reqacc,len(R)] )
                              
                            else:
                                FAIL.append([self._GENEIDtoSYMB[RES[0]],len(R),RES[1]])
                            
                        pbar.update(1)

                    else:

                        if len(R) == 1:
                            RESULT.append( [self._GENEIDtoSYMB[G[i]],float(self._GWAS[R[0]]),1] )
                        else:
                            TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"No SNPs"])
                           
                        pbar.update(1)
                else:
                    TOTALFAIL.append([self._GENEIDtoSYMB[G[i]],"Not in annotation"])
                    pbar.update(1)
            
            pbar.set_description(label+"(done          )")
            
        return RESULT,FAIL,TOTALFAIL
    
    def score(self,gene,parallel=1,unloadRef=False,method='auto',mode='',reqacc=1e-100,intlimit=1000000,nobar=False,threshold=True,autorescore=False):
        """
        Performs gene scoring for a given list of gene symbols
        
        Args:
        
            gene(list): gene symbols to score.
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
        G = []
        
        for i in range(0,len(gene)):
            if gene[i] in self._GENESYMB:
                G.append(self._GENESYMB[gene[i]])
            else:
                print("[WARNING]: "+gene[i]+" not in annotation -> ignoring")
        
        if parallel <= 1:
            R = self._scoremain(G,unloadRef,method,mode,reqacc,intlimit,'',0,nobar,threshold)
        else:
            R = [[],[],[]]
            S = np.array_split(G,parallel)
            
            result_objs = []
            pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
                
            for i in range(0,len(S)): 

                result = pool.apply_async(self._scoremain, (S[i],True,method,mode,reqacc,intlimit,'',i,nobar,threshold))
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
               
        if autorescore and len(R[1]) > 0:
            print("Rescoreing failed genes")
            R = self.rescore(R,method='ruben',mode='100d',reqacc=1e-100,intlimit=10000000,parallel=parallel,nobar=nobar)
        
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
        
    def rescore(self,RESULT,method='davies',mode='128b',reqacc=1e-100,intlimit=100000,parallel=1,nobar=False,threshold=True):
        """
        Function to re-score only the failed gene scorings of a previous scoring run with different scorer settings. 
       
       Args:
       
            RESULT(list): Return of one of the gene scorring methods
            parallel(int) : # of cores to use
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite')
            mode(string): Precision mode to use ('','128b','100d')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            threshold(bool): Threshold p-value to reqacc
            nobar(bool): Show progress bar
        
        Warning:
        
            Does NOT deep copy the input RESULT
        
        """
        GENES = np.array(RESULT[1])[:,0]
        G = []
        for i in range(0,len(GENES)):
            if GENES[i] in self._GENESYMB:
                G.append(self._GENESYMB[GENES[i]])
            else:
                print("[WARNING]: "+GENES[i]+" not in annotation -> ignoring")
        
        if parallel <= 1:
            RES = self._scoremain(G,True,method,mode,reqacc,intlimit)
        else:
            RES = [[],[],[]]
            S = np.array_split(G,parallel)
            
            pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
            
            result_objs = []
                
            for i in range(0,len(S)): 

                result = pool.apply_async(self._scoremain, (S[i],True,method,mode,reqacc,intlimit,'',i,nobar))
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
    
    
        RESULT[2].extend(RES[2])
        RESULT[1].clear()
        RESULT[1].extend(RES[1])
        
        return RESULT
    
    def score_chr(self,chrs,unloadRef=False,method='auto',mode='',reqacc=1e-100,intlimit=100000,parallel=1,nobar=False,threshold=True,autorescore=False):
        """
        Perform gene scoring for full chromosomes
        
        Args:
        
            chrs(list): List of chromosomes to score.
            unloadRef(bool): Keep only reference data for one chromosome in memory (True, False) per core
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite')
            mode(string): Precision mode to use ('','128b','100d')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            parallel(int) : # of cores to use
            nobar(bool): Show progress bar
            threshold(bool): Threshold p-value to reqacc
            autorescore(bool): Automatically try to re-score failed genes
        
        """
        
        S = np.array(chrs)
        
        RESULT = []
        FAIL = []
        TOTALFAIL = []
        
        if parallel==1:
            for c in S:
                RES = self._scoremain(self._CHR[str(c)][0],unloadRef,method,mode,reqacc,intlimit,label='[chr'+str(c)+'] ',nobar=nobar)
                RESULT.extend(RES[0])
                FAIL.extend(RES[1])
                TOTALFAIL.extend(RES[2])
        else: 
            pool = mp.Pool(max(1,min(parallel,mp.cpu_count())))
            
            result_objs = []
            for i in range(0,len(S)): 
                c = S[i]
                result = pool.apply_async(self._scoremain, (self._CHR[str(c)][0],True,method,mode,reqacc,intlimit,'[chr'+str(c)+'] ',i,nobar,threshold))
                result_objs.append(result)

            results = [result.get() for result in result_objs]    

            for r in results:
                RESULT.extend(r[0])
                FAIL.extend(r[1])
                TOTALFAIL.extend(r[2])
        
            pool.close()
            
      
        print(len(RESULT),"genes scored")
        if len(FAIL)>0:
            print(len(FAIL),"genes failed (try to .rescore with other settings)")
        if len(TOTALFAIL)>0:
            print(len(TOTALFAIL),"genes can not be scored (check annotation)")
      
        # Store in _SCORES:
        for X in RESULT:
            self._SCORES[X[0]] = float(X[1])
          
        
        if autorescore and len(FAIL) > 0:
            print("Trying to rescore failed genes")
            R = self.rescore((RESULT,FAIL,TOTALFAIL),method='ruben',mode='100d',reqacc=1e-100,intlimit=10000000,parallel=parallel,nobar=nobar)
            
            return R
     
        return RESULT, FAIL, TOTALFAIL
   

    def score_all(self,parallel=1,method='auto',mode='',reqacc=1e-100,intlimit=100000,nobar=False,threshold=True,autorescore=False):
        """
        Perform full gene scoring
        
        Args:
        
            parallel(int) : # of cores to use
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite')
            mode(string): Precision mode to use ('','128b','100d')
            reqacc(float): requested accuracy 
            intlimit(int) : Max # integration terms to use
            nobar(bool): Show progress bar
            threshold(bool): Threshold p-value to reqacc
            autorescore(bool): Automatically try to re-score failed genes
        
        """
        
        self._SCORES = {}
        
        return self.score_chr([i for i in range(1,23)],True,method,mode,reqacc,intlimit,parallel,nobar,threshold,autorescore)
      
    
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
    
    
    def _calcMultiGeneSNPcorr(self,cr,genes,REF,wAlleles=True):
        
        filtered = {}
        
        use = []
        RID = []
        pos = []
        
        for gene in genes:
            DATA = []
        
            G = self._GENEID[gene]
        
            P = list(REF[str(cr)][1].irange(G[1]-self._window,G[2]+self._window))
            
            DATA.extend(REF[str(cr)][0].get(P))
    
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
            C = np.corrcoef(use)
        else:
            C = np.ones((1,1))
            
        
        return C,np.array(RID),pos
    
                
    def plot_genesnps(self,G,show_correlation=False,tickspacing=10):
        """
        Plots the SNP p-values for a list of genes and the genotypic SNP-SNP correlation matrix
        
        Args:
        
            G(list): List of gene symbols
            show_correlation(bool): Plot the corresponding SNP-SNP correlation matrix 
            tickspacing(int): Spacing of ticks for correlation plot
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

            if len(self._GWAS_alleles) == 0:
                corr,RID,pos = self._calcMultiGeneSNPcorr(D[0],[self._GENESYMB[G]],REF,False)
            else:
                corr,RID,pos = self._calcMultiGeneSNPcorr(D[0],[self._GENESYMB[G]],REF,True)
      
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

            plt.bar(x,h1,color=plt.get_cmap("tab20")(5))    
            plt.ylabel("$-\log_{10}(p)$")
           
                
            for i in range(0,len(pos)-1):
                plt.axvline(x=pos[i], color='black', ls=':')
        
            plt.subplot(1, 2, 2)

            # Generate a custom diverging colormap

            cmap = sns.diverging_palette(230, 20, as_cmap=True)

            sns.heatmap(corr,cmap=cmap,square=True,vmin=-1,vmax=+1,xticklabels=tickspacing,yticklabels=tickspacing)
        
            for i in range(0,len(pos)-1):
                plt.axvline(x=pos[i], color='black', ls=':')
                plt.axhline(y=pos[i], color='black', ls=':')
        
        else:
            plt.bar(x,h1)    
            plt.ylabel("$-\log_{10}(p)$")
        
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
