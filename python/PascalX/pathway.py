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

from operator import itemgetter
from PascalX import hpstats, tools
import numpy as np
from scipy.stats import chi2
from abc import ABC

import random

class pathwayscorer(ABC):
    def __init__(self, genescorer, mergedist=100000,fuse=True):
        """
        Initialization:
        
        Args:
        
            genescorer(genescorer): The initialized genescorer to use to re-compute fused genes
            mergedist(int): Maximum inbetween distance of genes to merge
            fuse(bool): Fuse nearby genes to meta-genes
            
        """
        self._mergedist = mergedist
        self.set_genescorer(genescorer)
        self._fuse = fuse
        
    def set_genescorer(self,S):
        """
        Set the genescorer to use for re-computing p-values of fused genes
        
        Args:
        
            S(genescorer): The initialized genescorer to use
            
        """
        self._genescorer = S
        
    def load_modules(self,file,ncol=0,fcol=2):    
        """
        Load modules from tab separated file

        Args:
        
            file(string): path/filename
            ncol(int): Column with name of module
            fcol(int): Column with first gene (symbol) in module. Remaining genes have to follow tab (\t) separated
        
        """
        F = []
        
        f = open(file,'r')
        
        for line in f:
            
            L = line.rstrip('\n').split("\t")
            
            F.append([L[ncol],L[fcol:]])
            
        f.close()
        
        print(len(F),"modules loaded")
        
        return F
    
    def _genefusion_fuse(self,modules,chrs=None):
        FUSION_SET = []
        COMPUTE_SET = {}
        
        for M in modules:
            
            # Build up chr sets
            CHR_GENES = {}
            
            for G in M[1]:
                if G in self._genescorer._GENESYMB:
                    D = self._genescorer._GENEID[self._genescorer._GENESYMB[G]]

                    if chrs is None or D[0] in chrs:
                        if D[0] not in CHR_GENES:
                            CHR_GENES[D[0]] = [ D ]
                        else:
                            CHR_GENES[D[0]].append(D)
                #else:
                    #print("[WARNING]:",G,"not in annotation")
                    
            TOSCORE = []
               
            # Build up meta genes
            for C in CHR_GENES.keys():
                # Sort according to start position
                CHR_GENES[C] = sorted(CHR_GENES[C], key=itemgetter(1))
                   
                N = len(CHR_GENES[C])
                
                i = 0
                while i < N:
                    
                    meta = CHR_GENES[C][i][4]
                    spos = CHR_GENES[C][i][1]
                    epos = CHR_GENES[C][i][2]
                    
                    first = True
                    
                    while i < N-1:
                        i = i + 1
                        
                        if epos+self._mergedist >= CHR_GENES[C][i][1]:
                            if first:
                                meta = "METAGENE:"+meta+"_"+CHR_GENES[C][i][4]
                            else:
                                meta = meta +"_"+CHR_GENES[C][i][4]
                                
                            first = False
                            
                            epos = CHR_GENES[C][i][2]
                        else:
                            i = i - 1
                            break
                    
                    TOSCORE.append([C,spos,epos,'',meta])
                    i = i + 1
            
            # Score missing (meta)-genes
            F = []
            for G in TOSCORE:
                F.append(G[4])

                if not G[4] in self._genescorer._GENESYMB:
                    # Add to annotation
                    self._genescorer._GENESYMB[G[4]] = G[4]
                    self._genescorer._GENEID[G[4]] = G
                    self._genescorer._GENEIDtoSYMB[G[4]] = G[4]
                
                if G[0] not in COMPUTE_SET:
                    COMPUTE_SET[G[0]] = []
                   
                if not G[4] in self._genescorer._SCORES and not G[4] in COMPUTE_SET[G[0]]:     
                    # Store for each chr so that we process later more efficiently (I/O fileseek)
                    COMPUTE_SET[G[0]].append(G[4]) 
                    
                    # Debug:    
                    #if len(G[4].split("_")) > 10:
                    #    print("* [DEBUG]:",G[4])
                    #    print("         :",self._genescorer._GENEID[G[4]])
                    #    print("         :",M)
                
            FUSION_SET.append([M[0],F])
        
        return COMPUTE_SET, FUSION_SET
    
    def _genefusion(self,modules,method='auto',mode='auto',reqacc=1e-100,parallel=1,nobar=False,chrs=None):
        
        COMPUTE_SET, FUSION_SET = self._genefusion_fuse(modules,chrs)
        
        # Generate ordered list:
        SET = []
        for C in COMPUTE_SET:
            SET.extend(COMPUTE_SET[C])
        
        print("Scoring",len(SET),"missing (meta)-genes")
        
        # Compute missing (meta)-genes
        R = self._genescorer.score(SET,method=method,mode=mode,reqacc=reqacc,parallel=parallel,nobar=nobar,autorescore=True)
      
        #print(R)
        #print(FUSION_SET)
        
        return SET, FUSION_SET, R
    
    
    
    def _nogenefusion(self,modules):
        
        return [], modules, [[],[],[]]
    
    
    def get_sigpathways(self,RESULT,cutoff=1e-4):
        """
        Prints significant pathways in the result set 
        
        Args:
        
            RESULT(list) : Return of a pathwayscorer 
            cutoff(float) : Significance threshold to print pathways
            
        """
        p = []
        idx = []
        for i in range(0,len(RESULT[0])):
            if RESULT[0][i][3] < cutoff:
                idx.append(i)
                p.append(RESULT[0][i][3])
            
        # Sort and print
        so = np.argsort(p)
        indices = np.array(idx)[so]
        for i in indices:
            print(i,RESULT[0][i][0],"|",RESULT[0][i][3])

class chi2rank(pathwayscorer):
    """
    Pathway scoring via chi2 of ranked gene p-values. Nearby genes can be merged to form meta-genes
    
    Args:
    
        genescorer(genescorer): The initialized genescorer to use to re-compute fused genes
        mergedist(int): Maximum inbetween distance of genes to merge
        fuse(bool): Fuse nearby genes to meta-genes
        
    """
       
    def score(self,modules,method='saddle',mode='auto',reqacc=1e-100,parallel=1,nobar=False,genes_only=False,chrs_only=None):
        """
        Scores a set of pathways/modules
        
        Args:
        
            modules(list): List of modules to score
            samples(int): # of random gene sets to draw
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            nobar(bool): Show progress bar
            genes_only(bool): Compute only (fused)-genescores (accessible via genescorer method)
            chrs_only(list): Only consider genes on listed chromosomes. None for all.
        """
        
        # Compute fusion sets
        if self._fuse:
            COMPUTE_SET,FUSION_SET,R = self._genefusion(modules,method=method,mode=mode,reqacc=reqacc,parallel=parallel,nobar=nobar,chrs=chrs_only)
        else:
            COMPUTE_SET,FUSION_SET,R = self._nogenefusion(modules)
            
       

        if not genes_only: 
            # Build dictionary
            META_DIC = {}
            for m in self._genescorer._SCORES:
                if m[:9] == 'METAGENE:':
                    META_DIC[m] = self._genescorer._SCORES[m]
                    
            # Remove from ._SCORES to have almost same baseline for all modules
            for C in COMPUTE_SET:
                if C[:9] == 'METAGENE:' and C in self._genescorer._SCORES:
                    del self._genescorer._SCORES[C]
                  
            RESULT = []
            FAILS = R[1]

            # Score modules
            for F in FUSION_SET:

                # Note: Ranking inside the F loop because META-GENES have to be added on case by case basis
                # Rank gene scores
                L = list(self._genescorer._SCORES.keys())
                S = []

                for i in range(0,len(L)):
                    S.append(self._genescorer._SCORES[L[i]])

                # Add meta genes
                for i in range(0,len(F[1])):
                    if F[1][i][:9] == 'METAGENE:' and F[1][i] in META_DIC:
                        L.append(F[1][i])
                        S.append(META_DIC[F[1][i]])

                        # Remove metagene member genes from background gene list
                        mgenes = F[1][i][9:].split("_")

                        for g in mgenes:
                            if g in self._genescorer._SCORES:
                                I = L.index(g)
                                #print(I,"-", g)
                                del L[I]
                                del S[I]

                # Rank        
                ra = np.argsort(S)

                RANKS = {}
                for i in range(0,len(ra)):
                    RANKS[L[ra[i]]] = (i+1.)/(len(L)+1.) # +1: Ranking t start at 1

                # Calc chi2
                chi = np.zeros(len(F[1]))
                gpval = np.zeros(len(F[1]))
                fail = 0
                for i in range(0,len(F[1])):
                    if F[1][i] in RANKS:
                        chi[i] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(RANKS[F[1][i]])
                        gpval[i] = RANKS[F[1][i]]

                    else:
                        fail = fail + 1
                        gpval[i] = np.NaN

                        #print("[WARNING]: No gene score for",F[1][i])

                # Calc p-value
                df = len(chi)-fail

                if df > 0:
                    S = np.sum(chi)

                    p = hpstats.onemin_chi2_cdf(S,dof=df)


                    chi[chi == 0] = np.NaN
                    RESULT.append([F[0],F[1],gpval,p])
                else:
                    chi[chi == 0] = np.NaN
                    RESULT.append([F[0],F[1],gpval,np.NaN])

            # Cleanup
            for G in COMPUTE_SET:
                if G in self._genescorer._SCORES:
                    del self._genescorer._SCORES[G]
            
            
            # Return
            return [RESULT,FAILS,META_DIC]
        
        else:
            print("Only (fused)-gene scores computed")
        
        
        
class chi2perm(pathwayscorer):
    """
    Pathway scoring via testing summed inverse chi2 transformed gene p-values against equally size random samples of gene sets.
    
    Args:
    
        genescorer(genescorer): The initialized genescorer to use to re-compute fused genes
        mergedist(int): Maximum inbetween distance of genes to merge
        fuse(bool): Fuse nearby genes to meta-genes
        
    Note:
    
        Genes in the background gene sets are NOT fused.
    """
    
    def score(self,modules,samples=100000,method='saddle',mode='auto',reqacc=1e-100,parallel=1,nobar=False):
        """
        Scores a set of pathways/modules
        
        Args:
        
            modules(list): List of modules to score
            samples(int): # of random gene sets to draw
            method(string): Method to use to evaluate tail probability ('auto','davies','ruben','satterthwaite','pearson','saddle')
            mode(string): Precision mode to use ('','128b','100d','auto')
            reqacc(float): requested accuracy 
            nobar(bool): Show progress bar
            
        """
        # Compute fusion sets
        if self._fuse:
            COMPUTE_SET,FUSION_SET,R = self._genefusion(modules,method=method,mode=mode,reqacc=reqacc,parallel=parallel,nobar=nobar)
        else:
            COMPUTE_SET,FUSION_SET,R = self._nogenefusion(modules)
        
        # Build dictionary
        META_DIC = {}
        for m in R[0]:
            if m[0][:9] == 'METAGENE:':
                META_DIC[m[0]] = m[1]
       
        RESULT = []
        FAILS = R[1]
                
        # Compute chi2 values for all genes
        GENES = {}
        
        for G in self._genescorer._SCORES:
            GENES[G] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(self._genescorer._SCORES[G])
                
        G = list(GENES)
              
        for F in FUSION_SET:
            
            # Compute pathway score 
         
            # Calc chi2
            chi = np.zeros(len(F[1]))
            gpval = np.zeros(len(F[1]))
            
            fail = 0
            
            for i in range(0,len(F[1])):
                if F[1][i] in GENES:
                    chi[i] = GENES[F[1][i]]
                    gpval[i] = self._genescorer._SCORES[F[1][i]]
                    
                else:
                    fail = fail + 1
                    gpval[i] = np.NaN
                    
                    #print("[WARNING]: No gene score for",F[1][i])
            
            L = len(chi) - fail
            
            S = np.sum(chi)
            
            if L > 1:
                
                counter = 0.
                
                B = np.ones(samples)
                
                # Sample background
                for i in range(0,samples):
                    
                    # Draw len(L) random gene-sets
                    Rs = random.sample(G,L)
                    
                    b_score = np.zeros(L)
                    
                    for j in range(0,len(Rs)):
                        b_score[j] = GENES[Rs[j]]
                    
                    #print(b_score)
                    #print(np.sum(b_score))
                    
                    if np.sum(b_score) > S:
                        counter += 1.
                
                RESULT.append([F[0],F[1],gpval,(1+counter)/(1+samples)])
            
            else:
                if L == 1:
                    for j in range(0,len(F[1])):
                        if F[1][j] in self._genescorer._SCORES:
                            RESULT.append([F[0],F[1],chi,self._genescorer._SCORES[F[1][j]]])
                            break
                else:
                    RESULT.append([F[0],F[1],chi,np.NaN])
            
    
        return [RESULT,FAILS,META_DIC]