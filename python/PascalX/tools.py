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


import numpy as np
from scipy.stats import chi2
import gzip

import progressbar
import urllib
from urllib.request import urlretrieve
import re

#from numba import njit

#@njit
def normalInversionUpperTailApprox(p):
    lp = np.log(p)
    diff=1
    a1=1
    a=1
    while(diff>0.001):
        a=np.sqrt((-lp-np.log(np.sqrt(2*np.pi))-np.log(a1))*2)
        diff = np.abs(a-a1)
        a1=a;

    return a
    
def chiSquared1dfInverseCumulativeProbabilityUpperTail(p):
    p2=p/2.;
    if p2 < 1e-14:
        upper = normalInversionUpperTailApprox(p2);
        return upper**2;

    else:
        return chi2.ppf(1-p,1)

    
# Note: Its slow. Better to do via C lib
def read_vcf(filename,keepfilterfile=None,rsidOnly=True,qualityT=100):
    
    # Load filter info
    keep = set([])
    if keepfilterfile is not None:
        f = open(keepfilterfile,'r')
        for line in f:
            S = line.split("\t")[0]
            keep.add(S)
            
        f.close()
    
        print("Keeping",len(keep),"samples")
    
    # Detect compressed data
    if filename[-3:] == '.gz':
        f = gzip.open(filename,'r')
    else:
        f = open(filename,'r')
    
    sampleMap = {}
    
    dataMap = {}
    for i in range(0,23):
        dataMap[str(i)] = {}
        
    # Find header
    for line in f:
        # Detect infos and headers
        if line[:2] == "##":
            continue
            
        # Detect sample names
        if line[:2] == "#C":
            data = line.split("\t")
            tmp = data[9:]
            for i in range(0,len(tmp)):
                if (keepfilterfile is None) or (tmp[i] in keep):
                    sampleMap[i] = tmp[i]
                
            break
            
    sampleKeys = list(sampleMap.keys())
    
    print(len(sampleKeys),"unique samples found in")
    
    rsplit = re.compile("\||\/")
    
    # Main data import loop
    for line in f:
    
        # Data line
        data = line.split("\t")
        
        # Get GT pos
        tmp = data[8].split(":")
        GT = -1
        for i in range(0,len(tmp)):
            if tmp[i] == 'GT':
                GT = i
                break
        
        # Checks
        if (GT == -1) or (rsidOnly and data[2][:2] != 'rs') or (data[6] != 'PASS') or (int(data[5]) < qualityT):
            continue
    
        # Read genotype
        genotypes = data[9:]
        
        # Infer minor allele
        counter = np.zeros(len(data[4].split(","))+1,dtype='int')
        
        # Only read samples in sampleMap
        for i in sampleMap:
        
            geno = genotypes[i].split(":")[GT]
            #geno = rsplit(geno)
            geno = geno.replace("/","|").split("|")
            
            # Ignore half-calls
            if len(geno)==2 and geno[0] != "." and geno[1] != ".":
                counter[int(geno[0])] += 1
                counter[int(geno[1])] += 1
              
        minp = np.argmin(counter)
           
        gd = np.zeros(len(sampleKeys),dtype='B')
        
        for i in range(0,len(sampleKeys)):
           
            geno = genotypes[sampleKeys[i]].split(":")[GT]
            #geno = rsplit(geno)
            geno = geno.replace("/","|").split("|")
            
            # Ignore half-calls
            if len(geno)==2 and geno[0] != '.' and geno[1] != '.':

                if int(geno[0]) == minp:
                    gd[i] += 1

                if int(geno[1]) == minp:
                    gd[i] += 1
                

        dataMap[data[0]][int(data[1])] = [data[2],gd]
    
    f.close()
    
    return dataMap,sampleMap
    
    
class downloader:
    
    def __init__(self):
        self.pbar = None
        self.file = ''
 
    def show_progress(self,block_num, block_size, total_size):
        pbar = self.pbar
        file = self.file
        
        if pbar is None:
            widgets = [
            str(file)+' ', progressbar.Percentage(),
            ' ', progressbar.Bar(),
            ' ', progressbar.ETA(),
            ' ', progressbar.FileTransferSpeed(),
            ]
            pbar = progressbar.ProgressBar(widgets=widgets,maxval=total_size)
            pbar.start()   
        downloaded = block_num * block_size
        if downloaded < total_size:
            pbar.update(downloaded)
        else:
            pbar.finish()
            pbar = None

    # Helper function to download file        
    def download(self,url,path,filename):
        self.file = filename
        urlretrieve(url, path+filename, self.show_progress)