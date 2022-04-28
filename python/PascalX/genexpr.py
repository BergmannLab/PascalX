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

from PascalX.tools import downloader
from PascalX import genome
from PascalX import tools, hpstats
from PascalX import genescorer,pathway

import urllib
from urllib.request import urlopen
import io

import os
from os import path

import gzip
import zlib
import pickle

import numpy as np

import seaborn as sns


class genexpr:
    
    def __init__(self):
        
        self._GTEX_RNAseq_raw = 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
        self._GTEX_RNAseq = 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz' 
        self._GTEX_samples = 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
    
        self._GTEX_mapping = {}
        self._GTEX_tissues = {}
        
        self._GTEX_tissue_TPM = {}
        self._GTEX_raw_index = {}
        self._GTEX_file = ''
        
        self._ENSEMBL_gene_length = {}
        
        self._GENOME = None
        
        
    def get_GTEX_expr(self,filename):
        """
        Downloads GTEx v8 data and imports. 
        
        Args:
            
            filename(string): Filename to store downloaded GTEX data
            
            
        Note:
            The import may take several hours.
        
        """
        if not os.path.exists('GTEX'):
            os.mkdir('GTEX')
        
        D = downloader()
        
        # Download RNAseq data (reads)
        D.download('https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/'+self._GTEX_RNAseq_raw,'GTEX/',self._GTEX_RNAseq_raw)
        
        # Download RNAseq data (TPM)
        D.download('https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/'+self._GTEX_RNAseq,'GTEX/',self._GTEX_RNAseq)
        
        # Download sample info data
        D.download('https://storage.googleapis.com/gtex_analysis_v8/annotations/'+self._GTEX_samples,'GTEX/',self._GTEX_samples)
        
        print("GTEX data downloaded")

        # Download ensemble data
        cmd = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "start_position" /><Attribute name = "end_position" /></Dataset></Query>'
        url = 'http://www.ensembl.org/biomart/martservice?query='
        
        with urlopen(url+urllib.parse.quote(cmd)) as f:
          
            charset = f.info().get_content_charset()
            
            html = f.read().decode(charset)

            with open('GTEX/'+filename+'_ensembl_genelength.tsv', 'w') as f:
                f.write(html)
        
        print("Ensemble annotation downloaded")
        
        # Convert
        self._convert_GTEX_expr(filename,self._GTEX_RNAseq,self._GTEX_samples)
        
        
    def _create_GTEX_mapping(self,outfile,mapfile):
        self._GTEX_mapping = {}
        self._GTEX_tissues = {}
        
        with open('GTEX/'+mapfile,'r') as f:
                
                f.readline()

                for line in f:
                    line = line.rstrip().split('\t')
                      
                    self._GTEX_mapping[line[0]] = line[6]
                    
                    if line[6] not in self._GTEX_tissues:
                        self._GTEX_tissues[line[6]] = True
        
        # Store mapping for reload
        fp = gzip.open("GTEX/"+outfile+'_mapping.p.gz','wb')
        pickle.dump(self._GTEX_mapping, fp, protocol=pickle.HIGHEST_PROTOCOL)
        fp.close()              
        
        # Store tissues for reload
        fp = gzip.open("GTEX/"+outfile+'_tissues.p.gz','wb')
        pickle.dump(self._GTEX_tissues, fp, protocol=pickle.HIGHEST_PROTOCOL)
        fp.close()              
        
    
    def _convert_GTEX_expr(self,outfile,infile,mapfile):
        print("Starting GTEX data conversion...")
        
        # Load sample id -> sample type mapping
        self._create_GTEX_mapping(outfile,mapfile)
                
        # Load TPM for gene / tissue / sample
        with gzip.open('GTEX/'+infile,'rt') as f:        
            with gzip.open('GTEX/'+outfile+".tsv.gz",'wt') as fstore:
       
                f.readline()
                f.readline()

                # Get column ids
                header = f.readline().rstrip().split('\t')
                col_id = {}

                for i in range(2,len(header)):
                    col_id[header[i]] = i
                
                c = 0
                
                # Loop over data
                for line in f:
                    line = line.rstrip().split('\t')

                    # Set gene
                    eid = line[0][0:15]
                    
                    # Prepare data structure
                    GTEX_tissue_TPM = {}
                    for tissue in self._GTEX_tissues:
                        GTEX_tissue_TPM[tissue] = []
                        
                    # Read
                    for sample in col_id:
                        tpm = float(line[col_id[sample]])

                        GTEX_tissue_TPM[self._GTEX_mapping[sample]].append(tpm)
                     
                    # Write
                    for tissue in GTEX_tissue_TPM:
                        if len(GTEX_tissue_TPM[tissue]) > 0:
                            m = np.log2( 1 + np.median( np.array(GTEX_tissue_TPM[tissue]) ) )

                            fstore.write(eid+"\t"+tissue+"\t"+str(m)+"\n")
                    c = c + 1
        
        print("GTEX data converted [",c,"genes ]")
    
        self._convert_GTEX_expr_raw(outfile,infile,mapfile)
    
    
        
    def _convert_GTEX_expr_raw(self,outfile,infile,mapfile):
        
        print("Starting GTEX raw data conversion ...")

        
        # load gene length
        self._load_ensembl_gene_length('GTEX/'+outfile)
         
        INDEX = {}
        INDEX_g = {}
        
        # Load raw count for gene / tissue / sample
        with gzip.open('GTEX/'+infile,'rt') as f:
            with open('GTEX/'+outfile+"_raw.db",'wb') as fstore:
                with open('GTEX/'+outfile+"_tpm_raw.db","wb") as g:
                    
                    f.readline()
                    f.readline()

                    # Get column ids
                    header = f.readline().rstrip().split('\t')
                    col_id = {}
                    SRP = np.zeros(len(header)-2)
                    samples = ['']*(len(header)-2)
                    for i in range(2,len(header)):
                        col_id[header[i]] = i
                        samples[i-2] = header[i]
                    
                    self._GTEX_raw_samples = samples
                    
                    # Generate masks
                    self._generate_raw_masks()
        
                    c = 0
                    
                    # Loop over data
                    for line in f:
                        line = line.rstrip().split('\t')

                        # Set gene
                        eid = line[0][0:15]

                        if eid in self._ENSEMBL_gene_length:
                            # Only consider gene if in current ensemble annotation
                            c = c + 1
                            
                            # Write raw data
                            I = [0,0]
                            I[0] = fstore.tell()
                            fstore.write(zlib.compress(pickle.dumps(np.array(line[2:]).astype(float),protocol=pickle.HIGHEST_PROTOCOL)))
                            I[1] = fstore.tell()

                            INDEX[eid] = I

                            # Add to SRP
                            for i in range(0,len(samples)):
                                SRP[i] = SRP[i] + float(line[2+i])/self._ENSEMBL_gene_length[eid]*1000.
                                
                            tmp = {}
                            # Write raw masked RPK values over tissues 
                            for tissue in self._GTEX_raw_masks:
                                M = self._GTEX_raw_masks[tissue]
                                A = np.array(line[2:]).astype(float)[M]/self._ENSEMBL_gene_length[eid]*1000.
                                
                                tmp[tissue] = A
                            
                            J = [0,0]
                            J[0] = g.tell()  
                            g.write(zlib.compress(pickle.dumps(tmp,protocol=pickle.HIGHEST_PROTOCOL)))
                            J[1] = g.tell()
                            
                            INDEX_g[eid] = J
                    
                    
                    
                    # Store index 1
                    fp = gzip.open("GTEX/"+outfile+'_raw.idx.gz','wb')
                    pickle.dump(INDEX, fp, protocol=pickle.HIGHEST_PROTOCOL)
                    fp.close()

                    # Store index 2
                    fp = gzip.open("GTEX/"+outfile+'_tpm_raw.idx.gz','wb')
                    pickle.dump(INDEX_g, fp, protocol=pickle.HIGHEST_PROTOCOL)
                    fp.close()

                    # Store samples
                    fp = gzip.open("GTEX/"+outfile+'_raw_samples.p.gz','wb')
                    pickle.dump(samples, fp, protocol=pickle.HIGHEST_PROTOCOL)
                    fp.close()

                    # Store SRP
                    fp = gzip.open("GTEX/"+outfile+'_raw_SRP.p.gz','wb')
                    pickle.dump(SRP, fp, protocol=pickle.HIGHEST_PROTOCOL)
                    fp.close()
                    
            
        print("GTEX raw count data converted [",c," genes ]")
    
    def _load_ensembl_gene_length(self,filename):
        self._ENSEMBL_gene_length = {}
        
        with open(filename+'_ensembl_genelength.tsv',"r") as f:
            for line in f:
                line = line.rstrip().split('\t')
                
                # load length
                self._ENSEMBL_gene_length[line[0]] = int(line[2]) - int(line[1])
                  
                  
    def _load_genexpr_raw(self,genes):
        RET = {}
        with open(self._GTEX_file+"_raw.db","rb") as f:
            for g in genes:
                if g in self._GTEX_raw_index:
                    p = self._GTEX_raw_index[g]
                    f.seek(p[0])
                    data = f.read(p[1]-p[0])
                    RET[g] = pickle.loads( zlib.decompress(data) )
        
        return RET    
        
    def _generate_raw_masks(self):
        
        self._GTEX_raw_masks = {}
        A = np.array(self._GTEX_raw_samples)
        for sample in self._GTEX_mapping:
            tissue = self._GTEX_mapping[sample]

            if tissue not in self._GTEX_raw_masks:
                self._GTEX_raw_masks[tissue] = []

            pos = np.where(A==sample)[0]

            if len(pos) > 0:
                self._GTEX_raw_masks[tissue].append(pos[0])
    
    def load_expr(self,filename):
        """
        Loads the GTEx data for usage
        
        Args:
        
            filename(string): File to load
        """
        # Load mean TPM for gene / tissue 
        
        self._GTEX_file = filename
        self._GTEX_tissue_TPM = {}
        
        with gzip.open(filename+".tsv.gz",'rt') as f:
            c = 0
            for line in f:

                line = line.rstrip().split('\t')

                if line[0] not in self._GTEX_tissue_TPM:
                    self._GTEX_tissue_TPM[line[0]] = {}
                    c = c + 1
                self._GTEX_tissue_TPM[line[0]][line[1]] = float(line[2])
        
        print('Expression data for',c,'genes loaded')
        
        fp = gzip.open(filename+'_raw.idx.gz','rb')
        self._GTEX_raw_index  = pickle.load(fp)
        fp.close()
        
        fp = gzip.open(filename+'_raw_samples.p.gz','rb')
        self._GTEX_raw_samples  = pickle.load(fp)
        fp.close()
        
        fp = gzip.open(filename+'_raw_SRP.p.gz','rb')
        self._GTEX_raw_SRP  = pickle.load(fp)
        fp.close()
    
        print('Indices for raw expression data loaded')
      
        fp = gzip.open(filename+'_mapping.p.gz','rb')
        self._GTEX_mapping = pickle.load(fp)
        fp.close()
      
        print('Mapping for raw expression data loaded')
      
        fp = gzip.open(filename+'_tissues.p.gz','rb')
        self._GTEX_tissues = pickle.load(fp)
        fp.close()
    
        print("GTEX tissues loaded")
    
    
        self._load_ensembl_gene_length(filename)
        self._generate_raw_masks()
        
    
        self._GTEX_tissue_TPM_raw = {}
        
        fp = gzip.open(filename+'_tpm_raw.idx.gz','rb')
        GTEX_tpm_raw_index  = pickle.load(fp)
        fp.close()
        c = 0  
        with open(filename+"_tpm_raw.db","rb") as f:
          
            for g in GTEX_tpm_raw_index:
                c = c + 1
                p = GTEX_tpm_raw_index[g]
                f.seek(p[0])
                data = f.read(p[1]-p[0])
                self._GTEX_tissue_TPM_raw[g] = pickle.loads( zlib.decompress(data) )
        
        
        print('Raw count data for',c,'genes loaded')
        
        
        
    def load_genome(self,file,ccol=1,cid=0,csymb=5,cstx=2,cetx=3,cs=4,chrStart=0,splitchr='\t',NAgeneid='n/a',useNAgenes=False,header=False):
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
            * ``_SKIPPED`` (dict) - Genes (cid) which could not be imported

        Note:
           An unique gene id is automatically generated for n/a gene ids if ``useNAgenes=true``.

        Note:
           Identical gene ids occuring in more than one row are merged to a single gene, if on same chromosome and positional gap is < 1Mb. The strand is taken from the longer segment.    

        """
        GEN = genome.genome()
        GEN.load_genome(file,ccol,cid,csymb,cstx,cetx,cs,chrStart,splitchr,NAgeneid,useNAgenes,header)
        
        self._GENOME = GEN
        self._GENEID = GEN._GENEID
        self._GENESYMB = GEN._GENESYMB
        self._GENEIDtoSYMB = GEN._GENEIDtoSYMB
        self._CHR = GEN._CHR

        self._SKIPPED = GEN._SKIPPED
        
    
    def _chi2rank_fused(self,pathways):
        # Tissue ret init
        RET = {}
        FAILS = []
        RAW = []
        
        # Load raw data
        raw_data = {} 
        
        Nsamples = 0
        
        for pw in pathways:
            for gene in pw[1]:
                if gene not in raw_data:
                    # load
                    if gene in self._GENESYMB:
                        data = self._load_genexpr_raw([self._GENESYMB[gene]])
                        if self._GENESYMB[gene] in data:
                            raw_data[gene] = data[self._GENESYMB[gene]]
                            Nsamples = len(raw_data[gene])
        # Init
        GS = genescorer.chi2sum(genome=self._GENOME)
        PS = pathway.chi2rank(GS)
            
        _,FS = PS._genefusion_fuse(pathways)
        
            
        for i in range(0,len(FS)):
            
            #print("[DEBUG](FS):",FS[i])
    
            L = list(self._GTEX_tissue_TPM_raw.keys()) 
        
            METAGENES = {}
          
            SRP = np.copy(self._GTEX_raw_SRP)
            
            # Compute RPK for fused genes
            for g in FS[i][1]:
                
                if(g[:9]=='METAGENE:'):
               
                    data = np.zeros(Nsamples)
              
                    entities = g[9:].split('_')
                    
                    for e in entities:
                        if e in raw_data:
                            data = data + raw_data[e]

                            # Remove from background
                            if self._GENESYMB[e] in L:
                                # Background
                                #print(L.index(self._GENESYMB[e]),"-",self._GENESYMB[e],"-",e)
                                del L[L.index(self._GENESYMB[e])]

                                # SRP
                                E = self._GENEID[self._GENESYMB[e]]
                                length = E[2] - E[1]
                                SRP = SRP - raw_data[e]/(length/1000.) 
                            
                        
                    G = self._GENEID[self._GENESYMB[g]]
                    length = G[2]-G[1]
                 
                    # Divide by (meta-) gene length
                    data = data / (length/1000.)
                    
                   
                    # Add META gene to SRP
                    SRP = SRP + data
                    
                    # MASK 
                    METAGENES[g] = data                    
               
            # Adjust for SRP
            for mgene in METAGENES:
                tmp = {}
                for tissue in self._GTEX_raw_masks:
                    M = self._GTEX_raw_masks[tissue]
                    d = METAGENES[mgene][M]
                    
                    if len(d) > 0:
                        tmp[tissue] = np.median(d/SRP[M])
                    else:
                        tmp[tissue] = 0
                        
                METAGENES[mgene] = tmp  
                
                
            # Compute pathway score
            TISSUES = {}

            for gene in L:
                for tissue in self._GTEX_tissue_TPM_raw[gene]:

                    if tissue not in TISSUES:
                        TISSUES[tissue] = [[],[],[]]

                    M = self._GTEX_raw_masks[tissue]
                    TISSUES[tissue][0].append(gene)
                    if len(self._GTEX_tissue_TPM_raw[gene][tissue]) > 0:
                        TISSUES[tissue][1].append(np.median(self._GTEX_tissue_TPM_raw[gene][tissue]/SRP[M])) # Note: already masked on conversion
                    else: 
                        TISSUES[tissue][1].append(0)
                        
                    TISSUES[tissue][2].append(np.NaN)
            
            for metagene in METAGENES:
                for tissue in METAGENES[metagene]:
                    if tissue not in TISSUES:
                        TISSUES[tissue] = [[],[],[]]

                    TISSUES[tissue][0].append(metagene)
                    TISSUES[tissue][1].append(METAGENES[metagene][tissue])  
                    TISSUES[tissue][2].append(np.NaN)
            
            for tissue in TISSUES:
                TISSUES[tissue][0] = np.array(TISSUES[tissue][0])
                TISSUES[tissue][1] = np.array(TISSUES[tissue][1])
                TISSUES[tissue][2] = np.array(TISSUES[tissue][2])

                ra = np.argsort(TISSUES[tissue][1])[::-1]

                for n in range(0,len(ra)):
                    TISSUES[tissue][2][ra[n]] = (n+1.)/(len(ra)+1.) 
           
            R,F,W = self._calc_pw_enrichment(FS[i][1],TISSUES)
            
            RET[FS[i][0]] = R
            FAILS.append(F)
            RAW.append(W)
            
        return RET,FAILS,FS[i][1],RAW
        
    def _calc_pw_enrichment(self,pw,TISSUES):
        RET = {}
        
        FAILS = []
        
        RAW = {}
        
        # Calc pathway enrichment
        for tissue in TISSUES:
            if tissue not in RAW:
                RAW[tissue] = []
                
            A = np.zeros(len(pw))
            df = 0
            for i in range(0,len(pw)):
                
                if pw[i] in self._GENESYMB:
                    eid = self._GENESYMB[pw[i]]
                    pos = np.where(TISSUES[tissue][0]==eid)

                    if(len(pos[0]) > 0):

                        A[i] = tools.chiSquared1dfInverseCumulativeProbabilityUpperTail(TISSUES[tissue][2][pos[0][0]])
                        df = df + 1
                        RAW[tissue].append(TISSUES[tissue][2][pos[0][0]])
                    else:
                        if pw[i] not in FAILS:
                            FAILS.append(pw[i])
                            RAW[tissue].append(np.NaN)
                else:
                    if pw[i] not in FAILS:
                        FAILS.append(pw[i])
                        RAW[tissue].append(np.NaN)
                        
            if df > 0:
                S = np.sum(A)  
                p = hpstats.onemin_chi2_cdf(S,dof=df)
                RET[tissue] = p
            else:
                RET[tissue] = np.NaN
                
        return RET,FAILS,RAW
    
    
    def chi2rank(self, pathways, fuse=True):
        """
        Calculates tissue enrichment scores for set of genes using GTEx data
        
        Args:
        
            pathways(list): List of pathways to score [ ['name',['gene1','gene2',...]],...]
            fuse(bool): Fuse nearby genes
        """
        if not fuse:
            # Calc ranking
            TISSUES = {}

            for gene in self._GTEX_tissue_TPM:
                for tissue in self._GTEX_tissue_TPM[gene]:

                    if tissue not in TISSUES:
                        TISSUES[tissue] = [[],[],[]]

                    TISSUES[tissue][0].append(gene)
                    TISSUES[tissue][1].append(self._GTEX_tissue_TPM[gene][tissue])
                    TISSUES[tissue][2].append(np.NaN)

            for tissue in TISSUES:
                TISSUES[tissue][0] = np.array(TISSUES[tissue][0])
                TISSUES[tissue][1] = np.array(TISSUES[tissue][1])
                TISSUES[tissue][2] = np.array(TISSUES[tissue][2])

                ra = np.argsort(TISSUES[tissue][1])[::-1]

                for i in range(0,len(ra)):
                    TISSUES[tissue][2][ra[i]] = (i+1.)/(len(ra)+1.) 

            RET = {}        
            FAILS = []
            RAW = []
            
            for pw in pathways:
                R,F,W =  self._calc_pw_enrichment(pw[1],TISSUES)
                RET[pw[0]]  = R
                FAILS.append(F)
                RAW.append(W)
                
            return RET,FAILS,RAW
        
        else:
            return self._chi2rank_fused(pathways)
        
        
    def plot_genexpr(self,genes,tzscore=False,cbar_pos=(0., 0., 0.01, 0.5)):
        """
        Plots gene expression matrix for list of genes
        
        Args:
            
            genes(list) : list of genes
            tzscore(bool) : zscore over tissues per gene (true|false)
            cbar_pos(list): Position coordinates of color bar
        """   
        # Convert genes to ensemble ids
        G = []
        L = []
        
        for gene in genes:
            if gene in self._GENESYMB and self._GENESYMB[gene] in self._GTEX_tissue_TPM:
                G.append(self._GENESYMB[gene])
                L.append(gene)
                
        if len(G) > 1:
            
            if G[0] in self._GTEX_tissue_TPM:
                T = list(self._GTEX_tissue_TPM[G[0]].keys())
            else:
                print(G[0],"not found in expression data!")
                return
            
            # Generate array
            M = np.zeros((len(G),len(T)))

            for i in range(0,len(G)):
                for j in range(0,len(T)):
                    if T[j] in self._GTEX_tissue_TPM[G[i]]:
                        M[i,j] = self._GTEX_tissue_TPM[G[i]][T[j]]

                # Normalize over tissues
                if tzscore:
                    std = np.std(M[i])
                    if std != 0:
                        M[i] = (M[i] - np.mean(M[i]))/std
                    else:
                        M[i] = 0
                        
                    colors = "coolwarm"
                else:
                    colors = "Blues"
                  
            return sns.clustermap(M.T, cmap=colors,xticklabels=L,yticklabels=T,cbar_pos=cbar_pos)
            
        else:
            print("Not sufficient corresponding genes found in the annotation and/or expression data!")
        