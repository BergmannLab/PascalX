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
import urllib
from urllib.request import urlopen

import io 

class genome:
    """This class handles the genome annotation. It provides functionality for import of data from text files and automatic download of annotation data from ensembl.org.
    """
    def __init__(self):
        pass
    
    def gene_info(self,gene):
        """Prints the loaded information for given gene 
        
        Args:
            gene (string): Symbol of gene to query
        """
        if gene in self._GENESYMB:
            print(self._GENEID[self._GENESYMB[gene]])
        else:
            print("Gene",gene,"not in loaded annotation")
    
    
    def get_ensembl_annotation(self,filename,genetype='protein_coding',version='GRCh38'):
        """Gene annotation download function for ensembl.org BioMart data
        
        Args:
            filename(string): File to store downloaded annotation in
            genetype(string): Comma separated list of BioMart genetypes to download (protein_coding, pseudogene, ...)
            version(string): GRCh37 | CRCh38
            
        Example:
        ::
            from PascalX.genome import genome
            G = genome()
            G.get_ensemble_annotation('ensemble_hg38.txt')
        """
        
        genetype = genetype.replace(" ","")
        
        print("Downloading gene annotation from ensembl.org BioMart [",genetype,"] (",version,")")
        
        cmd = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"/><Filter name = "biotype" value = "'+genetype+'"/><Attribute name = "ensembl_gene_id" /><Attribute name = "chromosome_name" /><Attribute name = "transcript_start" /><Attribute name = "transcript_end" /><Attribute name = "strand" /><Attribute name = "external_gene_name" /><Attribute name = "band"/></Dataset></Query>'
          
        if version == 'GRCh38':
            url = 'http://www.ensembl.org/biomart/martservice?query='
        elif version == 'GRCh37':
            url = 'http://grch37.ensembl.org/biomart/martservice?query='
        else:
            print("Unknown requested version. Available options: GRCh37 or GRCh38")
            return
    
        with urlopen(url+urllib.parse.quote(cmd)) as f:
          
            charset = f.info().get_content_charset()
            
            html = f.read().decode(charset)

            with open(filename, 'w') as f:
                f.write(html)

        print("done")
    
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
        # Load genome
        self._GENEID = {}
        self._GENESYMB = {}
        self._GENEIDtoSYMB = {}
        self._CHR = {}
        self._BAND = {}
        self._SKIPPED = {}

        # Read first onlys not n/a genes
        with open(file,'r') as f:

            if header:
                f.readline()

            for line in f:

                line = line.rstrip().split(splitchr)

                # Continue if in skip list or N/A    
                if line[cid] in self._SKIPPED or line[cid] == NAgeneid:
                    continue

                # Continue if duplicate symbol with different id
                if line[cid] not in self._GENEID and line[csymb] in self._GENESYMB:
                    continue

                if line[cid] not in self._GENEID:
                    # Set
                    self._GENEID[line[cid]] = [line[ccol][chrStart:],int(line[cstx]),int(line[cetx]),str(line[cs]),line[csymb]]

                    if line[ccol][chrStart:] in self._CHR:

                        self._CHR[line[ccol][chrStart:]][0].append(line[cid])

                        # Collect min and max positions
                        if int(line[cstx]) < self._CHR[line[ccol][chrStart:]][1]:
                            self._CHR[line[ccol][chrStart:]][1] = int(line[cstx])

                        if int(line[cetx]) > self._CHR[line[ccol][chrStart:]][2]:
                            self._CHR[line[ccol][chrStart:]][2] = int(line[cstx])     

                    else:
                        self._CHR[line[ccol][chrStart:]] = [[line[cid]],int(line[cstx]),int(line[cetx]),0]

                    self._GENEIDtoSYMB[line[cid]] = line[csymb]
                    self._GENESYMB[line[csymb]] = line[cid]
                    
                    
                    # Store band info
                    if cb is not None:
                        pos = str(line[ccol])+str(line[cb])
                        if pos not in self._BAND:
                            self._BAND[pos] = []
                        
                        self._BAND[pos].append(line[csymb])
                        
                else:
                    # Update
                    if self._GENEID[line[cid]][0] == line[ccol][chrStart:] and self._GENEIDtoSYMB[line[cid]] == line[csymb] and (abs(self._GENEID[line[cid]][1] - int(line[cstx]) ) < 1e6  or abs(self._GENEID[line[cid]][2] - int(line[cetx]) ) < 1e6 ) :

                        # Update ending position
                        if self._GENEID[line[cid]][2] < int(line[cetx]):
                            self._GENEID[line[cid]][2] = int(line[cetx])

                        # Update start position
                        if self._GENEID[line[cid]][1] > int(line[cstx]):
                            self._GENEID[line[cid]][1] = int(line[cstx])

                        # Really ok ?
                        # Check strand
                        if self._GENEID[line[cid]][3] != str(line[cs]) and (int(line[cetx])-int(line[cstx]) > self._GENEID[line[cid]][2] -self._GENEID[line[cid]][1] ):
                            # Take strand from longer segment
                            self._GENEID[line[cid]][3] = str(line[cs])

                    else:
                        # Remove and SKIP
                        if self._GENEID[line[cid]][0] != line[ccol][chrStart:]:
                            self._SKIPPED[line[cid]] = [line[csymb],'chr jump']
                        else:
                            self._SKIPPED[line[cid]] = [line[csymb],'gap > 1e6']

                        self._CHR[self._GENEID[line[cid]][0]][0].remove(line[cid])
                        del self._GENEID[line[cid]]
                        del self._GENESYMB[self._GENEIDtoSYMB[line[cid]]]
                        del self._GENEIDtoSYMB[line[cid]]


        # Add n/a s                 
        if useNAgenes:
            NAMAP = {}  
            gid = 0

            # Add n/a genes
            with open(file) as f:
                if header:
                    f.readline()

                for line in f:
                    line = line.split()

                    if line[cid] in self._SKIPPED or (line[cid] != NAgeneid) or (line[csymb] in NAMAP and NAMAP[line[csymb]] not in self._GENEID ):
                        continue

                    if line[cid] == NAgeneid:
                        if line[csymb] in self._GENESYMB and not line[csymb] in NAMAP:
                            continue

                        if line[csymb] in NAMAP:
                            line[cid] = NAMAP[line[csymb]]
                        else:
                            line[cid] = 'X'+str(gid)
                            NAMAP[line[csymb]] = line[cid]
                            gid += 1

                    # Continue if duplicate symbol with different id
                    if line[cid] not in self._GENEID and line[csymb] in self._GENESYMB:
                        continue


                    if line[cid] not in self._GENEID:
                        # Set
                        self._GENEID[line[cid]] = [line[ccol][chrStart:],int(line[cstx]),int(line[cetx]),str(line[cs]),line[csymb]]

                        # ToDo: Move up
                        if line[ccol][chrStart:] in self._CHR:

                            self._CHR[line[ccol][chrStart:]][0].append(line[cid])

                            # Collect min and max positions
                            if int(line[cstx]) < self._CHR[line[ccol][chrStart:]][1]:
                                self._CHR[line[ccol][chrStart:]][1] = int(line[cstx])

                            if int(line[cetx]) > self._CHR[line[ccol][chrStart:]][2]:
                                self._CHR[line[ccol][chrStart:]][2] = int(line[cstx])     

                        else:
                            self._CHR[line[ccol][chrStart:]] = [[line[cid]],int(line[cstx]),int(line[cetx]),0]

                        self._GENEIDtoSYMB[line[cid]] = line[csymb]
                        self._GENESYMB[line[csymb]] = line[cid]

                    else:
                        # Update
                        if self._GENEID[line[cid]][0] == line[ccol][chrStart:] and self._GENEIDtoSYMB[line[cid]] == line[csymb] and (abs(self._GENEID[line[cid]][1] - int(line[cstx]) ) < 1e6  or abs(self._GENEID[line[cid]][2] - int(line[cetx])) < 1e6 ):

                            # Update ending position
                            if self._GENEID[line[cid]][2] < int(line[cetx]):
                                self._GENEID[line[cid]][2] = int(line[cetx])

                            # Update start position
                            if self._GENEID[line[cid]][1] > int(line[cstx]):
                                self._GENEID[line[cid]][1] = int(line[cstx])

                            # Really ok ?
                            # Check strand
                            if self._GENEID[line[cid]][3] != str(line[cs]) and (int(line[cetx])-int(line[cstx]) > self._GENEID[line[cid]][2] -self._GENEID[line[cid]][1] ):
                                # Take strand from longer segment
                                self._GENEID[line[cid]][3] = str(line[cs])

                        else:
                            # Remove and SKIP
                            if self._GENEID[line[cid]][0] != line[ccol][chrStart:]:
                                self._SKIPPED[line[cid]] = [line[csymb],'chr jump']
                            else:
                                self._SKIPPED[line[cid]] = [line[csymb],'gap > 1e6']

                            self._CHR[self._GENEID[line[cid]][0]][0].remove(line[cid])
                            del self._GENEID[line[cid]]
                            del self._GENESYMB[self._GENEIDtoSYMB[line[cid]]]
                            del self._GENEIDtoSYMB[line[cid]]




        # Add X,Y ?                

        # Calculate offsets (for plotting)
        last = 0
        for i in range(1,23):
            if str(i) in self._CHR:
                self._CHR[str(i)][3] = last
                last = last + self._CHR[str(i)][2] - self._CHR[str(i)][1] + 1000



        print(len(self._GENEID),"active genes")
        if len(self._SKIPPED) > 0:
            print(len(self._SKIPPED),"inconsistent genes removed (retrieve via ._SKIPPED)")



