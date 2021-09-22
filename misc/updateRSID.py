import sys
import gzip
import numpy as np
from fastnumbers import int
import re

def process_file(CHR,file,snpdb):
    c = 0
    t = 0
    with gzip.open(file[:-7]+".snpid.vcf.gz",'wt') as h:
        with gzip.open(snpdb,'rt') as g:
            for dbline in g:
                # Search for start
                if dbline[0]=='#' or dbline.split('\t')[0]!=CHR:
                    continue
                else:
                    break

            dbline = dbline.split('\t',5)
            CACHE = None
            with gzip.open(file,'rt') as f:
                for line in f:
                    found = False
                    
                    # Search for start
                    if line[0]=='#':
                        # Copy line to new file
                        h.write(line)
                        continue

                    line_split = line.split('\t',5)
                    
                    #print("READ:",line[0:5])
                    t += 1
                    if line_split[2][:2]!='rs':
                        ln = int(line_split[1])
                        
                        # Seek line position in db
                        while int(dbline[1])<ln and dbline[0]==CHR:
                            dbline = g.readline().split('\t',5)

                        # Read all SNPs at same position
                        if int(dbline[1])==ln and dbline[0]==CHR:
                            CACHE = [dbline]
                            while int(dbline[1])==ln and dbline[0]==CHR:
                                dbline = g.readline().split('\t',5)

                                if int(dbline[1])==ln and dbline[0]==CHR:
                                    CACHE.append(dbline)


                        #print(CACHE)    
                        if CACHE is not None and int(CACHE[0][1]) == ln:
                            R = line_split[3]
                            A = line_split[4]
                                          
                            for C in CACHE:
                                X = C[3]
                                Y = C[4].split(",")
                                          
                                if R==X and A in Y:
                                #if line_split[3] == C[3] and line_split[4] == C[4]:
                                    # Store with replace
                                    #print("***",C[0:5],line[0:5])
                                    
                                    h.write( re.sub(CHR+":(\w|:)+",C[2],line,1) )
                                    found = True
                                    c += 1
                                    break
                    
                    else:
                        # Store line
                        h.write(line)

                        c += 1
                        pass

        print("# SNPs in original ref data ( CHR",CHR,"):",t)
        print("# SNPs matched:",c)
     

if __name__ == "__main__":
    if len(sys.argv) == 4:
        process_file(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print("Arguments required: chr infile snpdb")
        