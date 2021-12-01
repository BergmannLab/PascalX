# PascalX

High precision gene and pathway scoring for GWAS summary statistics in python

## Features

- Internal quad and multi-precision arithmetic support for high precision gene scoring via exact CDF calculations (up to 100 digits)
- Fast random access to SNP reference panel genomic data with minimal memory footprint.
- Parallelization over chromosomes and/or genes
- Gene-wise coherence test between two GWAS

## Documentation

The full documentation can be found [here](https://bergmannlab.github.io/PascalX/index.html). For quickstart, continue below.

## Citation policy

If you make use of PascalX for your research, please cite PascalX via the doi: 10.5281/zenodo.4429922

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4429922.svg)](https://doi.org/10.5281/zenodo.4429922)

If you make use of the X-scorer (gene-wise coherence test between two GWAS), please cite the work 

*Krefl D., Bergmann S.*  
*Covariance of Interdependent Samples with Application to GWAS*  
[doi.org/10.1101/2021.05.16.21257289](https://doi.org/10.1101/2021.05.16.21257289)




## Installation 

### Debian/Ubuntu

Requirements:
  - Python3 with development headers
  - g++, make
  - [BOOST libraries](https://www.boost.org) 
  
  Install of requirements on Debian/Ubuntu:
  ``` bash
  sudo apt-get install python3 python3-dev python3-setuptools python3-pip g++ make libboost-all-dev
  ```

Export library path:
```bash
export LD_LIBRARY_PATH="/yourpath/PascalX/build/lib:$LD_LIBRARY_PATH"
```

Install PascalX via terminal:
```bash
git clone https://github.com/BergmannLab/PascalX.git
cd PascalX
make all
cd python
python3 setup.py install 
```

### Docker

Build the docker image 
```bash
git clone https://github.com/BergmannLab/PascalX.git
cd PascalX
docker build . -t pascalx:latest
```

Run the image in interactive mode with the host directory ```/your/workdir``` mounted as ```/data``` using the command
```bash
docker run --mount src=/your/workdir,target=/data,type=bind -p 8888:8888 -it pascalx bash
```
Jupyter notebook comes pre-installed and listens on port ```8888```.

## Usage

### Initialization:
**Import:**
```python
from PascalX import genescorer

Scorer = genescorer.chi2sum()
```

Options:

```window=50000```   : Window to add to gene start and end position

```varcutoff=0.99```: Variance to keep

```MAF=0.05```       : MAF threshold



**Set reference panel:**

1000 Genome Project reference data can be downloaded and converted via executing the script in the PascalX/misc folder as below (for GRCh37 replace 38 with 37).
```bash
bash get1KGGRCh38.sh pathtostore/ EUR 4
```
Note that the last number specified the # of cpu cores to utilize. 

The reference data has to be loaded by the gene scorer using the method
```python
Scorer.load_refpanel("path/filename")
```
The filename is required to not contain the ending ```.chr#....```. If the corresponding reference data has not been imported yet, PascalX will try to import the data from ```filename.chr#.tped.gz``` or ```filename.chr#.vcf.gz``` files in ```path/```. For .tped import the genotype information has to be supplied in gzip compressed 1-2-coded [plink](https://www.cog-genomics.org/plink/) tped files. The following plink options should do the job:

```bash
--recode 12 transpose
```

By default PascalX uses only one cpu core for the import. The number of cores to utilize can be set via the option ```parallel=```. For import of allele information raw .vcf files have to be used. In this case the plink step can be skipped. Note that later calls of ```load_refpanel``` will be fast as the converted reference data will be stored on disk and reloaded on the fly. 



**Load gene annotation:**
```python
Scorer.load_genome("path/filename",ccol=1,cid=5,csymb=7,cstx=3,cetx=4,cs=2,chrStart=3,NAgeneid='n/a',header=True)
```

```filename```: Text file containing gene information. Columns need to be separated by tabs (```\t``` delimiter)

```ccol```: Column containing chromosome 

```cid```: Column containing gene id

```csymb```: Column containing gene symbol

```cstx```: Column containing gene start position

```cetx```: Column containing gene end position

```cs```: Column containing strand

```chrStart```: # of first characters to skip in ccol 

```NAgeneid```: Identifier for not available gene id

```header```: First line is header ( True | False )


Note that PascalX can download human genome annotation automatically from [ensembl.org BioMart](https://www.ensembl.org/biomart/martview/). 

```python
from PascalX import genome

D = genome.genome()
D.get_ensembl_annotation(filename,genetype='protein_coding',version='GRCh38')
```

``` filename```: File to store the downloaded annotation 

``` genetype```: String of comma separated ensembl gene types to download

``` version ```: Annotation version to download ( GRCh37 | GRCh38 ) 

WARNING: 

PascalX matches genes with rsids via position overlap in the loaded reference panel. Both datasets should be based on the same annotation (for instance both hg19) for consistency.

**Load GWAS summary statistics to score:**
```python
Scorer.load_GWAS("path/filename",rscol=0,pcol=1,a1col=None,a2col=None)
```

```filename```: Text file containing variant ids and p-values with columns separated by tabs (```\t``` delimiter)

```rscol```: Column containing rsids

```pcol``` : Column containing p-values

```a1col```: Column containing alternate allele information (None for ignoring)

```a2col```: Column containing reference allele information (None for ignoring)


### Gene scoring:

**Exampe 1:** Score all genes in annotation

```python
R = Scorer.score_all()
```
**Example 2:** Score genes on chromosomes 21 and 22

```python
R = Scorer.score_chr(chrs=[21,22])
```

**Example 3:** Score genes WDR12 and FARP2

```python
R = Scorer.score(['WDR12','FARP2'])
```

**Options:**

```parallel=1```  : Number of cores to use

```nobar=False``` : Disable progress bar

**Return:**

```python
R = [R_SUCCESS,R_FAIL,R_TOTALFAIL]
```

```python
R_SUCCESS   = [ ['Symbol',p-value,Nsnps],...]
R_FAIL      = [ ['Symbol',[infos] ]     ,...]
R_TOTALFAIL = [ ['Symbol','Reason']     ,...] 
```

```R_SUCCESS```  : List of successfully scored genes

```R_FAIL```     : List of genes for which the scoring failed due to non-convergence of the scoring algorithm

```R_TOTALFAIL```: List of genes for which the scoring failed for other reasons


**Re-scoring:**

The genes in R_FAIL can be scored again with a manual choice of algorithm:

```python
R = Scorer.rescore(R,method='ruben',mode='128b',reqacc=1e-32,intlimit=100000)
```

```method= 'auto' | 'ruben' | 'davies' | 'satterthwaite'```: Algorithm to use 

```mode='' | '128b' | '100d'```: internal precision (double | quad | 100 digits)

```reqacc=float``` : Requested accuracy. 

```intlimit=integer```: Cutoff

NOTE:

```ruben``` and ```davies``` compute exactly up to requested precision (```reqacc```). ```satterthwaite``` is a second order approximation. ```auto``` tries to automatically select for given gene between davies and ruben to maximize throughput. ```auto``` is the default setting for the gene scorer. If it fails, it is recommended to rescore with ```ruben```.


TIP:

Ruben converges very slowly if the ratio between the largest and smallest eigenvalue is large. Try to reduce the ```varcutoff``` parameter in this case.

**Output:**

The scoring results ```R_SUCCESS``` can be saved in a a tab separated file via:

```python
scorer.save_scores('filename')
```


### Pathway scoring:

PascalX offers two different pathway scorers. 

**Initialization:**

Define a gene scorer as above and score or load scored genes for a GWAS. The scorers are then initialized as follows:

Rank based scoring

```python
Pscorer = pathway.chi2rank(scorer)
```

The rank scorer uniformizes the gene p-value distribution via ranking and aggregates p-values via inverse transform to chi^2 distributed random variables.


Monte-Carlo based scoring

```python
Pscorer = pathway.chi2perm(scorer)
```
Gene p-values are chi^2 inverse transformed and the sum of chi^2 values for a given pathway is compared against a randomly generated gene sets of equal size.


**Loading modules:**

Sets of modules to be scored can be loaded from a tab-separated file via the command

```python
M = P.load_modules("filename.tsv",ncol=0,fcol=2)
```

```ncol=0```: Column of name/id 

```fcol=2```: Column of first gene symbol

**Scoring:**

```python
RESULT = P.score(M)
```

```RESULT = [ [[name,[symbols],p-value],...], R_FAIL ]```

If the list R_FAIL is not empty, a gene re-scoring for the listed genes and re-scoring of affected pathways is recommended. Note that the genes and meta-genes with out gene score are removed from the pathway before pathway scoring.

### X scoring:

PascalX offers two different GWAS cross scorers.

**Coherence scorer:**

```python
from PascalX import xscorer

X = xscorer.zsum(leftTail=False)
X.load_genome('path/filename')
```
Note that the default initialization of the gene scoring above are used. ```leftTail=``` sets the side to test. False corresponds to anti-coherence and True to coherence. A gene annotation has to be loaded as for the standard Genescorer.

```python
X.load_GWAS('path/filenameA',name='GWAS A',rscol=0,pcol=1,bcol=2,a1col=3,a2col=4,header=False)
X.load_GWAS('path/filenameB',name='GWAS B',rscol=0,pcol=1,bcol=2,a1col=3,a2col=4,header=False)
```
In the GWAS data loading routine, we have to set in addition a name for each GWAS to be loaded via the ```name=``` argument, and it is necessary to specify the column with the raw betas ```bcol=```.

It is recommended to perform the scoring after filtering for matching alleles 

```python
X.matchAlleles('GWAS A','GWAS B')
```

This requires that allele information has been loaded for both GWAS, and for the reference panel (using .vcf import).

Note that it is also recommended to jointly QQ normalize p-values:

```python
X.jointlyRank('GWAS A','GWAS B')
```

The scoring is started via calling

```python
R = X.score_all(E_A='GWAS A',E_B='GWAS B')
```
The return ```R``` is as for the Genescorer class.

**Ratio scorer:**

As above, but with

```python
X = xscorer.rsum(leftTail=False)
```

NOTE: As the current cross scoring implementation consumes significantly more memory than the genescorer, it is recommended to keep ```parallel=1``` at the time being.


### Visualization:

```python
Scorer.plot_Manhattan(R[0])
```

**Options:**

```sigLine=pval```: Draw significance threshold line at p-value (0 for off)

```logsigThreshold=logpval```: Threshold above which to label genes (log p-value)

```labelSig=True```: Plot names of genes above logsigThreshold (True | False)

```labelList=[]```: List of gene symbols to label

```style='colorful'```:
![Colorful style](https://github.com/BergmannLab/PascalX/blob/main/misc/style_colorful.png)

```style='classic'```: 
![Classic style](https://github.com/BergmannLab/PascalX/blob/main/misc/style_classic.png)
