.. PascalX documentation master file, created by
   sphinx-quickstart on Wed Nov 18 11:20:06 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PascalX's documentation!
===================================

PascalX is a python3 library (`source <https://github.com/BergmannLab/PascalX>`_) for high precision gene and pathway scoring for GWAS summary statistics. Aggregation of SNP p-values to gene and pathway scores follows the `Pascal <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004714>`_ methodology, which is based on :math:`\chi^2` statistics. The cummulative distribution function of the weighted :math:`\chi^2` distribution is calculated exactly via a multi-precision C++ implementation of Ruben's and Davies algorithm. This allows to apply the Pascal methodology to modern UK BioBank scale GWAS. In addition, PascalX offers a novel coherence test between two different GWAS on the level of genes, based on the product-normal distribution, as described `here <https://doi.org/10.1101/2021.05.16.21257289>`_.    

**Highlights:**

    * Multi-precision weighted :math:`\chi^2` cdf calculation (up to 100 digits) 
    * Parallelization over chromosomes and/or genes
    * Fast random access to reference panel genomic data via own indexed SNP database
    * Gene-wise coherence test between two GWAS
    * Tissue enrichment test (experimental)
    
    
.. warning::
    PascalX is a research level tool. No warranty or guarantee whatsoever for its correct functionality is given. You should perform your own consistency checks on results PascalX implies.


**Why PascalX:**

In order to illustrate the impact of approximating the cdf of the weighted :math:`\chi^2` distribution at very small p-values via matching only the first two moments, we calcuated for an uniform random sample of weights and arguments :math:`\log_{10}` transformed p-values at :math:`N=10`, :math:`N=100` and :math:`N=1000` degrees of freedom via Ruben's algorithm at 100 digits precision and the Satterthwaite-Welch approximation. The resulting p-values are ordered, :math:`-\log_{10}` transformed, and plotted against each other (QQ-plots), see top row of the following figure. 

.. image:: SatterthwaiteVsRuben.png

In the second row, the ratio of p-values is plotted at various fixed arguments, sampled from different weights. The bullet point marks the mean and the error bar the range of values observed. We infer that the Satterthwaite-Welch approximation tends to overestimate the significance, in particular for small :math:`N`. Note that the error can go beyond several orders of magnitude.  




.. toctree::
   :maxdepth: 3
   :caption: Contents

   install
   usage
   PascalX
..


Index
==================

* :ref:`genindex`


License and citation policy
---------------------------

The PascalX package is an open-source package under AGPLv3. 

If you make use of PascalX for your research, please cite PascalX via the doi: 10.5281/zenodo.4429922

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4429922.svg
   :target: https://doi.org/10.5281/zenodo.4429922

If you make use of the X-scorer (gene-wise coherence test between GWAS), please cite the work:


| *Krefl D., Bergmann S.* 
| *Covariance of Interdependent Samples with Application to GWAS* 
| *doi:10.1101/2021.05.16.21257289*


Contact
-------

* For all technical issues (bug reports, etc.), please open a ticket on the `GitHub <https://github.com/BergmannLab/PascalX/issues>`_ page. 

* For general inquiries, please contact the main contributor: Daniel.Krefl(@unil.ch)
