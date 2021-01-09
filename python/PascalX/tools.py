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

import progressbar
import urllib
from urllib.request import urlretrieve

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