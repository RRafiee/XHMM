# Discovery of copy number variations (CNVs) from exome read depth using XHMM (eXome-Hidden Markov Model) 

Call copy number variation (CNV) from next-generation sequencing data, where exome capture was used (or targeted sequencing, more generally).

Input files:  

1) list of exome targets ("EXOME.interval_list", one column format of exome regions of one specific chromosome, e.g., 10:92895-94177)
2) exome sequencing reads ("DATA.RD.txt", matrix of all coverage means for all EXOME.interval_list regions of all samples, obtained from GATK pipeline)
3) Markov model parameters ("params.txt", all 9 model parameters)

Please see the following link for original source code: http://atgu.mgh.harvard.edu/xhmm/download.shtml

# Installing/using the XHMM R scripts

The XHMM R code is now available as an R library. There are still two options for using the scripts: 

# Newer (preferred) option

Install the xhmmScripts package via one of the following methods (you only need to do this once):

1) Compile it from the main XHMM source code downloaded from the Git repository above:
cd statgen-xhmm-*
make R
Then, run in R:
install.packages(list.files(path=".", pattern="xhmmScripts_.+\\.tar\\.gz"), repos=NULL, type="source")

2) Download and install the package from CRAN using the install.packages() command in R:
install.packages("xhmmScripts")

3) Download it from the xhmmScripts page at CRAN. Then, run in R:
install.packages(list.files(path=".", pattern="xhmmScripts_.+\\.tar\\.gz"), repos=NULL, type="source")

To use the R code, preface your R scripts with:
library(xhmmScripts)
