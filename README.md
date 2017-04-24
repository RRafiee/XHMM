# Discovery of copy number variations (CNVs) from exome read depth using XHMM (eXome-Hidden Markov Model) 

Call copy number variation (CNV) from next-generation sequencing data, where exome capture was used (or targeted sequencing, more generally).

Input files:  

1) list of exome targets ("EXOME.interval_list", one column format of exome regions of one specific chromosome, e.g., 10:92895-94177)
2) exome sequencing reads ("DATA.RD.txt", matrix of all coverage means for all EXOME.interval_list regions of all samples, obtained from GATK pipeline)
3) Markov model parameters ("params.txt", all 9 model parameters)

XHMM was developed assuming compilation on Linux with gcc version 4.4. Please see the following link for original source code
http://atgu.mgh.harvard.edu/xhmm/download.shtml
