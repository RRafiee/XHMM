#!/bin/bash -e
# 
# Dr Reza Rafiee, January 2017, Newcastle University, Northern Institiute for Cancer Research
# This pipeline provides an implementation of XHMM (eXome-Hidden Markov Model)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tput bold
echo "Reza Rafiee 2016-2017"
echo -e "Running XHMM pipeline ...\n"
tput sgr0

# Load settings for this run
source XHMMsettings.sh

# The XHMM C++ software suite was written to call copy number variation (CNV)
# from next-generation sequencing projects, where exome capture was used (or targeted sequencing, more generally).  
# XHMM uses principal component analysis (PCA) normalization and a hidden Markov model (HMM) to detect
# and genotype copy number variation (CNV) from normalized read-depth data from targeted sequencing experiments.
# XHMM was explicitly designed to be used with targeted exome sequencing at high coverage (at least 60x - 100x)
# using Illumina HiSeq (or similar) sequencing of at least ~50 samples.
# However, no part of XHMM explicitly requires these particular experimental conditions,
# just high coverage of genomic regions for many samples. Please see the following link for more details:
# http://atgu.mgh.harvard.edu/xhmm/tutorial.shtml
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ################################# Input #####################################
# 1) list of exome targets ("EXOME.interval_list", one column format of exome regions of one specific chromosome, e.g., 10:92895-94177)
# 2) exome sequencing reads ("DATA.RD.txt", matrix of all coverage means for all EXOME.interval_list regions of all samples, obtained from GATK pipeline)
# 3) Markov model parameters ("params.txt", all 9 model parameters)
# 
# It is assumed that the BAM files of "analysis-ready reads" being used have been generated using
# the "best practice" bwa + Picard + GATK pipeline, though variations on this pipeline,
# and even other pipelines to produce the BAM files for XHMM input, have been
# successfully employed by various users.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Base dir - should auto set to where this script resides
#BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

## System settings for launching java jobs
# On FMS cluster we need to use large pages have also set tmp dir to one
# provided by SoGE for each run

# Add in module for Java 1.8 (FMS cluster specific)
module add apps/java/jre-1.8.0_92

JAVA="/opt/software/java/jdk1.8.0_92/bin/java -XX:-UseLargePages -Djava.io.tmpdir=$TMPDIR"

## We need latest GCC libs for AVX hardware acceleration of pairHMM (FMS cluster
# specific)
module add compilers/gnu/4.9.3
## Latest version of R for plots (FMS cluster specific)

module add apps/R/3.3.1

# Optional 
# Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci):


# Till this stage, I already provided DATA.RD.txt (processig in R)


# Determine number of samples in master list
N=$(wc -l DATA.RD.txt | cut -d ' ' -f 1)
echo -e " - No of samples = $N-1\n"
tput sgr0

# main code starting from here
tput bold
echo -e "run GATK to calculate the per-target GC content and create a list of the targets with extreme GC content:\n"
#####################################################################################################################
# Optionally, run GATK to calculate the per-target GC content and create a list of the targets with extreme GC content:
java -Xmx3072m -jar /opt/software/bsu/bin/GenomeAnalysisTK-3.7.jar \
-T GCContentByInterval -L EXOME.interval_list \
-R /opt/databases/GATK_bundle/2.8/b37/human_g1k_v37_decoy.fasta \
-o DATA.locus_GC.txt
cat DATA.locus_GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' \
> extreme_gc_targets.txt
#####################################################################################################################
echo -e "run XHMM to filter .....:\n"
#####################################################################################################################
# Original, working without low_complexity_targets
# Filters samples and targets and then mean-centers the targets:
./xhmm --matrix -r ./DATA.RD.txt --centerData --centerType target \
-o ./DATA.filtered_centered.RD.txt \
--outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeTargets ./extreme_gc_targets.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

# I exclude the low_complexity_targets argument
#--excludeTargets ./extreme_gc_targets.txt --excludeTargets ./low_complexity_targets.txt \

#####################################################################################################################
echo -e "PCA on mean-centered data .....:\n"
# Runs PCA on mean-centered data:
./xhmm --PCA -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA
#####################################################################################################################
# Normalizes mean-centered data using PCA information:
./xhmm --normalize -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA \
--normalizeOutput ./DATA.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
####################################################################################################################
# Filters and z-score centers (by sample) the PCA-normalized data:
./xhmm --matrix -r ./DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
--maxSdTargetRD 30

##################################################################################################################### 
##################################################################################################################### 
#Filters original read-depth data to be the same as filtered, normalized data:
./xhmm --matrix -r ./DATA.RD.txt \
--excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o ./DATA.same_filtered.RD.txt
#####################################################################################################################
# ***********************************************************************
# Input CNV parameters file:
# ***********************************************************************
# 1e-08   6       70      -3      1       0       1       3       1
# ***********************************************************************
# translates into XHMM parameters of:
# ***********************************************************************
# Pr(start DEL) = Pr(start DUP) = 1e-08
# Mean number of targets in CNV [geometric distribution] = 6
# Mean distance between targets within CNV [exponential decay] = 70 KB
#
# DEL read depth distribution ~ N(mean=-3, var=1)
# DIP read depth distribution ~ N(mean=0, var=1)
# DUP read depth distribution ~ N(mean=3, var=1)
#***********************************************************************
#####################################################################################################################
# Discovers CNVs in normalized data:
./xhmm --discover -p ./params.txt \
-r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt \
-c ./DATA.xcnv -a ./DATA.aux_xcnv -s ./DATA
