#### VCF Stats ####

Pipeline to run basic statistics on vcf. 

Right now it executes the following functions:

    - Filters variants to keep only bi-allelic SNPS
    - Detects paralogs

To do :

    - Corrects genotypes based on callability by sample
        - from callable_all_sample.bed, set every badly called genotype as missing
    - Generates SFS
        - with easySFS
    - Statistics
        - pi by 100kb windows
        - Mean pi by chromosome
        - pi0/pi4 from gene annotation
        - thetaW (n snps) by 100kb windows and by chromosomes
        - Tajima's D

The minimum you need is a :

    - raw vcf, zipped and indexed (snpArcher output)
    - fai-like file giving the id and size of chr you want to analyse
    - list of bam paths
    - bed file summarizing callability (snpArcher output)
    - ngsparalogs locally installed and compiled
