#### VCF Stats ####

Pipeline to run basic statistics on vcf. 

Right now it executes the following functions:

    - Filters variants to keep only bi-allelic SNPS
    - Detects paralogs


To do :

    - Corrects genotypes based on callability by sample
    - Subsample vcf based on a given list of samples
    - Generates SFS
    - Measure Pi
    - ...


The minimum you need is a :

    - raw vcf, zipped and indexed
    - fai-like file giving the id and size of chr
    - file with list of samples to analyze (see format in toy_example/)
    - bed file summarizing callability (snpArcher output)
    - ngsparalogs locally installed and compiled
