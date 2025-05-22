# Mesure pi, fis, tajima's D, fst, pi0/pi4

    # see /projects/plantlp/02_VCF_PROCESSING/scripts/get_pi.sh for pi measures
rule pi:
    """ pi by sites """
    input: 
        vcf = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz"

rule pi_by_window:
    """ pi by 100kb windows """

rule mean_pi:
    """ mean pi by chromosomes """

rule sfs_preview:
    """ /projects/plantlp/02_VCF_PROCESSING/scripts/easySFS.sh """
rule sfs:
    """ sfs by chromosomes """
    # see /projects/plantlp/02_VCF_PROCESSING/scripts/get_sfs.sh

rule thetaW:
    """ Number of segregating sites """

    # n_snp/sum([1/i for i in range(1, n_indiv)])

rule Dtajima:
    """ Tajima's D statistic, per 100kb window or per chrom ?"""

    # Selon wikipedia ...
    # a1 <- sum(1 / n)
    # a2 <- sum(1 / n ^ 2)
    # b1 <- (n + 1) / (3 * (n - 1))
    # b2 <- 2 * (n ^ 2 + n + 3) / (9 * n * (n - 1))
    # c1 <- b1 - 1 / a1
    # c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1 ^ 2
    # e1 <- c1 / a1
    # e2 <- c2 / (a1 ^ 2 + a2)

    # D_obs <- (pi - S / a1) / sqrt(e1 * S + e2 * S * (S - 1))

    # Il faut d'abord récupérer pi et S... Listes avec les mêmes dimensions :
    # si par window :  chr_id   start   end pi & chr_id   start   end S
    # si par chr : chr_id   pi  & chr_id    S