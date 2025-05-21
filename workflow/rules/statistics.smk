# Mesure pi, fis, tajima's D, fst, pi0/pi4

    # see /projects/plantlp/02_VCF_PROCESSING/scripts/get_pi.sh for pi measures
rule pi:
    """ pi by sites """

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