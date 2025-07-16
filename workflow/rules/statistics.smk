# Mesure pi, fis, tajima's D, fst, pi0/pi4

    # see /projects/plantlp/02_VCF_PROCESSING/scripts/get_pi.sh for pi measures
rule pi:
    """ pi by sites """
    input:
        vcf = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz",
        vcf_idx = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz.tbi",
    output:
        sites_pi = "results/stats/pi/{prefix}.{chr_id}.sites.pi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        vcftools --gzvcf {input.vcf} --site-pi --out "results/stats/{wildcards.prefix}.{wildcards.chr_id}"
        """

rule pi_by_window:
    """ pi by 100kb windows """
    input:
        sites_pi = "results/stats/pi/{prefix}.{chr_id}.sites.pi",
        script = workflow.source_path("../scripts/pi_by_window.py")
    output:
        wdw_pi = "results/stats/pi/{prefix}.{chr_id}.wdw.pi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} \
            -i {input.sites_pi} \
            -w 100000 \
            -o {output.wdw_pi}
        """

rule :
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


rule sfs_preview:
    """ /projects/plantlp/02_VCF_PROCESSING/scripts/easySFS.sh """

rule sfs:
    """ sfs by chromosomes """
    # see /projects/plantlp/02_VCF_PROCESSING/scripts/get_sfs.sh