# Mesure pi, thetaW, tajima's D, pi0/pi4

rule filter_vcf:
    """ Filtre le vcf en enlevant les paralogs """
    input: 
        vcf= "results/callability/vcf/{prefix}.SNPS.NA.{chr_id}.vcf.gz",
        paralogs = "results/paralogs/bed/{prefix}.{chr_id}.paralogs.bed"
    output: 
        vcf = temp("results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf"),
        vcf_gz = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz",
        vcf_idx = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        bcftools view -T ^{input.paralogs} {input.vcf} -o {output.vcf}
        bgzip < {output.vcf} > {output.vcf_gz} && tabix -p vcf {output.vcf_gz}
        """

rule get_stats:
    """ pi, thetaW and tajimaD by 100kb windows """
    input:
        vcf = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz",
        vcf_idx = "results/paralogs/vcf/{prefix}.SNPS.NA.no_paralogs.{chr_id}.vcf.gz.tbi",
        pop_path = config["pop_path"]
    output:
        pi = "results/stats/{prefix}.{chr_id}_pi.txt",
        theta_W = "results/stats/{prefix}.{chr_id}_watterson_theta.txt",
        tajimaD = "results/stats/{prefix}.{chr_id}_tajima_d.txt"
    conda:
        "../envs/pixy.yml"
    shell:
        """
        pixy \
            --stats "pi" "watterson_theta" "tajima_d" \
            --vcf {input.vcf} \
            --populations {input.pop_path} \
            --window_size 100000 \
            --output_folder results/stats/ \
            --output_prefix {wildcards.prefix}.{wildcards.chr_id} \
            --bypass_invariant_check \
            --n_cores 4
        """
