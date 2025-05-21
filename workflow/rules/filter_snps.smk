include: "common.smk"

rule split_by_chr:
    """
        split vcf by chr
    """
    input:
        vcf = config["vcf_path"]
    output:
        splitted_vcf = temp("results/raw/vcf/{prefix}.raw.{chr_id}.vcf.gz"),
        splitted_vcf_idx = temp("results/raw/vcf/{prefix}.raw.{chr_id}.vcf.gz.tbi"),
        stats = "results/raw/stats/{prefix}.raw.{chr_id}.fai"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{chr_id}.log"
    shell:
        """
        bcftools view -Oz -r {wildcards.chr_id} -o {output.splitted_vcf} {input.vcf}
        tabix -p vcf {output.splitted_vcf}
        bcftools index -s {output.splitted_vcf} > {output.stats}
        """

rule split_bed:
    input:
        bed = config["bed_path"]
    output:
        splitted_bed = temp("results/raw/bed/{prefix}.raw.{chr_id}.callable.bed")
    log:
        "logs/{prefix}.{chr_id}.log"
    shell:
        """
        awk -v k="{wildcards.chr_id}" '$1 == k' "{input.bed}" > {output.splitted_bed}
        """

    #################################
    ### Filter on Bi-allelic SNPs ###
    #################################

rule filter_snps:
    """
        Remove indels and MNP
    """
    input:
        splitted_vcf = "results/raw/vcf/{prefix}.raw.{chr_id}.vcf.gz",
        splitted_vcf_idx = "results/raw/vcf/{prefix}.raw.{chr_id}.vcf.gz.tbi"
    output:
        vcf_snps = temp("results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf"), # need unziped vcf for later
        vcf_snps_gz = "results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf.gz",
        vcf_snps_idx = "results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf.gz.tbi",
        snps_stats = "results/snps/stats/{prefix}.SNPS.{chr_id}.stats"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{chr_id}.log"
    shell:
        """
        bcftools view -Oz -m2 -M2 -v snps -o {output.vcf_snps} {input.splitted_vcf}
        bgzip -c {output.vcf_snps} > {output.vcf_snps_gz}
        tabix -p vcf {output.vcf_snps_gz}
        bcftools index -s {output.vcf_snps_gz} > {output.snps_stats}
        """

#TODO: Un truc comme ça, j'ai récup d'un de mes autres snakemake donc surment des erreurs sur les paths et tout
rule correct_genotype:
    """
        Changes wrongly called genotype into NA
    """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf",
        splitted_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/correct_genotype.py")
    output:
        corrected_vcf = temp("results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf"),
        corrected_vcf_gz = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf.gz",
        corrected_vcf_idx = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} -i {input.vcf} -b {input.splitted_bed} -o {output.corrected_vcf}
        bgzip < {output.corrected_vcf} > {output.corrected_vcf_gz} && tabix -p vcf {output.corrected_vcf_gz}
        """