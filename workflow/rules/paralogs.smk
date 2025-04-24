# Detect paralogs with ngsparalogs
include:  "common.smk"

wildcard_constraints:
    chunk = "|".join(["%.2d" % i for i in range(int(config["chunks_count"]))])

rule get_snps_list:
    """
    List of all snp positions
    """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf.gz",
        vcf_idx = "results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf.gz.tbi"
    output:
        snps_list = "results/paralogs/{prefix}.{chr_id}.snps_list.bed"
    conda:
        "../envs/paralogs.yml"
    shell:
        """
        bcftools view --no-header -r {wildcards.chr_id} {input.vcf} | \
            awk -F '\t' '{{print $1 "\t" $2}}' > {output.snps_list}
        """

rule split_snps_list:
    input:
        snps_list = "results/paralogs/snps_list/{prefix}.{chr_id}.snps_list.bed",
        fai = "results/raw/stats/{prefix}.{chr_id}.fai"
    output:
        chunk = "results/paralogs/snps_list/{prefix}.{chr_id}.snps_list.{chunk}.bed"
    params:
        chunk = config["chunks_count"]
    shell:
        """
        split -n {params.chunk} -d --additional ".bed" \
            {input.snps_list} \
            results/paralogs/snps_list/{wildcards.prefix}.{wildcards.chr_id}.snps_list.
        """

rule ngs_paralog:
    """
    Main ngsparalog command
    """
    input:
        bam_list = config["bams_list"],
        chunk = "results/paralogs/snps_list/{prefix}.{chr_id}.snps_list.{chunk}.bed",
        ngsparalogs_path = config["ngsparalog_path"]
    output:
        paralogs_lr = "results/paralogs/{prefix}.{chr_id}.{chunk}.lr"
    conda:
        "../envs/paralogs.yml"
    log:
        "logs/paralogs/{prefix}.{chr_id}.{chunk}.log"
    params:
        minQ = config["minQ"],
        minind = config["minind"],
        mincov = config["mincov"]
    shell:
        """
        samtools mpileup -b {input.bam_list} -l {input.chunk} -q 0 -Q 0--ff UNMAP,DUP -d 0 | \
	        {input.ngsparalogs_path}/ngsParalog calcLR -infile - \
                -outfile {output.paralogs_lr} \
                -minQ {params.minQ} \
                -minind {params.minind} \
                -mincov {params.mincov}
        """

rule chi2:
    """
    Chi2 test for each site
    """
    input:
        paralogs_lr = "results/paralogs/{prefix}.{chr_id}.{chunk}.lr",
        script = workflow.source_path("../scripts/get_paralogs.py")
    output:
        paralogs_pval = "results/paralogs/{prefix}.{chr_id}.{chunk}.pval.tsv"
    conda:
        "../envs/paralogs.yml"
    shell:
        """
        python3 {input.script} -i {input.paralogs_lr} -o {output.paralogs_pval}
        """

rule merge_chunks:
    """
    """
    input:
        chunks = get_chunk_files
    output:
        merged_lr = "results/paralogs/lr/{prefix}.{chr_id}.pval.tsv"
    shell:
        "cat {input.chunks} > {output.merged_lr}"


rule get_paralogs:
    """
    Mark sites below 0.05 as sites belonging to a paralogous region
    """
    input:
        pval = "results/paralogs/lr/{prefix}.{chr_id}.pval.tsv"
    output:
        bed = "results/paralogs/bed/{prefix}.{chr_id}.paralogs.bed"
    shell:
        """
        cat {input.pval}| awk '$8 == "True" {{print $1 "\t" $2-1 "\t" $2}}' > {output.bed}
        """
