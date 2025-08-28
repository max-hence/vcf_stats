# Mesure pi, thetaW, tajima's D, pi0/pi4, sfs

rule remove_paralogs:
    """ Remove paralogs from classic vcf """

    input: 
        vcf= "results/callability/vcf/{prefix}.callability.{chr_id}.vcf.gz",
        paralogs = "results/paralogs/bed/{prefix}.{chr_id}.paralogs.bed"
    output: 
        vcf = temp("results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf"),
        vcf_gz = "results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf.gz",
        vcf_idx = "results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        bcftools view -T ^{input.paralogs} {input.vcf} -o {output.vcf}
        bgzip < {output.vcf} > {output.vcf_gz} && tabix -p vcf {output.vcf_gz}
        """

rule sfs_projection:
    """
        Run easySFS to run SFS projection for each sample size (2:2n)
    """
    input:
        vcf_snps = "results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf",
        pop_path = config["pop_path"],
        easySFS = config["easySFS_path"],
    output:
        preview = "results/sfs/{prefix}.{chr_id}.preview.txt",
    conda:
        "../envs/easySFS.yml"
    log:
        "logs/{prefix}.{chr_id}.log"
    shell:
        """
            python3 {input.easySFS} -i {input.vcf_snps} -p {input.pop_path} \
            --preview -v -a > {output.preview}
        """

rule get_best_params:
    """
        merge all results from easySFS run by chr and find best params based on two methods
        - one maximizes nbr of snps
        - one maximizes a likelihood value (including n_snps and n_indiv as params)
    """
    input:
        previews = get_previews("results/sfs/{prefix}.{chr_id}.preview.txt"),
        script = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample_size = "results/sfs/{prefix}.best_sample.txt"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.log"
    shell:
        """
            python3 {input.script} -i {input.previews} --method ml -o {output.best_sample_size}
        """

rule sfs:
    """
        easySFS with previously calculated sampling size
    """
    input:
        vcf = "results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf",
        best_sample = "results/sfs/{prefix}.best_sample.txt",
        pop_path = config["pop_path"],
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = directory("results/sfs/{prefix}.{chr_id}"),
        final_sfs = "results/sfs/{prefix}.{chr_id}.sfs"
    conda:
        "../envs/easySFS.yml"
    log:
        "logs/{prefix}.{chr_id}.log"
    shell:
        """
            sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 )))
            python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
            --dtype 'int' \
            --proj $sampling_size \
            --ploidy 2 \
            -v -a -f -o {output.sfs_dir}

            n_seq=$(sed -n '2 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk '{{print $2}}')
            echo $(sed -n '3 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk -v n_seq="$n_seq" '{{ for(i=2; i<=n_seq/2+1; i++) printf $i" "}}') \
                > {output.final_sfs}
        """

rule vcf2allsites:
    """ Prepare input for pixy """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf",
        bed = "results/raw/bed/{prefix}.{chr_id}.callable.bed",
        best_sample = "results/sfs/{prefix}.best_sample.txt",
        script = workflow.source_path("../scripts/vcf2allsites.py")
    output:
        vcf = temp("results/callability/vcf/{prefix}.allsites.{chr_id}.vcf"),
        vcf_gz = "results/callability/vcf/{prefix}.allsites.{chr_id}.vcf.gz",
        vcf_gz_idx = "results/callability/vcf/{prefix}.allsites.{chr_id}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        python3 {input.script} \
            -i {input.vcf} \
            --callability \
            -b {input.bed} \
            -t $sampling_size \
            -o {output.vcf}

        bgzip < {output.vcf} > {output.vcf_gz} && tabix -p vcf {output.vcf_gz}
        """


rule remove_paralogs_allsites:
    """ Remove paralogs from allsites vcf"""
    input: 
        vcf= "results/callability/vcf/{prefix}.allsites.{chr_id}.vcf.gz",
        paralogs = "results/paralogs/bed/{prefix}.{chr_id}.paralogs.bed"
    output: 
        vcf = temp("results/paralogs/vcf/{prefix}.no_paralogs.allsites.{chr_id}.vcf"),
        vcf_gz = "results/paralogs/vcf/{prefix}.no_paralogs.allsites.{chr_id}.vcf.gz",
        vcf_idx = "results/paralogs/vcf/{prefix}.no_paralogs.allsites.{chr_id}.vcf.gz.tbi"
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
        vcf = "results/paralogs/vcf/{prefix}.no_paralogs.allsites.{chr_id}.vcf.gz",
        vcf_idx = "results/paralogs/vcf/{prefix}.no_paralogs.allsites.{chr_id}.vcf.gz.tbi",
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

rule fis:
    #TODO: Utiliser bcftools view au lieu de ziper deziper à chaque fois ?
    """ fis by window """
    input:
        vcf_gz = "results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf.gz",
        vcf = "results/paralogs/vcf/{prefix}.callability.no_paralogs.{chr_id}.vcf",
        script = workflow.source_path("../scripts/fis_by_window.py")
    output:
        fis = "results/stats/{prefix}.{chr_id}.fis.tsv"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} \
            -i {input.vcf} \
            -s $(seq -s, 1 $(bcftools query -l {input.vcf_gz} | wc -l )) \
            -m 2 \
            -o {output.fis}
        """
