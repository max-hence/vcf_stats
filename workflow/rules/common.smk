from glob import glob

def check_fai_format(fai_path:str):
    with open(fai_path, "r") as fai:
        if len(fai_path.readline().split("\t")) != 2:
            raise (WorkflowError(" fai must contain only two tab-separated columns : chr_id\tchr_length"))


def get_chr_list(fai_path:str):
    """Return the list of scf that will analyzed

    Args:
        fai_path (str): tab separated table (fai format). 
                        Genome info, with only scf that will be analyzed
    Returns:
        list of scaffolds of interest
    """

    return [row.split('\t')[0] for row in open(fai_path, "r")]


def get_output():
    final_prefix = config["final_prefix"]
    chromosomes: list = get_chr_list(config["fai_path"])
    out: list  = []
    chunks = [i for i in range(int(config["chunks_count"]))]

    if final_prefix == "":
        raise (WorkflowError("'final_prefix' is not set in config."))

    # vcf filtered to keep only SNPs
    out.extend(expand("results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf.gz", prefix=final_prefix, chr_id=chromosomes))
    out.extend(expand("results/snps/vcf/{prefix}.SNPS.{chr_id}.vcf.gz.tbi", prefix=final_prefix, chr_id=chromosomes))

    # Paralogs
    if config["paralogs"]:
        out.extend(expand("results/paralogs/bed/{prefix}.{chr_id}.paralogs.bed", prefix=final_prefix, chr_id=chromosomes))
    
    return out


def get_previews(preview_path:str):
    chromosomes = get_chr_list(config["fai_path"])
    return expand(preview_path, prefix=config["final_prefix"], chr_id=chromosomes)


def get_chunk_files(wc):
    return expand(
        "results/paralogs/lr/{{prefix}}.{{chr_id}}.{chunk}.pval.tsv",
        chunk=["%.2d" % i for i in range(int(config["chunks_count"]))]
    )