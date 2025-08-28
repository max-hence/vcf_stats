#!/usr/bin/env python3
import sys
import gzip
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv


def all_sites_vcf(vcf_in:str, bed_path:str, threshold:int, vcf_out:str, callability:bool):
    """
    Add invariants sites, remove sites where less than <threshold> samples passed the callability filter
    optional : set each badly called genotype as missing
    """

    open_func = gzip.open if vcf_in.endswith(".gz") else open

    # callable bed file
    bed_df = read_csv(bed_path, delimiter="\t", header=None, names=["chrom", "chromStart", "chromEnd", "n_samples", "samples"], index_col=False)
    
    # check if only one chromosome :
    chrom = bed_df["chrom"].unique()
    if len(chrom) > 1: raise ValueError("This function can only analyze chromosomes one by one")
    else: chrom = str(chrom[0])

    chrom_length = bed_df["chromEnd"][bed_df.shape[0]-1] # total chromosome length
    interval_idx = 0
    bed_row = bed_df.iloc[interval_idx]
    variants = {}

    # Gather SNPs
    with open_func(vcf_in, "rt") as f:
        header = []
        samples = []
        for line in f:
            if line.startswith("##"):
                header.append(line.rstrip())
            elif line.startswith("#CHROM"):
                header.append(line.rstrip())
                samples = line.strip().split("\t")[9:]
            else:
                parts = line.strip().split("\t")
                variants[int(parts[1])] = {"info" :parts[0:9], "genotypes" : parts[9:]}

    # Écriture du VCF all-sites
    with open(vcf_out, "w") as out:
        for h in header: out.write(h + "\n")

        # Format de ligne invariant
        if callability:
            invariants = ["0|0"] * len(samples)
        else:
            invariants = "\t".join(["0|0"] * len(samples))

        for pos in range(1, chrom_length + 1):
            
            if pos < chrom_length:
                # Look in the right callable bed intervall
                if callability:
                    if pos == bed_row["chromEnd"]:
                        interval_idx += 1
                        bed_row = bed_df.iloc[interval_idx]

                    if int(bed_row["n_samples"]) >= threshold: # if callability threshold passed
                        missing_samples = find_missing_samples(bed_row, samples)
                        index_NA = [samples.index(sample) for sample in missing_samples] # get indices to modify the vcf
                        if pos in variants:
                            info = "\t".join(variants[pos]['info'])
                            out.write(f"{info}\t{set_missing_genotypes(variants[pos]['genotypes'], index_NA)}\n")
                        else:
                            info = "\t".join([chrom, str(pos), ".", "N", ".", ".", "PASS", ".", "GT"])
                            out.write(f"{info}\t{set_missing_genotypes(invariants, index_NA)}\n")
                else:
                    if pos in variants:
                        info = "\t".join(variants[pos]['info'])
                        str_variants = "\t".join(variants[pos]['genotypes'])
                        out.write(f"{info}\t{str_variants}\n")
                    else:
                        info = "\t".join([chrom, str(pos), ".", "N", ".", ".", "PASS", ".", "GT"])
                        out.write(f"{info}\t{invariants}\n")


def set_missing_genotypes(genotypes:list, index_NA:str):

    """for a given <row> apply NA (.|.) to every given <index_NA>"""

    return "\t".join([".|." if i in index_NA else genotype for i, genotype in enumerate(genotypes)])


def find_missing_samples(bed_row, all_samples:list):
    """Gives badly called samples for a given interval

    Args:
        bed_row (pandas array): bed entry
        all_samples (list): all sampling set

    Return:
        list of not called samples
    """
    return [sample for sample in all_samples if sample not in bed_row["samples"].split(",")]



def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', type=str, required=True,
        help="vcf file"
    )
    parser.add_argument('-b', '--bed', type=str, required=False,
        help="callable bed file"
    )
    parser.add_argument('--callability', action="store_true",
        help=""
    )
    parser.add_argument('-t', '--threshold', type=int, required=False,
        help="nbr of indivuals called for a site to be kept"
    )
    parser.add_argument('-o', '--output', type=str, required=True,
        help="path to AllSites vcf"
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_command_line()
    if args.callability:
        if args.threshold == None:
            raise ValueError("if --callability is set, you must precise a threshold (-t)")
    else:
        args.threshold = None, None
    all_sites_vcf(args.input, args.bed, args.threshold, args.output, args.callability)
