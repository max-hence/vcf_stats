# Correct vcf to mark missing sample genotypes in the bed file as missing data in the vcf

# from 
#   1	29571	.	A	T	.	PASS	.	GT	.|1	0|0	0|0	0|0	0|0

# and
# 1   1   30000   3   msp_0,msp_1,msp_2

# correct into
#   1	29571	.	A	T	.	PASS	.	GT	.|1	0|0	0|0	.|.	.|.

from pandas import read_csv
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def change_SNP_genotype(row:str, index_NA:str):

    """for a given <row> apply NA (.|.) to every given <index_NA>"""

    return "\t".join(row.split("\t")[0:9]) + "\t" + "\t".join([".|." if i in index_NA else genotype for i, genotype in enumerate(row.strip().split('\t')[9:])]) + "\n"


def correct_genotype(vcf_path:str, bed_path:str, correct_vcf_path:str):
    """Main function to assign NA to each badly called genotypes based on callable_all_sample.bed
    Works only for one chromosome !

    Args:
        vcf_path (str):
        bed_path (str):
        correct_vcf_path (str):
    """

    # Loading of big bed table (callable_all_sample.bed)
    bed_table = read_csv(bed_path, delimiter="\t", header=None, names=["chrom", "chromStart", "chromEnd", "n_samples", "samples"], index_col=False)
    total_intervals = bed_table.shape[0]
    print(f"BED TABLE LOADED : {total_intervals} rows", flush=True)

    interval_idx = 0
    new_idx = True
    bed_exceeded = False # True if SNPs is outside bed file (I don't get why it happens sometimes...)

    with open(correct_vcf_path, "w") as new_vcf:
        with open(vcf_path, "r") as vcf:
            for row in vcf:
                if row.startswith("##"): new_vcf.write(row); continue # skip comments
                if row.startswith("#CHROM"): # get sample ids
                    all_samples = row.strip().split("\t")[9:]
                    new_vcf.write(row)
                    continue

                pos = row.strip().split('\t')[1] #locus
                if not bed_exceeded:
                    # find the good interval on the bed file associated with snp position
                    # I do that instead of searching for the intervall for each SNP to gain time
                    while int(pos) >= int(bed_table.iloc[interval_idx]["chromEnd"]):
                        interval_idx += 1
                        new_idx = True
                        if interval_idx == total_intervals: 
                            bed_exceeded = True
                            break

                # Just to see the progress
                if interval_idx%10000 == 0 and new_idx and interval_idx != 0:
                    print(f"{interval_idx} rows", flush=True)
                    new_idx = False
                
                if not bed_exceeded:
                # Mark all genotype as NA if SNP outside bed file
                    missing_samples = find_missing_samples(bed_table.iloc[interval_idx], all_samples) # get not called samples

                    if len(missing_samples) == 0: new_vcf.write(row); continue # all samples called
                    else:
                        index_NA = [all_samples.index(sample) for sample in missing_samples] # get indices to modify the vcf
                else:
                    print(f"Warning : SNP at position {pos} is outside bed file limit. All genotypes will be markes as NA.",
                        flush=True)
                    index_NA = [i for i in range(len(all_samples))]

                new_row = change_SNP_genotype(row, index_NA)
                new_vcf.write(new_row)


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
    parser.add_argument('-i', '--input',type=str,
        help="Specifies the path to vcf file"
    )
    parser.add_argument('-b', '--bed', type=str,
        help="Specifies path to callability bed file"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Specifies the path to corrected vcf file"
    )

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = parse_command_line()
    correct_genotype(args.input, args.bed, args.output)
    