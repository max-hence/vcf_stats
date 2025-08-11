#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pandas import read_csv, DataFrame
from argparse import ArgumentParser, RawDescriptionHelpFormatter


def get_fis(vcf:str, samples:str, min_samples:int):
	"""
	Measure per site fis from A. Mackintosh script. 
	With correction of Weir & Cockerham (1984) for small sampling
	"""

	sample_idxs = [int(i) for i in samples.split(",")]
	overall_fis = []
	i = 0
	with open(vcf, "r") as file:

		for line in file: # for each SNPS

			line = line.strip('\n')

			if line.startswith("#"):
				pass
			else:
	
				if i % 100000 == 0:
					print(f"SNP n°{i}")

				line_list = line.split("\t")

				genotypes = [
					gt for i, field in enumerate(line_list[9:], start=1)
					if (gt := field.split(":")[0].replace('|', '/')) not in {"./.", "."}
					and (i in sample_idxs)
				]

				alleles = [allele for gt in genotypes for allele in gt.split("/")]
				n = len(genotypes)
				if n >= min_samples and len(set(alleles)) == 2: # check if enough called samples

					# make sur genotypes are in the format ./.
					alleles.sort()
					het_genotype_1 = alleles[0] + "/" + alleles[-1]
					het_genotype_2 = alleles[-1] + "/" + alleles[0]

					# observed Heterozygosity
					H_o = (genotypes.count(het_genotype_1) + genotypes.count(het_genotype_2)) / n

					p = alleles.count(alleles[0]) / (2 * n) # freq ancestral allele

					# He with Weir & Cockerham (1984) correction for small sampling
					b = n / (n - 1) * (p * (1 - p) - (2 * n - 1) / (4 * n) * H_o)
					c = H_o / 2

					fis = 1 - (c / (b + c))

					overall_fis.append({"POS": int(line_list[1]), "fis": fis})

				i += 1

		return DataFrame(overall_fis)


def fis_by_window(fis_df, window_size, outfile):
	""" Measure the mean by window """

	fis_df['window'] = (fis_df['POS'] - 1) // window_size

	window_fis = (
		fis_df.groupby(['window'], as_index=False)
			.agg(
				sum_fis=('fis', 'sum'), # sum of PI values
				count=('fis', 'count'), # number of positions
			)
		)

	window_fis['start'], window_fis['end'] = window_fis['window'] * window_size + 1, (window_fis['window'] + 1) * window_size   
	window_fis["mean_fis"] = round(window_fis["sum_fis"] / window_fis["count"], 8)

	window_fis.drop(columns=['window','sum_fis', 'count']).to_csv(outfile, index=None, sep="\t")


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help="Specifies the path to vcf file"
    )
    parser.add_argument('-s', '--samples', type=str,
        help="Comma delimited list of samples indices (1-based)"
    )
    parser.add_argument('-m', '--min_samples', type=int,
        help="Minimum number of genotypes for a SNP to contribute to the Fis estimate"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Path to fis by window"
    )

    args = parser.parse_args()

    return args

if __name__== "__main__":

	args = parse_command_line()
	fis_df = get_fis(args.input, args.samples, args.min_samples)
	fis_df.to_csv("hmmm.tsv", sep="\t", index=None)

	WINDOW = 100000
	fis_by_window(fis_df, WINDOW, args.output)