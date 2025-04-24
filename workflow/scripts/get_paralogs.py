# import libraries
from pandas import read_csv
from scipy.stats import chi2
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from argparse import ArgumentParser, RawDescriptionHelpFormatter


def chi2_test(infile:str, outfile:str):
    """Chi square test to check which site is reads covering a site derive form a single locus
    Prints out sites that as a corrected p-value below 0.05

    Args:
        infile (str): _description_
        outfile (str): csv file with p-value and sites labelled as paralogs
    """

    df_ngs = read_csv(infile, index_col=None, header=None,names=["chr", "pos", "loglikelihood", "alt_loglikelihood", "LR_mismapped"], sep="\t")
    df_ngs["chr"] = df_ngs["chr"].apply(strip_str)

    df_ngs["pval"] = 0.5*chi2.sf(df_ngs["LR_mismapped"], df=1) # sf means survival function = 1 - Cumulative density function = P(X > x)
    df_ngs["pval_adj"] = multipletests(df_ngs["pval"], method="bonferroni")[1]
    
    df_ngs["paralog"] = df_ngs["pval_adj"] < 0.05
    df_ngs.to_csv(outfile, index=None, sep="\t")


def manhattan_plot(infile:str, outfile:str):
    """Plot LRs from ngsParalog

    Args:
        infile (str): _description_
        outfile (str): _description_
    """
    df = read_csv(infile, index_col=None, sep="\t")

    # manhattan plot
    plt.figure(figsize=(18, 8))

    colors = df["paralog"].map({True: "#AB252D", False: "#B7B7B7"})
    plt.scatter(df["pos"], df['LR_mismapped'], color=colors, s=10)

    # show the graph
    plt.savefig(outfile)

def strip_str(string):
    return string.strip()

def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help="Input path"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Output path"
    )
    parser.add_argument('--pval',action="store_true",
        help="Estimate probabilities that sites are paralogs"
    )
    parser.add_argument('--plot',action="store_true",
        help="Save Manhattan plot to <output>.png or .svg"
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
 
    args = parse_command_line()

    if args.pval:
        chi2_test(args.input, args.output)

    elif args.plot:
        manhattan_plot(args.input, args.output)

    elif args.save:
        save_paralogs(args.input, args.output)

    else:
        raise ValueError("Choose btw --plot, --save et --pval")