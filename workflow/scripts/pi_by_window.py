from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import sys
import os

sys.path.append(os.path.abspath("/projects/plantlp/utils/scripts"))
from utils import get_files

def pi_by_wdw(infile:str, window_size:int, outfile:str):
    """Measures pi by windows and includes NAs when averaging

    Args:
        infile (str): pi by site file (from vcftools)
        window_size (int): int (default = 100000)
        outfile (str): pi by windows
    """
    sites_pi = read_csv(infile, sep="\t", index_col=None)
    sites_pi['window'] = (sites_pi['POS'] - 1) // window_size

    window_pi = (
    sites_pi.groupby(['CHROM', 'window'], as_index=False)
        .agg(
            sum_PI=('PI', 'sum'),         # sum of PI values
            count_valid=('PI', 'count'),  # number of non-NaN positions
            count_total=('POS', 'count')  # Total number of positions
        )
    )

    window_pi['start'], window_pi['end'] = window_pi['window'] * window_size + 1, (window_pi['window'] + 1) * window_size   
    window_pi["mean"] = round(window_pi["sum_PI"] / (window_size - window_pi["count_total"] + window_pi["count_valid"]), 6)

    window_pi.drop(columns=['window','sum_PI', 'count_valid', 'count_total']).to_csv(outfile, index=None, sep="\t")


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', type=str, required=True,
        help="Site-pi file"
    )
    parser.add_argument('-w', '--window', type=int, required=False,
        help="Window size"
    )
    parser.add_argument('-o', '--output', type=str, required=True,
        help="Window-pi file"
    )
    parser.add_argument('--centro', type=float, required=False,
        help="Relative centromeric position"
    )
    parser.add_argument('--plot', action="store_true",
        help="Plot pi by window"
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_command_line()
    if args.plot:
        plot_pi(args.input, args.output, args.centro)
    else:
        if args.window is None: args.window = 100000
        pi_by_wdw(args.input, args.window, args.output)