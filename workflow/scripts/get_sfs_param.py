from argparse import ArgumentParser, RawDescriptionHelpFormatter
from math import log, lgamma

def get_params(previews:list):
    """Parse inside preview.txt

    Args:
        previews (list): all preview files, by chr
    Returns:
        List of dictionnary
    """
    list_params = []
    for preview in previews:
        with open(preview, 'r') as file:
            for row in file:
                if row[0] == "(":
                    list_params.append({int(param[1:-1].split(", ")[0]) : int(param[1:-1].split(", ")[1]) for param in row.strip().split("\t")})

    return list_params


def get_sfs_param(previews: str, method: str, output_path:str):
    """ Get the highest sampling nbr that maximizes the nbr of segregating sites

    Args:
        previews (str): list of preview.txt for each chr
        output_path (str): 
    Returns:
        best_sample (str): best sampling size and nbr of segregating sites
    """

    all_params = get_params(previews)
    total_params = {n : 0 for n in all_params[0].keys()}
    for n in total_params: # sum all projections by chrom
        for chr_id in range(0, len(all_params)):
            total_params[n] += all_params[chr_id][n]

    if method == "max":
        best_sample = max(total_params, key=total_params.get)
    else: # method == "ml":
        best_sample = ml_sfs(total_params)
    
    with open(output_path, "w") as file:
        file.write("Sampling\tSNPS\n")
        file.write(f"{best_sample}\t{total_params[best_sample]}")


def ml_sfs(total_params:dict):
    """Maximul Likelihood function to measure what parameter sets are the best
    based on a Poisson distribution and a SFS at equilibrium

    Args:
        total_params (dict): _description_

    """
    # log(a x b) = log(a) + log(b)
    # log(1) = 0
    # donc log(n!) = sum([log(i) for i in [1:n])
    max_lnL = 0
    for n_ind, n_snp in total_params.items():
        if n_snp == 0: break # never happens except for a very small toy_example vcf
        sfs = neutral_sfs(n_ind, n_snp)
        lnL = lgamma(n_ind+1) + sum([cat*log(cat/n_ind)-lgamma(cat+1) for cat in sfs])

        if -lnL > max_lnL:
            max_lnL = -lnL
            best_param = n_ind

    return best_param


def neutral_sfs(n_ind:int, n_snp:int):
    thetaW = n_snp/sum([1/i for i in range(1, n_ind)])

    # unfolded neutral sfs
    return [thetaW/cat for cat in range(1, n_ind)]


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', nargs="+", required=True,
        help="Specifies the list of preview files"
    )
    parser.add_argument('-m', '--method', type=str, required=True,
        help=" max or ml to specify with technique to select"
    )
    parser.add_argument('-o', '--output', type=str, required=True,
        help="Specifies the output path"
    )

    args = parser.parse_args()

    return args


if __name__=="__main__":

    args = parse_command_line()

    if args.method not in ["ml", "max"]:
        raise ValueError("Arg method is either \"max\" or \"ml\"")

    get_sfs_param(args.input, args.method, args.output)