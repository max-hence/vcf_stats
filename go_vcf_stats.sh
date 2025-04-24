#! /bin/bash
#SBATCH --job-name=VCF_STATS
#SBATCH -o /projects/plantlp/02_VCF_PROCESSING/genouest_log/vcf_stats/%x.%j.out
#SBATCH -e /projects/plantlp/02_VCF_PROCESSING/genouest_log/vcf_stats/%x.%j.err
#SBATCH -p ecobio,genouest
#SBATCH --time=1-00
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

. /projects/plantlp/paths.sh
echo "Command :
$snp_scripts/go_vcf_stats.sh $@" >&2

echo
echo $(date) >&2
echo

usage() {
   echo "Usage: $0 [-i input] [-o output] [-p profile]" >&2
   exit 1
}

unlock=false
rerun=false
dryrun=false

while getopts "i:o:p:urn" opt; do
  case $opt in
    i)
      config="$OPTARG" # path/to/config.yaml
      ;;
    o)
      dir="$OPTARG" # path/to/analysis
      ;;
    p)
      profile="$OPTARG" # path/to/profile.yaml
      ;;
    u)
      unlock=true
      ;;
    r)
      rerun=true
      ;;
    n)
      dryrun=true
      ;;
    \?)
      usage
      ;;
  esac
done

. /local/env/envconda.sh
conda activate snparcher

if [ $unlock == true ]; then
  snakemake -s $vcf_scripts/vcf_stats/workflow/Snakefile \
    -d $dir \
    --unlock
  exit 1
fi

mkdir -p $dir/config
cp $config $dir/config/config.yml
cp $profile $vcf_scripts/vcf_stats/profiles/slurm/config.yaml

if [ $dryrun == true ]; then
  snakemake -s $vcf_scripts/vcf_stats/workflow/Snakefile \
    -d $dir \
    --dry-run
  exit 1
fi

if [ $rerun == true ]; then
  snakemake -s $vcf_scripts/vcf_stats/workflow/Snakefile \
  -d $dir \
  --workflow-profile $vcf_scripts/vcf_stats/profiles/slurm \
  --use-conda \
  --conda-frontend conda \
  --rerun-incomplete
else
  snakemake -s $vcf_scripts/vcf_stats/workflow/Snakefile \
    -d $dir \
    --workflow-profile $vcf_scripts/vcf_stats/profiles/slurm \
    --use-conda \
    --conda-frontend conda
fi