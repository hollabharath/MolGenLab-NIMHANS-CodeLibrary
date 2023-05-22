#!/bin/bash
N_THREADS=10
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# Function to display script usage
usage() {
  echo "Usage: bash $0 [OPTIONS]"
  echo "Options:"
  echo "  --ref_dir PATH_TO_REFERENCE     Path to the reference directory"
  echo "  --bim_prefix VALIDATION_BIM_PREFIX   Prefix of the validation BIM file"
  echo "  --sst_file SUM_STATS_FILE       Path to the summary statistics file"
  echo "  --n_gwas GWAS_SAMPLE_SIZE       Number of samples in the GWAS"
  echo "  --out_dir OUTPUT_DIR            Output directory with Prefix for Posterior Effects"
  echo "  --a VALUE                       Value for PARAM_A (optional, default: 1)"
  echo "  --b VALUE                       Value for PARAM_B (optional, default: 0.5)"
  echo "  --phi VALUE                     Value for PARAM_PHI (optional, default: auto)"
  echo "  --n_iter VALUE                  Value for MCMC_ITERATIONS (optional, default: 1000)"
  echo "  --n_burnin VALUE                Value for MCMC_BURNIN (optional, default: 500)"
  echo "  --thin VALUE                    Value for MCMC_THINNING_FACTOR (optional, default: 5)"
  echo "  --seed VALUE                    Value for SEED (optional, default: 1234)"
  echo
  echo "Example:"
  echo "bash $0 --ref_dir path_to_ref/ldblk_1kg_eur --bim_prefix validation --sst_file path_to_sumstats/sumstats.txt --n_gwas 1000 --out_dir path_to_output/prefix [--a 2 --b 0.8 --phi 0.3 --n_iter 2000 --n_burnin 1000 --thin 10 --seed 5678]"
}

# Set default values
a=1                 # "PARAM_A"
b=0.5               #"PARAM_B"
phi="auto"
n_iter=1000         #"MCMC_ITERATIONS"
n_burnin=500        #"MCMC_BURNIN"
thin=5              #"MCMC_THINNING_FACTOR"
seed=1234           #"SEED"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref_dir)
      ref_dir="${2}"
      shift 2
      ;;
    --bim_prefix)
      bim_prefix="${2}"
      shift 2
      ;;
    --sst_file)
      sst_file="${2}"
      shift 2
      ;;
    --n_gwas)
      n_gwas="${2}"
      shift 2
      ;;
    --out_dir)
      out_dir="${2}"
      shift 2
      ;;
    --a)
      a="${2}"
      shift 2
      ;;
    --b)
      b="${2}"
      shift 2
      ;;
    --phi)
      phi="${2}"
      shift 2
      ;;
    --n_iter)
      n_iter="${2}"
      shift 2
      ;;
    --n_burnin)
      n_burnin="${2}"
      shift 2
      ;;
    --thin)
      thin="${2}"
      shift 2
      ;;
    --seed)
      seed="${2}"
      shift 2
      ;;
    *)
      echo "Invalid argument: $1"
      usage
      exit 1
      ;;
  esac
done

#


# Check if required options are provided
if [[ -z $ref_dir || -z $bim_prefix || -z $sst_file || -z $n_gwas || -z $out_dir ]]; then
  echo "Error: Missing required options."
  usage
  exit 1
fi


# Specify the chromosome range
chromosomes=$(seq 1 22)

# Run PRScs.py command with specified options for each chromosome in parallel
echo "$chromosomes" | parallel -j0 python PRScs.py \
  --ref_dir "$ref_dir" \
  --bim_prefix "$bim_prefix" \
  --sst_file "$sst_file" \
  --n_gwas "$n_gwas" \
  --out_dir "$out_dir" \
  --a "$a" \
  --b "$b" \
  --n_iter "$n_iter" \
  --n_burnin "$n_burnin" \
  --thin "$thin" \
  --chrom {} \
  --seed "$seed"

# Combine the chromosome-specific summary statistics files into a single file
echo "Combining the chromosome-specific summary statistics files..."
for chr in $(seq 1 22); do
  cat "${out_dir}_pst_eff_a${a}_b${b}_phi${phi}_chr${chr}.txt" >> "${out_dir}_pst_eff_a${a}_b${b}_phi${phi}.txt"
done

# Remove duplicates from the score file
awk '!seen[$2]++' "${out_dir}_pst_eff_a${a}_b${b}_phi${phi}.txt" > "${out_dir}_pst_eff_a${a}_b${b}_phi${phi}_noduplicates.txt"

# Check for duplicates in the PLINK genotype file
echo ""
echo ""
echo "Checking for duplicates in the PLINK genotype file..."
echo ""
echo ""
plink2 \
  --bfile "$bim_prefix" \
  --rm-dup force-first list \
  --out "$bim_prefix"

# Check if duplicates were found
if [ -e "$bim_prefix".rmdup.list ]; then
  echo ""
  echo ""
  echo "Duplicates found in the PLINK genotype file. Removing duplicates..."
  echo ""
  echo ""
  plink2 \
    --bfile "$bim_prefix" \
    --exclude "$bim_prefix".rmdup.list \
    --make-bed \
    --out "$bim_prefix"_noduplicates

  bim_prefix="$bim_prefix"_noduplicates
fi

# Compute individual level polygenic scores using PLINK's --score function
echo "Computing individual level polygenic scores..."
plink \
  --bfile "$bim_prefix" \
  --score "${out_dir}_pst_eff_a${a}_b${b}_phi${phi}_noduplicates.txt" 2 4 6 sum center \
  --out "${out_dir}_a${a}_b${b}_phi${phi}_prs"

