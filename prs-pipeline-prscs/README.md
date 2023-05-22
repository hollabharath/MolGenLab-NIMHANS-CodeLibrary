# PRS Pipeline using PLINK and PRScs.py

This repository contains a pipeline for performing Polygenic Risk Score (PRS) analysis using PLINK and PRScs.py. PRS-CS is a Python-based command-line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors using GWAS summary statistics and an external LD reference panel.

## Requirements

Before running the pipeline, make sure you have the following prerequisites:

- [PLINK](https://www.cog-genomics.org/plink/) installed on your system.
- [PRScs.py](https://github.com/getian107/PRScs) installed on your system.
- External LD reference panel (can be downloaded from the [author's repository](https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference)). Extract the contents, for example: 

` tar -zxvf ldblk_ukbb_eur.tar.gz`
- Download and preprocess the summary statistics file. Example command:

` zcat pgc-bip2021-all.vcf.tsv.gz | awk 'NR==1 { print "SNP\tA1\tA2\tBETA\tP" } NR>1 && !/^##/ && !/^#/ { print $3, $4, $5, $6, $7 }' | column -t > pgc-bip2021-all.txt`
- GNU parallel installed on your system. If it's not already installed, you can typically install it via your package manager. For example, on Ubuntu, you can use the command: 
`sudo apt-get install parallel`
## Usage

```bash
bash PRScs_parallel_score.sh [OPTIONS]
```

**Options:**

- `--ref_dir PATH_TO_REFERENCE`: Path to the reference directory.
- `--bim_prefix VALIDATION_BIM_PREFIX`: Prefix of the validation BIM file.
- `--sst_file SUM_STATS_FILE`: Path to the summary statistics file.
- `--n_gwas GWAS_SAMPLE_SIZE`: Number of samples in the GWAS.
- `--out_dir OUTPUT_DIR`: Output directory with Prefix for Posterior Effects.
- `--a VALUE`: Value for PARAM_A (optional, default: 1).
- `--b VALUE`: Value for PARAM_B (optional, default: 0.5).
- `--phi VALUE`: Value for PARAM_PHI (optional, default: auto).
- `--n_iter VALUE`: Value for MCMC_ITERATIONS (optional, default: 1000).
- `--n_burnin VALUE`: Value for MCMC_BURNIN (optional, default: 500).
- `--thin VALUE`: Value for MCMC_THINNING_FACTOR (optional, default: 5).
- `--seed VALUE`: Value for SEED (optional, default: 1234).

## Example

```bash
bash PRScs_parallel_score.sh --ref_dir path_to_ref/ldblk_1kg_eur --bim_prefix validation --sst_file path_to_sumstats/sumstats.txt --n_gwas 1000 --out_dir path_to_output/prefix 
```

## Outputs
- Parallel execution of the PRScs.py command for multiple chromosomes, optimizing the efficiency of the analysis.
- Combines chromosome-specific posterior SNP effect sizes into a single file.
- [Check and Removes duplicates from the effect and genotype files.]
- Computes individual-level polygenic scores using PLINK's `--score` function.

