# PigGTEx-GFBLUP

Genomic prediction (GFBLUP) in pigs informed by PigGTEx QTL annotations.

> Note: The pipeline below is demonstrated with lncQTL (dataset data5, traits AGE/BF/TNB).
> To switch to other QTL types (eQTL/sQTL, etc.) or datasets, update the hard-coded PigGTEx
> input directory, RDS ID, and traits at the top of each script accordingly.

## Scripts

| Script | Function |
|---|---|
| `piggtex_d.R` | Matches PigGTEx significant QTLs to chip SNPs (exact + ±50 kb), outputs a per-tissue list of QTL SNPs |
| `qtl_piggtex_gs_50kb.R` | Performs 5-fold × 2-repeat cross-validation using GFBLUP with two random effects (QTL SNPs + background SNPs), evaluating prediction accuracy (cor) and bias |
| `result_summary.R` | Summarizes cross-validation results across trait × tissue combinations, reporting mean ± SD of cor and bias |
| `submit_qtl.sh` | SLURM array job for batch submission (one task per trait × tissue combination) |

## Usage

```bash
# 1. Generate QTL annotations (argument 1-5 selects the chip, see script comments)
Rscript piggtex_d.R 5

# 2. Submit cross-validation (first confirm the enabled Rscript line in submit_qtl.sh)
sbatch submit_qtl.sh

# 3. Summarize results
Rscript result_summary.R
```

Run the scripts in order.
