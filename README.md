# <img alt="yin-yang-two-colors" src="man/figures/CRAGR_v2.png" height="60"> ‎ ‎ ‎CRAGR_v2
CRAGR_v2 (**C**ell f**R**ee dn**A** fra**G**mentation in **R**) is a pipeline for [CRAGR](https://github.com/epifluidlab/cragr) that incorporates Irreproducible Discovery Rate (IDR) to hotspot calling enable consistency between replicates or sets of samples.

## Table of Contents

- [1. Installation](#installation)
- [2. Quick Start](#quick-start)
- [3. Parameters](#parameters)
- [4. Workflow Diagram](#workflow-diagram)
- [5. Citation](#citation)
- [6. Contact](#contact)
- [7. License](#license)

## Installation

1. Create and activate a `conda` environment with a Python version greater than 3.8 and an R version greater than 4.1.

```bash
conda create --name CRAGR python=3.11 r-base
conda activate CRAGR
```

2. Ensure that you have the proper dependencies (`snakemake`, `numpy`, `scipy`, `tabix`, `bgzip`, `samtools`, `bedtools`, `idr`) installed.
```bash
# For Snakemake:
pip install snakemake numpy scipy

# For Tabix, BGZIP, and SAMTools:
conda install bioconda::samtools

# For BEDTools:
conda install bioconda::bedtools

#For IDR:
Since IDR has such poor support, directions for this are TBD.
```

## Quick Start

1. Identify the path to the CRAGR_v2 pipeline by running the following in R.
```R
system.file("extdata/scripts/idr_pipeline", "idr_cragr.smk", package = "cragr")
```
2. Specify the relevant arguments for analysis in a `params.yaml` file. Check [here](#parameters) for more information.

3. Run the snakemake pipeline by running the following in bash.
```bash
snakemake -s {path_to_idr_crag_smk} --configfile params.yaml --cores 16
```

## Parameters
For a detailed understanding of the CRAG methods and analysis stages, we encourage you to read through our documentation here: [CRAGR](https://github.com/epifluidlab/cragr).

Our pipeline takes the following arguments in YAML format. An example `params.yaml` is linked [here](inst/extdata/scripts/idr_pipeline/params.yaml), for your guidance.

### Required Parameters
- `r_script`: (REQUIRED) This is a path to the `CRAGR.r` script.
  - You can find this path by running `system.file("extdata/scripts", "cragr.R", package = "cragr")` OR
  - Locate the source code directory and look at `int/ext/scripts/`.
- `samples`: (REQUIRED) This is a path to a file containing the paths to the relevant sample fragment files.
  - All of these fragment files must have an associated `tabix` index.
- `chroms`: (REQUIRED) This is a path to a file containing a list of chromosomes that hotspots and IFS scores should be retrieved for.
   - This file must be ordered, as this order determines the order of the sorting of the combined replicate fragment files.
- `genome`: (REQUIRED, OPTIONS=['hg19', 'hg38']) This is the name of the genome build to use for the CRAG pipeline. 
- `chrom_sizes`: (REQUIRED) This is the path to a file containing the chromosome sizes for the relevant `--genome`.
- `output_dir`: This is the path to a directory that all intermediate and output files will be stored in. 
   - If this directory does not exist, it will be created.

### Optional Parameters
- `seed`: (DEFAULT: None) This is an integer that sets the seed of the random splitting of the files for reproducibility.
- `exclude_regions`: (DEFAULT: None)  This is the path to a blacklist BED file, highlighting problematic areas of the genome to ignore.
  - You can learn more about the blacklist file purpose [here](https://www.nature.com/articles/s41598-019-45839-z).
- `high_mappability`: (DEFAULT: None)  This is a path to a file containing mappability scores in BED format. Restricts the analysis of the genome to high mappability regions.
- `gc_correct`: (DEFAULT: True)  This is a boolean parameter that determines whether or not to perform GC correction.
- `gc_correct_method`: (DEFAULT: standard, OPTIONS=['standard', 'caret']) This determines the method used in the GC correction.
- `gc_correct_N`: (DEFAULT: 1000000) This determines the maximal sample size for GC correction.
- `idr_threshold`: (DEFAULT: 1) This determines the FDR threshold for calling hotspots from the IDR value for each peak.
- `min_mapq`: (DEFAULT: 0) This is the minimum MAPQ score of a fragment to include in the analysis.
- `min_fraglen`: (DEFAULT: 0) This is the minimum fragment length to include in the analysis.
- `max_fraglen`: (DEFAULT: 1000) This is the maximum fragment length to include in the analysis.
- `window_size`: (DEFAULT: 200) This is the sliding window size to use for the CRAG IFS calculation.
- `step_size`: (DEFAULT: 20) The step size to use for the CRAG IFS calculation.

## Workflow Diagram

<p align="center">
  <img src="man/figures/CRAGR_v2_workflow.png" alt="CRAGR_v2 Workflow Diagram" height="800">
</p>

## Citation
Zhou X, Zheng H, Fu H, McKillip KL, Pinney SM, Liu Y. (2022) CRAG: De novo characterization of cell-free DNA fragmentation hotspots in plasma whole-genome sequencing. Genome Medicine [![Static Badge](https://img.shields.io/badge/DOI-10.1101/2020.07.16.201350-red?style=flat-square)](https://doi.org/10.1101/2020.07.16.201350)

## Contact

- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Kundan Baliga: kudosbeluga@gmail.com
- Yaping Liu: yaping@northwestern.edu

## License
This project falls under an MIT license. See the included `LICENSE` file for details.