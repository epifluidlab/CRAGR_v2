# <img alt="yin-yang-two-colors" src="man/figures/CRAGR_v2.png" height="60"> ‎ ‎ ‎CRAGR_v2
**CRAGR**_v2 (**C**ell f**R**ee dn**A** fra**G**mentation in **R**) is a pipeline for [CRAGR](https://github.com/epifluidlab/cragr) that incorporates Irreproducible Discovery Rate (IDR) to hotspot calling enable consistency between replicates or sets of samples.

## Table of Contents

- [1. Installation](#installation)
- [2. Quick Start](#quick-start)
- [3. Parameters](#parameters)
- [4. Workflow Diagram](#workflow-diagram)
- [5. Citation](#citation)
- [6. Contact](#contact)
- [7. License](#license)

## Installation

1. (`bash`) Create and activate a `conda` environment with a Python 3.9.

```bash
conda create --name idr_cragr python=3.9
conda activate idr_cragr
```

2. (`bash`) Ensure that you have the proper dependencies (`snakemake`, `tabix`, `bgzip`, `samtools`, `bedtools`, `idr`) installed.
```bash
# For Tabix, BGZIP, and SAMTools:
conda install bioconda::samtools

# For BEDTools:
conda install bioconda::bedtools

#For IDR:
conda install bioconda::idr 
# Open ~.conda/envs/idr_cragr/lib/python3.9/site-packages/idr/idr.py
# Replace numpy.int in lines 299 and 300 with numpy.int_

# For Snakemake:
pip install snakemake==7.32.4
```

3. (`R`) Install CRAGR.

```bash
install.packages("devtools")
devtools::install_github("epifluidlab/CRAGR_v2")
```

## Quick Start

1. (`R`) Identify the path to the CRAGR_v2 pipeline.
```R
system.file("extdata/scripts/idr_pipeline", "idr_cragr.smk", package = "cragr")
```
2. Specify the relevant arguments for analysis in a `params.yaml` file. Check [here](#parameters) for more information.

3. (`bash`) Run the pipeline.
```bash
snakemake -s {path_to_idr_crag_smk} --configfile params.yaml --cores 16
```

## Parameters
For a detailed understanding of the CRAG methods and analysis stages, we encourage you to read through our documentation here: [CRAGR](https://github.com/epifluidlab/cragr).

Our pipeline takes the following arguments in YAML format. An example `params.yaml` is linked [here](inst/extdata/scripts/idr_pipeline/params.yaml), for your guidance.

### Required Parameters
- `r_path`: (REQUIRED) This is a path to your Rscript executable. You can usually find this with `which Rscript`.
- `cragr_script`: (REQUIRED) This is a path to the `CRAGR.r` script.
  - You can find this path by running `system.file("extdata/scripts", "cragr.R", package = "cragr")` OR
  - Locate the source code directory and look at `int/ext/scripts/`.
- `samples`: (REQUIRED) This is a path to a file containing the paths to the relevant sample fragment files.
  - These files must start with an ID that is separated by a `.`. For example `{ID}.frag.bed.gz`
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
- `fdr_threshold`: (DEFAULT: 0.2) This determines the FDR cutoff for filtering peaks after Stage 2 analysis and before IDR.
- `merge_gap`: (DEFAULT: 200) If the distance between two hotspot intervals is less than this value, they will be merged into one larger hotspot.
- `idr_threshold`: (DEFAULT: 0.05) This determines the cutoff for filtering peaks from the IDR value for each peak.
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
- Yaping Liu: yaping@northwestern.edu

## License
This project falls under an MIT license. See the included `LICENSE` file for details.