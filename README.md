# <img alt="yin-yang-two-colors" src="man/figures/CRAGR_v2.png" height="60"> ‎ ‎ ‎CRAGR_v2
**CRAGR**_v2 (**C**ell f**R**ee dn**A** fra**G**mentation in **R**) is a pipeline for [CRAGR](https://github.com/epifluidlab/cragr) that incorporates Irreproducible Discovery Rate (IDR) to hotspot calling enable consistency between replicates or sets of samples.

## Table of Contents

- [1. Installation](#installation)
- [2. Quick Start](#quick-start)
- [3. Example Usage](#example-usage)
- [4. Parameters](#parameters)
- [5. Workflow Diagram](#workflow-diagram)
- [6. Citation](#citation)
- [7. Contact](#contact)
- [8. License](#license)

## Installation

1. (`bash`) Create and activate a `conda` environment with a Python 3.9.

```bash
conda create --name idr_cragr python=3.9
conda activate idr_cragr
```

2. (`bash`) Ensure that you have the proper dependencies (`snakemake`, `tabix`, `bgzip`, `samtools`, `bedtools`, `idr`, `libiconv`) installed.
```bash
# For Tabix, BGZIP, and SAMTools:
conda install bioconda::samtools

# For BEDTools:
conda install bioconda::bedtools

#For IDR:
conda install bioconda::idr 
# Open ~/.conda/envs/idr_cragr/lib/python3.9/site-packages/idr/idr.py
# Replace numpy.int in lines 299 and 300 with numpy.int_

# For Snakemake:
pip install snakemake==7.32.4
pip install pulp==2.7.0

# NOTE: For some HPC systems, libiconv.so may not be detected for the CRAGR installation process
# You may have to manually add libiconv to the path in this case, like so:
conda install conda-forge::libiconv
export $LD_LIBRARY_PATH=<PATH TO CONDA LIBRARY FOLDER>
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
snakemake -s {path_to_idr_crag_smk} --configfile params.yaml --cores {cores}
```

## Example Usage

For your ease of use, we have prepared an lightweight example usage [here](`inst/extdata/scripts/data/example`). You may run the pipeline:

```bash
# Enter the example usage directory.
cd ~/inst/extdata/scripts/data/example

# Modify the r_script parameter in example_params.yaml.

# Run the pipeline.
snakemake -s {path_to_idr_crag_smk} --configfile example_params.yaml --cores {cores}
```
> [!NOTE]  
> In order to remain lightweight, the sample fragment files are relatively low coverage. As a result, the pipeline will identify no peaks and fail at the IDR step (47%).



## Parameters
For a detailed understanding of the CRAG methods and analysis stages, we encourage you to read through our documentation here: [CRAGR](https://github.com/epifluidlab/cragr).

Our pipeline takes arguments in YAML format. A blank `params.yaml` is linked [here](inst/extdata/scripts/idr_pipeline/params.yaml), for your guidance.

### Required Parameters
- `r_path`: (REQUIRED) This is a path to your Rscript executable. You can usually find this with `which Rscript`.
- `cragr_script`: (REQUIRED) This is a full path to the `CRAGR.r` script.
  - You can find this by running `system.file("extdata/scripts", "cragr.R", package = "cragr")` and pasting the full output OR
  - Locate the source code directory and place the path to `inst/extdata/scripts/cragr.R`.
- `samples`: (REQUIRED) This is a path to a  `.txt` file containing a list of paths (new-line separated) to the relevant sample fragment files.
  - These files must have a consistent naming structure that starts with an ID separated by a `.`. For example `{ID}.frag.bed.gz`
  - All of these fragment files must have an associated `tabix` index.
- `chroms`: (REQUIRED) This is a path to a `.txt` file containing a list of chromosomes (new-line separated) that hotspots and IFS scores should be retrieved for.
   - This file must be ordered, as this order determines the order of the sorting of the combined replicate fragment files.
- `genome`: (REQUIRED, OPTIONS=['hg19', 'hg38']) This is the name of the genome build to use for the CRAG pipeline. 
- `chrom_sizes`: (REQUIRED) This is the path to a file containing the chromosome sizes for the relevant `--genome`.
  - See [`inst/extdata/scripts/data/`](inst/extdata/scripts/idr_pipeline/data) for .chrom.sizes files for `hg19` and `hg38`.
- `output_dir`: This is the path to a directory that all intermediate and output files will be stored in. 
   - If this directory does not exist, it will be created.

### Optional Parameters
- `split_method`: (DEFAULT: fragment_count, OPTIONS:['fragment_count', 'by_sample']) This is determines how the fragments are split into the replicates.
  - `fragment_count`: This will combine all of the fragment files and take `--subsample` fragments in each replicate with replacement. This inherently introduces overlaps to the replicates. To avoid very similar replicates, the total fragment count of all fragment files must be greater than `1.5*--subsample`.
  - `samples`: This will generate an equal split of the files into replicates. This ensures that the replicates have no data overlap.
- `seed`: (DEFAULT: None) This is an integer that sets the seed of the random splitting of the replicates for reproducibility.
- `subsample`: (DEFAULT: 200M) This argument sets the number of fragments to subsample from the combined fragments to create the two replicate fragment files.
- `exclude_regions`: (DEFAULT: None)  This is the path to a blacklist BED file, highlighting problematic areas of the genome to ignore.
  - You can learn more about the blacklist file purpose [here](https://www.nature.com/articles/s41598-019-45839-z).
  - See [`inst/extdata/scripts/data/`](inst/extdata/scripts/idr_pipeline/data) for .exclude.regions files for `hg19` and `hg38`.
- `high_mappability`: (DEFAULT: None)  This is a path to a file containing mappability scores in BED format. Restricts the analysis of the genome to high mappability regions.
- See [`inst/extdata/scripts/data/`](inst/extdata/scripts/idr_pipeline/data) for .high.mappability files for `hg19` and `hg38`.
- `gc_correct`: (DEFAULT: True)  This is a boolean parameter that determines whether or not to perform GC correction.
- `gc_correct_method`: (DEFAULT: standard, OPTIONS=['standard', 'caret']) This determines the method used in the GC correction.
- `gc_correct_N`: (DEFAULT: 1000000) This determines the maximal sample size for GC correction.
- `fdr_threshold`: (DEFAULT: 0.2) This determines the FDR cutoff for filtering peaks after Stage 2 analysis and before IDR.
- `merge_gap`: (DEFAULT: 200) If the distance between two hotspot intervals is less than this value, they will be merged into one larger hotspot.
- `idr_threshold`: (DEFAULT: 0.05) This determines the cutoff for filtering peaks from the IDR value for each peak.
- `min_mapq`: (DEFAULT: 30) This is the minimum MAPQ score of a fragment to include in the analysis.
- `min_fraglen`: (DEFAULT: 0) This is the minimum fragment length to include in the analysis.
- `max_fraglen`: (DEFAULT: 1000) This is the maximum fragment length to include in the analysis.
- `window_size`: (DEFAULT: 200) This is the sliding window size to use for the CRAG IFS calculation.
- `step_size`: (DEFAULT: 20) The step size to use for the CRAG IFS calculation.
- `total_fragment_min`  (DEFAULT: `1.5*--subsample`) This determines the minimum total fragments across your fragment files to allow the pipeline to run without throwing an error. Modify this threshold at your own risk.

## Workflow Diagram

<p align="center">
  <img src="man/figures/CRAGR_v2_workflow.png" alt="CRAGR_v2 Workflow Diagram" height="800">
</p>

## Citation
Zhou X, Zheng H, Fu H, McKillip KL, Pinney SM, Liu Y. (2022) CRAG: De novo characterization of cell-free DNA fragmentation hotspots in plasma whole-genome sequencing. Genome Medicine [![Static Badge](https://img.shields.io/badge/DOI-10.1101/2020.07.16.201350-red?style=flat-square)](https://doi.org/10.1101/2020.07.16.201350)

## Contact

- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Kundan Baliga: kundanbal2969@k12.ipsd.org
- Yaping Liu: yaping@northwestern.edu

## License
This project falls under an MIT license. See the included `LICENSE` file for details.
