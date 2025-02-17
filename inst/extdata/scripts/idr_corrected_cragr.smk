# idr_corrected_cragr.smk
#
# This pipeline is used to identify CRAG IFS hotspots that pass the Irreproducible Discovery Rate (IDR) threshold.
#
# The arguments to this script are the following:
#  - `--r_script` (required): The path to the crag.R script.
#  - `--samples` (required): A file containing the paths to the relevant fragment files. All fragment files must have a tabix index.`
#  - `--chroms` (required): A file containing the chromosome names to retrieve CRAG IFS hotspots for.
#  - `--genome` (required): The genome build to use.
#  - `--chrom_sizes` (required): A file containing the chromosome sizes.
#  - `--output_dir` (required): The path to the output directory.
#  - `--exclude_regions` (optional): A file containing the regions to exclude from the analysis. Default is None.
#  - `--high_mappability` (optional): A file containing mappability scores, in BED format. Default is None.`
#  - `--gc_correct` (optional): A flag to indicate whether to perform GC correction. Default is False.
#  - `--gc_correct_method` (optional): Methods used in GC correction. Should be either standard or caret. Default is "standard".
#  - `--gc_correct_N` (optional): Maximal sample size for GC correction model training. Default is 1000000.
#  - `--idr-threshold` (optional): The IDR threshold to use for filtering the CRAG IFS hotspots. Default is 1.
#  - `--min_mapq` (optional): The minimum mapping quality to use for filtering. Default is 0.
#  - `--min_fraglen` (optional): The minimum fragment size to use for filtering. Default is 0.
#  - `--max_fraglen` (optional): The maximum fragment size to use for filtering. Default is 1000.
#  - `--window_size` (optional): The sliding window size to use for the CRAG IFS hotspots. Default is 200.
#  - `--step_size` (optional): The step size to use for the CRAG IFS hotspots. Default is 20.
#
# The rough workflow of this pipeline is as follows:
# 1. Using the samples_file, randomly divide the samples into two groups, rep1 and rep2.
# 2. Combine the fragment files for each group, resort the fragments by chromosome (ordered by the order in --chroms), start, then stop. Tabix index the combined fragment files.
# 3. Using the output rep1.frag.gz and rep2.frag.gz files, run the ifs and peak functions from crag.R to identify CRAG IFS hotspots.
# 4. Using the rep1 and rep2 CRAG IFS peaks, run the idr function from https://github.com/nboley/idr to identify CRAG IFS hotspots that pass the IDR threshold.
# 5. Call the hotspots and merge nearby hotspots accordingly, resulting in an output hotspot BED file.
# 6. Calculate the IFS scores over the hotspots for all of the cfDNA samples, resulting in one output IFS score BED file. The columns of this file are the samples, the rows are the hotspots, and the value is the IFS score.
# 7. Output all of the relevant intermediate and final files to the output directory.