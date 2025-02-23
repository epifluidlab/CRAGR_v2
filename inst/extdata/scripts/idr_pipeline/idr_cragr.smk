import os
import sys
import glob
import gzip
import shutil
import random
import subprocess

# Required arguments.
r_script = config['r_script']
samples = config['samples']
samples_list = [line.strip() for line in open(samples)]
chroms = config['chroms']
chroms_list = [line.strip() for line in open(chroms)]
genome = config['genome']
chrom_sizes = config['chrom_sizes']
output_dir = config['output_dir']
os.makedirs(output_dir, exist_ok=True)

# Optional arguments.
seed = config.get('seed', None)
exclude_regions = config.get('exclude_regions', None)
high_mappability = config.get('high_mappability', None)
gc_correct = config.get('gc_correct', True)
gc_correct_method = config.get('gc_correct_method', 'standard')
gc_correct_N = config.get('gc_correct_N', 1000000)
idr_threshold = config.get('idr_threshold', 1)
min_mapq = config.get('min_mapq', 0)
min_fraglen = config.get('min_fraglen', 0)
max_fraglen = config.get('max_fraglen', 1000)
window_size = config.get('window_size', 200)
step_size = config.get('step_size', 20)

# Function to split the files into two replicates.
def split_files_to_reps(list_of_filenames, seed=None):
    if seed is not None:
        random.seed(seed)
    shuffled = list_of_filenames[:]
    random.shuffle(shuffled)
    mid = len(list_of_filenames) // 2
    rep1, rep2 = shuffled[:mid], shuffled[mid:]
    return rep1, rep2

# Function to combine the files into merged replicate file.
def combine_samples_to_reps(file_list, output_file, chroms_list):
    temp_output = output_file.replace(".sorted.gz.tbi", "")
    with open(temp_output, 'wb') as out_f:
        for file in file_list:
            with subprocess.Popen(["zcat", file], stdout=subprocess.PIPE) as proc:
                shutil.copyfileobj(proc.stdout, out_f)
    
    sorted_output = temp_output + ".sorted"
    with open(temp_output, 'r') as in_f, open(sorted_output, 'w') as out_f:
        lines = [line.strip().split('\t') for line in in_f if not line.startswith('#')]
        lines = [line for line in lines if line[0] in chroms_list]
        chrom_order = {chrom: i for i, chrom in enumerate(chroms_list)}
        lines.sort(key=lambda x: (chrom_order.get(x[0], float('inf')), int(x[1]), int(x[2])))
        for line in lines:
            out_f.write('\t'.join(line) + '\n')
    
    subprocess.run(["bgzip", "-f", sorted_output])
    subprocess.run(["tabix", "-pbed", sorted_output + ".gz"])

rule all:
    input:
        expand(output_dir + "/rep{rep}.frag.sorted.gz.tbi", rep=[1, 2])

# Split the samples into two replicate lists.
rule split_files:
    output:
        rep1 = output_dir + "/rep1_samples.txt",
        rep2 = output_dir + "/rep2_samples.txt"
    run:
        rep1_samples, rep2_samples = split_files_to_reps(samples_list, seed=seed)
        with open(output.rep1, 'w') as f1, open(output.rep2, 'w') as f2:
            f1.write("\n".join(rep1_samples) + "\n")
            f2.write("\n".join(rep2_samples) + "\n")

# Combine the samples into merged replicate files.
rule combine_samples:
    input:
        rep1_samples = output_dir + "/rep1_samples.txt",
        rep2_samples = output_dir + "/rep2_samples.txt"
    output:
        rep1_output = output_dir + "/rep1.frag.sorted.gz.tbi",
        rep2_output = output_dir + "/rep2.frag.sorted.gz.tbi"
    run:
        rep1_file_list = [line.strip() for line in open(input.rep1_samples)]
        rep2_file_list = [line.strip() for line in open(input.rep2_samples)]
        combine_samples_to_reps(rep1_file_list, output.rep1_output, chroms_list)
        combine_samples_to_reps(rep2_file_list, output.rep2_output, chroms_list)

# Run the CRAGR pipeline.
