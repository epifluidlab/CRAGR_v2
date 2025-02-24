import os
import sys
import glob
import gzip
import math
import shutil
import random
import subprocess

# Function to format sample names to IDS.
def format_sample(sample):
    return os.path.basename(sample).split('.')[0]

def get_basename(sample):
    return '.'.join(os.path.basename(sample).split('.')[1:])

def reverse_format_sample(id, basename):
    return f"{id}.{basename}"

def get_input_dir(sample):
    return os.path.dirname(sample)

# Required arguments.
r_path = config['r_path']
cragr_script = config['cragr_script']
samples = config['samples']
samples_list = [line.strip() for line in open(samples)]
ids_list = [format_sample(sample) for sample in samples_list]
basenames = [get_basename(sample) for sample in samples_list]
if len(set(basenames)) != 1:
    raise ValueError("The format {ID}.{BASENAME} has varied basenames in the fragment file input.")
basename = basenames[0]
input_dir = get_input_dir(samples_list[0])
chroms = config['chroms']
chroms_list = [line.strip() for line in open(chroms)]
genome = config['genome']
if genome =="hg19":
    genome = "GRCh37"
elif genome == "hg38":
    genome = "GRCh38"

chrom_sizes = config['chrom_sizes']
output_dir = config['output_dir']
os.makedirs(output_dir, exist_ok=True)
os.makedirs(output_dir + "/sample_hotspots", exist_ok=True)

# Optional arguments.
seed = config.get('seed', None)
exclude_regions = config.get('exclude_regions', None)
high_mappability = config.get('high_mappability', None)
gc_correct = config.get('gc_correct', True)
gc_correct_method = config.get('gc_correct_method', 'standard')
gc_correct_N = config.get('gc_correct_N', 1000000)
fdr_threshold = config.get('fdr_threshold', 0.2)
merge_gap = config.get('merge_gap', 200)
idr_threshold = config.get('idr_threshold', 0.05)
idr_threshold = -math.log10(idr_threshold)
min_mapq = config.get('min_mapq', 0)
min_fraglen = config.get('min_fraglen', 0)
max_fraglen = config.get('max_fraglen', 1000)
window_size = config.get('window_size', 200)
step_size = config.get('step_size', 20)
flank = (window_size - step_size) / 2

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
    temp_output = output_file.replace(".sorted.gz", "")
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
    subprocess.run(["rm", temp_output])
    subprocess.run(["tabix", "-pbed", sorted_output + ".gz"])

rule all:
    input:
        expand(output_dir + "/sample_hotspots/{sample}.signal.bedGraph.gz", sample=ids_list)

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
        rep1_output = output_dir + "/rep1.frag.sorted.gz",
        rep2_output = output_dir + "/rep2.frag.sorted.gz"
    run:
        rep1_file_list = [line.strip() for line in open(input.rep1_samples)]
        rep2_file_list = [line.strip() for line in open(input.rep2_samples)]
        combine_samples_to_reps(rep1_file_list, output.rep1_output, chroms_list)
        combine_samples_to_reps(rep2_file_list, output.rep2_output, chroms_list)

# Run the CRAGR IFS pipeline (Stage 1 Analysis).
rule run_ifs_per_chrom:
    input:
        frag_file = output_dir + "/rep{rep}.frag.sorted.gz"
    output:
        ifs_output = output_dir + "/rep{rep}.{chrom}.rawifs.bed.gz"
    params:
        gc_correct_flag = "--gc-correct" if gc_correct else "",
        gc_correct_method_param = f"--gc-correct-method={gc_correct_method}",
        gc_correct_n_param = f"--gc-correct-n={gc_correct_N}",
        mappability_param = f"--high-mappability={high_mappability}" if high_mappability else "",
        exclude_regions_param = f"--exclude-region={exclude_regions}" if exclude_regions else "",
        window_size_param = f"--window-size={window_size}",
        step_size_param = f"--step-size={step_size}",
        min_mapq_param = f"--min-mapq={min_mapq}",
        min_fraglen_param = f"--min-fraglen={min_fraglen}",
        max_fraglen_param = f"--max-fraglen={max_fraglen}"
    shell:
        """
        {r_path} {cragr_script} ifs \
        -i {input.frag_file} \
        -o {output.ifs_output} \
        --genome {genome} \
        --chrom {wildcards.chrom} \
        {params.gc_correct_flag} \
        {params.gc_correct_method_param} \
        {params.gc_correct_n_param} \
        {params.mappability_param} \
        {params.exclude_regions_param} \
        {params.window_size_param} \
        {params.step_size_param} \
        {params.min_mapq_param} \
        {params.min_fraglen_param} \
        {params.max_fraglen_param}
        """

rule concat_chrom_results:
    input:
        chrom_files = expand(output_dir + "/rep{{rep}}.{chrom}.rawifs.bed.gz", chrom=chroms_list)
    output:
        final_output = output_dir + "/rep{rep}.rawifs.bed.gz"
    run:
        temp_output = output.final_output.replace('.gz', '')
        with open(temp_output, 'w') as out_f:
            for chrom in chroms_list:
                chrom_file = f"{output_dir}/rep{wildcards.rep}.{chrom}.rawifs.bed.gz"
                with subprocess.Popen(["zcat", chrom_file], stdout=subprocess.PIPE, text=True) as proc:
                    for line in proc.stdout:
                        if not line.startswith('#'):
                            out_f.write(line)
        
        subprocess.run(["bgzip", "-f", temp_output])
        subprocess.run(["rm"] + glob.glob(f"{output_dir}/rep{wildcards.rep}.*.rawifs.bed.gz"))
        subprocess.run(["rm"] + glob.glob(f"{output_dir}/rep{wildcards.rep}.*.rawifs.bed.gz.tbi"))
        subprocess.run(["tabix", "-pbed", temp_output + ".gz"])

# Run the CRAGR PEAK pipeline (Stage 2 Analysis).
rule run_peak:
    input:
        frag_file = output_dir + "/rep{rep}.rawifs.bed.gz"
    output:
        ifs_output = output_dir + "/rep{rep}.ifs.bedGraph.gz"
    params:
        gc_correct_flag = "--gc-correct" if gc_correct else "",
        gc_correct_method_param = f"--gc-correct-method={gc_correct_method}",
        gc_correct_n_param = f"--gc-correct-n={gc_correct_N}",
        mappability_param = f"--high-mappability={high_mappability}" if high_mappability else "",
        window_size_param = f"--window-size={window_size}",
        step_size_param = f"--step-size={step_size}",
    shell:
        """
        export TMPDIR={output_dir}
        {r_path} {cragr_script} peak \
        -i {input.frag_file} \
        -o {output.ifs_output} \
        --genome {genome} \
        {params.gc_correct_flag} \
        {params.gc_correct_method_param} \
        {params.gc_correct_n_param} \
        {params.window_size_param} \
        {params.step_size_param} \
        """

# Run the IDR analysis.
rule idr_preprocess:
    input:
        ifs_file = output_dir + "/rep{rep}.ifs.bedGraph.gz",
    output:
        peak_format = output_dir + "/rep{rep}.hotspot.narrowpeak",
    shell:
        """
        zcat {input.ifs_file} | \
        awk -F'\t' -v OFS="\t" 'substr($1,1,1)!="#" && $17<={fdr_threshold} && $16!="."' | \
        bedtools slop -g {chrom_sizes} -i - -b {flank} -header | \
        bedtools merge -header -i - -d {merge_gap} -c 9,16,17 -o min | \
        perl -ne 'chomp;@f=split "\t";$f[3]=0-$f[3];$f[4]=0-log($f[4])/log(10);$f[5]=0-log($f[5])/log(10);print "$f[0]\t$f[1]\t$f[2]\t.\t0\t.\t$f[3]\t$f[4]\t$f[5]\t-1\n";' \
        > {output.peak_format}
        """

rule idr:
    input:
        rep1 = output_dir + "/rep1.hotspot.narrowpeak",
        rep2 = output_dir + "/rep2.hotspot.narrowpeak"
    output:
        idr_output = output_dir + "/hotspot.idr",
    shell:
        """
        idr --samples {input.rep1} {input.rep2} \
        --input-file-type narrowPeak \
        --rank p.value \
        --output-file {output.idr_output} \
        --output-file-type bed \
        --log-output-file {output.idr_output}.log \
        --plot \
        --peak-merge-method min
        """

rule filter_idr:
    input:
        idr_file = output_dir + "/hotspot.idr"
    output:
        idr_filtered = output_dir + "/hotspot.idr.filtered.bed.gz"
    shell:
        """
        awk '$8<={idr_threshold}' {input.idr_file} | \
        cut -f1-3,8 | \
        sort -k1,1 -k2,2n | \
        bgzip -c > {output.idr_filtered} && \
        tabix -p bed -b 2 -e 3 {output.idr_filtered}
        """

# Run the CRAGR Signal Method (Stage 4 Analysis).
rule run_signal:
    input:
        frag_file = input_dir + "/" + reverse_format_sample("{sample}", basename),
        idr_file = output_dir + "/hotspot.idr.filtered.bed.gz"
    output:
        signal_file = output_dir + "/{sample}.{chrom}.signal.bedGraph.gz"
    params:
        gc_correct_flag = "--gc-correct" if gc_correct else "",
        gc_correct_method_param = f"--gc-correct-method={gc_correct_method}",
        gc_correct_n_param = f"--gc-correct-n={gc_correct_N}",
        mappability_param = f"--high-mappability={high_mappability}" if high_mappability else "",
        exclude_regions_param = f"--exclude-region={exclude_regions}" if exclude_regions else "",
        window_size_param = f"--window-size={window_size}",
        min_mapq_param = f"--min-mapq={min_mapq}",
        min_fraglen_param = f"--min-fraglen={min_fraglen}",
        max_fraglen_param = f"--max-fraglen={max_fraglen}"
    shell:
        """
        {r_path} {cragr_script} signal \
        -i {input.frag_file} \
        --hotspot {input.idr_file} \
        -o {output.signal_file} \
        --genome {genome} \
        --chrom {wildcards.chrom} \
        {params.gc_correct_flag} \
        {params.gc_correct_method_param} \
        {params.gc_correct_n_param} \
        {params.mappability_param} \
        {params.exclude_regions_param} \
        {params.window_size_param} \
        {params.min_mapq_param} \
        {params.min_fraglen_param} \
        {params.max_fraglen_param}
        """
    
rule concat_signal_results:
    input:
        signal_files = expand(output_dir + "/{sample}.{chrom}.signal.bedGraph.gz", sample=ids_list, chrom=chroms_list),
        signal_files_index = expand(output_dir + "/{sample}.{chrom}.signal.bedGraph.gz.tbi", sample=ids_list, chrom=chroms_list)
    output:
        final_output = output_dir + "/sample_hotspots/{sample}.signal.bedGraph.gz"
    run:
        temp_output = output.final_output.replace('.gz', '')
        sorted_signal_files = sorted(input.signal_files, key=lambda f: chroms_list.index(f.split('.')[1]))
        sorted_signal_files_index = sorted(input.signal_files_index, key=lambda f: chroms_list.index(f.split('.')[1]))
        with open(temp_output, 'w') as temp_f:
            temp_f.write("#chrom\tstart\tend\tscore\tcov\tfrag_len\tgc\tscore_pre_gc\tz_score\n")
        with open(temp_output, 'a') as temp_f:
            for signal_file in sorted_signal_files:
                with gzip.open(signal_file, 'rt') as f:
                    for line in f:
                        if not line.startswith("#"):
                            temp_f.write(line)
        subprocess.run(["bgzip", temp_output])
        subprocess.run(["rm"] + sorted_signal_files)
        subprocess.run(["rm"] + sorted_signal_files_index)
        subprocess.run(["tabix", "-pbed", temp_output + ".gz"])
        subprocess.run(["rm", temp_output])
