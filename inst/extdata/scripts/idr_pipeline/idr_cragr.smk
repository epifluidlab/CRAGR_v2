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
hotspot_samples = config.get('hotspot_samples', None)
if hotspot_samples is not None:
    hotspot_samples_list = [line.strip() for line in open(hotspot_samples)]
    ids_list = [format_sample(sample) for sample in hotspot_samples_list]
else:
    ids_list = [format_sample(sample) for sample in samples_list]
basenames = [get_basename(sample) for sample in samples_list]
if len(set(basenames)) != 1:
    raise ValueError("The format {ID}.{BASENAME} has varied basenames in the fragment file input.")
basename = basenames[0]
input_dir = get_input_dir(samples_list[0])
chroms = config['chroms']
chroms_list = [line.strip() for line in open(chroms)]
chroms_string = ":".join(chroms_list)
genome = config['genome']
if genome =="hg19":
    genome = "GRCh37"
elif genome == "hg38":
    genome = "GRCh38"

chrom_sizes = config['chrom_sizes']
output_dir = config['output_dir']
if output_dir.endswith('/'):
    output_dir = output_dir[:-1]
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.join(output_dir, "sample_hotspots"), exist_ok=True)

# Optional arguments.
split_method = config.get('split_method', 'fragment_count')
if split_method not in ['fragment_count', 'by_sample']:
    raise ValueError("The split method must be either 'fragment_count' or 'by_sample'.")
seed = config.get('seed', None)
subsample = config.get('subsample', 200_000_000)
exclude_regions = config.get('exclude_regions', None)
high_mappability = config.get('high_mappability', None)
gc_correct = config.get('gc_correct', True)
gc_correct_method = config.get('gc_correct_method', 'standard')
gc_correct_N = config.get('gc_correct_N', 1000000)
fdr_threshold = config.get('fdr_threshold', 0.2)
merge_gap = config.get('merge_gap', 200)
idr_threshold = config.get('idr_threshold', 0.05)
idr_threshold = -math.log10(idr_threshold)
min_mapq = config.get('min_mapq', 30)
min_fraglen = config.get('min_fraglen', 0)
max_fraglen = config.get('max_fraglen', 1000)
window_size = config.get('window_size', 200)
step_size = config.get('step_size', 20)
flank = (window_size - step_size) / 2
total_fragment_min = config.get('total_fragment_min', 1.5*subsample)
threads = config.get('threads', 1)
n_cores = workflow.cores

# Configure resource settings
def get_resource_specs(rule_name):
    resource_specs = {
        'default': {
            'mem_mb': 4000,
            'time': '2:00:00'
        },
        'combine_samples': {
            'mem_mb': 16000,
            'time': '2:00:00'
        },
        'run_ifs_per_chrom': {
            'mem_mb': 2000,
            'time': '2:00:00'
        },
        'run_peak': {
            'mem_mb': 8000,
            'time': '2:00:00'
        },
        'run_signal': {
            'mem_mb': 8000,
            'time': '2:00:00'
        },
        'idr': {
            'mem_mb': 8000,
            'time': '1:00:00'
        }
    }
    return resource_specs.get(rule_name, resource_specs['default'])

# Function to split the files into two replicates.
def split_files_to_reps(list_of_filenames, split_method, seed=None):
    if split_method=="fragment_count":
        return [], []
    if seed is not None:
        random.seed(seed)
    shuffled = list_of_filenames[:]
    random.shuffle(shuffled)
    mid = len(list_of_filenames) // 2
    rep1, rep2 = shuffled[:mid], shuffled[mid:]
    return rep1, rep2

# Function to combine the files into merged replicate file.
def combine_samples_to_reps(file_list, output_file, chroms, chroms_list, split_method, samples_list, total_fragment_min, seed, subsample, min_fraglen, max_fraglen, min_mapq, threads):
    if split_method == "fragment_count":
        if seed is not None:
            random.seed(seed)
        combined_path = os.path.join(os.path.dirname(output_file), "combined_fragments")
        length_path = os.path.join(os.path.dirname(output_file), "combined_fragments_len")
        line_count = 0
        if not os.path.exists(combined_path):
            with open(combined_path, 'wb') as out_f, open(length_path, 'w') as len_f:
                for file in samples_list:
                    with subprocess.Popen(["zcat", file], stdout=subprocess.PIPE) as proc:
                        for line in proc.stdout:
                            line = line.decode('utf-8').strip()
                            if not line.startswith('#'):
                                fields = line.split('\t')
                                chrom = fields[0]
                                start = int(fields[1])
                                end = int(fields[2])
                                mapq = int(fields[4])
                                fraglen = end - start
                                if chrom in chroms_list and fraglen >= min_fraglen and fraglen <= max_fraglen and mapq >= min_mapq:
                                    out_f.write((line + '\n').encode('utf-8'))
                                    line_count += 1
                len_f.write(str(line_count))
        temp_output = output_file.replace(".sorted.gz", "")
        with open(length_path, 'r') as len_f:
            line_count = int(len_f.readline())
        if line_count < total_fragment_min:
            raise ValueError("The total number of fragments in the input files is less than the minimum number of fragments to run the pipeline. Check the documentation for more details.")
        selected_indices = set(random.sample(range(line_count), subsample))
        with open(combined_path, 'rb') as in_f, open(temp_output, 'wb') as out_f:
            for i, line in enumerate(in_f):
                if i in selected_indices:
                    out_f.write(line)
        sorted_output = temp_output + ".sorted.gz"
        print("Sorting the combined fragments.")
        sort_command = ["sort", "-V", "-k1,1", "-k2,2n", f"--parallel={threads}", temp_output]
        bgzip_command = ["bgzip", "-c", "-@", f"{threads}"]
        with open(sorted_output, "wb") as out_file:
            sort_process = subprocess.Popen(sort_command, stdout=subprocess.PIPE)
            subprocess.run(bgzip_command, stdin=sort_process.stdout, stdout=out_file)
        sort_process.wait()
        subprocess.run(["tabix", "-pbed", "-@", f"{threads}", sorted_output])


    else:
        temp_output = output_file.replace(".sorted.gz", "")
        line_count = 0
        with open(temp_output, 'wb') as out_f:
            for file in file_list:
                with subprocess.Popen(["zcat", file], stdout=subprocess.PIPE) as proc:
                    for line in proc.stdout:
                        line = line.decode('utf-8').strip()
                        if not line.startswith('#'):
                            fields = line.split('\t')
                            chrom = fields[0]
                            start = int(fields[1])
                            end = int(fields[2])
                            mapq = int(fields[4])
                            fraglen = end - start
                            if chrom in chroms_list and fraglen >= min_fraglen and fraglen <= max_fraglen and mapq >= min_mapq:
                                out_f.write((line + '\n').encode('utf-8'))
                                line_count += 1
        sorted_output = temp_output + ".sorted.gz"
        print("Sorting the combined fragments.")
        sort_command = ["sort", "-V", "-k1,1", "-k2,2n", f"--parallel={threads}", temp_output]
        bgzip_command = ["bgzip", "-c", "-@", f"{threads}"]
        with open(sorted_output, "wb") as out_file:
            sort_process = subprocess.Popen(sort_command, stdout=subprocess.PIPE)
            subprocess.run(bgzip_command, stdin=sort_process.stdout, stdout=out_file)
        sort_process.wait()
        subprocess.run(["tabix", "-pbed", "-@", f"{threads}", sorted_output])

rule all:
    input:
        expand(output_dir + "/sample_hotspots/{sample}.signal.bedGraph.gz", sample=ids_list),
        expand(output_dir + "/sample_hotspots/combined.signal.bedGraph.gz"),
        expand(output_dir + "/cleanup.done")

# Split the samples into two replicate lists.
rule split_files:
    output:
        rep1 = output_dir + "/rep1_samples.txt",
        rep2 = output_dir + "/rep2_samples.txt"
    resources:
        mem_mb=get_resource_specs('default')['mem_mb'],
        time=get_resource_specs('default')['time']
    run:
        rep1_samples, rep2_samples = split_files_to_reps(samples_list, split_method, seed=seed)
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
    threads: n_cores
    resources:
        mem_mb=get_resource_specs('combine_samples')['mem_mb'],
        time=get_resource_specs('combine_samples')['time']
    run:
        rep1_file_list = [line.strip() for line in open(input.rep1_samples)]
        rep2_file_list = [line.strip() for line in open(input.rep2_samples)]
        print("Combining samples into replicate 1 and replicate 2. Using threads: ", threads)
        combine_samples_to_reps(rep1_file_list, output.rep1_output, chroms, chroms_list, split_method, samples_list, total_fragment_min, seed, subsample, min_fraglen, max_fraglen, min_mapq, threads)
        combine_samples_to_reps(rep2_file_list, output.rep2_output, chroms, chroms_list, split_method, samples_list, total_fragment_min, seed+1, subsample, min_fraglen, max_fraglen, min_mapq, threads)

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
    threads: threads
    resources:
        mem_mb=get_resource_specs('run_ifs_per_chrom')['mem_mb'],
        time=get_resource_specs('run_ifs_per_chrom')['time']
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
        {params.max_fraglen_param} \
        -t {threads}
        """

rule concat_chrom_results:
    input:
        chrom_files = expand(output_dir + "/rep{{rep}}.{chrom}.rawifs.bed.gz", chrom=chroms_list)
    output:
        final_output = output_dir + "/rep{rep}.rawifs.bed.gz"
    threads: 1
    resources:
        mem_mb=get_resource_specs('default')['mem_mb'],
        time=get_resource_specs('default')['time']
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
        #subprocess.run(["rm"] + glob.glob(f"{output_dir}/rep{wildcards.rep}.*.rawifs.bed.gz"))
        #subprocess.run(["rm"] + glob.glob(f"{output_dir}/rep{wildcards.rep}.*.rawifs.bed.gz.tbi"))
        subprocess.run(["tabix", "-pbed", temp_output + ".gz"])

# Run the CRAGR PEAK pipeline (Stage 2 Analysis).
rule run_peak:
    input:
        frag_file = output_dir + "/rep{rep}.rawifs.bed.gz"
    output:
        ifs_output = output_dir + "/rep{rep}.ifs.bedGraph.gz"
    threads: threads
    resources:
        mem_mb=get_resource_specs('run_peak')['mem_mb'],
        time=get_resource_specs('run_peak')['time']
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
        -t {threads}
        """

# Run the IDR analysis.
rule idr_preprocess:
    input:
        ifs_file = output_dir + "/rep{rep}.ifs.bedGraph.gz",
    output:
        peak_format = output_dir + "/rep{rep}.hotspot.narrowpeak",
    threads: 1
    resources:
        mem_mb=get_resource_specs('default')['mem_mb'],
        time=get_resource_specs('default')['time']
    shell:
        """
        zcat {input.ifs_file} | \
        awk -F'\t' -v OFS="\t" 'substr($1,1,1)!="#" && $17<={fdr_threshold} && $17!="."' | \
        bedtools slop -g {chrom_sizes} -i - -b {flank} -header | \
        bedtools merge -header -i - -d {merge_gap} -c 9,16,17 -o min | \
        perl -ne 'chomp;@f=split "\t";$f[3]=0-$f[3];$f[4]=0-log($f[4])/log(10);$f[5]=0-log($f[5])/log(10);print "$f[0]\t$f[1]\t$f[2]\t.\t0\t.\t$f[3]\t$f[4]\t$f[5]\t-1\n";' \
        > {output.peak_format}

        line_count=$(wc -l < {output.peak_format})
        echo "# of peaks in {output.peak_format}: $line_count"
        if [ "$line_count" -lt 50000 ]; then
            echo "WARNING: {output.peak_format} has fewer than 50,000 peaks!"
        fi
        """

rule idr:
    input:
        rep1 = output_dir + "/rep1.hotspot.narrowpeak",
        rep2 = output_dir + "/rep2.hotspot.narrowpeak"
    output:
        idr_output = output_dir + "/hotspot.idr",
    threads: 1 
    resources:
        mem_mb=get_resource_specs('idr')['mem_mb'],
        time=get_resource_specs('idr')['time']
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
    threads: 1
    resources:
        mem_mb=get_resource_specs('default')['mem_mb'],
        time=get_resource_specs('default')['time']
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
        signal_file = output_dir + "/sample_hotspots/{sample}.signal.bedGraph.gz"
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
    threads: threads
    resources:
        mem_mb=get_resource_specs('run_signal')['mem_mb'],
        time=get_resource_specs('run_signal')['time']
    shell:
        """
        {r_path} {cragr_script} signal \
        -i {input.frag_file} \
        --hotspot {input.idr_file} \
        -o {output.signal_file} \
        --genome {genome} \
        --chrom {chroms_string} \
        {params.gc_correct_flag} \
        {params.gc_correct_method_param} \
        {params.gc_correct_n_param} \
        {params.mappability_param} \
        {params.exclude_regions_param} \
        {params.window_size_param} \
        {params.min_mapq_param} \
        {params.min_fraglen_param} \
        {params.max_fraglen_param} \
        -t {threads}
        """

# Combine the signal files into a single file.
def combine_samples(list_of_filenames, output_dir):
    import pandas as pd
    list_of_filenames = list(list_of_filenames)
    dfs = []
    for file in list_of_filenames:
        id = format_sample(file)
        df = pd.read_csv(file, sep='\t', header=None, comment='#', usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'signal'])
        df.rename(columns={'signal': id}, inplace=True)
        df.replace('.', float('nan'), inplace=True)
        df[id] = df[id].astype(float)
        dfs.append(df)
    combined_df = dfs[0]
    for df in dfs[1:]:
        combined_df = pd.merge(combined_df, df, on=['chrom', 'start', 'end'], how='outer')
    combined_df.to_csv(output_dir + "/sample_hotspots/combined.signal.bedGraph.gz", sep='\t', index=False, compression='gzip')


rule combine_signal:
    input:
        signal_files = expand(output_dir + "/sample_hotspots/{sample}.signal.bedGraph.gz", sample=ids_list)
    output:
        combined_signal = output_dir + "/sample_hotspots/combined.signal.bedGraph.gz"
    threads: 1
    resources:
        mem_mb=16000,
        time='2:00:00'
    run:
        combine_samples(input.signal_files, output_dir)

rule cleanup:
    input:
        expand(output_dir + "/sample_hotspots/{sample}.signal.bedGraph.gz", sample=ids_list),
        output_dir + "/sample_hotspots/combined.signal.bedGraph.gz"
    output:
        output_dir + "/cleanup.done"
    threads: 1
    resources:
        mem_mb=2000,
        time='0:30:00'
    shell:
        """
        # Remove intermediate files from rep1 and rep2
        rm -f {output_dir}/rep*.frag
        rm -f {output_dir}/rep*.*.rawifs.bed.gz
        rm -f {output_dir}/rep*.*.rawifs.bed.gz.tbi
        
        # Remove intermediate bed and bedGraph files
        # rm -f {output_dir}/rep*.hotspot.narrowpeak
        
        # Remove other temporary files
        rm -f {output_dir}/combined_fragments*
        rm -f {output_dir}/*.log
        
        # Keep the final IDR results and signal files
        echo "Cleanup completed. All intermediate files have been removed."
        touch {output_dir}/cleanup.done
        """