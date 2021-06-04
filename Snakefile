import pandas as pd
import math
import os

# configfile: "config.yaml"

outdir = config["outdir"]

snp_file = config["snp_file"]

random_snps = config["random_snps"]

batch_size = config["batch_size"]

snp_den = config["snp_den"]

if random_snps:
    n_batch = config["n_batches"]
else:
    snp_info = pd.read_table(snp_file, delim_whitespace=True, header=None)
    n_batch = math.ceil(len(snp_info)/batch_size)

batchid = list(range(1, n_batch + 1))

elai = config["elai_source"]
source_genotypes = config["source_genotypes"]
test_genotype = config["test_genotype"]
elai_params = list(config["elai_params"].keys())
elai_ext = ["admix.txt", "em.txt", "log.txt", "ps21.txt", "snpinfo.txt"]

if random_snps:
    finish_input = expand(os.path.join(outdir, "batch_{batchid}/output_{elai_params}/elai_results.html"), batchid = batchid, elai_params = elai_params)
else:
    finish_input = expand(os.path.join(outdir, "final_results/output_{elai_params}/elai_results.html"), elai_params = elai_params)

rule finish:
    input:
        #expand(os.path.join(outdir, "final_results/output_{elai_params}/{files}"), elai_params = elai_params, files = ["overall_admixture.tsv", "local_dosage.tsv"])
        finish_input

rule split_snps:
    input: snp_file
    output: expand(os.path.join(outdir, "batch_{batchid}/snp_pos"), batchid = batchid)
    params:
        random_snps = random_snps,
        batch = batchid,
        batch_size = batch_size,
        snp_den = snp_den,
        logname = "split_snps",
        logdir = os.path.join(outdir, "log")
    threads: config['split_snps_cores']
    script: "script/split_snps.R"

# for i in batchid:
#     os.makedirs(os.path.join(outdir, "batch_{}".format(i)), exist_ok=True)
#     start = 1000*(i-1)
#     end = 1000*i
#     snp_batch = snp_info[start:end]
#     snp_batch.to_csv(os.path.join(outdir, "batch_{}/snp_pos".format(i)), sep='\t', header=False, index=False)

n_sources = len(source_genotypes)
source_files = ""
for i in range(0, n_sources):
    file = source_genotypes[i]
    source_files+= "-g {file} -p 1{i} ".format(file=file, i=i)


rule elai:
    input:
        source_genotypes = source_genotypes,
        test_file = test_genotype,
        snp_file = os.path.join(outdir, "batch_{batchid}/snp_pos"),
    output:
        expand(os.path.join(outdir, "batch_{batchid}/output_{elai_params}/elai_r.{elai_ext}"), batchid = "{batchid}", elai_params = "{elai_params}", elai_ext = elai_ext)
    params:
        source_files = source_files,
        options = lambda wildcards: config["elai_params"][wildcards.elai_params],
        # out_dir = lambda wildcards, output: os.path.basename(os.path.dirname(output[0])),
        workdir = lambda wildcards, output: os.path.dirname(output[0]),
        logname = "elai_{batchid}_{elai_params}",
        logdir = os.path.join(outdir, "log")
    resources:
        mem_gb = config["elai_mem_gb"]
    shell:
        """
        cd {params.workdir}
        echo "$({elai} \\
        {params.source_files} \\
        -g {input.test_file} -p 1 \\
        -pos {input.snp_file} \\
        -o elai_r \\
        -C {n_sources} {params.options})"

        mv output/* ./
        rm -rf output
        """

if random_snps:
    aggregate_outdir = os.path.join(outdir, "batch_{batchid}")
else:
    aggregate_outdir = os.path.join(outdir, "final_results")

rule read_elai:
    input: expand(os.path.join(outdir, "batch_{batchid}/output_{elai_params}/elai_r.{elai_ext}"), batchid = "{batchid}" if random_snps else batchid, elai_params = "{elai_params}", elai_ext = ["admix.txt", "ps21.txt"])
    output: os.path.join(aggregate_outdir, "output_{elai_params}/elai_results.html")
    params:
        elai_options = lambda wildcards: config["elai_params"][wildcards.elai_params],
        adm_file = os.path.join(aggregate_outdir, "output_{elai_params}/overall_admixture.tsv"),
        dosage_file = os.path.join(aggregate_outdir, "output_{elai_params}/local_dosage.tsv"),
        true_dosage_file = config["true_inference"],
        genome_file = config["genome"],
        logname = "read_elai_{batchid}_{elai_params}" if random_snps else "read_elai_{elai_params}",
        logdir = os.path.join(outdir, "log")
    conda: "conda_rmarkdown.yaml"
    script: "script/read_elai.Rmd"

rule test:
    shell: "echo {source_files}"

