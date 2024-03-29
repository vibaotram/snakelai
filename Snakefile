import pandas as pd
import math
import os
import vcf
import re


# configfile: "config.yaml"


def unique_chrom(vcf):

    # intilize a null list
    unique_list = []

    contigs = vcf.contigs
    for c in contigs:
        if re.match("Contig\\d+", c):
            pass
        else:
            unique_list.append(c)
    return unique_list


outdir = config["outdir"]



SIM_SOURCE = config["simulate_source"]
TEST_FILE = config["test_file"]


# elai = config["elai_source"]

######################
#### PREPARE SOURCE GENOTYPES

###### SPLIT SOURCE GENOTYPES BY CHROMOSOME

source_vcf = config['source_vcf']
source_name = os.path.splitext(os.path.splitext(os.path.basename(source_vcf))[0])[0]
vcftools_module = config['vcftools']

chrom_id= config['chromosome']
vcf_reader = vcf.Reader(open(source_vcf, 'r'))
source_chromosome = unique_chrom(vcf_reader)
if chrom_id == '@':
    chromosome = source_chromosome
# elif not all(c in source_chromosome for c in chromosome):
#     raise ValueError('Specified chromosomes do not present in the vcf file')
else:
    chrom_id = [int(i) - 1 for i in chrom_id]
    chromosome = [source_chromosome[i] for i in chrom_id]
    # pass

snp_selection = list(config['snp_selection'].keys())

elai_params = list(config['elai_params'].keys())

rule finish:
    input:
        #expand(os.path.join(outdir, "final_results/output_{elai_params}/{files}"), elai_params = elai_params, files = ["overall_admixture.tsv", "local_dosage.tsv"])
        expand(os.path.join(outdir, "elai_results", "{chromosome}", "{snp_selection}", "{elai_params}", "elai_results.html"), chromosome=chromosome, snp_selection=snp_selection, elai_params=elai_params)

rule split_source_chrom:
    input: source_vcf
    output: os.path.join(outdir, 'source', '{chromosome}', source_name+'_{chromosome}.recode.vcf')
    params:
        logname = "split_source_chrom_{chromosome}",
        logdir = os.path.join(outdir, "log")
    envmodules: vcftools_module
    # conda: "conda/conda_vcftools.yaml"
    # singularity: "shub://vibaotram/singularity-container:guppy4.0.14gpu-conda-api"
    shell:
        """
        cd $(dirname {output})
        vcftools --vcf {input} --chr {wildcards.chromosome} --out {source_name}_{wildcards.chromosome} --recode
        """

###### CREATE ELAI SOURCE FILE

nb_groups = config['nb_groups']
k = range(1, nb_groups+1)
nb_genotypes = config['nb_genotypes']

rule simulate_source:
    input: rules.split_source_chrom.output
    output: expand(os.path.join(outdir, 'source', '{chromosome}', 'group_{k}'), chromosome = "{chromosome}", k = k)
    params:
        # nb_snps = nb_snps,
        # chrom_length = chrom_length,
        # nb_batches = nb_batches,
        output = lambda wildcards, output: ','.join(output),
        nb_groups = nb_groups,
        nb_genotypes = nb_genotypes,
        script = "script/simulate_source.R",
        logname = "simulate_source_{chromosome}",
        logdir = os.path.join(outdir, "log")
    threads: config['simulate_cores']
    resources:
        mem_gb = config["simulate_mem_gb"]
    singularity: config['singularity']
    shell:
        """
        Rscript {params.script} -o {params.output} -i {input} -g {params.nb_groups} -n {params.nb_genotypes} -t {threads}
        """



######################
#### PREPARE TEST FILE AND SNP FILE

test_vcf = config['test_vcf']
test_id = config['test_id']

vcftools_indv = []

if test_id == "@":
    vcftools_indv = ''
elif TEST_FILE == "" and not all(id in vcf_reader.samples for id in test_id):
    raise ValueError("Specified genotype IDs do not present in the test_vcf file")
else:
    for id in test_id:
        vcftools_indv.append("--indv {id}".format(id=id))



rule split_test_chrom:
    input: test_vcf
    output: os.path.join(outdir, 'test', '{chromosome}', 'test_{chromosome}.recode.vcf')
    params:
        indv = vcftools_indv,
        logname = "split_test_chrom_{chromosome}",
        logdir = os.path.join(outdir, "log")
    envmodules: vcftools_module
    # conda: "conda/conda_vcftools.yaml"
    # singularity: "shub://vibaotram/singularity-container:guppy4.0.14gpu-conda-api"
    shell:
        """
        cd $(dirname {output})
        vcftools --vcf {input} --chr {wildcards.chromosome} {params.indv} --out test_{wildcards.chromosome} --recode
        """

rule test_file:
    input: rules.split_test_chrom.output
    output:
        test_file = os.path.join(outdir, 'test', '{chromosome}', 'test_{chromosome}.geno'),
        snp_file = os.path.join(outdir, 'pos', '{chromosome}', 'all_snp.pos')
    threads: 1 #config['split_snps_cores']
    params:
        # test_id = test_id,
        script = "script/create_test_file.R",
        logname = "test_file_{chromosome}",
        logdir = os.path.join(outdir, "log")
    # conda: "conda/conda_rmarkdown.yaml"
    singularity: config['singularity']
    shell:
        """
        Rscript {params.script} {output.test_file} {output.snp_file} {input}
        """

######################
#### SELECT SNPS

# snp_file = config["snp_file"]
#
# random_snps = config["random_snps"]
#
# batch_size = config["batch_size"]
#
# snp_den = config["snp_den"]
#
# if random_snps:
#     n_batch = config["n_batches"]
# else:
#     snp_info = pd.read_table(snp_file, delim_whitespace=True, header=None)
#     n_batch = math.ceil(len(snp_info)/batch_size)
#
# batchid = list(range(1, n_batch + 1))

rule select_snps:
    input: rules.test_file.output.snp_file if TEST_FILE == "" else config["snp_file"]
    output:
        os.path.join(outdir, 'pos', '{chromosome}', "{snp_selection}", 'select_snps.done')
        # select_snps_output
        # expand(os.path.join(outdir, "batch_{batchid}/snp_pos"), batchid = batchid)
    params:
        nb_snps = lambda wildcards: config["snp_selection"][wildcards.snp_selection]["nb_snps"],
        window_size = lambda wildcards: config["snp_selection"][wildcards.snp_selection]["window_size"],
        script = "script/split_snps.R",
        logname = "select_snps_{chromosome}_{snp_selection}",
        logdir = os.path.join(outdir, "log")
    threads: config['split_snps_cores']
    # conda: "conda/conda_rmarkdown.yaml"
    singularity: config['singularity']
    shell:
        """
        Rscript {params.script} -i {input} -o {output} -n {params.nb_snps} -w {params.window_size} -t {threads}
        """

elai = config['elai']
# source_genotypes = config["source_genotypes"]
# test_genotype = config["test_genotype"]
elai_ext = ["admix.txt", "em.txt", "log.txt", "ps21.txt", "snpinfo.txt"]

# source_genotypes = rules.simulate_source.output
# n_sources = len(source_genotypes)
# source_files = ""
# for i in range(0, n_sources):
#     file = os.path.dirname(source_genotypes[i])
#     source_files+= "-g {file} -p 1{i} ".format(file=file, i=i)




rule elai:
    input:
        source_files = rules.simulate_source.output if SIM_SOURCE else config['source_files'],
        test_file = rules.test_file.output.test_file if TEST_FILE == "" else TEST_FILE,
        snp_file = rules.select_snps.output,
    output:
        os.path.join(outdir, "elai_results", "{chromosome}", "{snp_selection}", "{elai_params}", "elai_results.html")
        # expand(os.path.join(outdir, "batch_{batchid}/output_{elai_params}/elai_r.{elai_ext}"), batchid = "{batchid}", elai_params = "{elai_params}", elai_ext = elai_ext)
        # expand(os.path.join(outdir, "elai", "{chromosome}", "{snp_selection}", "{elai_params}", "elai_r.{elai_ext}"), chromosome= "{chromosome}", snp_selection="{snp_selection}", elai_params = "{elai_params}", elai_ext=elai_ext)
    params:
        elai = elai,
        options = lambda wildcards: config["elai_params"][wildcards.elai_params],
        # out_dir = lambda wildcards, output: os.path.basename(os.path.dirname(output[0])),
        workdir = os.path.join(outdir, "elai", "{chromosome}", "{snp_selection}", "{elai_params}"),
        logname = "elai_{chromosome}_{snp_selection}_{elai_params}",
        logdir = os.path.join(outdir, "log"),
        mem_gb = config["elai_mem_gb"]
    resources:
        mem_gb = config["elai_mem_gb"]
    # singularity: "shub://vibaotram/singularity-container:guppy4.0.14gpu-conda-api"
    conda: "conda/conda_rmarkdown.yaml"
    script: "script/elai.Rmd"

merged_snps = config["merged_snps"]
merged_params = config["merged_params"]

rule merge_elai:
    input:
        require = rules.finish.input,
        # files = expand(rules.elai.output, chromosome = chromosome, snp_selection = merged_snps, elai_params = merged_params)
    output: os.path.join(outdir, "Merged_ELAI_Results.html")
    params:
        merged_snps = merged_snps,
        merged_params = merged_params,
        genome_fa = config["genome"]
    conda: "conda/conda_rmarkdown.yaml"
    script: "script/merge_elai.Rmd"

# rule read_elai:
#     input: expand(rules.elai.output, chromosome = "{chromosome}", snp_selection="{snp_selection}", elai_params = "{elai_params}", elai_ext = ["admix.txt", "ps21.txt"])
#     output: os.path.join(outdir, "elai_results", "{chromosome}", "{snp_selection}", "{elai_params}", "elai_results.html")
#     params:
#         elai_options = lambda wildcards: config["elai_params"][wildcards.elai_params],
#         adm_file = lambda wildcards, output: os.path.join(os.path.dirname(output[0]), "overall_admixture.tsv"),
#         dosage_file = lambda wildcards, output: os.path.join(os.path.dirname(output[0]), "local_dosage.tsv"),
#         # true_dosage_file = config["true_inference"],
#         # genome_file = config["genome"],
#         logname = "read_elai_{chromosome}_{snp_selection}_{elai_params}",
#         logdir = os.path.join(outdir, "log")
#     conda: "conda/conda_rmarkdown.yaml"
#     script: "script/read_elai.Rmd"
