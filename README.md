# Snakemake workflow to run ELAI with different parameters and SNP sets, and scripts to get average of multiple runs

*Scripts here are adapted on IFB cluster. Branch 'itrop' is for scripts adapted on iTrop.*

## Run ELAI on simulated genotypes for testing (1 chromosome only)
### Input
- source files, test file, and snp file in ELAI input format containing all SNPs
- define SNP sets

Here is my [config file](./example/config_sim_data.yaml).

### Steps
- generate or spit SNP sets as predefined
- run ELAI for each set of parameters and snp sets
- evaluate the results by [custom scripts](https://github.com/vibaotram/vn_robusta_ELAI/tree/master/validate_elai)

Here is my [script](./example/run_validate_elai.sh) to run Snakemake.

## Run ELAI on real data (multiple chromosomes possibly)
### Input
- vcf file containing all the source individuals
- vcf file containing all the test individuals
- define SNP sets

Here is my [config file](./example/config_real_data.yaml).

### Steps
- spit the source and test genotypes by chromosomes
- simulate source populations using the source individuals and generate source files for ELAI
- generate test file for ELAI from the vcf file
- generate or spit SNP sets as predefined
- run ELAI for each set of parameters and snp sets
- rerun ELAI by `snakemake --forcerun elai`
- average results from ELAI runs by [custom script](./script/average_elai.R)
- merge results from runs of different snp sets by [custom script](./script/merge_elai.R)

Here is my [script](./example/run_validate_elai.sh) to run Snakemake.

## Some notes
- the workflow runs on singularity container: shub://vibaotram/singularity-container:myr_4-0-2_rstudio_1.3.sif
- may need to change R lib path in the [script files](./script), and install all the packages before running to prevent failure of parallel jobs
