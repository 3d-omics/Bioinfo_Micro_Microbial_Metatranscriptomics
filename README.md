# Snakemake workflow: `Bioinfo_Macro_Microbial_Metatranscriptomics`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `<description>`


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.

##Set up required softwares

## Usage
  ```
  #Clone the git repository in your terminal
  git clone git@github.com:3d-omics/Bioinfo_Macro_Microbial_Metatranscriptomics.git
  #Change directory to the one you cloned in the previous step
  cd Bioinfo_Macro_Microbial_Metatranscriptomics
  #Activate conda environment where you have snakemake
  conda activte Snakemake
  #run the pipeline with the test data, it will download all the necesary software through conda. It should take less than 5 minutes.
  snakemake --use-conda --jobs 8 all
  ``` 

- Run it with your own data:
  - Edit `config/samples.tsv` and add your samples and where are they located. Here is an example of the tsv table filled with the information
    
    <img width="828" alt="image" src="https://github.com/3d-omics/Bioinfo_Macro_Microbial_Metatranscriptomics/assets/103645443/1afac14d-8bc5-47da-b154-6e0c75f0e255">


  - Edit `config/features.yml` with information regarding the reference and the mags you are
    using like in this example.
    
    <img width="476" alt="image" src="https://github.com/3d-omics/Bioinfo_Macro_Microbial_Metatranscriptomics/assets/103645443/99af1063-9a25-4ef8-adab-1a7af1df4e64">


  - Edit `config/params.yml` to change the execution of the steps like in this example
    
    <img width="548" alt="image" src="https://github.com/3d-omics/Bioinfo_Macro_Microbial_Metatranscriptomics/assets/103645443/33187236-3b46-40ac-997a-2d1b49466c60">


## Features
- FASTQ processing with `fastp`
- taxonomic classification of preprocessed reads with `kraken2`
- detection of ribosomes with `ribodetector`
- allignment of reads to the host with `star`
- allignment of non host reads to mag catalogue with `bowtie2`
- SAM/BAM/CRAM processing with `samtools`
- coverage eveluation with `coverm`
- Reports with `multiqc` and `FastQC`
  
## DAG

![image](https://github.com/3d-omics/Bioinfo_Macro_Microbial_Metatranscriptomics/assets/103645443/59fb81c6-47a8-4307-9056-1eab5f180303)



## References

- [`fastp`](https://github.com/OpenGene/fastp)
- [`STAR`](https://github.com/alexdobin/STAR)
- [`samtools`](https://github.com/samtools/samtools)
- [`FastQC`](https://github.com/s-andrews/FastQC)
- [`multiqc`](https://github.com/ewels/MultiQC)

