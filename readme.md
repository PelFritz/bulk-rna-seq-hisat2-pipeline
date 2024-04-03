### Bulk RNA-seq nextflow pipeline
This is a bulk-RNA-seq processing pipeline built using nextflow. To run this pipeline, one requires
nextflow and conda or docker. Please visit the nextflow, docker and anaconda websites for installation instructions. This pipeline uses hisat2, trimmomatic, samtools, featureCounts, fastqc and ncbi sratoolkit. This pipeline assumes that you reference genome, accessions.txt and .gff3 files are downloaded in the current working directory. That is the directory containing the main.nf and module.nf files. 

Pipeline can be run in 2 ways:
using docker:

#### Run with docker
Build a docker image using the dockerfile inside the docker folder. Use the following steps.

`
cd docker
`

`
docker build --no-cache -t bulk-rna-seq .
`

`
cd ..
`

To run the pipeline in normal mode which stores all intermediate files (this can quickly fill up your storage if you have several accessions).

`
nextflow run main.nf --species Arabidopsis --annot-format gff3 --accession_list accessions.txt --with-docker bulk-rna-seq --optimized false
`

Else run the pipeline in optimized mode which saves only the fastqc reports and read counts files.

`
nextflow run main.nf --species Arabidopsis --annot-format gff3 --accession_list accessions.txt--with-docker bulk-rna-seq --optimized true
`

#### Run with conda 

In case you wish to run with a conda environment and nor docker, then use the following commands.

`
conda env create --file docker/env.yml
`
The above command will create a conda environment using the env.yml file that has environment requirements for the pipeline. The environment created is called RNASeq. Once this is created, you can run nextflow with conda using the command:

For not optimized:

`
nextflow run main.nf --species Arabidopsis --annot-format gff3 --accession_list accessions.txt --with-conda RNASeq --optimized false
`

For optimized:


`
nextflow run main.nf --species Arabidopsis --annot-format gff3 --accession_list accessions.txt --with-conda RNASeq --optimized true
`

The accession.txt file is simply a csv file with two columns: the sra accession and the librarylayout. An example accession.txt can be found with the repository. These accession belong to *Arabidopsis thaliana*.

NB:

--species should match the beginning of your genome and gff3 files e.g **Arabidopsis_thaliana.TAIR10.dna.toplevel.fa** and **Arabidopsis_thaliana.TAIR10.58.gff3**