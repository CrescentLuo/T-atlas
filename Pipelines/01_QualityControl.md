# Computational Analysis Pipeline

## List of softwares

| Software      | Version             |
| :-----------: | :-----------------: |
| fastqc        |                     |
| cutadapt      |                     |
| trim          |                     |


## Steps

## Key contents of quanlity control

1. Basic checking of raw reads (see table S1)
2. Decide extraction methods
3. Quality control of raw reads
4. Remove adapters
5. Trim barcodes - addititonal

## Conmands


1. Extract fastq file from SRA files.

    ```bash
    # For paired-end RNA-seq
    for i in `ls *.sra`
    do
        fastq-dump --split-files --defline-qual + --defline-seq @\$ac-\$si/\$ri $i &
    done
    ```
2. Quality control

    ```bash
    for i in `ls *.fastq`
    do
        fastqc $i
    done
    ```
3. Cut adapter 

## Supplemental materials

1. [Adapters of Illumunia RNA-seq data][1]
2. [Document of trimmomatic][2]


[1]: https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-09.pdf
[2]: http://www.usadellab.org/cms/index.php?page=trimmomatic