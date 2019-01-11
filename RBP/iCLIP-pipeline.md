# iCLIP-seq pipeline

## Usage

Pipeline for iCLIP-seq data processing 

## Data

* RNA-seq: fastq file
* Reference: reference genome

## Dependency

* bamToBed
* bowtie / tophat2 / hisat2
* genomeCoverageBed

## Usage：

-----------------------------------
For short reads from iCLIP-seq, they usually have barcode at 5' end, you need remove and record these barc

## iCLIP pipeline

1. Remove adapter

    **materials:**
    * iCLIP-seq Fastq files

    **command:**
    ```bash
    cutadapt -b AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC --match-read-wildcards -m 20 -o sample1.rmAdapt.fastq sample1.fastq &> sample1.rmAdapt.log
    ```

2. Trim barcode

    **materials:**
    * Fastq files without adapters from last steps

    **command:**
    ```bash
    barcode_5.py sample1.rmAdapt.fastq
    ```

    **results:**
    * sample1.rmAdapt.rmBc.fastq

3. Mapping

    **materials:**
    * Fastq files trimmed barcoded

    **command:**
    ```bash
    bowtie -a -v 1 -m 1 --sam -p 10 $index_path/hg19_all sample1.rmAdapt.rmBc.fastq sample1.bowtie.uniq.sam &> sample1.bowtie.uniq.log
    ```

4. Remove duplicates

    **materials:**
    * bam file and barcode file

    ```bash
    rm_pcr.py sample1.rmAdapt.bc sample1.bowtie.uniq.bam &
    ```
    **results:**
    

5. Call peak

    * call peak
    * cluster.py
    * gencode V19
    * Gene Yeo's paper

6. Statisical test
7. 