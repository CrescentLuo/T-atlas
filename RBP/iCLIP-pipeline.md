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

1. Remove adapter

    Use cutadapt remove 3' adapters

2. Remove barcodes 

    This step will strip the barcode from sequencing reads and record them. In the result file you will find the barcode of every read.

    ```bash
    ./barcode.py flag_remove_adapter.fastq &
    ```

3. Mapping to genome

    ```bash
    nohup tophat2 -g 1 -a 6 --microexon-search -m 2 --bowtie1 -p 10 \
    -G /picb/rnomics2/xiaoou/data/new/human/knownGene.gtf \
    -o /picb/rnomics2/xiaoou/data/iCLIP/new2/flag_tophat \
    /picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/seqs/index/hg19_all flag_rb.fastq,../new/flag_rb.fastq >! nohup_flag_tophat.out &
    ```

4. Remove duplicates
    According to the barcodes striped from last step, remove duplicates with same barcode and mapped to the same location.

    ```bash
    ./rm_pcr.py flag_barcode_new.txt flag_tophat.bam &
    ```

5. Peak calling
    
    ```bash
    bamToBed -split -bed12 -i flag_tophat2_combine.bam > flag_tophat2_combine.bed

    perl -alne 'if($F[5] eq "+"){$s=$F[1]-1;$e=$F[1]}else{$s=$F[2];$e=$F[2]+1};$"="\t";print "$F[0]\t$s\t$e\t@F[3..5]"' flag_tophat2_combine.bed > flag_tophat2_combine_tag.bed

    genomeCoverageBed -bg -trackline -trackopts name=flag_tag -i flag_tophat2_combine_tag.bed -g /picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/chrom.size >! flag_tophat2_combine_tag.bedgraph

    sortBed -i flag_tophat2_combine_tag.bed|bgzip >! flag_tags.bed.gz

    tabix -p bed flag_tags.bed.gz

    nohup ./cluster.py -p 20 -r wgEncodeGen
    codeCompV19.bed -c flag_tags.bed.gz -o flag_all_cluster.txt &

    ./summary.py /picb/rnomics2/xiaoou/data/new/human/mart_export.txt flag_all_cluster.txt flag_all_cluster_summary.txt

    ./cluster.py -r mir302.bed -c flag_tags.bed.gz -o flag_mir302.txt

    ./cluster.py -f 1 -r pri_miRNA.bed -c flag_tags.bed.gz -o flag_miRNA.txt
    ```


-----------------------------------

##Note

Gene annotation file **refFlat.txt** is in the format ([Gene Predictions and RefSeq Genes with Gene Names](https://genome.ucsc.edu/FAQ/FAQformat.html#format9)) below (see details in [the example file](https://github.com/YangLab/TERate/blob/master/example/refFlat.txt)).

| Field       | Description                   |
| :---------: | :---------------------------- |
| geneName    | Name of gene                  |
| isoformName | Name of isoform               |
| chrom       | Reference sequence            |
| strand      | + or - for strand             |
| txStart     | Transcription start position  |
| txEnd       | Transcription end position    |
| cdsStart    | Coding region start           |
| cdsEnd      | Coding region end             |
| exonCount   | Number of exons               |
| exonStarts  | Exon start positions          |
| exonEnds    | Exon end positions            |
