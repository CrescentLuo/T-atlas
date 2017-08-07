#!/usr/bin/env bash
# Author: Crescent Luo <crescentluozheng@gmail.com>
# Description: Convert BAM files to bigWig files

function help {
    echo "bamTobw.sh -- convert stranded sequencing BAM file to bigWig file"
    echo "Usage: bamTobw.sh -b <bamlist> [-s] [-d] [-p]"
    echo "-b <bamlist> -- file contains bam files (one file per line)"
    echo "-s -- if set, related files will be scaled to HPB"
    echo "-d -- if set, bam file will be divide into strand plus and strand minus"
    echo "-p -- if set, bam file were generated from PE sequencing"
    exit 0
}

function bam_divide {
    local bam_file=$1
    local name=$2.bam
    local flag=$3
    echo "Create $name" | tee -a bamTobw.log
    if (( $flag )); then
        samtools view -f 16 -h -b -o $name $bam_file
    else
        samtools view -F 16 -h -b -o $name $bam_file
    fi
    samtools index $name
}

function bam_divide_pe {
    local bam_file=$1
    local temp=$2.temp.bam
    local name=$2.bam
    local name83=$2.flag83.bam
    local name163=$2.flag163.bam
    local name99=$2.flag99.bam
    local name147=$2.flag147.bam
    local flag=$3
    echo "Create $name" | tee -a bamTobw.log
    if (( $flag )); then
        samtools view -h -b -f83 -o $name83 $bam_file 
        samtools view -h -b -f163 -o $name163 $bamfile
        samtools merge $temp $name83 $name163
        samtools sort $temp -o $name
        rm $name83
        rm $name163
        rm $temp 
    else
        samtools view -h -b -f99 -o $name99 $bam_file 
        samtools view -h -b -f147 -o $name147 $bam_file 
        samtools merge $temp $name99 $name147
        samtools sort $temp -o $name 
        rm $name99 
        rm $name147 
        rm $temp
    fi 
    samtools index $name 
}

function bam_to_bedgraph {
    local bam_file=$1.bam
    local name=$2.bedgraph
    local chrom_sizes=$3
    local scale_flag=$4
    local flag=$5
    local count=$6
    # determine read length using a perl script borrowed from Shanshan Zhu
    local read_length=$(samtools view $bam_file | perl -lane 'print scalar(split //,$F[9]) and last if $F[5]=~/^[\dM]*$/;')
    local ratio=1
    if (( $scale_flag )); then
        ratio=`echo "scale=8;r=1000000000/$count/$read_length;if(length(r)==scale(r)) print 0;print r" | bc`
    fi
    if !(( $flag )); then
        ratio=-$ratio
    fi
    echo "Convert $bam_file to $name" | tee -a bamTobw.log
    echo "read counts: $count" | tee -a bamTobw.log
    echo "read length: $read_length" | tee -a bamTobw.log
    echo "scale ratio: $ratio" | tee -a bamTobw.log
    # round signal value with a simple perl script
    genomeCoverageBed -bg -split -ibam $bam_file -g $chrom_sizes -scale $ratio |perl -alne '$"="\t"; $F[-1]=int($F[-1]+0.5); print "@F"'> $name
}

function bedgraph_to_bw {
    local bedgraph_file=$1.bedgraph
    local name=$1.bw
    local chrom_sizes=$2
    echo "Convert $bedgraph_file to $name" | tee -a bamTobw.log
    echo "bedGraphToBigWig $bedgraph_file $chrom_sizes $name" | tee -a bamTobw.log
    bedGraphToBigWig $bedgraph_file $chrom_sizes $name
}

scale_flag=0
stranded_flag=0
paired_flag=0

while getopts ":b:sdp" optname; do
    case $optname in
        b)
            bam_list=$OPTARG;;
        s)
            scale_flag=1;;
        d)
            stranded_flag=1;;
        p)
            paired_flag=1;;
        :)
            help;;
        ?)
            help;;
    esac
done

if [[ !(-e "$bam_list") ]]; then
    help
fi

echo "Start bamTobw.sh at "`date` | tee -a bamTobw.log
echo "Parameters you input:" | tee -a bamTobw.log
echo "bam list file: $bam_list" | tee -a bamTobw.log
if (( $scale_flag )); then
    echo "convert to HPB: Yes" | tee -a bamTobw.log
else
    echo "convert to HPB: No" | tee -a bamTobw.log
fi
if (( $stranded_flag )); then
    echo "stranded: Yes" | tee -a bamTobw.log
else
    echo "stranded: No" | tee -a bamTobw.log
fi

if (( $paired_flag )); then
    echo "paired_end: Yes" | tee -a bamTobw.log
else
    echo "paired_end: No" | tee -a bamTobw.log
fi

while read line; do
    echo "Deal with $line at "`date` | tee -a bamTobw.log
    echo "Index $line" | tee -a bamTobw.log
    samtools index $line
    prefix=${line%%.bam}
    samtools idxstats $line | perl -alne 'print "$F[0]\t$F[1]" if $F[0]!~/\*/' > chrom_sizes.tmp
    chrom_sizes=chrom_sizes.tmp
    # count total reads using samtools idxstats, it reads infomation from heads of BAM files
    count=$(samtools idxstats $line | perl -ane '$a+=$F[2];END{print "$a"}')
    if (( $stranded_flag )); then
        plus_name=$prefix'_plusS'
        minus_name=$prefix'_minusS'
        echo "Divide bam into strand plus bam and strand minus bam" | tee -a bamTobw.log
        if (( $paired_flag )); then
            bam_divide_pe $line $plus_name 1
            bam_divide_pe $line $minus_name 0
        else
            bam_divide $line $plus_name 1
            bam_divide $line $minus_name 0
        fi
        echo "Convert bam to bedgraph" | tee -a bamTobw.log
        bam_to_bedgraph $plus_name $plus_name $chrom_sizes $scale_flag 1 $count
        bam_to_bedgraph $minus_name $minus_name $chrom_sizes $scale_flag 0 $count
        echo "Convert bedgraph to bw" | tee -a bamTobw.log
        bedgraph_to_bw $plus_name $chrom_sizes
        bedgraph_to_bw $minus_name $chrom_sizes
    else
        line=$prefix
        echo "Convert bam to bedgraph" | tee -a bamTobw.log
        bam_to_bedgraph $line $prefix $chrom_sizes $scale_flag 1 $count
        echo "Convert bedgraph to bw" | tee -a bamTobw.log
        bedgraph_to_bw $prefix $chrom_sizes
    fi
    rm chrom_sizes.tmp
done < $bam_list
echo "End bamTobw.sh at "`date` | tee -a bamTobw.log
