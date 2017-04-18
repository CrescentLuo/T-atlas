#!/usr/bin/env perl
#
# Name:        /picb/extprog/biopipeline/bin/colCount.pl
# Author:      Shanshan Zhu <sszhu1007@gmail.com>
# License:     GPL
# Created:     Sun Jun. 03, 2012 18:33:36
# Last Change: Shanshan Zhu <sszhu1007@gmail.com> Mon May. 12, 2014 11:06:37
#
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $version='1';
my $command=join"\t",@ARGV;
print "Command:\t$0\t$command\n"; 
local $|=1;
####################### USAGE;
### colDel.pl
### perl -i -ne 'chop;@F=split /\t/; print join("\t",@F[0..$#F-1])."\n" ' mergePES.txt
sub usage {
    print <<"END_USAGE";
Usage: perl $0
        --input         FILENAMES 
        --noheader      default (true)
        --conditions    default (none), eg '2=~/SNP135/ and 4!~/SNP/ or 6=="Yes"'
        --columns       required, support ',' and '..'
        --combinations  default (none)
        --scombinations default (none)
        --transform     default 'no', eg 'uc', '2|uc','4|uc,5|lc'
        --ctransform    default 'no', transform for selected all cols, eg 'ifBiger 1,'
        --symbol        default, '*' for NA
        --test          print col infor
        --help
        colCount.pl -i \$mergePES 
        colCount.pl -i \$mergePESFH -com '78-79-80-81-82' 
        colCount.pl -i \$mergePES -com '78-79-80-81-82' > \$count
        colCount.pl -i \$mergePES3 -col '2' -tr 'uc' >>~/count.txt
END_USAGE
    exit;
}

####################### get options to overlap default values;
my ($input,$noheader,$conditions,$columns,$combs,$scombs,$tr,$ctr,$help,$test);
my $symbol="*";
GetOptions (
			'input=s'		=>	\$input,
            'noheader'      =>  \$noheader,
            'conditions=s'  =>  \$conditions,
            'columns=s'     =>  \$columns,
            'combinations=s'   =>  \$combs,
            'scombinations=s'     =>  \$scombs,
            'transform=s'     =>  \$tr,
            'ctransform=s'     =>  \$ctr,
            'symbol'        =>  \$symbol,
            'test'        =>  \$test,
			'help'			=>	\$help,
			) or usage();
#usage() if $help or !$input or (! $test and !$columns and ! $combs);
usage() if $help or !$input;
#$columns='58..59,66..71,73..76';
$test=1 if !$columns and ! $combs;

$conditions=~s/(\d+)(?==|!|>|<)/\$F[$1]/g if $conditions;
warn "Used conditions: $conditions" if defined $conditions;

my @cols;
eval '@cols=('.$columns.')' if $columns;
#print $columns,@cols;
my ($max)=sort{$b<=>$a} @cols;
$max=0 unless $max;
my (%tcols,@tcols);
if($tr){
    if($tr=~/\|/){
        for(split /,/,$tr){
            my($col,$var)=split /\|/;
            $tcols{$col}=$var;
        }
        push @tcols,sort keys %tcols;
    }else{
        map{$tcols{$_}=$tr} @cols;
    }
}
my (@combs,@scombs,%uniq,%comb2tag2count);
if($combs){
    for (split /,/,$combs){   ##### get combs
        next if $uniq{$_}++;
        my @vals=split /-/;
        push @combs,\@vals;
    }
}
if ($scombs){                 ##### sum combinations
    for (split /,/,$scombs){
        next if $uniq{$_}++;
        my @vals=split /-/;
        push @scombs,\@vals;
    }
}

open my $in,$input or die;
my (@titles,%count);
my $total;
while(<$in>){
    next if /^$/;
    my @F=map {chomp($_);$_} split "\t",$_; 
    @titles=@F and next if !$noheader and $.==1 ;                     #### get headers
    @titles=map {"column $_"} 0..$#F if $noheader and $.==1 ;
    if ($test){                                                       #### get information for columns
        print "$_\t$titles[$_]\n" for 0..$#titles;
        exit;
    }
    eval '$_='."$tr ".'$_' if @tcols<1 and $tr;                       #### only one transform
    @F=map {chomp($_);$_} split "\t",$_; 
    die "Too many columns are selected!!" if $#F<$max;
    #map{eval '$F[$_]='."$tcols{$_} \"$F[$_]\"";} grep {$F[$_] and $tcols{$_}} @tcols if $tr ;
    map{eval '$F[$_]='."$tcols{$_} ".'$F[$_]';} grep {$F[$_] and $tcols{$_}} @tcols if $tr ;
    map{eval '$F[$_]='."$ctr ".'$F[$_]';} grep {$F[$_]} @cols if $ctr ;
    if ($conditions){
        my $tt=0;                                                         #### check conditions
        eval '$tt=1 if '.$conditions;
        next unless $tt;
    }
    $total++;
    map{my $tt=($F[$_] and length($F[$_])>0)?$F[$_]:"*";$count{$_}{$tt}++;} @cols;                                  #### count for columns
    for(@combs,@scombs){                                                                               #### count for combinations
        my $comb=join("\t",map {$titles[$_]} @$_);
        my $tag=join("\t",map {($F[$_] and length($F[$_])>0)?$F[$_]:"*"} @$_);
        $comb2tag2count{$comb}{$tag}++;
    }
}
#print @cols;
for my $id(@cols){
    print "-->$titles[$id]:\tTerm\tCount($total)\tPercent\n";
    print "\t$_\t$count{$id}{$_}\t".$count{$id}{$_}/$total."\n" for sort keys %{$count{$id}} 
}
#for my $comb(sort keys %comb2tag2count){
for (@combs){
    my $comb=join("\t",map {$titles[$_]} @$_);
    my (%temp,@scom,@scom_key);
    map {$temp{$_}++} split /\t/,$comb;
    for (@scombs){
        my $check=grep {!$temp{$_}} map{$titles[$_]} @$_;       #### not all in comb
        next if $check;
        push @scom,[map {$titles[$_]} @$_];
        push @scom_key,join("\t",map {$titles[$_]} @$_);
    }
    my $scomb_str="";                                                 #### sum comb
    $scomb_str="\t".join("\t",map {"Percent_in_".(join "-",@$_)} @scom) if @scom;
    print "-->\t$comb\tCount($total)$scomb_str\n";
    #print "\t$_\t$comb2tag2count{$comb}{$_}\n" for sort keys %{$comb2tag2count{$comb}} 
    my @names=split /\t/,$comb;                                        #### get sum ids
    my (%name2id,@tags_ids); 
    map {$name2id{$names[$_]}=$_} 0..$#names;
    for (@scom){
        push @tags_ids,[map {$name2id{$_}}@$_];
    }
    for my $tag(sort keys %{$comb2tag2count{$comb}}){
        my @tags=split /\t/,$tag;                                      #### get sum tags
        my @stags;
        for(@tags_ids){
            push @stags,join("\t",map {$tags[$_]} @$_);
        }
        if(@stags){
            print "\t$tag\t$comb2tag2count{$comb}{$tag}\t".join("\t",map {$comb2tag2count{$comb}{$tag}/$comb2tag2count{$scom_key[$_]}{$stags[$_]}} 0..$#stags)."\n";
        }else{
            print "\t$tag\t$comb2tag2count{$comb}{$tag}\n";
        }
    }
}
sub ifBiger{
    my ($cut,$val)=@_;
    if ($val>=$cut){return 1}
    else{return 0}
}
sub Biger($){
    my $val=shift;
    if ($val>=1){return 1}
    else{return 0}
}
sub getFiles{
    my $input=$_[0];
    my (@files,%files);
    for (split /,/,$input){#### get file names based on $methy_files;
        if(/\|/){
            my $bas=basename($_);
            my $dir=dirname($_);
            for my $f(split /\|/,$bas){
                map {push @files,$_ if ! $files{$_}++} glob "$dir/$f";     ### updated by sszhu1007@gmail.com 12.10.07 12:54:13 ###
            }
        }else{ map {push @files,$_ if ! $files{$_}++} glob $_; }
    }
    #@files=sort keys %files;
    return \@files;
}
=tt
perldoc -feval
cd /picb/rnomics1/sszhu/RNA_seq/editing/merge_split_reads_16/sort_vcf/SNP135_5HPB519/merge_sites/multiple_tissues
mergeH='/picb/rnomics1/sszhu/RNA_seq/editing/merge_split_reads_16/sort_vcf/SNP135_5HPB519/merge_sites/multiple_tissues/mergeH.txt'
head $mergeH | colCount.pl -in - -col '58..59,66..71,73..75' 
head $mergeH | colCount.pl -in - -col '58..75' | head
head $mergeH | colCount.pl -in - -no -col '58..75' -tr 'uc' | head -n 20
head $mergeH | colCount.pl -in - -no -col '58..75' -tr '58|uc,59|lc' | head -n 20
colCount.pl -i tt.txt -com '29-30-36-37' -scom '29-36'    ### updated by sszhu1007@gmail.com 13.03.09 00:44:59 ###
colCount.pl -i tt.txt -com '29-30-36-37' -con '29=~/1/ and 30=~/SNP/' -scom '29-36,29-37'    ### updated by sszhu1007@gmail.com 13.03.09 01:48:48 ###
colCount.pl -i tt.txt -col '29,30'
# vim:fdm=marker

