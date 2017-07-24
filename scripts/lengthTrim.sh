#!/bin/bash
perl -alne 'print length $_ if $. % 4 ==2 ' $1 |sort |uniq -c |sort -k2 -nr | awk '{$3=c+=$1}1' > $1.temp
awk 'FNR==NR{sum+=$1;next}; {print $0,($3/sum)}' $1.temp{,} > ${1%.*}.lenSum
rm $1.temp