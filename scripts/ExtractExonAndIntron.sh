#!/bin/bash
awk -F'[,\t]' '{$1=$1; print $0}' $1
