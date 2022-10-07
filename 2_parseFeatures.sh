#!/bin/bash

#### first sort the output bed file 
BED = .path/to/bedfile ### this is the bedfile that is the output of "get_context_from_vcfs.R"
bedtools sort -i $BED > ${BED/.bed/sorted.bed}

#### find the closest genebody 

A=.path/to/bedfile ### this is the sorted bed file
B=./features/genebody.sorted.bed ### bed file with genomic locations of gene bodies
OUT=${A/.bed/.GENEBODY.bed}

bedtools closest -a $A -b $B -d > $OUT

### then merge overlap
BED=$OUT #### may want to run seperately, input is output BED from closest 

bedtools merge -i $BED  -s -c 4,5,6,13 -o distinct,distinct,distinct,min > ${BED/.bed/.merged.bed}

#### find closest simple repeat 

B=./features/simpleRepeats.hg38.sorted.bed #### bed file with genomic locations of simple repeats
OUT=${A/.bed/.SIMPLEREPEAT.bed}

bedtools closest -a $A -b $B -d  > $OUT

### and merge overlap

BED = $OUT 
bedtools merge -i $BED -d -1 -s -c 4,5,6,10 -o distinct,distinct,distinct,min > ${BED/.bed/.merged.bed}

##### call transcribed or untranscribed 

FEATURE_BED=features/genebody.sorted.bed

bedtools intersect -a $A -b $FEATURE_BED  -wa -loj > ${BED/.bed/.TSB.bed}

bedtools merge -i $BED -d -1 -s -c 4,5,6,10,12 -o distinct,distinct,distinct,distinct,distinct > ${BED/.bed/.merged.bed}

##### intersect with repliseq to get replication timing 

B=./features/repliseq.sorted.bed #### bed files with replication timing per genomic region 
OUT=${A/.bed/.REPLISEQ.bed}

bedtools intersect -a $A -b $B -loj > $OUT

BED = $OUT
bedtools merge -i $BED -d -1 -s -c 4,5,6,10 -o distinct,distinct,distinct,median > ${BED/.bed/.merged.bed}





