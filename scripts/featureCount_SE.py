#!/usr/bin/env python
"""

Date : November 20, 2017

Author : Heather Landry Drexler

Count reads from a bam file that map to gene features: genes, exons, and introns. The length of each feature
is also calculated as well as RPKM calculations for the whole gene, only  exons, and only introns.

This script is design for single end Illumina sequencing reads.

use : python featureCount_SE.py     iMerge (input file with 
                                    iBAM (input BAM file with only unique alignments) [1]
                                    oBAM (output BAM file containing only non duplicated reads) [2]


"""
import sys, pysam, os, numpy, re

# open files from command line arguments
iMerge = open(sys.argv[1], 'r')
iBAM = pysam.Samfile(sys.argv[2], 'rb')
oTest = open(sys.argv[3], 'w')

# save the gene name in a set if it's not already in it
gene_name = set()

# calculate million mapped reads for RPKM analysis
mmr = float(iBAM.mapped)/float(1000000)

# loop through a file with feature coordinates sorted by gene name
# record information about each feature
for line in iMerge:
    chrom = line.split('\t')[0]                     # chromosome
    start = int(line.split('\t')[1])                # start coordinate of feature
    end = int(line.split('\t')[2])                  # end coordinate of feature
    gene = line.split('\t')[3]                      # gene name
    feature = line.split('\t')[4].split('_')[0]     # feature: gene, exon, or intron
    strand = line.split('\t')[5][0]                 # strand of gene with intron
    feature_length = end-start                      # length of feature
    read_length = 80

    # check if gene name is in the set, if not save it 
    # start variables for counting the number of reads in features
    if gene not in gene_name:
        
        # print header if this is the first line
        if len(gene_name)==0:
            oTest.write('\t'.join(['gene','gene_count','exon_count','intron_count','gene_RPKM','exon_RPKM','intron_RPKM','exon_intron_ratio'])+'\n')
            
        elif gene_length <= 10 or exon_length <= 10 or intron_length <= 10:
            print(str(gene)+' thrown out')
            oTest.write('\t'.join([str(gene_output),str(gene_count),str(exon_count),str(intron_count),'NULL','NULL','NULL','NULL'])+'\n')
            
        else:
            gene_RPKM=float(gene_count)/float(gene_length)*1000 / float(mmr)
            exon_RPKM=float(exon_count)/float(exon_length)*1000 / float(mmr)
            intron_RPKM=float(intron_count)/float(intron_length)*1000 / float(mmr)
            
            if intron_RPKM==0:
                ratio="Inf"
            else:
                ratio=exon_RPKM/intron_RPKM
            oTest.write('\t'.join([str(gene_output),str(gene_count),str(exon_count),str(intron_count),str(gene_RPKM),str(exon_RPKM),str(intron_RPKM),str(ratio)])+'\n')

        # set variables for the next gene
        found = set()
        gene_name.add(gene)
        gene_count=0
        exon_count=0
        intron_count=0
        gene_length=0
        exon_length=0
        intron_length=0
        
    # count reads that span the gene
    if feature=='gene' and feature_length>10:
        gene_length+=feature_length
        for read in iBAM.fetch(chrom, int(start)+3, int(end)-3):
            if (strand=='+' and read.is_reverse
                or (strand=='-' and not read.is_reverse)):
                if (start-read_length)<read.aend<(end+read_length) and (start-read_length)<read.pos<(end+read_length):
                    key = str(read.qname)+"_gene_"+str(gene)
                    if key not in found:
                        found.add(key)
                        gene_count+=1
     
    # count the number of reads that span exons in this gene
    # the code counts per exon, but does not double count over multiple exons in the same gene
    if feature=='exon' and feature_length>10:
        exon_length+=feature_length
        for read in iBAM.fetch(chrom, int(start)+3, int(end)-3):
            if (strand=='+' and read.is_reverse
                or (strand=='-' and not read.is_reverse)):
                if (start<read.aend<end or start<read.pos<end):
                    key = str(read.qname)+"_exon_"+str(gene)
                    if key not in found:
                        found.add(key)
                        exon_count+=1

    # count the number of reads that span introns in this gene
    # the code counts per intron, but does not double count over multiple introns in the same gene
    if feature=='intron' and feature_length>10:
        intron_length+=feature_length
        for read in iBAM.fetch(chrom, int(start)+3, int(end)-3):
            if (strand=='+' and read.is_reverse
                or (strand=='-' and not read.is_reverse)):
                if (start<read.aend<end or start<read.pos<end):
                    key = str(read.qname)+"_intron_"+str(gene)
                    if key not in found:
                        found.add(key)
                        intron_count+=1   

    # store gene name before reading the next one
    gene_output=gene

iMerge.close()
iBAM.close()
oTest.close()
